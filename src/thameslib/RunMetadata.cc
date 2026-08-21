/**
@file RunMetadata.cc
@brief Implementation of the runmeta:: namespace declared in RunMetadata.h.

Design notes:
  - All state is kept in one file-scope Metadata struct plus an initialized_
    flag. THAMES runs one hydration per process, so a global is appropriate
    and avoids threading a metadata pointer through every CSV writer.
  - JSON is emitted by hand (no vendored JSON serializer). The schema is
    small and stable; escaping is limited to `"` and `\\` in the few
    free-text fields we ever store.
  - simparams.json is parsed with simdjson (already used everywhere else in
    THAMES-Hydration for JSON reading — see ChemicalSystem.cc).
  - SHA-256 uses picosha2 (single-header MIT library vendored alongside).
*/

#include "RunMetadata.h"
#include "picosha2.h"
#include "../version.h"

// Fallback identity macros for editor / partial-build contexts where CMake
// has not yet regenerated version.h from version.h.in. A real build always
// overrides these via configure_file — see backend/thames-hydration/CMakeLists.txt.
#ifndef THAMES_VERSION_STRING
#  define THAMES_VERSION_STRING "unknown"
#endif
#ifndef THAMES_GIT_HASH
#  define THAMES_GIT_HASH "unknown"
#endif
#ifndef THAMES_BUILD_DATE
#  define THAMES_BUILD_DATE "unknown"
#endif
#ifndef THAMES_COMPILER_ID
#  define THAMES_COMPILER_ID "unknown"
#endif
#ifndef THAMES_COMPILER_VERSION
#  define THAMES_COMPILER_VERSION "unknown"
#endif
#ifndef THAMES_CXX_FLAGS
#  define THAMES_CXX_FLAGS "unknown"
#endif
#ifndef THAMES_CMAKE_VERSION
#  define THAMES_CMAKE_VERSION "unknown"
#endif
#ifndef THAMES_BUILD_TYPE
#  define THAMES_BUILD_TYPE "unknown"
#endif

#include <chrono>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>
#include <unistd.h>

#if defined(_WIN32)
#  include <winsock2.h>          // gethostname on Windows
#else
#  include <sys/utsname.h>
#endif

#include "../GEMS3K-standalone/simdjson/simdjson.h"

namespace fs = std::filesystem;

namespace runmeta {

namespace {

struct Metadata {
  // Backend identity (compile-time via version.h).
  std::string thames_backend_version = THAMES_VERSION_STRING;
  std::string thames_backend_git_hash = THAMES_GIT_HASH;
  std::string thames_backend_build_date = THAMES_BUILD_DATE;

  // Build tooling — added schema_version 2. Diagnostic-only; used when a
  // cross-host numerical divergence investigation needs to rule the
  // compiler / optimization / flag set in or out.
  std::string build_compiler_id = THAMES_COMPILER_ID;
  std::string build_compiler_version = THAMES_COMPILER_VERSION;
  std::string build_cxx_flags = THAMES_CXX_FLAGS;
  std::string build_cmake_version = THAMES_CMAKE_VERSION;
  std::string build_type = THAMES_BUILD_TYPE;

  // UI identity (from simparams.json ui_context block; empty when standalone).
  std::string thames_ui_version;
  std::string python_version;

  // Platform / host (either from ui_context or captured at runtime).
  std::string platform_system;
  std::string platform_release;
  std::string platform_machine;
  std::string hostname;
  bool hostname_from_ui = false;

  // DCH fingerprint.
  std::string dch_path;
  std::string dch_sha256;
  long long   dch_size_bytes = 0;

  // Run identity.
  std::string operation_name;
  std::string run_started_utc;
  std::string run_finished_utc;
  std::string exit_reason = "in_progress";
  int         exit_code = 0;
  std::string diagnostics;

  // Output folder we own.
  std::string output_folder;
};

Metadata g;
bool g_initialized = false;
bool g_finalized = false;   // set true by first finalize call; subsequent
                            // calls become no-ops so specific diagnostic
                            // reasons (from catch blocks in thames.cc)
                            // aren't overwritten by generic fallback
                            // reasons (from deleteDynAllocMem).

std::string iso8601UtcNow() {
  const auto now = std::chrono::system_clock::now();
  const auto tt = std::chrono::system_clock::to_time_t(now);
  std::tm tm{};
#if defined(_WIN32)
  gmtime_s(&tm, &tt);
#else
  gmtime_r(&tt, &tm);
#endif
  char buf[32];
  std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
  return buf;
}

std::string captureHostname() {
  char buf[256];
  buf[0] = '\0';
  if (::gethostname(buf, sizeof(buf) - 1) == 0) {
    buf[sizeof(buf) - 1] = '\0';
    return std::string(buf);
  }
  return "unknown";
}

// uname() gives us system/release/machine on POSIX. On Windows we set a
// short static fallback — the UI can (and normally does) supply richer
// platform info via ui_context, so this branch is best-effort.
void capturePlatformFallback(Metadata &m) {
#if defined(_WIN32)
  m.platform_system = "Windows";
  m.platform_release = "";
  m.platform_machine = "";
#else
  struct utsname u{};
  if (::uname(&u) == 0) {
    m.platform_system = u.sysname;
    m.platform_release = u.release;
    m.platform_machine = u.machine;
  }
#endif
}

// Trim trailing whitespace + CR from a string. Windows-authored input.in
// files reach us as CRLF; getline() on POSIX keeps the trailing \r, which
// would otherwise turn "thames-dat.lst" into "thames-dat.lst\r" — a
// filename that doesn't exist.
std::string rstrip(std::string s) {
  while (!s.empty() && (s.back() == '\r' || s.back() == '\n' ||
                        s.back() == ' ' || s.back() == '\t')) {
    s.pop_back();
  }
  return s;
}

// Extract the first quoted filename from the GEMS input list. The file
// format is:
//   -t "thames-dch.dat" "thames-ipm.dat" "thames-dbr.dat"
// (all on one line, or one entry per line). Matches the parse thames.cc
// does inline in prepOutputFolder.
std::string extractDchName(const std::string &gemInputList) {
  const std::string path = rstrip(gemInputList);
  std::ifstream in(path);
  if (!in.is_open()) return "";
  std::string buff;
  in >> buff;   // discard flag ("-t" or similar)
  buff.clear();
  in >> buff;   // dch filename, possibly quoted
  if (!buff.empty() && (buff.front() == '"' || buff.front() == '\'')) {
    // Strip both surrounding quotes if present.
    buff = buff.substr(1);
    if (!buff.empty() && (buff.back() == '"' || buff.back() == '\'')) {
      buff.pop_back();
    }
  }
  return rstrip(buff);
}

// Compute SHA-256 hex + size of a file. Returns empty hash on failure.
void hashDchFile(const std::string &path, std::string &hexOut,
                 long long &sizeOut) {
  hexOut.clear();
  sizeOut = 0;
  std::ifstream in(path, std::ios::binary);
  if (!in.is_open()) return;
  std::vector<unsigned char> buf(std::istreambuf_iterator<char>(in), {});
  sizeOut = static_cast<long long>(buf.size());
  hexOut = picosha2::hash256_hex_string(buf);
}

// Read `ui_context` block from simparams.json if present. All fields are
// optional; missing fields leave the struct field untouched.
void readUiContext(const std::string &simParamsPath, Metadata &m) {
  simdjson::dom::parser parser;
  simdjson::dom::element root;
  auto err = parser.load(simParamsPath).get(root);
  if (err) return;

  simdjson::dom::element ui;
  if (root["ui_context"].get(ui) != simdjson::SUCCESS) return;

  auto tryString = [&](const char *key, std::string &dst) {
    std::string_view v;
    if (ui[key].get(v) == simdjson::SUCCESS) dst = std::string(v);
  };
  tryString("ui_version", m.thames_ui_version);
  tryString("python_version", m.python_version);
  tryString("platform_system", m.platform_system);
  tryString("platform_release", m.platform_release);
  tryString("platform_machine", m.platform_machine);
  tryString("operation_name", m.operation_name);

  std::string_view host;
  if (ui["hostname"].get(host) == simdjson::SUCCESS) {
    m.hostname = std::string(host);
    m.hostname_from_ui = true;   // UI-supplied honors the privacy toggle
  }
}

// JSON string escape for the small subset of characters we might emit.
std::string jsonEscape(const std::string &in) {
  std::string out;
  out.reserve(in.size() + 8);
  for (char c : in) {
    switch (c) {
      case '"':  out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n";  break;
      case '\r': out += "\\r";  break;
      case '\t': out += "\\t";  break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char esc[8];
          std::snprintf(esc, sizeof(esc), "\\u%04x", c);
          out += esc;
        } else {
          out += c;
        }
    }
  }
  return out;
}

std::string sidecarPath() {
  return (fs::path(g.output_folder) / "run_metadata.json").string();
}

// Serialize the metadata block to the sidecar. Writes to a temp file then
// renames — cheap atomicity on same filesystem, avoids torn reads if the UI
// polls the file while the backend rewrites it.
void writeSidecar() {
  const std::string finalPath = sidecarPath();
  const std::string tmpPath = finalPath + ".tmp";

  std::ofstream ofs(tmpPath);
  if (!ofs) {
    std::clog << "WARNING: RunMetadata could not write " << tmpPath
              << std::endl;
    return;
  }

  ofs << "{\n"
      << "  \"schema_version\": 2,\n"
      << "  \"operation_name\": \"" << jsonEscape(g.operation_name) << "\",\n"
      << "  \"thames_backend\": {\n"
      << "    \"version\": \"" << jsonEscape(g.thames_backend_version) << "\",\n"
      << "    \"git_hash\": \"" << jsonEscape(g.thames_backend_git_hash) << "\",\n"
      << "    \"build_date\": \"" << jsonEscape(g.thames_backend_build_date) << "\"\n"
      << "  },\n"
      << "  \"build\": {\n"
      << "    \"compiler_id\": \"" << jsonEscape(g.build_compiler_id) << "\",\n"
      << "    \"compiler_version\": \"" << jsonEscape(g.build_compiler_version) << "\",\n"
      << "    \"cxx_flags\": \"" << jsonEscape(g.build_cxx_flags) << "\",\n"
      << "    \"cmake_version\": \"" << jsonEscape(g.build_cmake_version) << "\",\n"
      << "    \"build_type\": \"" << jsonEscape(g.build_type) << "\"\n"
      << "  },\n"
      << "  \"thames_ui\": {\n"
      << "    \"version\": \"" << jsonEscape(g.thames_ui_version) << "\",\n"
      << "    \"python_version\": \"" << jsonEscape(g.python_version) << "\"\n"
      << "  },\n"
      << "  \"platform\": {\n"
      << "    \"system\": \"" << jsonEscape(g.platform_system) << "\",\n"
      << "    \"release\": \"" << jsonEscape(g.platform_release) << "\",\n"
      << "    \"machine\": \"" << jsonEscape(g.platform_machine) << "\",\n"
      << "    \"hostname\": \"" << jsonEscape(g.hostname) << "\"\n"
      << "  },\n"
      << "  \"gems_dch\": {\n"
      << "    \"path\": \"" << jsonEscape(g.dch_path) << "\",\n"
      << "    \"sha256\": \"" << jsonEscape(g.dch_sha256) << "\",\n"
      << "    \"size_bytes\": " << g.dch_size_bytes << "\n"
      << "  },\n"
      << "  \"run\": {\n"
      << "    \"started_utc\": \"" << jsonEscape(g.run_started_utc) << "\",\n"
      << "    \"finished_utc\": \"" << jsonEscape(g.run_finished_utc) << "\",\n"
      << "    \"exit_code\": " << g.exit_code << ",\n"
      << "    \"exit_reason\": \"" << jsonEscape(g.exit_reason) << "\",\n"
      << "    \"success\": " << ((g.exit_code == 0 && g.exit_reason != "in_progress") ? "true" : "false") << ",\n"
      << "    \"diagnostics\": \"" << jsonEscape(g.diagnostics) << "\"\n"
      << "  }\n"
      << "}\n";
  ofs.close();

  std::error_code ec;
  fs::rename(tmpPath, finalPath, ec);
  if (ec) {
    // Rename failed (e.g. cross-device on some Windows setups). Fall back
    // to copy+remove; the sidecar still lands, just not atomically.
    fs::copy_file(tmpPath, finalPath,
                  fs::copy_options::overwrite_existing, ec);
    fs::remove(tmpPath);
  }
}

}  // namespace

void initialize(const std::string &outputFolder,
                const std::string &gemInputList,
                const std::string &simParamsPath) {
  // Strip trailing CR/whitespace from every argument. thames.cc gets
  // gemInputList and simParamsPath via getline(cin, ...) which keeps the
  // trailing \r on Windows-authored input.in files (CRLF line endings).
  // Without this, opening any of these paths on POSIX fails silently.
  const std::string outputFolderClean = rstrip(outputFolder);
  const std::string gemInputListClean = rstrip(gemInputList);
  const std::string simParamsPathClean = rstrip(simParamsPath);

  g.output_folder = outputFolderClean;

  // Safety net for the elastic-only path, which bypasses prepOutputFolder
  // and would otherwise silently fail to write the sidecar if OutputFolder
  // doesn't already exist. mkdir-p semantics: no-op if it already exists.
  std::error_code ec;
  fs::create_directories(outputFolderClean, ec);

  // Runtime platform / hostname capture. UI values overwrite in
  // readUiContext below if simparams supplies them.
  capturePlatformFallback(g);
  g.hostname = captureHostname();
  g.hostname_from_ui = false;

  // UI context (may overwrite platform_* and hostname).
  readUiContext(simParamsPathClean, g);

  // DCH fingerprint.
  g.dch_path = extractDchName(gemInputListClean);
  if (!g.dch_path.empty()) {
    hashDchFile(g.dch_path, g.dch_sha256, g.dch_size_bytes);
  }

  g.run_started_utc = iso8601UtcNow();
  g.exit_reason = "in_progress";
  g.exit_code = 0;

  g_initialized = true;
  writeSidecar();

  std::clog << "\nRunMetadata initialized:\n"
            << "  backend   = " << g.thames_backend_version
            << " (" << g.thames_backend_git_hash << ")\n"
            << "  ui        = " << (g.thames_ui_version.empty() ? "(standalone)" : g.thames_ui_version) << "\n"
            << "  platform  = " << g.platform_system;
  if (!g.platform_release.empty()) std::clog << " " << g.platform_release;
  std::clog << "\n"
            << "  host      = " << g.hostname << "\n"
            << "  dch       = " << g.dch_path << " (" << g.dch_size_bytes
            << " bytes, sha256=" << (g.dch_sha256.empty() ? "unknown"
                                                          : g.dch_sha256.substr(0, 8))
            << "...)\n"
            << "  sidecar   = " << sidecarPath() << std::endl;
}

void finalize(int exitCode,
              const std::string &exitReason,
              const std::string &diagnostics) {
  if (!g_initialized) return;
  if (g_finalized) {
    // Preserve the first (most specific) finalize reason. A common pattern:
    // Controller::doCycle calls finalize with the true termination cause,
    // then deleteDynAllocMem calls finalize again with a generic fallback
    // reason. First-call-wins keeps the specific one.
    std::clog << "\nRunMetadata: ignoring redundant finalize (already "
              << "finalized as \"" << g.exit_reason << "\")" << std::endl;
    return;
  }
  g.exit_code = exitCode;
  g.exit_reason = exitReason.empty() ? "unknown" : exitReason;
  g.diagnostics = diagnostics;
  g.run_finished_utc = iso8601UtcNow();
  g_finalized = true;
  writeSidecar();
  std::clog << "\nRunMetadata finalized: " << g.exit_reason
            << " (exit_code=" << g.exit_code << ")" << std::endl;
}

std::string csvCommentLine() {
  std::ostringstream oss;
  oss << "# thames_backend=" << (g_initialized ? g.thames_backend_version : std::string("unknown"))
      << " build=" << (g_initialized ? g.thames_backend_git_hash : std::string("unknown"))
      << " platform=" << (g_initialized && !g.platform_system.empty()
                          ? g.platform_system : std::string("unknown"))
      << " ui=" << (g_initialized && !g.thames_ui_version.empty()
                    ? g.thames_ui_version : std::string("standalone"))
      << " run_started=" << (g_initialized ? g.run_started_utc : iso8601UtcNow())
      << " dch_sha256=" << (g_initialized && !g.dch_sha256.empty()
                            ? g.dch_sha256.substr(0, 16) : std::string("unknown"));
  return oss.str();
}

bool isInitialized() { return g_initialized; }

}  // namespace runmeta
