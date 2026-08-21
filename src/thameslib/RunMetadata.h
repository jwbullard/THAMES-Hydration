/**
@file RunMetadata.h
@brief Provenance stamp for every THAMES-Hydration output.

Populates a single per-run metadata block that identifies:
  - the exact backend build (version + short git hash + build date)
  - the UI that launched it (version + platform + hostname), if supplied
  - the platform running the backend
  - a SHA-256 fingerprint of the GEMS3K DCH thermodynamic database
  - run start/finish times and exit status

The block is serialized to `<OutputFolder>/run_metadata.json` at run start
(with `exit_status = "in_progress"`) and rewritten on completion with the
real exit status. Every CSV writer in Controller.cc also emits a one-line
`#`-prefixed comment header carrying a shorthand of the same identity so
individual files remain traceable after separation from their folder.

`ui_context` fields (thames_ui_version, python_version, platform_*,
hostname, operation_name) are optionally read from the `ui_context` block
of simparams.json — the UI writes them there at operation launch. When
running the backend standalone (no UI) these fields stay empty.

@sa docs/NIST-diagnostic.md
*/

#ifndef SRC_THAMESLIB_RUNMETADATA_H_
#define SRC_THAMESLIB_RUNMETADATA_H_

#include <string>

namespace runmeta {

/// One-shot initializer, called from thames.cc after prepOutputFolder.
///
/// @param outputFolder Directory that will hold Result files (typically "Result").
/// @param gemInputList Path to the GEMS input list file (e.g. "thames-dat.lst");
///                     the first quoted entry is the DCH filename to hash.
/// @param simParamsPath Path to simparams.json (may contain a "ui_context" block).
///
/// Captures backend identity (version + git hash + build date), platform
/// info via uname/gethostname, hashes the DCH file, records the start
/// timestamp, and writes the sidecar with exit_status="in_progress".
void initialize(const std::string &outputFolder,
                const std::string &gemInputList,
                const std::string &simParamsPath);

/// Called from Controller after each hydration completes (or fails).
/// Updates the sidecar with finish timestamp and the exit reason. Also
/// used by the failure path in Controller::doCycle. Safe to call more
/// than once — last write wins.
///
/// @param exitCode 0 on success, non-zero on failure (mirrors old exit_status.json).
/// @param exitReason Short human-readable phrase (e.g. "FINAL_TIME", "IC_DEPLETION").
/// @param diagnostics Long-form context; may be empty.
void finalize(int exitCode,
              const std::string &exitReason,
              const std::string &diagnostics);

/// Compact identity string for CSV comment headers. One line, no trailing
/// newline. Format:
///   # thames_backend=X.Y.Z build=abc12345 platform=darwin ui=1.0.0-alpha.2 \
///     run_started=2026-08-21T14:33:19Z dch_sha256=e3b0c442
std::string csvCommentLine();

/// True if initialize() has been called successfully. Guards optional
/// consumers (e.g. csvCommentLine returns "# thames_backend=unknown ..." if
/// called too early — never crashes).
bool isInitialized();

}  // namespace runmeta

#endif  // SRC_THAMESLIB_RUNMETADATA_H_
