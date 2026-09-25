// Standalone unit test for the periodic Euclidean distance transform.
// Compile with build_and_run.sh — no THAMES dependencies needed.

#include "DistanceTransform.h"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

int gFailed = 0;

void expect(bool cond, const char *tag) {
  std::printf("%-62s %s\n", tag, cond ? "PASS" : "FAIL");
  if (!cond) ++gFailed;
}

void expectNear(double actual, double expected, double tol, const char *tag) {
  const double diff = std::fabs(actual - expected);
  const bool ok = diff <= tol;
  std::printf("%-62s actual=%.6f expected=%.6f %s\n", tag, actual, expected,
              ok ? "PASS" : "FAIL");
  if (!ok) ++gFailed;
}

inline int idx(int ix, int iy, int iz, int nx, int ny) {
  return ((iz * ny) + iy) * nx + ix;
}

// Brute-force squared distance to the nearest background voxel, minimum image
// convention. Slow, but obviously correct — the point of comparison.
double bruteForce(int ix, int iy, int iz, int nx, int ny, int nz,
                  const std::vector<char> &fg) {
  double best = -1.0;
  for (int kz = 0; kz < nz; ++kz) {
    for (int ky = 0; ky < ny; ++ky) {
      for (int kx = 0; kx < nx; ++kx) {
        if (fg[idx(kx, ky, kz, nx, ny)])
          continue;
        int dx = std::abs(kx - ix);
        if (dx > nx - dx) dx = nx - dx;
        int dy = std::abs(ky - iy);
        if (dy > ny - dy) dy = ny - dy;
        int dz = std::abs(kz - iz);
        if (dz > nz - dz) dz = nz - dz;
        const double d2 = static_cast<double>(dx) * dx +
                          static_cast<double>(dy) * dy +
                          static_cast<double>(dz) * dz;
        if (best < 0.0 || d2 < best)
          best = d2;
      }
    }
  }
  return best;
}

// --- Test group 1: background voxels measure zero ---
void testBackgroundIsZero() {
  const int n = 6;
  std::vector<char> fg(n * n * n, 1);
  fg[idx(2, 3, 4, n, n)] = 0;

  const std::vector<double> d = edt::squaredDistanceToBackground(n, n, n, fg);
  expectNear(d[idx(2, 3, 4, n, n)], 0.0, 1e-12, "background voxel: zero");
  expectNear(d[idx(3, 3, 4, n, n)], 1.0, 1e-12, "its face neighbour: 1");
  expectNear(d[idx(3, 4, 4, n, n)], 2.0, 1e-12, "its edge neighbour: 2");
  expectNear(d[idx(3, 4, 5, n, n)], 3.0, 1e-12, "its corner neighbour: 3");
}

// --- Test group 2: distance wraps across the seam ---
void testPeriodicity() {
  const int n = 8;
  std::vector<char> fg(n * n * n, 1);

  // A background plane at ix == 0. The voxel at ix == n-1 is one step away
  // going the short way round, not n-1 steps.
  for (int iz = 0; iz < n; ++iz)
    for (int iy = 0; iy < n; ++iy)
      fg[idx(0, iy, iz, n, n)] = 0;

  const std::vector<double> d = edt::squaredDistanceToBackground(n, n, n, fg);
  expectNear(d[idx(n - 1, 3, 3, n, n)], 1.0, 1e-12,
             "across the seam: 1, not (n-1)^2");
  expectNear(d[idx(n / 2, 3, 3, n, n)], 16.0, 1e-12,
             "mid-slab: half the box away");
}

// --- Test group 3: agrees with brute force on a scattered set ---
void testAgainstBruteForce() {
  const int nx = 7, ny = 6, nz = 5;
  std::vector<char> fg(nx * ny * nz, 1);

  // A handful of background voxels in no particular pattern, including ones
  // on the faces so the wrap gets exercised from both sides.
  const int seeds[][3] = {{0, 0, 0}, {6, 5, 4}, {3, 2, 1},
                          {1, 4, 3}, {5, 0, 2}, {2, 5, 0}};
  for (const auto &s : seeds)
    fg[idx(s[0], s[1], s[2], nx, ny)] = 0;

  const std::vector<double> d =
      edt::squaredDistanceToBackground(nx, ny, nz, fg);

  double worst = 0.0;
  for (int iz = 0; iz < nz; ++iz) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int ix = 0; ix < nx; ++ix) {
        const double want = bruteForce(ix, iy, iz, nx, ny, nz, fg);
        const double got = d[idx(ix, iy, iz, nx, ny)];
        const double diff = std::fabs(got - want);
        if (diff > worst)
          worst = diff;
      }
    }
  }
  expectNear(worst, 0.0, 1e-9, "scattered background: exact everywhere");
}

// --- Test group 4: a pore geometry, which is the real use ---
void testPoreDepthOrdering() {
  const int n = 12;
  std::vector<char> fg(n * n * n, 0); // solid everywhere

  // A large spherical pore and a small one. The center of the large pore must
  // rank deeper than every voxel of the small one, which is exactly what the
  // 7^3 window could not resolve.
  const int bigC = 3, smallC = 9;
  for (int iz = 0; iz < n; ++iz)
    for (int iy = 0; iy < n; ++iy)
      for (int ix = 0; ix < n; ++ix) {
        const double rb = (ix - bigC) * (ix - bigC) + (iy - bigC) * (iy - bigC) +
                          (iz - bigC) * (iz - bigC);
        const double rs = (ix - smallC) * (ix - smallC) +
                          (iy - smallC) * (iy - smallC) +
                          (iz - smallC) * (iz - smallC);
        if (rb <= 9.0 || rs <= 1.0)
          fg[idx(ix, iy, iz, n, n)] = 1;
      }

  const std::vector<double> d = edt::squaredDistanceToBackground(n, n, n, fg);

  const double bigDepth = d[idx(bigC, bigC, bigC, n, n)];
  const double smallDepth = d[idx(smallC, smallC, smallC, n, n)];
  expect(bigDepth > smallDepth, "large pore center ranks deeper than small");
  // The small "sphere" is the 7-voxel plus shape, so the nearest solid to its
  // center is the face diagonal, not the voxel two steps along an axis.
  expectNear(smallDepth, 2.0, 1e-12, "small pore center: face diagonal away");

  // The deepest voxel in the whole box is the large pore's center: that is
  // where the first cavity nucleates.
  int argmax = 0;
  for (int i = 1; i < n * n * n; ++i)
    if (d[i] > d[argmax])
      argmax = i;
  expect(argmax == idx(bigC, bigC, bigC, n, n),
         "deepest voxel is the large pore center");
}

// --- Test group 5: degenerate inputs ---
void testDegenerate() {
  const int n = 4;

  {
    const std::vector<char> fg(n * n * n, 0);
    const std::vector<double> d = edt::squaredDistanceToBackground(n, n, n, fg);
    expect(d.size() == static_cast<std::size_t>(n * n * n),
           "all background: right size");
    expectNear(d[5], 0.0, 1e-12, "all background: all zero");
  }

  {
    // No background anywhere: nothing to measure to, so everything comes back
    // beyond any real distance rather than garbage.
    const std::vector<char> fg(n * n * n, 1);
    const std::vector<double> d = edt::squaredDistanceToBackground(n, n, n, fg);
    expect(d[0] > 3.0 * n * n, "no background: beyond any real distance");
  }

  {
    const std::vector<char> fg(n * n * n - 1, 1);
    const std::vector<double> d = edt::squaredDistanceToBackground(n, n, n, fg);
    expect(d.empty(), "size mismatch: rejected safely");
  }
}

} // namespace

int main() {
  std::printf("=== EDT: background voxels ===\n");
  testBackgroundIsZero();
  std::printf("\n=== EDT: periodicity ===\n");
  testPeriodicity();
  std::printf("\n=== EDT: against brute force ===\n");
  testAgainstBruteForce();
  std::printf("\n=== EDT: pore depth ordering ===\n");
  testPoreDepthOrdering();
  std::printf("\n=== EDT: degenerate inputs ===\n");
  testDegenerate();

  std::printf("\n%s (%d failure%s)\n", gFailed == 0 ? "ALL PASS" : "FAILURES",
              gFailed, gFailed == 1 ? "" : "s");
  return (gFailed == 0) ? 0 : 1;
}
