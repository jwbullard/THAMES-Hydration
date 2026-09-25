// Standalone unit test for percolation assessment.
// Compile with build_and_run.sh — no THAMES dependencies needed.

#include "Percolation.h"

#include <cstdio>
#include <vector>

namespace {

int gFailed = 0;

void expect(bool cond, const char *tag) {
  std::printf("%-62s %s\n", tag, cond ? "PASS" : "FAIL");
  if (!cond) ++gFailed;
}

void expectNear(double actual, double expected, double tol, const char *tag) {
  const double diff = (actual > expected) ? actual - expected : expected - actual;
  const bool ok = diff <= tol;
  std::printf("%-62s actual=%.6f expected=%.6f %s\n", tag, actual, expected,
              ok ? "PASS" : "FAIL");
  if (!ok) ++gFailed;
}

inline int idx(int ix, int iy, int iz, int nx, int ny) {
  return ((iz * ny) + iy) * nx + ix;
}

// Everything bonds on contact: plain geometric connectivity, which is what
// capillary percolation uses.
std::vector<char> allBonding(int n) { return std::vector<char>(n, 1); }

// --- Test group 1: a bar spanning one axis only ---
void testBarSpansOneAxis() {
  const int n = 8;
  const int numSites = n * n * n;
  std::vector<char> participates(numSites, 0);

  // A single line of voxels along x at (y, z) = (3, 3).
  for (int ix = 0; ix < n; ++ix)
    participates[idx(ix, 3, 3, n, n)] = 1;

  const percolation::Result r =
      percolation::assess(n, n, n, participates, allBonding(numSites),
                          std::vector<int>());

  expect(r.numParticipating == n, "bar: participating count");
  expect(r.direction[0].spans, "bar: spans x");
  expect(!r.direction[1].spans, "bar: does not span y");
  expect(!r.direction[2].spans, "bar: does not span z");
  expect(!r.spansAllDirections(), "bar: not set (needs all three)");
  expectNear(r.direction[0].connectedFraction, 1.0, 1e-12,
             "bar: x connected fraction");
  expectNear(r.minConnectedFraction(), 0.0, 1e-12,
             "bar: min connected fraction is zero");
}

// --- Test group 2: touching grains must not connect ---
void testTouchingGrainsDoNotConnect() {
  const int n = 8;
  const int numSites = n * n * n;
  std::vector<char> participates(numSites, 0);
  std::vector<char> bonds(numSites, 0); // nothing grew in place
  std::vector<int> particleId(numSites, 1); // 1 == no particle

  // Two half-bars along x that touch at ix == 3 | 4, belonging to DIFFERENT
  // original grains. Physically these are two flocculated particles in
  // contact, which must not carry load across the contact.
  for (int ix = 0; ix < n; ++ix) {
    const int site = idx(ix, 3, 3, n, n);
    participates[site] = 1;
    particleId[site] = (ix < n / 2) ? 7 : 9;
  }

  percolation::Result r =
      percolation::assess(n, n, n, participates, bonds, particleId);
  expect(!r.direction[0].spans, "two grains touching: does NOT span x");
  expectNear(r.direction[0].connectedFraction, 0.0, 1e-12,
             "two grains touching: zero connected fraction");

  // Same geometry, one grain: now it is a single particle interior and spans.
  for (int ix = 0; ix < n; ++ix)
    particleId[idx(ix, 3, 3, n, n)] = 7;
  r = percolation::assess(n, n, n, participates, bonds, particleId);
  expect(r.direction[0].spans, "one grain: spans x");

  // Two grains again, but a hydrate voxel that grew in place bridges them.
  for (int ix = 0; ix < n; ++ix)
    particleId[idx(ix, 3, 3, n, n)] = (ix < n / 2) ? 7 : 9;
  bonds[idx(n / 2, 3, 3, n, n)] = 1; // this voxel adheres to anything
  r = percolation::assess(n, n, n, participates, bonds, particleId);
  expect(r.direction[0].spans, "bridged by in-place growth: spans x");
}

// --- Test group 3: periodicity handled per axis ---
void testPeriodicity() {
  const int n = 8;
  const int numSites = n * n * n;

  // (a) A bar along x that is one voxel SHORT of spanning. Under a fully
  // periodic labeling its two ends would meet across the seam and it would
  // falsely percolate; with the spanning axis non-periodic it must not.
  {
    std::vector<char> participates(numSites, 0);
    for (int ix = 0; ix < n - 1; ++ix)
      participates[idx(ix, 3, 3, n, n)] = 1;
    const percolation::Result r =
        percolation::assess(n, n, n, participates, allBonding(numSites),
                            std::vector<int>());
    expect(!r.direction[0].spans,
           "gap at far face: does NOT span x (seam must not close it)");
  }

  // (b) A path that spans x while wrapping across the y seam. The transverse
  // axes stay periodic, so this must count: the microstructure is an RVE.
  {
    std::vector<char> participates(numSites, 0);
    for (int ix = 0; ix < n / 2; ++ix)
      participates[idx(ix, 0, 3, n, n)] = 1;            // first leg at y = 0
    participates[idx(n / 2 - 1, n - 1, 3, n, n)] = 1;   // step across the seam
    for (int ix = n / 2 - 1; ix < n; ++ix)
      participates[idx(ix, n - 1, 3, n, n)] = 1;        // second leg at y = n-1
    const percolation::Result r =
        percolation::assess(n, n, n, participates, allBonding(numSites),
                            std::vector<int>());
    expect(r.direction[0].spans,
           "wrapping in y while spanning x: DOES span x");
  }
}

// --- Test group 4: degenerate inputs ---
void testDegenerate() {
  const int n = 4;
  const int numSites = n * n * n;

  {
    const std::vector<char> participates(numSites, 0);
    const percolation::Result r =
        percolation::assess(n, n, n, participates, allBonding(numSites),
                            std::vector<int>());
    expect(r.numParticipating == 0, "empty set: nothing participates");
    expect(!r.spansAllDirections(), "empty set: does not span");
    expectNear(r.minConnectedFraction(), 0.0, 1e-12,
               "empty set: zero connected fraction");
  }

  {
    const std::vector<char> participates(numSites, 1);
    const percolation::Result r =
        percolation::assess(n, n, n, participates, allBonding(numSites),
                            std::vector<int>());
    expect(r.spansAllDirections(), "full box: spans all three");
    expectNear(r.minConnectedFraction(), 1.0, 1e-12,
               "full box: connected fraction 1");
  }

  {
    // Mismatched array size must be rejected rather than read out of bounds.
    const std::vector<char> participates(numSites - 1, 1);
    const percolation::Result r =
        percolation::assess(n, n, n, participates, allBonding(numSites),
                            std::vector<int>());
    expect(r.numParticipating == 0 && !r.spansAllDirections(),
           "size mismatch: rejected safely");
  }
}

// --- Test group 5: the capillary-percolation use, and its loss ---
void testCapillaryDepercolation() {
  const int n = 8;
  const int numSites = n * n * n;

  // A slab of pore voxels spanning every direction, as at early age.
  std::vector<char> participates(numSites, 0);
  for (int iz = 0; iz < n; ++iz)
    for (int iy = 0; iy < n; ++iy)
      for (int ix = 0; ix < n; ++ix)
        if (iy < 2)
          participates[idx(ix, iy, iz, n, n)] = 1;

  percolation::Result r =
      percolation::assess(n, n, n, participates, allBonding(numSites),
                          std::vector<int>());
  expect(r.direction[0].spans && r.direction[2].spans,
         "capillary slab: spans x and z");

  // Products seal one plane of it: the path along x is cut everywhere.
  for (int iz = 0; iz < n; ++iz)
    for (int iy = 0; iy < 2; ++iy)
      participates[idx(n / 2, iy, iz, n, n)] = 0;

  r = percolation::assess(n, n, n, participates, allBonding(numSites),
                          std::vector<int>());
  expect(!r.direction[0].spans, "after sealing a plane: x depercolates");
  expect(r.direction[2].spans, "after sealing a plane: z still spans");
}

} // namespace

int main() {
  std::printf("=== Percolation: bar spans one axis ===\n");
  testBarSpansOneAxis();
  std::printf("\n=== Percolation: touching grains ===\n");
  testTouchingGrainsDoNotConnect();
  std::printf("\n=== Percolation: periodicity ===\n");
  testPeriodicity();
  std::printf("\n=== Percolation: degenerate inputs ===\n");
  testDegenerate();
  std::printf("\n=== Percolation: capillary depercolation ===\n");
  testCapillaryDepercolation();

  std::printf("\n%s (%d failure%s)\n", gFailed == 0 ? "ALL PASS" : "FAILURES",
              gFailed, gFailed == 1 ? "" : "s");
  return (gFailed == 0) ? 0 : 1;
}
