// Standalone unit test for TransportStats aggregation math.
// Compile with build_and_run.sh — no THAMES dependencies needed.

#include "TransportStats.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

namespace {

int gFailed = 0;

void expect(bool cond, const char *tag) {
  std::printf("%-55s %s\n", tag, cond ? "PASS" : "FAIL");
  if (!cond) ++gFailed;
}

void expectNear(double actual, double expected, double tol,
                const char *tag) {
  const bool ok = std::fabs(actual - expected) <= tol;
  std::printf("%-55s actual=%.6e expected=%.6e tol=%.2e %s\n",
              tag, actual, expected, tol, ok ? "PASS" : "FAIL");
  if (!ok) ++gFailed;
}

// --- Test group 1: monodisperse δ distribution ---
void testMonodisperse() {
  std::printf("\n--- monodisperse δ = 3e-6 m, 100 sites, K = 5 ---\n");
  const int N = 100;
  std::vector<double> deltas(N, 3.0e-6);
  std::vector<int> pids(N, 42);
  auto stats = xport::aggregateShellDistribution(deltas, pids, 5);
  expect(stats.numSitesReached == 100, "all sites reached");
  expect(stats.bins.size() == 5, "K bins produced");
  // All bins should have the same representative δ (monodisperse).
  for (size_t i = 0; i < stats.bins.size(); ++i) {
    char tag[64];
    std::snprintf(tag, sizeof(tag), "bin %zu deltaRep", i);
    expectNear(stats.bins[i].deltaRep, 3.0e-6, 1e-12, tag);
    std::snprintf(tag, sizeof(tag), "bin %zu siteFraction", i);
    expectNear(stats.bins[i].siteFraction, 0.2, 1e-9, tag);
    std::snprintf(tag, sizeof(tag), "bin %zu dominantShellPhaseId", i);
    expect(stats.bins[i].dominantShellPhaseId == 42, tag);
  }
  expectNear(stats.deltaHarmonicRaw, 3.0e-6, 1e-12,
             "monodisperse harmonic mean = δ");
  expectNear(stats.deltaArithmeticRaw, 3.0e-6, 1e-12,
             "monodisperse arithmetic mean = δ");
}

// --- Test group 2: bimodal δ distribution (kinetic + diffusion limited) ---
void testBimodal() {
  std::printf("\n--- bimodal 50/50, δ ∈ {1e-6, 100e-6}, K = 2 ---\n");
  std::vector<double> deltas;
  std::vector<int> pids;
  for (int i = 0; i < 50; ++i) { deltas.push_back(1.0e-6); pids.push_back(10); }
  for (int i = 0; i < 50; ++i) { deltas.push_back(100.0e-6); pids.push_back(20); }
  auto stats = xport::aggregateShellDistribution(deltas, pids, 2);
  expect(stats.numSitesReached == 100, "100 sites reached");
  expect(stats.bins.size() == 2, "K=2 bins");
  // Bin 0 should be the thin-shell bucket, bin 1 the thick.
  expectNear(stats.bins[0].deltaRep, 1.0e-6, 1e-12, "bin 0 median (thin)");
  expectNear(stats.bins[1].deltaRep, 100.0e-6, 1e-12, "bin 1 median (thick)");
  expectNear(stats.bins[0].siteFraction, 0.5, 1e-9, "bin 0 fraction 50%");
  expectNear(stats.bins[1].siteFraction, 0.5, 1e-9, "bin 1 fraction 50%");
  expect(stats.bins[0].dominantShellPhaseId == 10,
         "bin 0 dominant phase = 10");
  expect(stats.bins[1].dominantShellPhaseId == 20,
         "bin 1 dominant phase = 20");
  // Harmonic mean = 2 / (1/1e-6 + 1/100e-6) = 2 / (1.01e6) ≈ 1.980e-6.
  expectNear(stats.deltaHarmonicRaw, 2.0 / (1.0e6 + 1.0e4), 1e-12,
             "bimodal harmonic mean");
  // Arithmetic mean = (1e-6 + 100e-6) / 2 = 50.5e-6.
  expectNear(stats.deltaArithmeticRaw, 50.5e-6, 1e-12,
             "bimodal arithmetic mean");
  // Ratio arith/harmonic = 50.5e-6 / 1.98e-6 ≈ 25.5 → pathological.
  expect(xport::distributionIsPathological(stats, 5.0),
         "bimodal distribution flagged pathological");
}

// --- Test group 3: unreached sites (infinite δ) filtered out ---
void testUnreachedFiltered() {
  std::printf("\n--- 60 finite + 40 infinite, K = 3 ---\n");
  std::vector<double> deltas;
  std::vector<int> pids;
  for (int i = 0; i < 20; ++i) { deltas.push_back(2.0e-6); pids.push_back(1); }
  for (int i = 0; i < 20; ++i) { deltas.push_back(4.0e-6); pids.push_back(2); }
  for (int i = 0; i < 20; ++i) { deltas.push_back(6.0e-6); pids.push_back(3); }
  for (int i = 0; i < 40; ++i) {
    deltas.push_back(std::numeric_limits<double>::infinity());
    pids.push_back(-1);
  }
  auto stats = xport::aggregateShellDistribution(deltas, pids, 3);
  expect(stats.numSitesTotal == 100, "100 sites total (incl unreached)");
  expect(stats.numSitesReached == 60, "60 sites reached");
  expect(stats.bins.size() == 3, "K=3 bins");
  // Bins should represent the three finite δ values, one per bin.
  expectNear(stats.bins[0].deltaRep, 2.0e-6, 1e-12, "bin 0 = thinnest");
  expectNear(stats.bins[1].deltaRep, 4.0e-6, 1e-12, "bin 1 = medium");
  expectNear(stats.bins[2].deltaRep, 6.0e-6, 1e-12, "bin 2 = thickest");
  // Each bin has 1/3 of the 60 reached sites.
  for (int b = 0; b < 3; ++b) {
    char tag[64];
    std::snprintf(tag, sizeof(tag), "bin %d fraction 1/3", b);
    expectNear(stats.bins[b].siteFraction, 1.0 / 3.0, 1e-9, tag);
  }
}

// --- Test group 4: K = 1 collapses to a single scalar ---
void testKEqualsOne() {
  std::printf("\n--- K = 1 collapse, δ ∈ {1e-6, 2e-6, 3e-6, ...} ---\n");
  std::vector<double> deltas;
  std::vector<int> pids;
  for (int i = 1; i <= 10; ++i) {
    deltas.push_back(i * 1.0e-6);
    pids.push_back(i);
  }
  auto stats = xport::aggregateShellDistribution(deltas, pids, 1);
  expect(stats.bins.size() == 1, "single bin");
  expectNear(stats.bins[0].siteFraction, 1.0, 1e-12,
             "single-bin fraction = 1");
  // Median of 1e-6..10e-6 = 6e-6 (idx = 10/2 = 5 → deltas[5]).
  expectNear(stats.bins[0].deltaRep, 6.0e-6, 1e-12,
             "K=1 rep = median");
  // Arithmetic mean = 5.5e-6.
  expectNear(stats.deltaArithmeticRaw, 5.5e-6, 1e-12,
             "K=1 arithmetic mean");
}

// --- Test group 5: empty input yields a safe empty stats ---
void testEmpty() {
  std::printf("\n--- empty input ---\n");
  std::vector<double> deltas;
  std::vector<int> pids;
  auto stats = xport::aggregateShellDistribution(deltas, pids, 5);
  expect(stats.numSitesTotal == 0, "no total sites");
  expect(stats.numSitesReached == 0, "no reached sites");
  expect(stats.bins.size() == 1, "one placeholder empty bin");
  expect(stats.bins[0].siteFraction == 0.0, "empty bin has zero fraction");
}

// --- Test group 6: zero-shell site short-circuits the harmonic mean ---
void testZeroShellSite() {
  std::printf("\n--- one zero-shell site among many ---\n");
  std::vector<double> deltas = {0.0, 1.0e-6, 2.0e-6, 3.0e-6};
  std::vector<int> pids = {0, 1, 2, 3};
  auto stats = xport::aggregateShellDistribution(deltas, pids, 2);
  expect(stats.numSitesReached == 4, "4 reached sites");
  // Physical: any site with δ=0 provides infinite conductance
  // (parallel resistor short-circuit) → aggregate harmonic → 0.
  expectNear(stats.deltaHarmonicRaw, 0.0, 1e-30,
             "zero-shell → harmonic mean 0");
}

}  // namespace

int main() {
  std::printf("=== test_transport_stats ===\n");

  testMonodisperse();
  testBimodal();
  testUnreachedFiltered();
  testKEqualsOne();
  testEmpty();
  testZeroShellSite();

  std::printf("\n%d probe(s) failed\n", gFailed);
  return gFailed == 0 ? 0 : 1;
}
