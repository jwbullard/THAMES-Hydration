/**
@file test_nucleation_rate.cc
@brief Standalone verification of the CNT free functions in NucleationRate.

Reproduces the sanity numbers from Cell 5 of the Python prototype at
~/Research/THAMES-Tests-2026/Scripts/NucleationCNT-Prototype.ipynb using the
calibrated portlandite defaults (gamma = 0.044 J/m^2, theta = 180 deg,
A0 = 1e30 /(m^3 s)) with T = 298.15 K, v_m = 33.08e-6 m^3/mol,
dt = 36 s (= 0.01 h), 200^3 lattice at 1 um resolution, 50% porosity.

Expected log10 voxel-per-cycle for the four probe saturation ratios,
recomputed by hand from the same textbook CNT (matches the prototype's
Cell 5 verbatim; supersedes an earlier looser sanity table where the
S=2 target of "10^-55" was confused with the exp-argument):
  S = 2.0   -> log10(N) ~ -43.3
  S = 4.5   -> log10(N) ~ -0.4   (calibrated onset region; N ~ O(1) family)
  S = 10.0  -> log10(N) ~  5.85
  S = 50.0  -> log10(N) ~  8.5

Tolerance: log10 within +/- 0.5 of expected. Any drift beyond that indicates
a real formula bug.

Build & run: ./build_and_run.sh in this directory. Zero external deps
beyond the free-function pair NucleationRate.{h,cc} and NucleationParameters.h.
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "NucleationParameters.h"
#include "NucleationRate.h"

int main() {
  NucleationParameters np;
  np.gamma = 0.044;
  np.theta_deg = 180;
  np.A0 = 1.0e30;

  const double v_m = 33.08e-6;         // m^3/mol (33.08 cm^3/mol, portlandite)
  const double T_K = 298.15;
  const double dt_s = 36.0;             // 0.01 h in seconds
  const double V_voxel = 1.0e-18;       // 1 um edge
  const double V_electrolyte = 4.0e-12; // 50% of 200^3 * V_voxel

  const double f180 = cnt::fTheta(180);
  std::printf("Sanity check: f(180 deg) = %.15f (expect exactly 1.0)\n", f180);
  if (std::fabs(f180 - 1.0) > 1e-12) {
    std::printf("FAIL: f(180) must equal 1 for the homogeneous limit\n");
    return 1;
  }

  struct Probe {
    double S;
    double expected_log10_N;
  } probes[] = {
      {2.0,  -43.3},
      {4.5,   -0.4},   // calibrated onset; N ~ O(1) family
      {10.0,   5.85},
      {50.0,   8.5}
  };
  const double log10_tol = 0.5;

  std::printf("\n");
  std::printf("%6s %10s %12s %10s %14s %14s %14s\n",
              "S", "ln(S)", "r*[m]", "f(theta)",
              "J[1/(m^3 s)]", "N/cycle", "PASS?");

  int fails = 0;
  for (const Probe &p : probes) {
    const double ln_S = std::log(p.S);
    const double r_star = cnt::criticalRadius(np.gamma, v_m, T_K, ln_S);
    const double J = cnt::nucleationRate(np, v_m, T_K, p.S);
    const double N = cnt::voxelsPerCycle(np, v_m, T_K, p.S,
                                          dt_s, V_electrolyte, V_voxel);
    const double log10_N = (N > 0.0) ? std::log10(N) : -1e99;
    const bool ok = std::fabs(log10_N - p.expected_log10_N) <= log10_tol;
    if (!ok) fails++;
    std::printf("%6.2f %10.4f %12.3e %10.6f %14.3e %14.3e   %s (expect log10~%.1f, got %.2f)\n",
                p.S, ln_S, r_star, cnt::fTheta(np.theta_deg),
                J, N,
                ok ? "PASS" : "FAIL",
                p.expected_log10_N, log10_N);
  }

  std::printf("\n%d probe(s) failed\n", fails);
  return fails == 0 ? 0 : 1;
}
