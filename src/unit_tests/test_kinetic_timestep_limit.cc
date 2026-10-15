// Standalone unit test for KineticTimestepLimit.h: the one-voxel floor on the
// kinetics-based timestep limit. Header only; no THAMES dependencies.
// Compile with build_and_run.sh.
//
// Expected values are computed by hand in the comments, not by calling the
// code under test. Round inputs: N = 1e6 sites (100^3), V_m = 5e-5 m^3/mol,
// rho_0 = 1.6 g solid per cm^3 of microstructure (about w/c 0.32), f = 0.05.

#include "KineticTimestepLimit.h"

#include <cmath>
#include <cstdio>

namespace {

int gFailed = 0;

void expectNear(double actual, double expected, double relTol,
                const char *tag) {
  const double scale = std::fabs(expected) > 0.0 ? std::fabs(expected) : 1.0;
  const bool ok = std::fabs(actual - expected) <= relTol * scale;
  std::printf("%-62s actual=%.9g expected=%.9g %s\n", tag, actual, expected,
              ok ? "PASS" : "FAIL");
  if (!ok)
    ++gFailed;
}

void expectTrue(bool cond, const char *tag) {
  std::printf("%-62s %s\n", tag, cond ? "PASS" : "FAIL");
  if (!cond)
    ++gFailed;
}

const int N = 1000000;
const double VM = 5.0e-5;
const double RHO0 = 1.6;
const double F = 0.05;
const double TOL = 1.0e-12;

// 100 / (1e6 * 1e6 * 5e-5 * 1.6) = 100 / 8e7 = 1.25e-6 mol/100 g
const double NVOX = 1.25e-6;

// --- molesPerVoxel ---
void testMolesPerVoxel() {
  std::printf("\n--- molesPerVoxel ---\n");
  expectNear(kinetic_dt::molesPerVoxel(N, VM, RHO0), NVOX, TOL,
             "100 / (N 1e6 V_m rho_0)");

  // Frame identity: all N voxels of this DC would occupy the whole
  // microstructure, whose volume per 100 g of solid is 100 / rho_0 cm^3.
  // N * n_vox * V_m [m^3] * 1e6 [cm^3/m^3] = 100 / 1.6 = 62.5 cm^3.
  const double totalCm3 = static_cast<double>(N) *
                          kinetic_dt::molesPerVoxel(N, VM, RHO0) * VM * 1.0e6;
  expectNear(totalCm3, 62.5, TOL,
             "N voxels fill the microstructure: 100/rho_0 cm^3");

  // Doubling the lattice halves a voxel's worth.
  expectNear(kinetic_dt::molesPerVoxel(2 * N, VM, RHO0), NVOX / 2.0, TOL,
             "2N sites -> half the moles per voxel");

  expectNear(kinetic_dt::molesPerVoxel(0, VM, RHO0), 0.0, TOL,
             "no sites -> 0 (floor disabled)");
  expectNear(kinetic_dt::molesPerVoxel(N, 0.0, RHO0), 0.0, TOL,
             "zero molar volume -> 0 (floor disabled)");
  expectNear(kinetic_dt::molesPerVoxel(N, VM, -1.0), 0.0, TOL,
             "nonpositive solid mass -> 0 (floor disabled)");
}

// --- timestepLimit ---
void testTimestepLimit() {
  std::printf("\n--- timestepLimit ---\n");
  // Well-populated phase: the floor is irrelevant.
  // 0.05 * 1e-2 / 1e-3 = 0.5 h
  expectNear(kinetic_dt::timestepLimit(1.0e-2, NVOX, 1.0e-3, F), 0.5, TOL,
             "large phase unchanged: f n / r");
  expectNear(kinetic_dt::timestepLimit(-1.0e-2, NVOX, -1.0e-3, F), 0.5, TOL,
             "signs of n and r do not matter");

  // Exactly one voxel: floored and unfloored agree.
  // 0.05 * 1.25e-6 / 6e-8 = 6.25e-8 / 6e-8 = 25/24 h
  expectNear(kinetic_dt::timestepLimit(NVOX, NVOX, 6.0e-8, F), 25.0 / 24.0,
             TOL, "one voxel: floor and amount coincide");

  // Sub-voxel remainder (the 2026-10-14 Anhydrite case, ~0.008 voxel).
  // Unfloored: 0.05 * 1e-8 / 6e-8 = 0.05 / 6 = 1/120 h
  // Floored:   0.05 * 1.25e-6 / 6e-8 = 25/24 h
  expectNear(kinetic_dt::timestepLimit(1.0e-8, 0.0, 6.0e-8, F), 1.0 / 120.0,
             TOL, "sub-voxel, floor disabled: f n / r");
  expectNear(kinetic_dt::timestepLimit(1.0e-8, NVOX, 6.0e-8, F), 25.0 / 24.0,
             TOL, "sub-voxel, floored: f n_vox / r");

  // Precipitation from a trace amount: unfloored dt is below the
  // controller's 1e-5 h floor (stepTimeTHR_); floored dt is not.
  // Unfloored: 0.05 * 1e-10 / 1e-6 = 5e-6 h
  // Floored:   0.05 * 1.25e-6 / 1e-6 = 0.0625 h
  expectNear(kinetic_dt::timestepLimit(1.0e-10, 0.0, 1.0e-6, F), 5.0e-6,
             TOL, "trace precipitating, floor disabled");
  expectNear(kinetic_dt::timestepLimit(1.0e-10, NVOX, 1.0e-6, F), 0.0625,
             TOL, "trace precipitating, floored");
}

// --- the geometric decay the floor removes ---
// Dissolve a sub-voxel remainder at constant rate r, each step sized by the
// limit and clamped to what is left (as calculateKineticStep's depletion
// clamp does).
void testDepletionSequence() {
  std::printf("\n--- depletion sequence at constant rate ---\n");
  const double n0 = 1.0e-8;
  const double r = 6.0e-8;

  // Unfloored: each step removes f of the remainder, so after k steps
  // n_k = n0 (1 - f)^k and t_k = (n0 - n_k) / r < n0 / r = 1/6 h.
  double n = n0, t = 0.0;
  for (int k = 0; k < 200; ++k) {
    const double dt = kinetic_dt::timestepLimit(n, 0.0, r, F);
    n -= std::fmin(r * dt, n);
    t += dt;
  }
  expectNear(n, n0 * std::pow(1.0 - F, 200), 1.0e-9,
             "unfloored: n_200 = n0 (1 - f)^200, still > 0");
  expectTrue(n > 0.0 && t < n0 / r,
             "unfloored: 200 steps never reach t = n0 / r = 1/6 h");

  // Floored: the first step asks for r dt = 0.05 n_vox = 6.25e-8 > n0, the
  // clamp takes exactly n0, and the phase is gone after one step.
  n = n0;
  int steps = 0;
  while (n > 0.0 && steps < 10) {
    const double dt = kinetic_dt::timestepLimit(n, NVOX, r, F);
    n -= std::fmin(r * dt, n);
    ++steps;
  }
  expectTrue(n == 0.0 && steps == 1,
             "floored: remainder cleared in one step");
}

} // namespace

int main() {
  testMolesPerVoxel();
  testTimestepLimit();
  testDepletionSequence();
  std::printf("\n%s: %d failure(s)\n", gFailed ? "FAILED" : "ALL PASSED",
              gFailed);
  return gFailed ? 1 : 0;
}
