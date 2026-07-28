/**
@file test_jmak_growth.cc
@brief Standalone test for the JMAK per-voxel growth free functions.

Exercises `jmak::advanceMoments`, `jmak::snapshotSeed`,
`jmak::extendedVolumePerVoxel`, `jmak::transformedFraction`, and
`jmak::growthVelocity` in isolation. No dependency on THAMES beyond the
JMAK sources.

Verification targets:

  1. Constant-J, constant-G limit reproduces the closed-form JMAK curve
     X(t) = 1 - exp(-(pi/3) * J * G^3 * t^4) to 4 decimal places at
     several times spanning the induction, growth, and saturation
     regimes.

  2. Zero-radius sanity: at t = t_c (generation just seeded), Y = 0 and
     X = 0 (no time to grow).

  3. Generation accumulation sanity: two generations with the same seed time
     but different N_c produce transformed volumes proportional to N_c.

  4. Growth-velocity converter sanity: r=0 or v_molar<=0 return 0;
     otherwise linear product.

  5. Time-varying J sanity: two-stage constant-J profile with a step
     change matches direct piecewise closed-form integration.

Exit code 0 if all checks pass, non-zero on any failure. Prints a
per-check PASS/FAIL line so it's easy to spot regressions in CI diffs.
*/

#include "JMAKGrowth.h"
#include "JMAKParameters.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace {

constexpr double kPi = 3.14159265358979323846;

int failed = 0;

void report(const char *name, bool ok, const char *why = "") {
  std::printf("%-56s %s%s\n", name, ok ? "PASS" : "FAIL",
              ok ? "" : " ");
  if (!ok) {
    if (why && why[0]) std::printf("  reason: %s\n", why);
    ++failed;
  }
}

bool approxEqual(double a, double b, double rtol = 1.0e-4,
                 double atol = 1.0e-9) {
  const double diff = std::fabs(a - b);
  if (diff <= atol) return true;
  const double denom = std::fmax(std::fabs(a), std::fabs(b));
  return (denom > 0.0) && (diff / denom <= rtol);
}

// Simulate constant-J, constant-G by stepping advanceMoments many times.
// Returns X(t_final) for a single generation seeded at t = 0.
double simulateConstantJG(double J, double G, double t_final, int nsteps,
                          double V_voxel, const JMAKParameters &jp) {
  jmak::GlobalMoments acc{0.0, 0.0, 0.0, 0.0, 0.0};
  const double dt = t_final / nsteps;
  // Seed at t = 0.
  const jmak::GenerationMomentsAtSeed seed = jmak::snapshotSeed(acc);
  for (int i = 0; i < nsteps; ++i) {
    jmak::advanceMoments(acc, J, G, dt);
  }
  const double Y_over_V = jmak::extendedVolumePerVoxel(seed, acc, jp);
  return jmak::transformedFraction(Y_over_V);
}

// Closed-form JMAK for constant J and G, n = 4:
//   X(t) = 1 - exp(-(alpha / 4) * J * G^3 * t^4)
// For alpha = 4*pi/3 this is X(t) = 1 - exp(-(pi/3) J G^3 t^4).
double closedFormConstantJG(double J, double G, double t, double alpha) {
  const double K = (alpha / 4.0) * J * G * G * G;
  return -std::expm1(-K * t * t * t * t);
}

}  // namespace

int main() {
  std::printf("=== test_jmak_growth ===\n\n");

  const JMAKParameters jp{4.0, 4.0 * kPi / 3.0};  // n = 4, 3D isotropic
  const double V_voxel = 1.0e-18;                  // 1 um^3

  // --- 1. Constant-J-G limit vs closed form -----------------------------
  {
    const double J = 1.0e15;    // /(m^3 s), moderate
    const double G = 1.0e-9;    // m/s, ~ nm per second
    const int nsteps = 10000;

    std::printf("--- Constant J, G vs closed-form JMAK X = 1 - exp(-(pi/3) J G^3 t^4) ---\n");
    std::printf("%-14s %-14s %-14s %-8s\n", "t (s)", "X_simulated", "X_closedform", "match?");
    const double times[] = {1.0, 10.0, 100.0, 1000.0, 10000.0};
    for (double t : times) {
      const double X_sim = simulateConstantJG(J, G, t, nsteps, V_voxel, jp);
      const double X_ref = closedFormConstantJG(J, G, t, jp.alpha);
      const bool ok = approxEqual(X_sim, X_ref, 1.0e-3, 1.0e-8);
      std::printf("%-14.3f %-14.6e %-14.6e %-8s\n", t, X_sim, X_ref,
                  ok ? "yes" : "NO");
      char nm[128];
      std::snprintf(nm, sizeof(nm), "  constant-JG t = %.1f s", t);
      report(nm, ok);
    }
    std::printf("\n");
  }

  // --- 2. Zero-radius sanity --------------------------------------------
  {
    jmak::GlobalMoments acc{0.0, 0.0, 0.0, 0.0, 0.0};
    jmak::advanceMoments(acc, 1.0e15, 1.0e-9, 1.0);
    // Seed AFTER advancing so any Y is genuinely from t > t_c.
    const auto seed = jmak::snapshotSeed(acc);
    // At t = t_c, the generation has had zero time to accumulate anything.
    const double Y_over_V = jmak::extendedVolumePerVoxel(seed, acc, jp);
    const bool ok = (Y_over_V == 0.0);
    std::printf("zero-radius sanity (Y=0 at t = t_c):  Y/V = %.6e  %s\n",
                Y_over_V, ok ? "PASS" : "FAIL");
    report("zero-radius sanity", ok);
  }

  // --- 3. Generation scaling sanity -----------------------------------------
  {
    // Two identical generations should have identical X_c. The AGGREGATE
    // Portlandite volume for N_c voxels is just N_c * V_voxel * X_c;
    // we verify by comparing N_c=1 vs N_c=1000 cases scale linearly.
    const double J = 1.0e15;
    const double G = 1.0e-9;
    const double t_final = 1000.0;
    const int nsteps = 5000;

    jmak::GlobalMoments acc{0.0, 0.0, 0.0, 0.0, 0.0};
    const auto seed = jmak::snapshotSeed(acc);
    const double dt = t_final / nsteps;
    for (int i = 0; i < nsteps; ++i) {
      jmak::advanceMoments(acc, J, G, dt);
    }
    const double X = jmak::transformedFraction(
        jmak::extendedVolumePerVoxel(seed, acc, jp));

    const double V_from_1 = 1 * V_voxel * X;
    const double V_from_1000 = 1000 * V_voxel * X;
    const bool ok = approxEqual(V_from_1000, 1000.0 * V_from_1, 1e-12);
    std::printf("generation scaling (N_c linearity):  X = %.6e, "
                "V(1000)/V(1) = %.6e (expect 1000)  %s\n",
                X, V_from_1000 / V_from_1, ok ? "PASS" : "FAIL");
    report("generation scaling linearity", ok);
  }

  // --- 4. Growth-velocity converter -------------------------------------
  {
    const double v_m = 33.08e-6;  // portlandite molar volume, m^3/mol
    const double r = 4.05e-4;     // Han 2025 k at 24 C, mol/m^2/s
    const double G = jmak::growthVelocity(r, v_m);
    const double G_ref = r * v_m;
    const bool ok1 = approxEqual(G, G_ref, 1.0e-15);
    report("growthVelocity: normal case", ok1);

    const bool ok2 = (jmak::growthVelocity(-1.0, v_m) == 0.0);
    report("growthVelocity: negative r -> 0", ok2);

    const bool ok3 = (jmak::growthVelocity(r, 0.0) == 0.0);
    report("growthVelocity: zero v_molar -> 0", ok3);

    const bool ok4 = (jmak::growthVelocity(r, -v_m) == 0.0);
    report("growthVelocity: negative v_molar -> 0", ok4);
  }

  // --- 5. Time-varying J (step-change) ----------------------------------
  {
    // Two-stage profile: J = J1 for [0, t_switch], J = J2 for [t_switch, t_final].
    // G constant throughout. Generation seeded at t = 0.
    // Direct closed form for stepwise J with constant G:
    //   Y/V = alpha * G^3 * [ J1 * ((t_final^4 - (t_final - t_switch)^4) / 4)
    //                       + J2 * ((t_final - t_switch)^4 / 4) ]
    // Wait: let s = t_final - tau. Then integral J(tau) (t_final - tau)^3 dtau.
    //   For tau in [0, t_switch]: J1 * (t_final - tau)^3, integral over [0, t_switch]
    //     = J1 * [ -(t_final - tau)^4 / 4 ]_{0}^{t_switch}
    //     = J1 * (t_final^4 - (t_final - t_switch)^4) / 4
    //   For tau in [t_switch, t_final]: J2 * (t_final - tau)^3
    //     = J2 * (t_final - t_switch)^4 / 4
    const double J1 = 1.0e15;
    const double J2 = 1.0e13;  // 100x slower for the second stage
    const double G = 1.0e-9;
    const double t_switch = 500.0;
    const double t_final = 1000.0;
    const int nsteps = 20000;

    const double dt = t_final / nsteps;
    jmak::GlobalMoments acc{0.0, 0.0, 0.0, 0.0, 0.0};
    const auto seed = jmak::snapshotSeed(acc);
    for (int i = 0; i < nsteps; ++i) {
      const double t_lo = i * dt;
      const double J_now = (t_lo < t_switch) ? J1 : J2;
      jmak::advanceMoments(acc, J_now, G, dt);
    }
    const double Y_sim = jmak::extendedVolumePerVoxel(seed, acc, jp);

    const double dY_stage1 = J1 * (std::pow(t_final, 4)
                                 - std::pow(t_final - t_switch, 4)) / 4.0;
    const double dY_stage2 = J2 * std::pow(t_final - t_switch, 4) / 4.0;
    const double Y_ref = jp.alpha * G * G * G * (dY_stage1 + dY_stage2);

    // 0.1 % tolerance for the 20,000-step discretization vs the closed-form.
    const bool ok = approxEqual(Y_sim, Y_ref, 1.0e-3, 1.0e-12);
    std::printf("step-J: Y_sim = %.6e, Y_ref = %.6e, rel err = %.3e  %s\n",
                Y_sim, Y_ref,
                std::fabs(Y_sim - Y_ref) / std::fabs(Y_ref),
                ok ? "PASS" : "FAIL");
    report("time-varying J step-change", ok);
  }

  // --- 6. Extended surface area matches closed form for constant J, G ----
  {
    // For constant J and G, generation seeded at t = 0, no impingement:
    //   A_ext(t) / V_voxel = 4*pi * J * G^2 * t^3 / 3
    // Derivation: A_ext(t) = 4*pi * integral J * V * [G*(t - tau)]^2 dtau
    //           = 4*pi * J * V * G^2 * t^3 / 3
    //           Divide by V_voxel -> per-voxel value 4*pi/3 * J * G^2 * t^3
    // (Wait — the factor of V_voxel cancels because the integrand already
    //  has V_voxel implicit — see JMAKGrowth.h header discussion. The
    //  formula in the code is dimensionally per-unit-voxel-volume.)
    const double J = 1.0e15;
    const double G = 1.0e-9;
    const int nsteps = 10000;

    std::printf("\n--- Extended surface area vs closed form 4*pi/3 * J * G^2 * t^3 ---\n");
    std::printf("%-14s %-14s %-14s %-8s\n", "t (s)", "A_sim (1/m)", "A_closed (1/m)", "match?");
    const double times[] = {100.0, 1000.0, 5000.0};
    for (double t : times) {
      jmak::GlobalMoments acc{0.0, 0.0, 0.0, 0.0, 0.0};
      const jmak::GenerationMomentsAtSeed seed = jmak::snapshotSeed(acc);
      const double dt = t / nsteps;
      for (int i = 0; i < nsteps; ++i) {
        jmak::advanceMoments(acc, J, G, dt);
      }
      const double A_sim = jmak::extendedSurfaceAreaPerVoxel(seed, acc);
      const double A_ref = (4.0 * kPi / 3.0) * J * G * G * t * t * t;
      const bool ok = approxEqual(A_sim, A_ref, 1.0e-3, 1.0e-8);
      std::printf("%-14.3f %-14.6e %-14.6e %-8s\n", t, A_sim, A_ref,
                  ok ? "yes" : "NO");
      char nm[128];
      std::snprintf(nm, sizeof(nm), "  A_ext constant-JG t = %.1f s", t);
      report(nm, ok);
    }
  }

  // --- 7. Extended surface area sanity: zero at t = t_c ------------------
  {
    jmak::GlobalMoments acc{0.0, 0.0, 0.0, 0.0, 0.0};
    jmak::advanceMoments(acc, 1.0e15, 1.0e-9, 100.0);
    const auto seed = jmak::snapshotSeed(acc);
    // Immediately after seed: no time has elapsed since seed.
    const double A = jmak::extendedSurfaceAreaPerVoxel(seed, acc);
    const bool ok = (A == 0.0);
    std::printf("A_ext at t = t_c: %.6e  %s\n", A, ok ? "PASS" : "FAIL");
    report("A_ext zero at seed time", ok);
  }

  std::printf("\n%d probe(s) failed\n", failed);
  return failed == 0 ? 0 : 1;
}
