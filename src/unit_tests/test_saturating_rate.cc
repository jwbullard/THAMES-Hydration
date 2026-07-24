/**
@file test_saturating_rate.cc
@brief Standalone verification of the saturating-rate free functions.

Verifies `sat::dissolutionRate` and `sat::precipitationRate` against
values computed from Han et al. 2025 CEJ Eq. (7) with Han's published
Portlandite calibration at 24 °C:
    k = 4.05e-4 mol/(m^2 s), B = 0.74, n = 1.9.

Also spot-checks arrheniusScale for a small temperature shift.

Build & run: ./build_and_run.sh in this directory. Zero external deps
beyond `<cmath>` and the sat:: free-function pair SaturatingRate.{h,cc}.
*/

#include <cmath>
#include <cstdio>

#include "PhysicalConstants.h"
#include "SaturatingRate.h"

namespace {

// Reference rates computed by hand from Han Eq. (7) with the Portlandite
// calibration. Each expected value is documented with the derivation so
// a maintainer can verify by inspection.
struct DissolutionProbe {
  double Omega;
  double expected_rate;      // mol / m^2 / s
  const char *derivation;    // how the expected value was obtained
};

// Wide-tolerance comparison: allow 1% relative difference. The formula
// itself is smooth and analytic, so agreement should be much tighter
// than 1%; this leaves room for future double-precision quirks.
bool close(double got, double expected, double tol_rel = 0.01) {
  if (expected == 0.0) return std::fabs(got) < 1e-15;
  return std::fabs(got - expected) / std::fabs(expected) < tol_rel;
}

}  // namespace

int main() {
  const double k = 4.05e-4;  // Han Table 2, 24 C
  const double B = 0.74;
  const double n = 1.9;

  int fails = 0;

  std::printf("=== Portlandite saturating-rate check (Han 2025, 24 C) ===\n");
  std::printf("k = %g mol/m^2/s, B = %g, n = %g\n\n", k, B, n);

  //
  // Section 1 — dissolutionRate reference values at six probe Omega.
  //
  // Reference values were locked in by computing
  //   r = k * (1 - exp[-(-B ln Omega)^n])
  // in double precision with the calibration constants above. Each row
  // was hand-verified against textbook algebra (ln, pow, exp) to confirm
  // the formula returns the expected physical behavior (monotonic
  // decrease as Omega -> 1, saturation at k as Omega -> 0).
  //
  // These serve as a regression test: any change to the formula, the
  // free-function API, or the numerical guards must reproduce these
  // values to within 1% relative.
  //
  DissolutionProbe probes[] = {
      {0.05,  4.006700e-04, "far from equilibrium; saturating toward k"},
      {0.10,  3.791776e-04, "clearly undersaturated"},
      {0.30,  2.235657e-04, "midrange; where the etch-pit -> step-flow transition begins"},
      {0.50,  9.929315e-05, "close to half of the calibrated range"},
      {0.70,  3.098454e-05, "approaching equilibrium; rate ~1% of k"},
      {0.90,  3.165051e-06, "near equilibrium; rate ~0.1% of k"},
  };

  std::printf("%-8s %-16s %-16s %-6s\n", "Omega", "computed", "expected", "PASS?");
  for (const DissolutionProbe &p : probes) {
    const double got = sat::dissolutionRate(k, B, n, p.Omega);
    const bool ok = close(got, p.expected_rate);
    if (!ok) fails++;
    std::printf("%-8.3f %-16.6e %-16.6e %-6s\n",
                p.Omega, got, p.expected_rate, ok ? "PASS" : "FAIL");
  }
  std::printf("\n");

  //
  // Section 2 — boundary conditions.
  //
  std::printf("=== Boundary conditions ===\n");

  const double at_eq = sat::dissolutionRate(k, B, n, 1.0);
  const bool eq_ok = (at_eq == 0.0);
  if (!eq_ok) fails++;
  std::printf("dissolutionRate at Omega=1.0 = %g  %s (must be exactly 0)\n",
              at_eq, eq_ok ? "PASS" : "FAIL");

  const double at_zero = sat::dissolutionRate(k, B, n, 0.0);
  const bool zero_ok = close(at_zero, k);
  if (!zero_ok) fails++;
  std::printf("dissolutionRate at Omega=0.0 = %g  %s (must saturate at k = %g)\n",
              at_zero, zero_ok ? "PASS" : "FAIL", k);

  const double far_below = sat::dissolutionRate(k, B, n, 1e-12);
  const bool far_ok = close(far_below, k);
  if (!far_ok) fails++;
  std::printf("dissolutionRate at Omega=1e-12 = %g  %s (must approach k)\n",
              far_below, far_ok ? "PASS" : "FAIL");

  // dissolutionRate must return 0 if called incorrectly with Omega > 1.
  const double wrong_diss = sat::dissolutionRate(k, B, n, 2.0);
  const bool wrong_diss_ok = (wrong_diss == 0.0);
  if (!wrong_diss_ok) fails++;
  std::printf("dissolutionRate at Omega=2.0 (misuse) = %g  %s (returns 0)\n",
              wrong_diss, wrong_diss_ok ? "PASS" : "FAIL");

  std::printf("\n");

  //
  // Section 3 — precipitationRate at Omega > 1.
  //
  // r_precip(Omega) = k * (1 - exp[-(B ln Omega)^n])
  // Same shape mirrored: saturates at k for large Omega, zero at Omega = 1.
  //
  std::printf("=== Precipitation branch (Omega > 1) ===\n");

  const double pr_eq = sat::precipitationRate(k, B, n, 1.0);
  const bool pr_eq_ok = (pr_eq == 0.0);
  if (!pr_eq_ok) fails++;
  std::printf("precipitationRate at Omega=1.0 = %g  %s (must be exactly 0)\n",
              pr_eq, pr_eq_ok ? "PASS" : "FAIL");

  const double pr_2 = sat::precipitationRate(k, B, n, 2.0);
  const double pr_5 = sat::precipitationRate(k, B, n, 5.0);
  const double pr_far = sat::precipitationRate(k, B, n, 1e6);
  const bool pr_monotonic = (pr_2 < pr_5) && (pr_5 < pr_far);
  const bool pr_saturate = close(pr_far, k);
  if (!pr_monotonic) fails++;
  if (!pr_saturate) fails++;
  std::printf("precipitationRate(2)=%g precipitationRate(5)=%g precipitationRate(1e6)=%g\n",
              pr_2, pr_5, pr_far);
  std::printf("  monotonic increase: %s   saturates at k: %s\n",
              pr_monotonic ? "PASS" : "FAIL",
              pr_saturate ? "PASS" : "FAIL");

  const double wrong_precip = sat::precipitationRate(k, B, n, 0.5);
  const bool wrong_precip_ok = (wrong_precip == 0.0);
  if (!wrong_precip_ok) fails++;
  std::printf("precipitationRate at Omega=0.5 (misuse) = %g  %s (returns 0)\n",
              wrong_precip, wrong_precip_ok ? "PASS" : "FAIL");

  std::printf("\n");

  //
  // Section 4 — arrheniusScale for the Han et al. Portlandite dH.
  //
  //   dH_a = 13.9 kJ/mol = 13900 J/mol; refT = 297.15 K.
  //   Expected factor at T = 313 K:
  //     exp[(13900 / 8.314) * (1/297.15 - 1/313)]
  //     = exp[1671.9 * 1.7017e-4]
  //     = exp[0.28455]
  //     = 1.3292
  //
  std::printf("=== arrheniusScale ===\n");
  const double dHact = 13900.0;    // J / mol
  const double T_ref = 297.15;     // K (24 C, Han reference)
  const double T = 313.0;          // K (40 C)
  const double k_scaled = sat::arrheniusScale(k, T_ref, T, dHact);
  const double expected_factor = 1.3292;
  const double expected_k = k * expected_factor;
  const bool arr_ok = close(k_scaled, expected_k);
  if (!arr_ok) fails++;
  std::printf("k(24 C) = %g, k(40 C) = %g, factor = %.4f  %s (expected ~1.3292)\n",
              k, k_scaled, k_scaled / k, arr_ok ? "PASS" : "FAIL");

  //
  // Report.
  //
  std::printf("\n%d probe(s) failed\n", fails);
  return fails == 0 ? 0 : 1;
}
