/**
@file JMAKGrowth.cc
@brief Implementations of the JMAK free functions declared in JMAKGrowth.h.

See JMAKGrowth.h for the physical model, moment-decomposition
derivation, and unit conventions.
*/

#include "JMAKGrowth.h"

#include <cmath>

namespace jmak {

void advanceMoments(GlobalMoments &acc, double J, double G, double dt) {
  // End-of-interval convention: r_acc advances first, then M_k use the
  // updated r_acc. See file header for rationale (equivalent to
  // right-endpoint numerical integration; midpoint would be marginally
  // more accurate but the discrepancy is second-order in dt).
  acc.r_acc += G * dt;
  const double Jdt = J * dt;
  acc.M0 += Jdt;
  acc.M1 += Jdt * acc.r_acc;
  acc.M2 += Jdt * acc.r_acc * acc.r_acc;
  acc.M3 += Jdt * acc.r_acc * acc.r_acc * acc.r_acc;
}

GenerationMomentsAtSeed snapshotSeed(const GlobalMoments &acc) {
  GenerationMomentsAtSeed s;
  s.M0_c = acc.M0;
  s.M1_c = acc.M1;
  s.M2_c = acc.M2;
  s.M3_c = acc.M3;
  return s;
}

double extendedVolumePerVoxel(const GenerationMomentsAtSeed &seed,
                              const GlobalMoments &acc,
                              const JMAKParameters &jp) {
  // Difference-of-moments form (see file header). Requires n = 4 in the
  // continuous-nucleation convention (cubic growth kernel). alpha is
  // the morphology coefficient (4*pi/3 for 3D isotropic spheres).
  const double g = acc.r_acc;
  const double dM0 = acc.M0 - seed.M0_c;
  const double dM1 = acc.M1 - seed.M1_c;
  const double dM2 = acc.M2 - seed.M2_c;
  const double dM3 = acc.M3 - seed.M3_c;
  const double bracket = g * g * g * dM0
                       - 3.0 * g * g * dM1
                       + 3.0 * g * dM2
                       - dM3;
  // Guard against negative values from floating-point roundoff near
  // t = t_c (all dMk small and nearly cancel). Physically Y >= 0 always.
  const double y = jp.alpha * bracket;
  return (y > 0.0) ? y : 0.0;
}

double transformedFraction(double Y_over_V) {
  // X = 1 - exp(-Y/V), computed as -expm1(-Y/V) for numerical stability
  // when Y is small (naive 1-exp gets catastrophic cancellation).
  if (Y_over_V <= 0.0) return 0.0;
  return -std::expm1(-Y_over_V);
}

double growthVelocity(double r_surface_normal, double v_molar) {
  if (r_surface_normal <= 0.0) return 0.0;
  if (v_molar <= 0.0) return 0.0;
  return r_surface_normal * v_molar;
}

}  // namespace jmak
