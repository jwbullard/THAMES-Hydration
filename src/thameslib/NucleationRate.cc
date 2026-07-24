/**
@file NucleationRate.cc
@brief Implementations of the CNT rate free functions declared in NucleationRate.h.

Formulas follow Cell 5 of the Python prototype notebook verbatim (see
~/Research/THAMES-Tests-2026/Scripts/NucleationCNT-Prototype.ipynb).

Physical constant: k_B = 1.380649e-23 J/K (exact SI-2019 value).
*/

#include "NucleationRate.h"

#include <cmath>

#include "PhysicalConstants.h"

namespace {
// Convert molar volume [m^3/mol] to volume per molecule [m^3/molecule].
// The CNT formulas (with k_B, not R) require per-molecule volume, but GEMS
// and the JSON schema report molar volume; this hides the conversion inside
// the free functions so callers pass what they naturally have.
inline double molecularVolume(double V_molar) {
  return V_molar / AVOGADROCONSTANT;
}
}  // namespace

namespace cnt {

double fTheta(int theta_deg) {
  const double theta_rad = static_cast<double>(theta_deg) * Pi / 180.0;
  const double c = std::cos(theta_rad);
  return 0.25 * (2.0 + c) * (1.0 - c) * (1.0 - c);
}

double criticalRadius(double gamma, double v_m, double T_K, double ln_S) {
  const double v_molec = molecularVolume(v_m);
  return 2.0 * gamma * v_molec / (BOLTZMANNCONSTANT * T_K * ln_S);
}

double deltaGStarHom(double gamma, double v_m, double T_K, double ln_S) {
  const double v_molec = molecularVolume(v_m);
  const double numer = 16.0 * Pi * gamma * gamma * gamma * v_molec * v_molec;
  const double denom = 3.0 * std::pow(BOLTZMANNCONSTANT * T_K * ln_S, 2.0);
  return numer / denom;
}

double nucleationRate(const NucleationParameters &np,
                      double v_m, double T_K, double S) {
  if (S <= 1.0) return 0.0;
  const double ln_S = std::log(S);
  const double dG = deltaGStarHom(np.gamma, v_m, T_K, ln_S) * fTheta(np.theta_deg);
  return np.A0 * std::exp(-dG / (BOLTZMANNCONSTANT * T_K));
}

double voxelsPerCycle(const NucleationParameters &np,
                      double v_m, double T_K, double S,
                      double dt_seconds,
                      double V_electrolyte, double V_voxel) {
  if (S <= 1.0) return 0.0;
  const double ln_S = std::log(S);
  const double J = nucleationRate(np, v_m, T_K, S);
  const double r_star = criticalRadius(np.gamma, v_m, T_K, ln_S);
  const double V_crit = (4.0 / 3.0) * Pi * r_star * r_star * r_star;
  return J * V_electrolyte * dt_seconds * V_crit / V_voxel;
}

}  // namespace cnt
