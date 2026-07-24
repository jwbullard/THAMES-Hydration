/**
@file SaturatingRate.cc
@brief Implementations of the sat:: free functions declared in SaturatingRate.h.

Numerical notes:
- At Omega = 1 exactly, `log(1) = 0`, so `(±B ln Omega)^n = 0^n = 0`,
  `exp(-0) = 1`, and `1 - 1 = 0`. Both dissolution and precipitation
  branches return 0 at equilibrium naturally, no special-case needed.
- As Omega -> 0 (dissolution far from equilibrium), `-B ln Omega`
  grows without bound; `pow(large, n)` is finite in IEEE doubles up
  to the pow overflow threshold, at which point `exp(-inf) = 0` and
  the rate saturates at k. No numerical guard needed for typical
  cement-chemistry Omegas.
- As Omega -> +inf (precipitation far from equilibrium), similarly
  saturates at k.
*/

#include "SaturatingRate.h"

#include <cmath>

#include "PhysicalConstants.h"

namespace sat {

double dissolutionRate(double k, double B, double n, double Omega) {
  // Guard against non-physical inputs.
  if (Omega >= 1.0) return 0.0;
  if (Omega <= 0.0) return k;   // Undefined at Omega = 0; return the asymptote.

  const double arg = -B * std::log(Omega);   // > 0 for Omega < 1, B > 0
  const double powered = std::pow(arg, n);
  return k * (1.0 - std::exp(-powered));
}

double precipitationRate(double k, double B, double n, double Omega) {
  // Guard against non-physical inputs.
  if (Omega <= 1.0) return 0.0;

  const double arg = B * std::log(Omega);    // > 0 for Omega > 1, B > 0
  const double powered = std::pow(arg, n);
  return k * (1.0 - std::exp(-powered));
}

double arrheniusScale(double k_ref, double T_ref, double T, double dHact) {
  return k_ref *
         std::exp((dHact / GASCONSTANT) * ((1.0 / T_ref) - (1.0 / T)));
}

}  // namespace sat
