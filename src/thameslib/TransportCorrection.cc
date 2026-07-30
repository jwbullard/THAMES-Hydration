/**
@file TransportCorrection.cc
@brief Phase-2 stub implementations. Phase 3 fills in Newton solver
       and per-shell D_eff mapping.
*/

#include "TransportCorrection.h"

#include "TransportParameters.h"

namespace xport {

double kineticRate(double k, double area, double driving_force) {
  return k * area * driving_force;
}

double diffusionRate(double dEff, double area, double deltaC,
                     double delta) {
  if (delta <= 0.0) return 0.0;
  return dEff * area * deltaC / delta;
}

double solveSurfaceConcentrationLinear(double k, double C_eq, double dEff,
                                       double delta, double C_bulk) {
  // Solve  k · (1 - C_surf / C_eq)  =  dEff · (C_surf - C_bulk) / δ
  // Rearrange:  k - k · C_surf / C_eq  =  dEff · C_surf / δ  -  dEff · C_bulk / δ
  //             k + dEff · C_bulk / δ  =  C_surf · (k / C_eq + dEff / δ)
  //             C_surf  =  (k + dEff · C_bulk / δ)  /  (k / C_eq + dEff / δ)
  if (delta <= 0.0) return C_bulk;      // no shell, C_surf = C_bulk
  if (C_eq <= 0.0) return C_bulk;       // degenerate
  const double num = k + dEff * C_bulk / delta;
  const double den = k / C_eq + dEff / delta;
  if (den <= 0.0) return C_bulk;
  return num / den;
}

double pickDEff(int /*shellPhaseId*/, const TransportParameters &params) {
  // Phase 2 stub: return the block's global dEff regardless of shell
  // composition. Phase 3 will introduce a per-shell-phase map so that
  // e.g. Ca+2 through C-S-H and Ca+2 through AFm get different values.
  return params.dEff;
}

}  // namespace xport
