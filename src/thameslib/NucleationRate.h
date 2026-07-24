/**
@file NucleationRate.h
@brief Classical Nucleation Theory (CNT) rate formulas as free functions.

Namespace-scoped free functions with no hidden state. Callers are
StandardKineticModel and PozzolanicModel via their `computeNucleationVoxels`
methods. See `docs/CNT_ARCHITECTURE.md` for how these functions fit into the
per-cycle placement pipeline.

Convention: S = IAP/Ksp (no log). ChemicalSystem::getMicroPhaseSI already
returns S in this linear form (it converts internally from the GEMS log10
convention). Callers pass S directly.

Physical units are SI throughout:
  gamma [J/m^2]
  v_m   [m^3/mol] (molar volume from GEMS; converted to per-molecule
                   internally by dividing by Avogadro's number)
  T     [K]
  dt    [s]
  V_electrolyte [m^3]
  V_voxel       [m^3]

Verified against the calibrated portlandite defaults from
~/Research/THAMES-Tests-2026/Scripts/NucleationCNT-Prototype.ipynb by the
standalone binary at src/unit_tests/test_nucleation_rate.cc.
*/

#ifndef SRC_THAMESLIB_NUCLEATIONRATE_H_
#define SRC_THAMESLIB_NUCLEATIONRATE_H_

#include "NucleationParameters.h"

namespace cnt {

/// Turnbull-Vonnegut heterogeneity factor: (2 + cos theta)(1 - cos theta)^2 / 4.
/// f(180) = 1 exactly (homogeneous limit).
double fTheta(int theta_deg);

/// Critical nucleus radius r* = 2 gamma v_m / (k_B T ln S), meters.
/// Requires ln_S > 0; caller must guard against undersaturation.
double criticalRadius(double gamma, double v_m, double T_K, double ln_S);

/// Homogeneous barrier deltaG*_hom = 16 pi gamma^3 v_m^2 / [3 (k_B T ln S)^2], Joules.
double deltaGStarHom(double gamma, double v_m, double T_K, double ln_S);

/// Full nucleation rate J = A0 * exp(-f(theta) * deltaG*_hom / (k_B T)), [1/(m^3 s)].
/// Returns 0.0 if S <= 1 (undersaturated / at equilibrium).
double nucleationRate(const NucleationParameters &np,
                      double v_m, double T_K, double S);

/// Fractional voxels expected to nucleate in one cycle:
///   N = J * V_electrolyte * dt * V_crit / V_voxel
/// Returns a double so the caller can accumulate sub-voxel amounts.
/// Returns 0.0 if S <= 1.
double voxelsPerCycle(const NucleationParameters &np,
                      double v_m, double T_K, double S,
                      double dt_seconds,
                      double V_electrolyte, double V_voxel);

}  // namespace cnt

#endif  // SRC_THAMESLIB_NUCLEATIONRATE_H_
