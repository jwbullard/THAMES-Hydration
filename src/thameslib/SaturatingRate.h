/**
@file SaturatingRate.h
@brief Saturating dissolution/precipitation rate law (Bullard 2015 / Han 2025).

Namespace-scoped free functions with no hidden state. Callers use them
from `SaturatingRateModel::calculateKineticStep`. Unit-tested by
`src/unit_tests/test_saturating_rate.cc` against the Portlandite
values Han et al. published at 24 °C.

Origin: Bullard et al. 2015 CCR Eq. (2) introduced the form
  r = k · A_eff · (1 − exp[−((ln K − ln Q)/C₁)^r_exp])
for alite dissolution. Han et al. 2025 CEJ Eq. (7) re-derived the
algebraically equivalent form
  r = k · (1 − exp[−(−B ln Ω)^n])
with `B ≡ 1/C₁`, `n ≡ r_exp`, `Ω ≡ Q/K`. We use Han's notation because
the current Portlandite calibration (k, B, n) uses it.

Convention: S = IAP/Ksp (linear ratio, no log). Matches
`ChemicalSystem::getMicroPhaseSI` which returns the same linear form
(GEMS internally uses log10 but converts before storing).

Rate saturates at k in the far-from-equilibrium limit; drops smoothly
to zero at equilibrium (Ω = 1). See `docs/CNT_ARCHITECTURE.md` for the
integration architecture; SaturatingRateModel is a 4th concrete
KineticModel subclass alongside Standard, Pozzolanic, ParrotKilloh.

Physical units are SI:
  k     [mol/(m^2 · s)]
  B, n  [dimensionless]
  Omega [dimensionless]
  T     [K]
  Ea    [J/mol]
*/

#ifndef SRC_THAMESLIB_SATURATINGRATE_H_
#define SRC_THAMESLIB_SATURATINGRATE_H_

namespace sat {

/// Dissolution rate at Omega <= 1 (undersaturated / at equilibrium).
/// Returns a positive value in mol / m^2 / s. Callers apply the sign
/// convention for their DC-moles bookkeeping.
/// Returns 0 exactly at Omega = 1 (equilibrium).
/// Callers should NOT call this with Omega > 1 — use precipitationRate.
double dissolutionRate(double k, double B, double n, double Omega);

/// Precipitation rate at Omega >= 1 (supersaturated / at equilibrium).
/// Uses the mirrored form
///     r = k * (1 - exp[-(B ln Omega)^n])
/// which is numerically well-behaved for Omega > 1 (positive log
/// keeps the exponentiation base positive).
/// Returns a positive value; callers apply the sign for growth.
/// Returns 0 exactly at Omega = 1 (equilibrium).
/// Callers should NOT call this with Omega < 1 — use dissolutionRate.
double precipitationRate(double k, double B, double n, double Omega);

/// Temperature-scaled rate constant. Returns k_ref * exp[Ea/R * (1/T_ref - 1/T)],
/// i.e. the Arrhenius factor by which the rate at temperature T
/// differs from the rate at temperature T_ref, times the given
/// reference rate constant. Both T values in Kelvin; Ea in J/mol.
double arrheniusScale(double k_ref, double T_ref, double T,
                      double dHact);

}  // namespace sat

#endif  // SRC_THAMESLIB_SATURATINGRATE_H_
