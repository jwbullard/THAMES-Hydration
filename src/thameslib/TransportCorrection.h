/**
@file TransportCorrection.h
@brief Series-resistance rate math for shell-diffusion-limited kinetics.

Free functions in namespace `xport`. Consumed by the concrete kinetic
models (Standard / SaturatingRate / Pozzolanic) to compute an effective
rate that combines the surface-reaction kinetic law with a Fickian
diffusion resistance through the product shell.

Phase 2 lands this header with declarations only + stubbed bodies;
Phase 3 fills in the Newton solver, the D_eff picker, and the diagnostic
r_kin / r_diff exposure. Keeping the API surface stable across the two
phases so the rate models can wire against this header once and not
change again when Phase 3 fills in the bodies.

No THAMES dependencies. Header-only from the model perspective; only
this file's translation unit pulls in <cmath>.
*/

#ifndef SRC_THAMESLIB_TRANSPORTCORRECTION_H_
#define SRC_THAMESLIB_TRANSPORTCORRECTION_H_

#include "TransportParameters.h"
#include "TransportStats.h"

namespace xport {

/**
@brief Kinetic-only surface reaction rate for logging / calibration.

Returns k · area · driving_force. Trivial wrapper — its purpose is
to make the "r_kinetic" contribution to the series resistance visible
in verbose logs, so calibration can distinguish "wrong k" from
"wrong D_eff" when the effective rate deviates from measurement.

Units are the caller's responsibility (typically mol / s in SI).
*/
double kineticRate(double k, double area, double driving_force);

/**
@brief Fickian diffusion rate through a shell of thickness delta.

Returns dEff · area · deltaC / delta. Complementary to kineticRate;
exposed for the same calibration-transparency reason.

Units are the caller's responsibility (mol / s if deltaC is in
mol/m^3, area is m^2, delta is m, dEff is m^2/s).
*/
double diffusionRate(double dEff, double area, double deltaC, double delta);

/**
@brief Solve the steady-state flux-balance equation for the reactant-
       surface concentration under a linear-Fick driving-force closure.

Assumes the kinetic law is f(Ω) = 1 - Ω (dissolution) or f(Ω) = Ω - 1
(precipitation) near equilibrium. Returns C_surf that satisfies

    k · (1 - C_surf / C_eq)  =  dEff · (C_surf - C_bulk) / delta

which is the closed form for the linear regime.
*/
double solveSurfaceConcentrationLinear(double k, double C_eq, double dEff,
                                       double delta, double C_bulk);

/**
@brief Solve the steady-state flux-balance equation for arbitrary
       driving-force f(Ω) via Brent's root-finder.

Solves for C_surf that satisfies

    k · f(C_surf / C_eq)  =  dEff · (C_surf - C_bulk) / delta

The caller supplies a driving-force functor. For dissolution, f is
positive at Ω < 1 and zero at Ω = 1 (Standard's (1-Ω^p)^q, SR's
saturating form). For precipitation the sign convention is caller's
choice — the function only cares that the flux balance holds.

Brent bracketing:
  - Dissolution (C_bulk < C_eq): C_surf is bracketed by [C_bulk, C_eq].
  - Precipitation (C_bulk > C_eq): bracket is [C_eq, C_bulk].
  - Equilibrium (C_bulk == C_eq): C_surf = C_bulk = C_eq; return early.

Converges to relative tolerance 1e-8 within ~30 iterations for typical
kinetic laws. Returns C_bulk (no-op fallback) if the bracket is
degenerate or the function values don't have opposite signs at the
endpoints (shouldn't happen for well-behaved f but be defensive).

Cost: ~5-10 f() evaluations per call for typical shell parameters.
Suitable for K bins × N phases per cycle without dominating wall time.
*/
double solveSurfaceConcentration(double k, double C_eq, double dEff,
                                 double delta, double C_bulk,
                                 double (*f_driving)(double omega));

/**
@brief Pick D_eff for a specific shell composition.

Looks up the phase pair (reactant, shell) in the transport
parameters and returns the corresponding effective diffusivity.
Phase 2 stub returns the block's global dEff regardless of shell
composition. Phase 3 will introduce a per-shell-phase D_eff map so
different shell materials (C-S-H vs AFm vs calcite around the same
reactant) get different values.

@param shellPhaseId  microstructure phase id of the dominant shell
                     traversed in the walk (from ShellBin)
@param params        the phase's transport parameters
@return effective diffusivity for that (reactant, shell) pair
*/
double pickDEff(int shellPhaseId, const TransportParameters &params);

}  // namespace xport

#endif  // SRC_THAMESLIB_TRANSPORTCORRECTION_H_
