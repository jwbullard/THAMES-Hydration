/**
@file TransportCorrection.h
@brief Series-resistance rate math for shell-diffusion-limited kinetics.

Free functions in namespace `xport`. Consumed by the concrete kinetic
models (Standard / SaturatingRate / Pozzolanic) to compute an effective
rate that combines the surface-reaction kinetic law with a Fickian
diffusion resistance through the product shell.

Framework landed 2026-07-30 (Session 55). Design conversation preserved
verbatim at `docs/transport_kinetics_brainstorm.md`. Opt-in per phase
via a `transport` sub-block under `kinetic_data`; absent that block,
kinetic models fall back to their prior no-shell rate calculation.

Current state of the pieces exposed by this header:
  - `shellCorrectionFactor` — the production path. Closed-form
    linear-rate correction summed over the per-phase K-bin δ
    histogram. All three kinetic models call this today.
  - `solveSurfaceConcentration` — Brent iteration on an arbitrary
    driving-force functor. Implemented and unit-tested (see
    `test_transport_correction`) but NOT called from any kinetic
    model as of 2026-09-12. Reserved for a future switch when we
    want the true SR nonlinearity at the flux balance instead of
    the linear closed-form. Tracked in `docs/POST_ALPHA_TODOS.md`.
  - `pickDEff` — currently returns the block's global `dEff`
    regardless of shell composition. A per-shell-phase D_eff map
    (Ca+2 through C-S-H vs Ca+2 through AFm, etc.) is a deferred
    refinement; the plumbing to select on `bin.dominantShellPhaseId`
    is already in `shellCorrectionFactor`.

C_eq caveat for readers coming in cold. The kinetic models compute
C_eq at the call site from the reactant's own saturation index via
C_eq = C_bulk / SI^(1/stoich). Numerically physical values (mM to
hundreds of mM for typical species) depend on the equilibrium
constant being physically correct. Pre-Session-56, ln K(C3S) was
off by ~26 units, S(Alite) was ~1e-14 in a real paste, and C_eq
came out ~1e14 mol/m³ — an artifact, not a formulation error.
The S56 / S57 ln K corrections restored physical C_eq values;
future clinker phases migrating from PK to SR/Standard need the
same K audit before their transport blocks will behave sensibly.

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

NOT CURRENTLY WIRED: all three kinetic models use
`shellCorrectionFactor` (linear closed-form) as of 2026-09-12.
Wiring this true-nonlinear path is a POST_ALPHA refinement — see
"Shell-diffusion long-duration validation run" in POST_ALPHA_TODOS
and its "Refinements that may be worth landing first" list.

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

CURRENT IMPLEMENTATION: returns the block's global `dEff`
regardless of shell composition (shellPhaseId argument unused).
The per-shell-phase D_eff map that would let e.g. Ca+2 through
C-S-H and Ca+2 through AFm get different values is a deferred
refinement — the calling loop in `shellCorrectionFactor` already
passes `bin.dominantShellPhaseId` through, so wiring a real map
here is the only remaining step.

@param shellPhaseId  microstructure phase id of the dominant shell
                     traversed in the walk (from ShellBin); currently
                     ignored, but preserved in the API so the future
                     per-shell-phase lookup does not need a call-site
                     change.
@param params        the phase's transport parameters
@return effective diffusivity for that (reactant, shell) pair
*/
double pickDEff(int shellPhaseId, const TransportParameters &params);

/**
@brief Compute the bin-weighted series-resistance rate correction
       factor for a phase with a shell.

Consumes the K-bin ShellStats produced by
Lattice::computeShellStats. For each bin with representative
thickness δ_bin and dominant shell composition c_bin, computes the
per-bin Damköhler number

    Da_bin = k * δ_bin / (D_eff(c_bin) * C_eq)

and the per-bin correction factor 1 / (1 + Da_bin). Weights by
bin.siteFraction and sums.

Returned factor multiplies the caller's kinetic-only rate to give
the shell-corrected rate:

    r_effective = r_kinetic_at_bulk_omega * factor

Derivation of the closed form (for the linear rate law
`r = k · (1 − C_surf/C_eq)`, dissolution): steady-state flux
balance across the shell reads

    k · (1 − C_surf/C_eq)  =  D_eff · (C_surf − C_bulk) / δ

Solving for C_surf and substituting back into `r`, then dividing
by the no-shell rate `r_bulk = k · (1 − C_bulk/C_eq)`, collapses
to `1 / (1 + Da)` with `Da = k · δ / (D_eff · C_eq)`. Summing
that per-bin correction weighted by `siteFraction` gives the
returned factor.

For the linear driving-force this is EXACT (per-bin steady-state
flux balance). For nonlinear f (Standard's (1-Ω^p)^q, SR's
saturating form), it is a first-order approximation that is exact
near equilibrium and degrades gracefully far from equilibrium.
Wiring `solveSurfaceConcentration` per bin for exact nonlinear
behavior is deferred (see that function's doc); the API of this
function will not change when that switch happens.

C_eq numerical sanity: the caller derives C_eq from the reactant's
saturation index (`C_eq = C_bulk / SI^(1/stoich)`). Values only
land in a physical range (mM to hundreds of mM) when the phase's
equilibrium constant is itself physical — see the file-level
docstring's note about the S56 / S57 ln K corrections.

Guards: returns 1.0 (no correction) if ShellStats is empty, if any
bin has degenerate δ_bin ≤ 0 or C_eq ≤ 0, or if D_eff comes back
non-positive from pickDEff.

@param k        kinetic rate constant [mol/m^2/s]
@param C_eq     equilibrium concentration of the limiting DC
                [mol/m^3, or any consistent unit]
@param stats    per-phase K-bin ShellStats
@param params   the phase's transport parameters (for pickDEff)
@return correction factor in (0, 1]
*/
double shellCorrectionFactor(double k, double C_eq,
                             const ShellStats &stats,
                             const TransportParameters &params);

}  // namespace xport

#endif  // SRC_THAMESLIB_TRANSPORTCORRECTION_H_
