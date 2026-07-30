/**
@struct TransportParameters
@brief Shell-diffusion parameters for one phase, used by Phase-3 of
       the mass-transport-kinetics extension.

Populated by KineticController::parseTransportBlock from a `transport`
sub-block inside a phase's `kinetic_data`. Consumed at simulation time
by Standard / SaturatingRate / Pozzolanic when computing the effective
dissolution or precipitation rate through a semipermeable product
shell.

The full plan (approved 2026-07-30) is at
~/.claude/plans/effervescent-forging-diffie.md. Phase 2 lands this
struct plus the parser plus xport-namespace scaffolding; Phase 3 wires
the rate models to actually consume it.

Fields:
- `dEff`: effective diffusivity of the limiting dependent component (DC)
  through the product shell around this phase [m^2/s]. Literature order:
  Ca+2 through C-S-H ~1e-13, SO4-2 through AFm ~1e-14, Ca+2 through
  calcite shell ~1e-12. Users will calibrate.
- `normalRadiusVoxels`: radius of the ball used to estimate the outward
  normal at each surface site (Bullard's method — see TransportStats.h).
  Voxel units. Too small collapses to nearest-neighbor and is noisy;
  too large lets non-local topography bias the direction. Default 2.5.
- `numShellBins`: K, the number of equal-frequency bins in the per-phase
  δ histogram. K=1 collapses to a single scalar (debugging / byte-parity);
  K=5 default (Phase-3 evaluation); K → N_sites approaches per-site
  rate summation (Option 4 fallback).
- `limitingDCName`: name of the DC whose diffusion through the shell
  rate-limits the reaction. Not an IC name — the DC. Ca+2 for
  Alite/Belite/Portlandite, SO4-2 for ettringite/gypsum, SiO2@ for
  silica fume, etc. Empty string ("") means the transport block is
  parsed but disabled (rate model falls back to no shell correction).
- `stoich`: stoichiometry of the limiting DC in the solid's dissolution
  reaction (e.g., 1 for Ca in Portlandite, 6 for Ca in ettringite,
  3 for SO4 in ettringite). Used to convert bulk SI back to C_eq.

Absence of a `transport` block in the JSON leaves the phase's
KineticData::transport optional empty; the rate model then uses its
prior (no-shell) rate calculation. Backwards-compatible with every
existing config.
*/

#ifndef SRC_THAMESLIB_TRANSPORTPARAMETERS_H_
#define SRC_THAMESLIB_TRANSPORTPARAMETERS_H_

#include <string>

struct TransportParameters {
  double dEff = 0.0;               /**< effective D of limiting DC [m^2/s] */
  double normalRadiusVoxels = 2.5; /**< ball radius for normal estimate */
  int numShellBins = 5;            /**< K bins for δ histogram; 1 = scalar */
  int maxWalkSteps = 50;           /**< walk cap; beyond this shell → ∞ */
  std::string limitingDCName;      /**< e.g. "Ca+2"; empty = disabled */
  double stoich = 1.0;             /**< stoichiometric coefficient */
};

#endif  // SRC_THAMESLIB_TRANSPORTPARAMETERS_H_
