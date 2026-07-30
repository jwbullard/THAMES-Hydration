/**
@file TransportStats.h
@brief POD types for shell-thickness distribution and normal-vector
       estimation used by mass-transport kinetics.

Shell-thickness δ is computed at each dissolution-surface voxel by
(1) estimating the outward normal via a porosity-weighted centroid
inside a ball of radius r (Bullard's method — see plan phase 1) and
(2) walking outward one voxel at a time through the product shell
until reaching an electrolyte voxel or a step cap. The per-site δ
distribution is aggregated per phase into a K-bin equal-frequency
histogram so the nonlinear rate law can be evaluated bin-wise in
Phase 3 without collapsing to a single scalar (see
`docs/POST_ALPHA_TODOS.md` and `docs/transport_kinetics_brainstorm.md`
for the physics rationale).

No THAMES dependencies. Header-only, testable in isolation.
*/

#ifndef SRC_THAMESLIB_TRANSPORTSTATS_H_
#define SRC_THAMESLIB_TRANSPORTSTATS_H_

#include <cstdint>
#include <limits>
#include <vector>

namespace xport {

/**
@brief Minimal 3D vector for lattice-space directions. Unit vectors
       stored in double so we don't lose precision when snapping to
       the nearest face-neighbor.
*/
struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

/**
@brief Result of a single shell walk from one surface site.

`deltaVoxels` counts product-phase voxels traversed before reaching
electrolyte (0 if the surface site is directly adjacent to electrolyte).
`reachedElectrolyte` is false if the walk terminated at the step cap
without hitting electrolyte — those sites are treated as having δ
above the maximum resolvable shell thickness and are excluded from the
harmonic-mean aggregate (they contribute zero to the diffusion rate).
`shellPhaseIds` records the microstructure-phase id of each product
voxel visited, in walk order, so Phase 3 can pick a per-bin D_eff from
the dominant shell composition.
*/
struct ShellHit {
  int deltaVoxels = 0;
  bool reachedElectrolyte = false;
  std::vector<int> shellPhaseIds;  // one entry per traversed voxel
};

/**
@brief One K-bin of the per-phase δ histogram.

`deltaRep` is the bin representative δ in physical units (m), taken
as the median of the sites that fall in the bin (robust to outliers
within a bin). `siteFraction` is the fraction of the phase's total
dissolution sites that fall in this bin. `dominantShellPhaseId` is
the mode of the intermediate-phase histogram across the sites in the
bin — the shell composition Phase 3 will use for D_eff lookup.
*/
struct ShellBin {
  double deltaRep = 0.0;              // representative shell thickness [m]
  double siteFraction = 0.0;          // in [0, 1]
  int dominantShellPhaseId = -1;      // -1 if no shell traversed (δ=0)
};

/**
@brief Per-phase shell-thickness statistics for one cycle.

`bins` is the K-bin equal-frequency histogram consumed by Phase 3's
rate calculation. `deltaHarmonicRaw` and `deltaArithmeticRaw` are
computed from the full unbinned distribution and retained for
diagnostics (calibration debugging, visual verification that the bins
are capturing the distribution). `numSitesTotal` includes both
reached-electrolyte sites AND unreached (trapped-by-thick-shell) sites;
`numSitesReached` is the count that contributed to the aggregates.
Difference = the trapped-by-thick-shell tail, which is worth logging
because it flags a scenario where per-site (Option 4) may be needed.
*/
struct ShellStats {
  std::vector<ShellBin> bins;         // length K
  double deltaHarmonicRaw = 0.0;      // [m] from reached sites only
  double deltaArithmeticRaw = 0.0;    // [m] from reached sites only
  int numSitesTotal = 0;
  int numSitesReached = 0;
};

/**
@brief Aggregate a raw per-site δ distribution into a K-bin
       equal-frequency histogram.

Inputs:
  deltasMeters  — one entry per dissolution surface site of a phase;
                  δ in meters. Non-negative. May include +infinity
                  sentinels for sites whose walk did not reach
                  electrolyte within the step cap; those are dropped
                  from the aggregate.
  shellPhaseIds — same length as deltasMeters; per-site dominant
                  shell phase id (mode of the phases traversed by
                  the walk). Ignored for +inf-δ sites.
  K             — number of histogram bins. K=1 collapses to a single
                  scalar (median δ + arithmetic-mode shell composition).

Returns a ShellStats populated with the K bins, the diagnostic
harmonic and arithmetic means over reached sites, and the reached /
total site counts.

Purely functional; deterministic on identical input. No THAMES types.
*/
ShellStats aggregateShellDistribution(
    const std::vector<double> &deltasMeters,
    const std::vector<int> &shellPhaseIds,
    int K);

/**
@brief Return true if the finite-δ sites in the distribution are
       inadequately captured by K bins.

Rule of thumb used by the fallback-to-Option-4 detector: if the
ratio of arithmetic mean to harmonic mean exceeds a threshold, the
distribution has a long-tailed or bimodal shape that fixed-K
equal-frequency binning cannot approximate to acceptable tolerance.
Returns false for degenerate distributions (all δ = 0, no sites, etc).
*/
bool distributionIsPathological(const ShellStats &stats,
                                double arithToHarmonicRatioThreshold = 5.0);

}  // namespace xport

#endif  // SRC_THAMESLIB_TRANSPORTSTATS_H_
