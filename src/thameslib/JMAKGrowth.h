/**
@file JMAKGrowth.h
@brief Johnson-Mehl-Avrami-Kolmogorov per-voxel growth machinery.

Pure-math free functions for the JMAK model applied at the lattice-voxel
scale, plus the moment-decomposition bookkeeping that makes the
per-generation update O(1) per generation per cycle regardless of the length of
the nucleation history. Same design pattern as `NucleationRate.{h,cc}`
does for CNT — no state, no `this`, no external dependencies beyond
`<cmath>` and `PhysicalConstants.h`. Callers are `KineticController`
(generation storage + drain), the concrete kinetic models (growth-velocity
provider), and the standalone unit test.

## Physical model

JMAK inside a single voxel that has been seeded (received its first
nucleation event) at time t_c:

    Y_c(t) / V_voxel  =  alpha * integral_{t_c}^{t} J(tau) * r(tau, t)^m dtau
    r(tau, t)         =  integral_{tau}^{t} G(s) ds   (linear growth radius)
    X_c(t)            =  1 - exp(-Y_c(t) / V_voxel)

where m = n - 1 for continuous nucleation. For n = 4 (3D continuous)
m = 3 and alpha = 4*pi/3.

`Y_c(t)` is the "extended volume" inside the voxel, summed over all
nucleation events that have occurred inside since t_c. `X_c(t)` is
JMAK's Poisson-impingement correction giving the actually-transformed
fraction (impingement rescales the extended volume via 1 - exp).

## Moment decomposition (integer n = 4)

Expand `[r_acc(t) - r_acc(tau)]^3` and separate the tau-dependent piece
from the t-dependent piece. Define global moment accumulators:

    M_k(t)  =  integral_{0}^{t} J(tau) * r_acc(tau)^k dtau,  k = 0, 1, 2, 3

Then, per-generation:

    Y_c(t) / V_voxel  =  alpha * [ r_acc(t)^3 * (M_0(t) - M_0(t_c))
                                 - 3 * r_acc(t)^2 * (M_1(t) - M_1(t_c))
                                 + 3 * r_acc(t) * (M_2(t) - M_2(t_c))
                                 - (M_3(t) - M_3(t_c)) ]

Each generation stores M_k at its seed time (4 scalars). The global
accumulator holds current M_k and r_acc (5 scalars per phase).
Per-cycle work is O(1) per generation per phase, which is the point of the
whole decomposition.

## Growth velocity from surface-normal rate

The rate laws (StandardKineticModel, SaturatingRateModel, PozzolanicModel)
return a molar precipitation rate per unit surface area (mol/m^2/s).
Converting to linear growth velocity for JMAK requires the molar
volume:

    G(m/s) = r_surface_normal(mol/m^2/s) * v_molar(m^3/mol)

Physical constants and units are SI throughout.

## Verified by

`src/unit_tests/test_jmak_growth.cc` — reproduces the closed-form
JMAK curve X(t) = 1 - exp(-(pi/3) J G^3 t^4) for constant J and G, plus
sanity checks that time-varying J via moment decomposition matches
direct numerical integration.
*/

#ifndef SRC_THAMESLIB_JMAKGROWTH_H_
#define SRC_THAMESLIB_JMAKGROWTH_H_

#include "JMAKParameters.h"

namespace jmak {

/**
@brief Global moment accumulator for one nucleating phase.

Advanced each cycle by adding contributions from the current
nucleation rate J and growth velocity G. Generations snapshot the current
`M_k` at seed time to compute their per-generation extended volume.

Units:
  r_acc [m]                      -- time-integrated growth rate; equals
                                    the radius that a hypothetical
                                    nucleus born at t=0 would have
                                    reached by time t. A nucleus born at
                                    tau has radius r_acc(t) - r_acc(tau)
                                    at time t (JMAK's r(tau, t)).
  M_k   [ (1/(m^3 s)) * m^k * s ] = [ m^(k-3) * m^0 ]
                                 -- physical meaning of M_k is
                                    "number-density weighted k-th
                                     moment of growth radius"

Naming note: JMAK/Kolmogorov convention reserves the symbol G for
the growth rate (m/s), not for anything time-integrated. The
integrated quantity is a radius (m), and we call it r_acc here so
the distinction is unambiguous in the code.
*/
struct GlobalMoments {
  double r_acc;  /**< integral_0^t G(s) ds     [m]    (integrated growth
                      rate; radius of a hypothetical nucleus born at
                      t=0 by time t) */
  double M0;     /**< integral_0^t J(tau) dtau [1/m^3]           */
  double M1;     /**< integral_0^t J r_acc dtau                  */
  double M2;     /**< integral_0^t J r_acc^2 dtau                */
  double M3;     /**< integral_0^t J r_acc^3 dtau                */
};

/**
@brief Generation snapshot of the global moments at a generation's seed time.

Only stored ONCE per generation at creation; the difference against the
current global moments gives the per-generation extended volume.
*/
struct GenerationMomentsAtSeed {
  double M0_c;
  double M1_c;
  double M2_c;
  double M3_c;
};

/**
@brief Advance the global moments by one cycle at current (J, G).

    r_acc(t + dt)  =  r_acc(t) + G * dt
    M_k(t + dt)    =  M_k(t) + J * r_acc(t+dt)^k * dt

Using r_acc after the increment (midpoint of the interval would be
more accurate but negligibly so for the timesteps THAMES uses; the
end-of-interval convention is cleanest and matches the discrete-cycle
semantics elsewhere in KineticController).

@param acc   global accumulator to update (in/out)
@param J     nucleation rate this cycle [1/(m^3 s)]
@param G     linear growth velocity this cycle [m/s]
@param dt    cycle duration [s]
*/
void advanceMoments(GlobalMoments &acc, double J, double G, double dt);

/**
@brief Snapshot the global moments as a generation's seed-time state.

Called once at generation creation. `GenerationMomentsAtSeed` values are
immutable thereafter.
*/
GenerationMomentsAtSeed snapshotSeed(const GlobalMoments &acc);

/**
@brief Per-voxel extended volume for one generation at the current time.

    Y_c / V_voxel  =  alpha * [ r_acc^3 (M0 - M0_c)
                              - 3 r_acc^2 (M1 - M1_c)
                              + 3 r_acc (M2 - M2_c)
                              - (M3 - M3_c) ]

Returns 0.0 if the generation has no accumulated extended volume yet
(e.g., zero elapsed time since seed).

@param seed   moments at generation creation
@param acc    current global moments
@param jp     JMAK parameters (alpha coefficient; n unused here — this
              function assumes n = 4 via the moment decomposition)
@return Y_c / V_voxel, dimensionless (extended volume per voxel volume)
*/
double extendedVolumePerVoxel(const GenerationMomentsAtSeed &seed,
                              const GlobalMoments &acc,
                              const JMAKParameters &jp);

/**
@brief Per-voxel *extended* surface area for one generation at the current time.

    A_ext / V_voxel  =  4*pi * [ r_acc^2 (M0 - M0_c)
                               - 2 * r_acc (M1 - M1_c)
                               + (M2 - M2_c) ]

Derived from A_ext = integral J(tau) * 4*pi * [r(tau,t)]^2 dtau; expanding
the square and separating the tau-dependent factors from the t-dependent
factors, as with the volume decomposition but requiring only three moments
(M_0, M_1, M_2). Returns [1/m] (m^2 per m^3 of voxel volume).

Uses the same integer-exponent identity as `extendedVolumePerVoxel`; the
4*pi coefficient is fixed here because it comes from the sphere surface
area, not from the morphology-dependent volume coefficient alpha. If the
generalization to non-spherical crystal habit is done later, this function
would need a shape-dependent surface prefactor too.

Guards against small negative values from floating-point roundoff at
t = t_c (where all differences are near zero).

@param seed   moments at generation creation
@param acc    current global moments
@return A_ext / V_voxel [1/m]
*/
double extendedSurfaceAreaPerVoxel(const GenerationMomentsAtSeed &seed,
                                   const GlobalMoments &acc);

/**
@brief JMAK Poisson-impingement correction: X = 1 - exp(-Y/V).

Uses `-expm1` for numerical stability at small Y (X = 1 - exp(-Y) is
approximately Y for small Y and would suffer catastrophic cancellation
at exp(-tiny) ~ 1).

@param Y_over_V  extended volume divided by voxel volume
@return transformed fraction X, in [0, 1]
*/
double transformedFraction(double Y_over_V);

/**
@brief Convert a surface-normal molar rate to a linear growth velocity.

    G = r * v_molar

Both inputs are in SI. `r` is what the rate laws return (positive for
precipitation, mol/m^2/s), and `v_molar` is the molar volume of the
solid phase (m^3/mol, from GEMS via ChemicalSystem::getDCMolarVolume).

Returns 0.0 if `r` is non-positive (dissolution regime, or below
equilibrium; no growth).

@param r_surface_normal  surface-normal precipitation rate [mol/m^2/s]
@param v_molar           molar volume of the solid phase [m^3/mol]
@return linear growth velocity [m/s], always non-negative
*/
double growthVelocity(double r_surface_normal, double v_molar);

}  // namespace jmak

#endif  // SRC_THAMESLIB_JMAKGROWTH_H_
