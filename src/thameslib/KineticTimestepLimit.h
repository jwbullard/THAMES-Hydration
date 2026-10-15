/**
@file KineticTimestepLimit.h
@brief The kinetics-based timestep limit, with its one-voxel floor.

KineticController::computeKineticsBasedMaxTimestep limits each step so that
no kinetically controlled DC changes by more than a fraction f
(maxRelativeChange, default 0.05) of its amount:

    dt = f * max(|n|, n_vox) / |r|

with n the DC moles and r the molar rate, both in the per-100-g-solid frame
of DCMoles_ (Lattice::normalizePhaseMasses), and n_vox the moles of that DC
that would fill one voxel, in the same frame:

    n_vox = 100 / (N * 1e6 * V_m * rho_0)

N is the number of lattice sites, V_m the DC molar volume in m^3/mol, and
rho_0 the initial solid mass in g per cm^3 of microstructure
(Lattice::getInitSolidMass). One voxel is 1/N of each cm^3 of
microstructure, i.e. 1e-6/N m^3, so it holds 1e-6/(N V_m) mol per cm^3;
multiplying by 100/rho_0 converts to mol per 100 g of solid. This is the
conversion the CNT placement blocks in KineticController.cc already use,
with the molar mass cancelled.

Why the floor. Without it, a phase whose rate does not fall in proportion to
its remaining amount gets a dt sized to remove f of what remains, so the
remainder falls by (1 - f) each step and dt with it. Simulated time then
converges on a limit it never reaches (observed 2026-10-14: Anhydrite near
214 h, Gypsum and Bassanite near 93.7 h, cem151 w/c 0.32 sealed). With the
floor, once the remainder is below one voxel the next step requests more
than is left, and the DC-depletion clamp in calculateKineticStep dissolves
exactly the remainder through commitSolidICTransfer, conserving mass. The
departure from the kinetic law is bounded by one voxel's worth of each phase
finishing within one step. The same floor serves a phase precipitating from a
trace amount, whose dt the unfloored rule would also throttle.

Header only, with no THAMES dependencies, so it can be unit tested on its
own (unit_tests/test_kinetic_timestep_limit.cc).
*/

#ifndef SRC_THAMESLIB_KINETICTIMESTEPLIMIT_H_
#define SRC_THAMESLIB_KINETICTIMESTEPLIMIT_H_

#include <algorithm>
#include <cmath>

namespace kinetic_dt {

/**
@brief Moles of a DC that fill one voxel, per 100 g of initial solid.

@param numSites number of lattice sites N
@param vMolar DC molar volume [m^3/mol]
@param initSolidMass initial solid mass [g per cm^3 of microstructure]
@return n_vox [mol/100 g], or 0 (no floor) if any input is not positive
*/
inline double molesPerVoxel(int numSites, double vMolar,
                            double initSolidMass) {
  if (numSites <= 0 || vMolar <= 0.0 || initSolidMass <= 0.0) {
    return 0.0;
  }
  return 100.0 /
         (static_cast<double>(numSites) * 1.0e6 * vMolar * initSolidMass);
}

/**
@brief Largest timestep that changes a DC by at most a fraction of its
amount, with that amount floored at one voxel's worth.

@param dcMoles current DC moles [mol/100 g]
@param molesPerVoxel one voxel's worth of the DC [mol/100 g]; 0 disables
the floor
@param molarRate current molar rate [mol/100 g/h], nonzero
@param maxRelativeChange allowed fractional change per step
@return dt [h]
*/
inline double timestepLimit(double dcMoles, double molesPerVoxel,
                            double molarRate, double maxRelativeChange) {
  const double refMoles = std::max(std::fabs(dcMoles), molesPerVoxel);
  return maxRelativeChange * refMoles / std::fabs(molarRate);
}

} // namespace kinetic_dt

#endif // SRC_THAMESLIB_KINETICTIMESTEPLIMIT_H_
