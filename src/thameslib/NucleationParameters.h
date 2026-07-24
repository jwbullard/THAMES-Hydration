/**
@struct NucleationParameters
@brief Classical Nucleation Theory (CNT) parameters for one phase.

Populated by KineticController::parseNucleationBlock from a `nucleation`
sub-block inside a phase's `kinetic_data`. Read at simulation time by
StandardKineticModel or PozzolanicModel when computing nucleation rate.

Contact angle convention: theta_deg = 180 means homogeneous nucleation
(placement over electrolyte voxels). theta_deg in [1, 179] means
heterogeneous nucleation (placement over substrate-interface voxels).
This is a placement-branch flag, not just a barrier modulator: physically
theta = 180 IS the homogeneous limit (Turnbull-Vonnegut f(theta) = 1),
but the placement location differs.
*/

#ifndef SRC_THAMESLIB_NUCLEATIONPARAMETERS_H_
#define SRC_THAMESLIB_NUCLEATIONPARAMETERS_H_

struct NucleationParameters {
  double gamma;   /**< Interfacial energy [J/m^2] */
  int theta_deg;  /**< Contact angle [degrees]; 180 = homogeneous placement */
  double A0;      /**< Kinetic prefactor [1/(m^3 s)] */
};

#endif  // SRC_THAMESLIB_NUCLEATIONPARAMETERS_H_
