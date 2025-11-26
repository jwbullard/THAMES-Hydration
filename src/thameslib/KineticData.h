/**
@struct KineticData
@brief Stores data about each phase possible in the system for ease of parsing
the input files.

In THAMES, phases are identified either thermodynamically--- in the
GEM data repository--- or microstructurally.  A microstructural phase can
be one, or a combination of more than on, thermodynamically defined phase.

The structure is defined to make it easier to parse the input file for the
kinetic model. It is not used elsewhere in the code.  In fact,
the same members are identified as class variables in the `KineticModel` class.

Most of the members have self-evident meanings:
    - `name` is the name of the microstructure phase
    - `microPhaseId` is the integer id for the phase in the microstructure
    - `GEMPhaseId` is the integer id of the same phase in the GEM Chemical
System Definition (CSD)
    - `DCId` is the integer id of the GEM dependent component making up the
phase
    - `type` is a string specifying whether the phase is under kinetic control
        or thermodynamic control
    - `surfaceAreaMultiplier` is a dimensionless multiplication factor for
surface area to account for unresolved internal structure, roughness, etc.
kg
    - `scaledMass` is the number of grams of the phase per 100 grams of solid
    - `k1` is the Parrot and Killoh <i>K</i><sub>1</sub> parameter for the phase
    - `k2` is the Parrot and Killoh <i>K</i><sub>2</sub> parameter for the phase
    - `k3` is the Parrot and Killoh <i>K</i><sub>3</sub> parameter for the phase
    - `n1` is the Parrot and Killoh <i>N</i><sub>1</sub> parameter for the phase
    - `n3` is the Parrot and Killoh <i>N</i><sub>3</sub> parameter for the phase
    - `dissolutionRateConst` is a generic rate constant for dissolution or
growth (flux units)
    - `diffusionRateConstEarly` is a generic rate constant for early-age
diffusion (flux units)
    - `diffusionRateConstLate` is a generic rate constant for later-age
diffusion (flux units)
    - `dissolvedUnits` is the number of DC units produced by unit dissolution
reaction
    - `siexp` is the exponent on the saturation index in the rate equation
    - `dfexp` is the exponent on the driving force term in the rate equation
    - `dorexp` is the exponent on the degree of reaction in the diffusion rate
equations
    - `ohexp` is the exponent on the hydroxyl ion activity in the rate equation
    - `sio2' is the mole fraction of SiO<sub>2</sub> in the component
    - `al2o3' is the mole fraction of Al<sub>2</sub>O<sub>3</sub> in the
component
    - `cao' is the mole fraction of CaO in the component
    - `activationEnergy` is the activation energy [J/mol/K]
    - `critDOR` is the critical degree of reaction used in the equation for
        calculating the influence of w/c ratio.
*/

#ifndef SRC_THAMESLIB_KINETICDATA_H_
#define SRC_THAMESLIB_KINETICDATA_H_

#include "global.h"

struct KineticData {
  std::string name;       /**< Name of the microstructure phase */
  int microPhaseId;       /**< Integer id of the microstructure phase */
  int GEMPhaseId;         /**< Integer id of the phase in the GEM CSD */
  int DCId;               /**< Integer id of the DC making up the phase */
  std::string type;       /**< Specifies kinetic or thermodynamic control */
  double scaledMass; /**< Mass percent on a total solids basis */
  double surfaceAreaMultiplier; /**< How much to multiply the surface area to
                                   account  for internal structure, roughness,
                                   etc. */
  double temperature;           /**< Temperature [K] */
  double reftemperature;        /**< Reference temperature [K] */
  double activationEnergy;      /**< Activation energy [J/mol/ */
  double k1;        /**< Parrot and Killoh <i>K</i><sub>1</sub> parameter */
  double k2;        /**< Parrot and Killoh <i>K</i><sub>2</sub> parameter */
  double k3;        /**< Parrot and Killoh <i>K</i><sub>3</sub> parameter */
  double n1;        /**< Parrot and Killoh <i>N</i><sub>1</sub> parameter */
  double n3;        /**< Parrot and Killoh <i>N</i><sub>3</sub> parameter */
  double dorHcoeff; /**< Lothenbach-Kulik H coefficient to compute critDOR */
  double critDOR;   /**< Critical degree of reaction (hydration) for w/c effect */
  double dissolutionRateConst;    /**< Generic rate constant [mol/m2/s] */
  double diffusionRateConstEarly; /**< Generic rate constant [mol/m2/s] */
  double diffusionRateConstLate;  /**< Generic rate constant [mol/m2/s] */
  double dissolvedUnits; /**< Number of DC units produced by unit dissolution */
  double siexp;          /**< Exponent on saturation index [unitless] */
  double dfexp;          /**< Exponent on driving force [unitless] */
  double dorexp;         /**< Exponent on degree of reaction [unitless] */
  double ohexp;          /**< Exponent on OH ion activity [unitless] */
  double loi;            /**< Loss on ignition [unitless] */
  double sio2;           /**< Mole fraction SiO2 in material */
  double al2o3;          /**< Mole fraction Al2O3 in material */
  double cao;            /**< Mole fraction CaO in material */
  double rh_;            /**< relative humidity */
  double rhFactor_;      /**< relative humidity factor, i.e. the correction
                              of the hydration rate taking into account the
                              ambient relative humidity */

  // for new C3S_dissKM <=> eq. 2 from Bullard et al.,Cem. Concr. Res. 74(2015)26:
  double kC3S; /**< dissolution rate constant [mol/m2/s] */
  double C1;   /**< unitless constant - determined by fitting experimental 
                    data from Nicoleau et al., Cem. Concr. Res. 47(2013)14 */
  double r;    /**< unitless constant - determined by fitting experimental 
                    data from Nicoleau et al., Cem. Concr. Res. 47(2013)14 */

};
#endif // SRC_THAMESLIB_KINETICDATA_H_
