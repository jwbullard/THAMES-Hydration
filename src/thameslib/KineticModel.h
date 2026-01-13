/**
@file KineticModel.h
@brief Declaration of the KineticModel class.

@section Introduction

In THAMES, the `KineticModel` class can be perceived as the engine that
calculates the kinetic changes in the system during a given time increment.  The
primary kinetic aspect that is calculated is the extent of dissolution of
mineral phases in the original cement.  This is just the base class.  It is not
used.

*/

#ifndef SRC_THAMESLIB_KINETICMODEL_H_
#define SRC_THAMESLIB_KINETICMODEL_H_

#include "global.h"
#include "Exceptions.h"
#include "ChemicalSystem.h"
#include "KineticData.h"
#include "Lattice.h"

/**
@class KineticModel
@brief Base class for all kinetic models

THAMES allows some flexibility in defining different types of kinetic models.
*/

class KineticModel {

protected:
  std::string modelName_; /**< The kinetic model name in the kinetic model */
  int numPhases_;         /**< Total number of phases in the kinetic model */
  ChemicalSystem *chemSys_;  /**< Pointer to the ChemicalSystem object for
                                  this simulation */
  Lattice *lattice_;         /**< Pointer to the lattice object holding the
                                  microstructure */
  double initSolidMass_;     /**< initial total mass of solids controlled by
                                  this model [g] */
  double temperature_;       /**< Temperature [K] */
  double refT_;              /**< Reference temperature for PK model [K] */
  double sulfateAttackTime_; /**< Time at which sulfate attack simulation
                                  starts [h] */
  double leachTime_;         /**< Time at which leaching simulation starts
                                  [h] */

  std::string name_;      /**< Name of phase controlled by this kinetic model */
  int microPhaseId_;      /**< Microstructure id controlled by this model */
  int DCId_;              /**< List of DC ids from the ChemicalSystem object */
  int GEMPhaseId_;        /**< List of phase ids from the ChemicalSystem object */

  std::vector<std::string> ICName_;      /**< Names of ICs */
  std::vector<std::string> DCName_;      /**< Names of DCs */
  double scaledMass_;             /**< Phase mass percent, total solids basis */
  double initScaledMass_;         /**< Initial phase scaled mass */
  double activationEnergy_;       /**< Apparent activation energy for the reaction
                                       [J/mol/K] */
  double specificSurfaceArea_;    /**< Specific surface area (m2/kg) */
  double refSpecificSurfaceArea_; /**< Reference specific surface area (m2/kg) */
  double ssaFactor_;              /**< Reference specific surface area (m2/kg) */
  double degreeOfReaction_;       /**< Degree of reaction of this component (mass
                                       basis) */
  double lossOnIgnition_;         /**< Loss on ignition of this component (ignited
                                       mass basis) */
  bool verbose_;                  /**< Flag for verbose output */
  bool warning_;                  /**< Flag for warnining output */

public:
  /**
  @brief Default constructor.

  This constructor is not used in THAMES.  It just establishes default values
  for all the member variables.

  @note NOT USED.
  */
  KineticModel();

  /**
  @brief Destructor does nothing.
  */
  virtual ~KineticModel() {}

  /**
  @brief Get the ChemicalSystem object for the simulation used by the kinetic
  model.

  @note NOT USED.

  @return a pointer to the ChemicalSystem object
  */
  ChemicalSystem *getChemSys() const { return chemSys_; }

  /**
  @brief Set the multiplicative adjustment to clinker phase rate constant due to
  pozzolans

  @param pfk is the specific surface area [m<sup>2</sup>/kg]
  */
  virtual void setPfk(double pfk) {
    // Most components do not have a pfk variable, generically return
    return;
  }

  /**
  @brief Get the type of kinetic model

  @return a string indicating the model type
  */
  virtual std::string getType() const { return (GenericType); }

  /**
  @brief Set the specific surface area

  @note NOT USED.

  @param sval is the specific surface area [m<sup>2</sup>/kg]
  */
  void setSpecificSurfaceArea(double sval) { specificSurfaceArea_ = sval; }

  /**
  @brief Get the specific surface area.

  @note NOT USED.

  @return the specific surface area fineness [m<sup>2</sup>/kg]
  */
  double getSpecificSurfaceArea() const { return specificSurfaceArea_; }

  /**
  @brief Set the reference specific surface area.

  The value set in the Parrot and Killoh model is 385 m<sup>2</sup>/kg, and
  there is no particular reason to change it.

  @note NOT USED.

  @param rsval is the reference specific surface area [m<sup>2</sup>/kg]
  */
  void setRefSpecificSurfaceArea(const double rsval) {
    refSpecificSurfaceArea_ = rsval;
  }

  /**
  @brief Get the reference specific surface area.

  @note NOT USED.

  @return the reference specific surface area [m<sup>2</sup>/kg]
  */
  double getRefSpecificSurfaceArea() const { return refSpecificSurfaceArea_; }

  /**
  @brief Set the ratio of the true specific surface area to the model reference
  value.

  @note NOT USED.

  @param sfact is the ratio of the actual specific surface area to the reference
  value
  */
  void setSsaFactor(const double sfact) { ssaFactor_ = sfact; }

  /**
  @brief Get the ratio of the true specific surface area to the model reference
  value.

  @note NOT USED.

  @return the ratio of the actual specific surface area to the reference value
  */
  double getSsaFactor() const { return ssaFactor_; }

  /**
  @brief Set the degree of reaction of this component (mass basis)

  @note NOT USED.

  @param dor is the degree of reaction to set
  */
  virtual void setDegreeOfReaction(const double dor) {
    degreeOfReaction_ = dor >= 0.0 ? dor : 0.0;
  }

  /**
  @brief Get the degree of reaction of this component (mass basis)

  @note NOT USED.

  @return the degree of reaction of this component
  */
  virtual double getDegreeOfReaction() const { return degreeOfReaction_; }

  /**
  @brief Set the loss on ignition of this component (ignited mass basis)

  @note NOT USED.

  @param loi is the degree of reaction to set
  */
  virtual void setLossOnIgnition(const double loi) { lossOnIgnition_ = loi; }

  /**
  @brief Get the loss on ignition of this component (ignited mass basis)

  @note NOT USED.

  @return the loss on ignition of this component
  */
  virtual double getLossOnIgnition() const { return lossOnIgnition_; }

  /**
  @brief Set the SiO2 content of pozzolanic materials

  @note This is a pure virtual function

  @param sio2
  */
  // virtual void setSio2(const double sio2) { return; }

  /**
  @brief Get the SiO2 content of pozzolanic materials

  @return the SiO2 content (mass percent)
  */
  virtual double getSio2() const { return 0.0; }

  /**
  @brief Compute normalized initial microstructure phase masses

  Given the initial masses of all phases in the microstructure,
  this method scales them to 100 grams of solid.  In the process,
  this method also sets the initial moles of water in the
  chemical system definition.
  */

  // void getPhaseMasses(void);

  /**
  @brief Get the microstructure id in the KineticModel.

  @return the list of all microstructure ids.
  */
  int getMicroPhaseId() const { return microPhaseId_; }

  /**
  @brief Set the DC index for liquid water

  @note NOT USED.

  @param waterid is the index value to use
  */
  // void setWaterId(const int waterid) { waterId_ = waterid; }

  /**
  @brief Get the DC index for liquid water

  @return the DC index for liquid water
  */
  // int getWaterId() const { return waterId_; }

  /**
  @brief Set the number of ICs

  @note NOT USED.

  @param icnum is the number of ICs to specify
  */
  // void setICNum(const int icnum) { ICNum_ = icnum; }

  /**
  @brief Get the number of ICs

  @return the number of ICs
  */
  // int getICNum() const { return ICNum_; }

  /**
  @brief Set the number of DCs

  @note NOT USED.

  @param dcnum is number of DCs to specify
  */
  // void setDCNum(const int dcnum) { DCNum_ = dcnum; }

  /**
  @brief Get the number of DCs

  @return the number of DCs
  */
  // int getDCNum() const { return DCNum_; }

  int getDCId() const { return DCId_; }

  /**
  @brief Set the number of GEM phases

  @note NOT USED.

  @param gpnum is number of GEM phases to specify
  */
  // void setGEMPhaseNum(const int gpnum) { GEMPhaseNum_ = gpnum; }

  /**
  @brief Get the number of GEM phases

  @return the number of GEM phases
  */
  // int getGEMPhaseNum() const { return GEMPhaseNum_; }

  /**
  @brief Set the IC names

  @note NOT USED.

  @param icname is the list of IC names
  */
  // void setICName(std::vector<std::string> icname) { ICName_ = icname; }

  /**
  @brief Get the IC names

  @return the list of IC names
  */
  // std::vector<std::string> getICName() const { return ICName_; }

  /**
  @brief Set the DC names

  @note NOT USED.

  @param dcname is the list of DC names
  */
  // void setDCName(std::vector<std::string> dcname) { DCName_ = dcname; }

  /**
  @brief Get the DC names

  @return the list of DC names
  */
  // std::vector<std::string> getDCName() const { return DCName_; }

  /**
  @brief Set the total number of phases in the kinetic model.

  @note NOT USED.

  @param numphases is the total number of phases in the kinetic model
  */
  // void setNumPhases(const unsigned int numphases) { numPhases_ = numphases; }

  /**
  @brief Get the total number of phases in the kinetic model.

  @note NOT USED.

  @return the total number of phases in the kinetic model
  */
  // int getNumPhases() const { return numPhases_; }

  /**
  @brief Set the simulation time at which to begin external sulfate attack.

  @param sattacktime is the simulation time to begin sulfate attack [hours]
  */
  // void setSulfateAttackTime(double sattacktime) { sulfateAttackTime_ =
  // sattacktime; }

  /**
  @brief Get the simulation time at which to begin external sulfate attack.

  @note NOT USED.

  @return the simulation time to begin sulfate attack [hours]
  */
  // double getSulfateAttackTime(void) const { return sulfateAttackTime_; }

  /**
  @brief Set the simulation time at which to begin leaching.

  @param leachtime is the simulation time to begin leaching [hours]
  */
  // void setLeachTime(double leachtime) { leachTime_ = leachtime; }

  /**
  @brief Get the simulation time at which to begin leaching.

  @note NOT USED.

  @return the simulation time to begin leaching [hours]
  */
  // double getLeachTime(void) const { return leachTime_; }

  /**
  @brief Get the list of phase names used by the kinetic model.

  @note NOT USED.

  @return the vector of names of phases in the kinetic model
  */
  std::string getName() const { return name_; }

  /**
  @brief Get the list of activation energies for the phases in the kinetic
  model.

  @note NOT USED.

  @return the vector of activation energies [J/mol/K]
  */
  // double getActivationEnergy() const { return activationEnergy_; }

  /**
  @brief Set the absolute temperature.

  @note NOT USED.

  @param tval is the absolute temperature [K]
  */
  // void setTemperature(double tval) { temperature_ = tval; }

  /**
  @brief Get the absolute temperature.

  @note NOT USED.

  @return the absolute temperature [K]
  */
  // double getTemperature() const { return temperature_; }

  /**
  @brief Set the model reference temperature.

  @note NOT USED.

  @param rtval is the reference temperature [K]
  */
  // void setRefT(double rtval) { refT_ = rtval; }

  /**
  @brief Get the model reference temperature.

  @note NOT USED.

  @return the reference temperature [K]
  */
  // double getRefT() const { return refT_; }

  /**
  @brief Get the scaled mass of the phase in the kinetic model.

  The scaled mass of a phase is its mass percent on a total solids basis.

  @note NOT USED.

  @return the vector of scaled masses [percent solids]
  */
  // double getScaledMass() const { return scaledMass_; }

  /**
  @brief Set the <i>initial</i> mass of the phase in the kinetic model.

  The scaled mass of a phase is its mass in grams per unit system volume

  @param initscaledmass is the value to set
  */
  // void setInitScaledMass(const double initscaledmass) {
  //   if (initscaledmass < 0.0) {
  //     initScaledMass_ = 0.0;
  //   } else {
  //     initScaledMass_ = initscaledmass;
  //   }
  //   return;
  // }

  /**
  @brief Get the <i>initial</i> mass of the phase in the kinetic model.

  The scaled mass of a phase is its mass in grams per unit system volume

  @return the initial scaled mass
  */
  // double getInitScaledMass() const { return initScaledMass_; }

  /**
  @brief Set the <i>initial</i> scaled moles of the phase in the kinetic model.

  The scaled moles of a phase is its moles per unit system volume

  @param initscaledmoles is the value to set
  */
  // void setInitScaledMoles(const double initscaledmoles) {
  //   if (initscaledmoles < 0.0) {
  //     initScaledMoles_ = 0.0;
  //   } else {
  //     initScaledMoles_ = initscaledmoles;
  //   }
  //   return;
  // }

  /**
  @brief Get the <i>initial</i> scaled moles of the phase in the kinetic model.

  The scaled moles of a phase is its moles per unit system volume

  @return the initial scaled moles
  */
  // double getInitScaledMoles() const { return initScaledMoles_; }

  /**
  @brief Master method for implementing one kinetic time step.

  In a given time step, a certain number of moles of kinetically controlled
  phases may dissolve or precipitate.
  This function determines the number of moles of each phase to change,
  based on the time interval being simulated.

  This is now a pure virtual function.

  @remark This method is very long and several parts are hard-coded when they
  should be made more general.

  @todo Split this method into more convenient chunks
  @todo Make the methods more general, less hardwiring of parameters
  @todo Make the local variable names more descriptive

  @param timestep is the time interval to simulate [hours]
  @param temperature is the absolute temperature during this step [K]
  @param rh is the internal relative humidity
  @param scaledMass is C-style array of the normalized mass of each
  microstructure phase [g/100 g]
  @param massDissolved is the C-style array of dissolved mass of each
  microstructure phase [g/100g]
  @param cyc is the cycle number (iteration of main loop)
  @param totalDOR is the total degree of reaction [dimensionless]
  */
  virtual void calculateKineticStep(const double timestep, double &scaledMass,
                                    double &massDissolved, int cyc,
                                    double totalDOR, bool doTweak,
                                    bool &doNotModif) = 0;

  /**
  @brief Estimate the initial dissolution rate for this phase.

  This method estimates the dissolution rate at the start of the simulation
  (when DOR is near zero) for use in adaptive timestep initialization.
  The rate is returned as mass fraction dissolved per hour.

  @return the estimated dissolution rate [1/hour], or 0 if not applicable
  */
  virtual double estimateInitialDissolutionRate() const { return 0.0; }

  /**
  @brief Get the initial scaled mass of this phase.

  @return the initial scaled mass [g/100g solids]
  */
  double getInitScaledMass() const { return initScaledMass_; }

  /**
  @brief Set up the number of moles of dependent components in the kinetic
  phases.

  This method loops over the <i>kinetically</i> controlled phases in the kinetic
  model, gets the DC stoichiometry of each phase, and determines the number of
  moles of each DC component based on the number of moles of the kinetically
  controlled phases.
  */
  void setKineticDCMoles();

  /**
  @brief Set the number of moles of dependent components to zero.

  This method loops over the <i>kinetically</i> controlled phases in the kinetic
  model, and sets the number of moles of each DC component of that phase to
  zero.
  */
  void zeroKineticDCMoles();

  /**
  @brief Set the verbose flag

  @param isverbose is true if verbose output should be produced
  */
  // void setVerbose(const bool isverbose) { verbose_ = isverbose; }

  /**
  @brief Get the verbose flag

  @return the verbose flag
  */
  // bool getVerbose() const { return verbose_; }

  /**
  @brief Set the warning flag

  @param iswarning is true if verbose output should be produced
  */
  // void setWarning(const bool iswarning) { warning_ = iswarning; }

  /**
  @brief Get the warning flag

  @return the warning flag
  */
  // bool getWarning() const { return warning_; }

  /**
  @brief Get the name of this kinetic model

  @return the name of this kinetic model
  */
  std::string getModelName(void) { return modelName_; }
}; // End of KineticModel class

#endif // SRC_THAMESLIB_KINETICMODEL_H_
