/**
@file C3SKineticModel.h
@brief Declaration of the StandardKineticModel class.

@section Introduction
This class implements a dissolution equation like Bullard et al. for
C3S ([1] eq. 2).

@section References

    1. Jeffrey W. Bullard, George W. Scherer, Jeffrey J. Thomas,
Time dependent driving forces and the kinetics of tricalcium silicate hydration,
Cement and Concrete Research 74(2015)26–34.

*/

#ifndef SRC_THAMESLIB_C3SKINETICMODEL_H_
#define SRC_THAMESLIB_C3SKINETICMODEL_H_

#include "global.h"
#include "Exceptions.h"
#include "ChemicalSystem.h"
#include "KineticController.h"
#include "KineticData.h"
#include "KineticModel.h"
#include "Lattice.h"

/**
@class C3SKineticModel
@brief Handles the kinetic model of C3S

*/

class C3SKineticModel : public KineticModel {

protected:
  double
      surfaceAreaMultiplier_;   /**< Dimensionless factor to multiply the
                                     calculated surface area to account for
                                     unresolved roughness, etc. */
  double dissolutionRateConst_; /**< Rate constant for dissolution (mol/m2/h) */
  /**
  @brief Number of dissolved DC units per unit dissolution reaction
  */
  double dissolvedUnits_;
  double siexp_; /**< Exponent on saturation index (unitless) */
  double dfexp_; /**< Exponent on driving force (unitless) */

  double rh_;        /**< relative humidity */
  double rhFactor_;  /**< relative humidity factor, i.e. the correction of the
                          hydration rate taking into account the ambient relative
                          humidity */
  double arrhenius_; /**< arrhenius factor */

  double kC3S_;
  double C1_;
  double r_;

public:
  /**
  @brief Default constructor.

  This constructor is not used in THAMES.  It just establishes default values
  for all the member variables.

  @note NOT USED.
  */
  C3SKineticModel();

  /**
  @brief Overloaded constructor.

  This constructor is the one invoked by THAMES.  It can only be called once the
  various other objects for the simulation are allocated and constructed.

  @param cs is a pointer to the ChemicalSystem object for the simulation
  @param lattice is a pointer to the Lattice object holding the microstructure
  @param kineticData is the collection of kinetic parameters already stored
  @param verbose is true if verbose output should be produced
  @param warning is false if suppressing warning output
  */
  C3SKineticModel(ChemicalSystem *cs, Lattice *lattice,
                       struct KineticData &kineticData, const bool verbose,
                       const bool warning);

  /**
  @brief Get the type of kinetic model

  @return a string indicating the model type
  */
  std::string getType() const { return (C3SType); }

  /**
  @brief Set the surface area multiplier

  This is a dimensionless multiplication factor to the surface
  area to account for unresolved internal roughness, etc.

  @param sam is the rate constant value to use
  */
  void setSurfaceAreaMultiplier(const double sam) {
    surfaceAreaMultiplier_ = max(sam, 0.0);
  }

  /**
  @brief Set the dissolution rate constant

  @param rc is the rate constant value to use
  */
  void setDissolutionRateConst(const double rc) {
    dissolutionRateConst_ = max(rc, 0.0);
  }

  /**
  @brief Get the dissolution rate constant

  @note NOT USED.

  @return the dissolution rate constant
  */
  double getDissolutionRateConst() const { return dissolutionRateConst_; }

  /**
  @brief Set the number of dissolved DC units per unit dissolution

  @note NOT USED.

  @param dissolvedUnits is the value to set
  */
  void setDissolvedUnits(const double dissolvedUnits) {
    dissolvedUnits_ = max(dissolvedUnits, 1.0);
  }

  /**
  @brief Get the number of dissolved DC units per unit dissolution

  @note NOT USED.

  @return the number of dissolved DC units
  */
  double getDissolvedUnits() const { return dissolvedUnits_; }

  /**
  @brief Set the exponent on the saturation index

  @note NOT USED.

  @param siexp is the exponent value to use
  */
  void setSiexp(const double siexp) { siexp_ = max(siexp, 0.0); }

  /**
  @brief Get the exponent on the saturation index

  @note NOT USED.

  @return the exponent on the saturation index
  */
  double getSiexp() const { return siexp_; }

  /**
  @brief Set the exponent on the driving force

  @note NOT USED.

  @param dfexp is the exponent value to use
  */
  void setDfexp(const double dfexp) { dfexp_ = max(dfexp, 0.0); }

  /**
  @brief Get the exponent on the driving force

  @note NOT USED.

  @return the exponent on the driving force
  */
  double getDfexp() const { return dfexp_; }

  /**
  @brief Master method for implementing one kinetic time step.

  Overloaded from base class to handle pozzolanic materials

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
                                    bool &doNotModif);

}; // End of C3SKineticModel class

#endif // SRC_THAMESLIB_C3SKINETICMODEL_H_
