/**
@file PozzolanicModel.h
@brief Declaration of the PozzolanicModel class.

@section Introduction
Extends the power-law form of StandardKineticModel with the empirical
factors Dove and Crerar identified for silica dissolution — OH-activity
dependence (ohexp), water-activity squared, LOI correction, sio2 mass
fraction, and Ca/Na/K Langmuir adsorption on the base rate constant.
Also adds an early-vs-late diffusion-rate switch for products whose
kinetics change regime as the surrounding gel densifies.

Used for silica fume, fly-ash, metakaolin, and other reactive SCMs.

@section References

    -# PM Dove, N Han, AF Wallace, JJ De Yoreo, Kinetics of amorphous silica
dissolution and the paradox of the silica polymorphs, Proceedings of the
National Academy of Sciences USA, 105 (2008) 9903–9908.
    -# PM Dove, DA Crerar, Kinetics of quartz dissolution in electrolyte
solutions using a hydrothermal mixed flow reactor, Geochimica et
Cosmochimica Acta, 54 (1990) 955–969.
*/

#ifndef SRC_THAMESLIB_POZZOLANICMODEL_H_
#define SRC_THAMESLIB_POZZOLANICMODEL_H_

#include <optional>

#include "ChemicalSystem.h"
#include "Exceptions.h"
#include "KineticController.h"
#include "KineticData.h"
#include "KineticModel.h"
#include "Lattice.h"
#include "NucleationParameters.h"
#include "TransportCorrection.h"
#include "TransportParameters.h"
#include "global.h"

/**
@class PozzolanicModel
@brief Kinetic model for reactive SCMs (silica fume, fly ash, metakaolin).
Extends the power-law rate with Dove/Crerar factors for silica/glass
dissolution kinetics. CNT scaffold present for parity with
StandardKineticModel; rarely fires in practice because pozzolanic phases
dissolve rather than nucleate.
*/

class PozzolanicModel : public KineticModel {

protected:
  double surfaceAreaMultiplier_;   /**< Dimensionless factor to multiply the
                                        calculated surface area to account for
                                        unresolved internal porosity, roughness,
                                        etc. */
  double dissolutionRateConst_;    /**< Rate constand for dissolution
                                        (mol/m2/h) */
  double diffusionRateConstEarly_; /**< Rate constant for early-age diffusion
                                        (mol/m2/h) */
  double diffusionRateConstLate_;  /**< Rate constant for later-age diffusion
                                        (mol/m2/h) */
  /**
  @brief Number of dissolved DC units per unit dissolution reaction
  */
  double dissolvedUnits_;
  double siexp_;  /**< Exponent on saturation index (unitless) */
  double dfexp_;  /**< Exponent on driving force (unitless) */
  double dorexp_; /**< Exponent on degree of reaction (unitless) */
  double ohexp_;  /**< Exponent on OH ion activity (unitless) */
  double sio2_;   /**< Mass fraction of SiO2 (unitless) */
  double al2o3_;  /**< Mass fraction of Al2O3 (unitless) */
  double cao_;    /**< Mass fraction of CaO (unitless) */

  double rh_;        /**< relative humidity */
  double rhFactor_;  /**< relative humidity factor, i.e. the correction of
                          the hydration rate taking into account the ambient
                          relative humidity */
  double arrhenius_; /**< arrhenius factor */

  std::optional<NucleationParameters> nucleation_;
      /**< CNT parameters for this phase; empty = CNT disabled */
  double nucleationAccumulator_ = 0.0;
      /**< fractional-voxel accumulator drained when it crosses 1.0 */

  std::optional<TransportParameters> transport_;
      /**< Shell-diffusion parameters (Phase 3 of transport-kinetics
           plan). Empty = no shell correction. */
  int limitingDCId_ = -1;
      /**< DC id of the limiting species; see StandardKineticModel.h
           for semantics. */

public:
  /**
  @brief Default constructor.

  This constructor is not used in THAMES.  It just establishes default values
  for all the member variables.

  @note NOT USED.
  */
  PozzolanicModel();

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
  PozzolanicModel(ChemicalSystem *cs, Lattice *lattice,
                  struct KineticData &kineticData, const bool verbose,
                  const bool warning);

  /**
  @brief Get the type of kinetic model

  @return a string indicating the model type
  */
  std::string getType() const { return (PozzolanicType); }

  /**
  @brief Set the surface area multiplier

  This is a dimensionless multiplication factor to the surface
  area to account for unresolved internal porosity, roughness, etc.

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
  @brief Set the early-age diffusion rate constant

  @note NOT USED.

  @param rc is the early-age diffusion rate constant value to use
  */
  void setDiffusionRateConstEarly(const double rc) {
    diffusionRateConstEarly_ = max(rc, 0.0);
  }

  /**
  @brief Get the dissolution rate constant

  @note NOT USED.

  @return the early-age diffusion rate constant
  */
  double getDiffusionRateConstEarly() const { return diffusionRateConstEarly_; }

  /**
  @brief Set the later-age diffusion rate constant

  @note NOT USED.

  @param rc is the later-age diffusion rate constant value to use
  */
  void setDiffusionRateConstLate(const double rc) {
    diffusionRateConstLate_ = max(rc, 0.0);
  }

  /**
  @brief Get the later-age diffusion rate constant

  @note NOT USED.

  @return the later-age diffusion rate constant
  */
  double getDiffusionRateConstLate() const { return diffusionRateConstLate_; }

  /**
  @brief Set the exponent on the degree of reaction in the diffusion rate
  equation

  @note NOT USED.

  @param dorexp is the exponent value to use
  */
  void setDorexp(const double dorexp) { dorexp_ = max(dorexp, 0.0); }

  /**
  @brief Get the exponent on the degree of reaction in the diffusion rate
  equation

  @note NOT USED.

  @return the exponent on the degree of reaction
  */
  double getDorexp() const { return dorexp_; }

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
  @brief Set the exponent on the hydroxyl ion activity

  @note NOT USED.

  @param ohexp is the exponent value to use
  */
  void setOhexp(const double ohexp) { ohexp_ = max(ohexp, 0.0); }

  /**
  @brief Get the exponent on the hydroxyl ion activity

  @note NOT USED.

  @return the exponent on the hydroxyl ion activity
  */
  double getOhexp() const { return ohexp_; }

  /**
  @brief Set the mass fraction of SiO2

  @note NOT USED.

  @param sio2 is the mass fraction of SiO2
  */
  void setSio2(const double sio2) {
    sio2_ = max(sio2, 0.0);
    if (sio2_ > 1.0)
      sio2_ = 1.0;
  }

  /**
  @brief Get the SiO2 mass fraction

  @note NOT USED.

  @return the SiO2 mass fraction
  */
  double getSio2() const { return sio2_; }

  /**
  @brief Set the mass fraction of Al2O3

  @note NOT USED.

  @param al2o3 is the mass fraction of Al2O3
  */
  void setAl2o3(const double al2o3) {
    al2o3_ = max(al2o3, 0.0);
    if (al2o3_ > 1.0)
      al2o3_ = 1.0;
  }

  /**
  @brief Get the Al2O3 mass fraction

  @note NOT USED.

  @return the Al2O3 mass fraction
  */
  double getAl2o3() const { return al2o3_; }

  /**
  @brief Set the mass fraction of CaO

  @note NOT USED.

  @param cao is the mass fraction of CaO
  */
  void setCao(const double cao) {
    cao_ = max(cao, 0.0);
    if (cao_ > 1.0)
      cao_ = 1.0;
  }

  /**
  @brief Get the CaO mass fraction

  @note NOT USED.

  @return the CaO mass fraction
  */
  double getCao() const { return cao_; }

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

  /**
  @brief Estimate the initial dissolution rate for this pozzolanic phase.

  Uses the dissolution rate constant with conservative assumptions
  about surface area, OH- activity, and saturation index.

  @return the estimated dissolution rate [1/hour]
  */
  double estimateInitialDissolutionRate() const override;

  /**
  @brief Estimate the initial dissolution rate using actual saturation index.

  Uses the actual SI from pre-equilibration to compute the driving force
  term (1 - SI^siexp)^dfexp rather than assuming SI ≈ 0.

  @param saturationIndex the saturation index from GEMS equilibration
  @return the estimated dissolution rate [1/hour]
  */
  double estimateInitialDissolutionRate(double saturationIndex) const override;

  /**
  @brief Get the current molar rate of dissolution/precipitation.

  Calculates the current rate based on the Pozzolanic model rate equation,
  which includes OH- activity dependence in addition to the Standard model terms.

  @param scaledMass the current scaled mass of this phase [g/100g solids]
  @return the current molar rate [mol/hour per 100g solids]
          Positive = dissolution, Negative = precipitation
  */
  double getCurrentMolarRate(double scaledMass) const override;

  /**
  @brief Compute the fractional number of voxels expected to nucleate this cycle.

  CNT rate for this Pozzolanic phase. Guarded on:
  - `nucleation_.has_value()` — CNT is configured for this phase
  - `S > 1.0` — thermodynamic driver exists
  - `lattice_->getGrowthInterfaceSize()[microPhaseId_] == 0` — the Session-50
    same-phase-interface gating rule that prevents CNT from double-counting
    when Pozzolanic's empirical acceleration parameters already handle growth

  Result is fractional; caller accumulates and places voxels when the
  accumulator crosses 1.0.
  */
  double computeNucleationVoxels(double dt_hours) const override;

  /**
  @brief Report whether this phase has a CNT nucleation block configured.
  */
  bool hasNucleation() const override { return nucleation_.has_value(); }

  /**
  @brief Add fractional voxels to this phase's CNT accumulator.
  */
  void accumulateNucleation(double dN) override;

  /**
  @brief Return floor(accumulator) and subtract that many whole voxels.
  */
  int drainNucleationInteger() override;

  /**
  @brief CNT nucleation rate J at supersaturation S.
  Delegates to `cnt::nucleationRate` when a nucleation block is configured.
  */
  double getNucleationRate(double S) const override;

  /**
  @brief Linear growth velocity at supersaturation S.

  Returns the Pozzolanic-model precipitation rate at S > 1 multiplied
  by the DC molar volume to convert mol/m^2/s -> m/s. Returns 0.0 for
  S <= 1 (dissolution or equilibrium; no growth on the precipitation
  branch). Uses the same rate expression as Standard-Eq-6 but with the
  Pozzolanic-specific parameters — since Pozzolanic phases in practice
  rarely nucleate (they dissolve as sources of reactive silica), this
  override is present mostly for parity with Standard/SaturatingRate.
  */
  double getGrowthVelocity(double S) const override;

}; // End of PozzolanicModel class

#endif // SRC_THAMESLIB_POZZOLANICMODEL_H_
