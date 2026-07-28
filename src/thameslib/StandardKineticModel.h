/**
@file StandardKineticModel.h
@brief Declaration of the StandardKineticModel class.

@section Introduction
General power-law rate model:
    r = k * (1 - Omega^siexp)^dfexp
where Omega = IAP/Ksp and k is dissolutionRateConst. Used for phases whose
kinetics do NOT require the OH-activity / water-activity / LOI / sio2
factors that PozzolanicModel adds — clinker sulfates (Anhydrite, Bassanite,
Gypsum), and (with CNT) Portlandite for what-if exploration.

Han et al. 2025 CEJ found that this Eq. (6) power-law form diverges at
high Omega for portlandite dissolution; their Eq. (7) saturating form is
the follow-on `SaturatingRateModel` work (see
`memory/project_saturating_rate_model_next.md`).

CNT scaffold (nucleation-controlled bootstrapping of zero-mass phases)
lives on this class via `nucleation_` + accumulator + four virtual
overrides; details in `docs/CNT_ARCHITECTURE.md`.
*/

#ifndef SRC_THAMESLIB_STANDARDKINETICMODEL_H_
#define SRC_THAMESLIB_STANDARDKINETICMODEL_H_

#include <optional>

#include "global.h"
#include "Exceptions.h"
#include "ChemicalSystem.h"
#include "KineticController.h"
#include "KineticData.h"
#include "KineticModel.h"
#include "Lattice.h"
#include "NucleationParameters.h"

/**
@class StandardKineticModel
@brief Power-law kinetic model for non-pozzolanic phases (Anhydrite,
Bassanite, Gypsum, and Portlandite when CNT is enabled).
*/

class StandardKineticModel : public KineticModel {

protected:
  double
      surfaceAreaMultiplier_;   /**< Dimensionless factor to multiply the
                                     calculated surface area to account for
                                     unresolved internal porosity, roughness, etc. */
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

  std::optional<NucleationParameters> nucleation_;
      /**< CNT parameters for this phase; empty = CNT disabled */
  double nucleationAccumulator_ = 0.0;
      /**< fractional-voxel accumulator drained when it crosses 1.0 */

public:
  /**
  @brief Default constructor.

  This constructor is not used in THAMES.  It just establishes default values
  for all the member variables.

  @note NOT USED.
  */
  StandardKineticModel();

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
  StandardKineticModel(ChemicalSystem *cs, Lattice *lattice,
                       struct KineticData &kineticData, const bool verbose,
                       const bool warning);

  /**
  @brief Get the type of kinetic model

  @return a string indicating the model type
  */
  std::string getType() const { return (StandardType); }

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

  /**
  @brief Estimate the initial dissolution rate for this phase.

  Uses the dissolution rate constant with conservative assumptions
  about surface area and saturation index to estimate the initial rate.

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

  Calculates the current rate based on the Standard model rate equation:
    rate = k * rh * area * f(SI) * arrhenius
  where f(SI) is the driving force term based on saturation index.

  @param scaledMass the current scaled mass of this phase [g/100g solids]
  @return the current molar rate [mol/hour per 100g solids]
          Positive = dissolution, Negative = precipitation
  */
  double getCurrentMolarRate(double scaledMass) const override;

  /**
  @brief Compute the fractional number of voxels expected to nucleate this cycle.

  Classical Nucleation Theory (CNT) rate for this phase. Guarded: returns 0.0
  when kineticData_.nucleation is empty (CNT not configured for this phase)
  or when S <= 1.0 (undersaturated / at equilibrium). Otherwise pulls
  gamma / theta / A0 from kineticData_.nucleation, v_m from chemSys_,
  T from temperature_, V_voxel and V_electrolyte from lattice_, and
  delegates to cnt::voxelsPerCycle in NucleationRate.h.

  Result is fractional; the caller accumulates and places voxels when the
  accumulator crosses 1.0. Defined in Step 4 of the CNT integration plan
  but NOT YET CALLED from calculateKineticStep — that wiring happens in
  Step 5, along with the batched random placement and the
  AdaptiveTimeController N_cap integration.

  @param dt_hours proposed cycle timestep [hours]
  @return fractional voxel count for this phase this cycle
  */
  double computeNucleationVoxels(double dt_hours) const override;

  /**
  @brief Report whether this phase has a CNT nucleation block configured.
  */
  bool hasNucleation() const override { return nucleation_.has_value(); }

  /**
  @brief Add fractional voxels to this phase's CNT accumulator (Session-50 spec).
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

  Returns the Standard-Eq-6 precipitation rate (dissolutionRateConst_
  times the S-driving-force term) at S > 1, multiplied by the DC
  molar volume to convert mol/m^2/s -> m/s. Returns 0.0 for S <= 1
  (dissolution or equilibrium; no growth on the precipitation branch).
  */
  double getGrowthVelocity(double S) const override;

}; // End of StandardKineticModel class

#endif // SRC_THAMESLIB_STANDARDKINETICMODEL_H_
