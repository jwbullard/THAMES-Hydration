/**
@file SaturatingRateModel.h
@brief Declaration of SaturatingRateModel — the 4th concrete KineticModel.

Implements Bullard 2015 CCR Eq. (2) / Han 2025 CEJ Eq. (7), the
saturating rate law:

  r = k * (1 - exp[-(-B ln Omega)^n])           (dissolution branch)
  r = k * (1 - exp[-(B ln Omega)^n])            (precipitation branch, mirrored)

The dissolution and precipitation branches carry independent
parameters (rateConstant, B, n) — near equilibrium they mirror each
other by microscopic reversibility; far from equilibrium the rate laws
can differ. The precipitation block is optional; if absent, the model
uses the dissolution parameters as the working symmetric default and
logs one line noting the fallback.

Rate saturates at k in the far-from-equilibrium limit rather than
diverging as `StandardKineticModel`'s Eq. 6 power law does — this is
the mathematical property that motivated adding a new class. Used for
Portlandite and (eventually) Alite, both of which have published
Han-style dissolution calibrations.

CNT support (`nucleation_`, accumulator, four virtuals) is symmetric
with StandardKineticModel / PozzolanicModel. See `docs/CNT_ARCHITECTURE.md`
for the CNT integration architecture.

@section References
    -# J.W. Bullard, G.W. Scherer, J.J. Thomas, Time dependent driving
       forces and the kinetics of tricalcium silicate hydration,
       Cement and Concrete Research, 74 (2015) 26-34.
    -# Y. Han, M.G. Ucak-Astarlioglu, J.F. Burroughs, J.W. Bullard,
       Calcium hydroxide dissolution rates: Dependence on temperature
       and saturation, Chemical Engineering Journal, 515 (2025) 162484.
*/

#ifndef SRC_THAMESLIB_SATURATINGRATEMODEL_H_
#define SRC_THAMESLIB_SATURATINGRATEMODEL_H_

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
#include "SaturatingRateParameters.h"
#include "global.h"

/**
@class SaturatingRateModel
@brief Saturating dissolution/precipitation rate law (Bullard 2015 / Han
2025 Eq. 7). Alternative to StandardKineticModel with a bounded
far-from-equilibrium rate. Currently a scaffold; rate calculation
lands in Step S3 of the SaturatingRateModel plan.
*/
class SaturatingRateModel : public KineticModel {

protected:
  double surfaceAreaMultiplier_;
      /**< Dimensionless factor multiplying computed surface area. */
  double dissolvedUnits_;
      /**< Number of DC units produced per unit dissolution reaction. */
  double rh_;
      /**< Relative humidity. */
  double rhFactor_;
      /**< RH correction factor applied to the rate. */
  double arrhenius_;
      /**< Arrhenius scale factor precomputed at construction from
           activationEnergy, refT, and temperature (see StandardKineticModel
           for the analogous computation). */

  std::optional<SaturatingRateParameters> dissolution_;
      /**< Bullard/Han Eq. 7 params for the dissolution direction.
           Required at construction; if absent the constructor throws. */
  std::optional<SaturatingRateParameters> precipitation_;
      /**< Bullard/Han Eq. 7 params for the precipitation direction.
           Optional; empty means fall back to dissolution values via
           microscopic reversibility. */

  std::optional<NucleationParameters> nucleation_;
      /**< CNT parameters for this phase; empty = CNT disabled */
  double nucleationAccumulator_ = 0.0;
      /**< Fractional-voxel accumulator drained when it crosses 1.0 */

  std::optional<TransportParameters> transport_;
      /**< Shell-diffusion parameters (Phase 3 of transport-kinetics
           plan). Empty = no shell correction; rate law uses the
           uncorrected saturating form. */
  int limitingDCId_ = -1;
      /**< DC id of the limiting species. See StandardKineticModel.h
           for the resolution/fallback semantics. */

  bool precipitationFallbackLogged_ = false;
      /**< First-time-only latch: logs one line when the model falls
           back to dissolution_ parameters for the precipitation branch
           because the input JSON omitted the precipitation block. */

public:
  /**
  @brief Default constructor. NOT USED.
  */
  SaturatingRateModel();

  /**
  @brief Constructor invoked by KineticController::makeModel.
  Requires kineticData.saturatingDissolution to be populated; throws
  DataException if not. Precipitation block is optional.
  */
  SaturatingRateModel(ChemicalSystem *cs, Lattice *lattice,
                      struct KineticData &kineticData, const bool verbose,
                      const bool warning);

  /**
  @brief Get the type of kinetic model.
  */
  std::string getType() const { return (SaturatingRateType); }

  /**
  @brief Master method for implementing one kinetic time step.
  Currently a placeholder — returns immediately. Filled in during
  Step S3 (see docs/SATURATING_RATE.md).
  */
  virtual void calculateKineticStep(const double timestep, double &scaledMass,
                                    double &massDissolved, int cyc,
                                    double totalDOR, bool doTweak,
                                    bool &doNotModif);

  /**
  @brief Estimate the initial dissolution rate for this phase.
  Currently returns 0 — placeholder until S3.
  */
  double estimateInitialDissolutionRate() const override;

  /**
  @brief Estimate the initial dissolution rate using an actual SI value.
  Currently returns 0 — placeholder until S3.
  */
  double estimateInitialDissolutionRate(double saturationIndex) const override;

  /**
  @brief Get the current molar rate of dissolution/precipitation.
  Currently returns 0 — placeholder until S3.
  */
  double getCurrentMolarRate(double scaledMass) const override;

  /**
  @brief Compute the fractional number of voxels expected to nucleate.
  Symmetric with StandardKineticModel's implementation; guards on
  nucleation_ presence and S > 1; delegates to cnt::voxelsPerCycle.
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

  Returns `cnt::nucleationRate(*nucleation_, v_m, T_K, S)` when a
  nucleation block is configured, else 0. Consumed by
  KineticController's JMAK-per-voxel machinery.
  */
  double getNucleationRate(double S) const override;

  /**
  @brief Linear growth velocity at supersaturation S.

  Returns `sat::precipitationRate(...) * v_molar_DC` for S > 1
  (precipitation regime). Returns 0.0 for S <= 1 (dissolution or
  equilibrium; no growth on the precipitation branch). Uses the
  precipitation block's rateConstant/B/n if present, else falls back
  to the dissolution block by microscopic reversibility (the same
  fallback as calculateKineticStep).
  */
  double getGrowthVelocity(double S) const override;

}; // End of SaturatingRateModel class

#endif // SRC_THAMESLIB_SATURATINGRATEMODEL_H_
