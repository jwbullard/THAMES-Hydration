/**
@file SaturatingRateModel.cc
@brief Method definitions for SaturatingRateModel.

Rate calculation implements Bullard 2015 CCR Eq. (2) / Han 2025 CEJ
Eq. (7):
    r_diss = k * (1 - exp[-(-B ln Omega)^n])           for Omega < 1
    r_prec = k * (1 - exp[-( B ln Omega)^n])           for Omega > 1
Dispatched through the sat:: free functions defined in SaturatingRate.h
(exercised directly by test_saturating_rate). This class's role is to
wrap those calls with the surface-area, RH, Arrhenius, and
dissolvedUnits scaling that the KineticController expects, using the
same bookkeeping pattern as StandardKineticModel::calculateKineticStep.

CNT hooks (computeNucleationVoxels, accumulateNucleation,
drainNucleationInteger, hasNucleation) are structurally identical to
StandardKineticModel's.
*/

#include "SaturatingRateModel.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

#include "NucleationRate.h"
#include "SaturatingRate.h"

using std::cout;
using std::endl;

SaturatingRateModel::SaturatingRateModel() {
  // Default constructor is NOT USED; provides safe defaults only.
  surfaceAreaMultiplier_ = 1.0;
  dissolvedUnits_ = 1.0;
  rh_ = 1.0;
  rhFactor_ = 1.0;
  arrhenius_ = 1.0;
  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;
}

SaturatingRateModel::SaturatingRateModel(ChemicalSystem *cs, Lattice *lattice,
                                         struct KineticData &kineticData,
                                         const bool verbose,
                                         const bool warning) {
#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  chemSys_ = cs;
  lattice_ = lattice;

  // A SaturatingRate phase MUST carry a dissolution block. The parser
  // enforces this; guard here in case a caller ever bypasses the parser.
  if (!kineticData.saturatingDissolution.has_value()) {
    throw DataException(
        "SaturatingRateModel", "SaturatingRateModel",
        "SaturatingRate phase constructed without a dissolution block");
  }

  dissolution_ = kineticData.saturatingDissolution;
  precipitation_ = kineticData.saturatingPrecipitation;  // optional; may be empty

  // CNT block if configured for this phase.
  nucleation_ = kineticData.nucleation;

  surfaceAreaMultiplier_ = kineticData.surfaceAreaMultiplier;
  dissolvedUnits_ = kineticData.dissolvedUnits;

  initSolidMass_ = 100.0;

  temperature_ = kineticData.temperature;
  refT_ = kineticData.reftemperature;

  modelName_ = "SaturatingRateModel";
  name_ = kineticData.name;
  microPhaseId_ = kineticData.microPhaseId;
  DCId_ = kineticData.DCId;
  GEMPhaseId_ = kineticData.GEMPhaseId;
  activationEnergy_ = kineticData.activationEnergy;
  scaledMass_ = kineticData.scaledMass;
  initScaledMass_ = kineticData.scaledMass;

  double critporediam = lattice_->getLargestSaturatedPore(); // in nm
  critporediam *= 1.0e-9;                                    // in m
  rh_ = exp(-6.23527e-7 / critporediam / temperature_);
  rh_ = rh_ > 0.55 ? rh_ : 0.551;
  rhFactor_ = rh_;

  arrhenius_ = exp((activationEnergy_ / GASCONSTANT) *
                   ((1.0 / refT_) - (1.0 / temperature_)));

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;
}

void SaturatingRateModel::calculateKineticStep(const double timestep,
                                               double &scaledMass,
                                               double &massDissolved, int cyc,
                                               double /*totalDOR*/,
                                               bool doTweak,
                                               bool &doNotModif) {
  double dissrate = 0.0;

  try {
    scaledMass_ = scaledMass;

    // CNT bypass (guard G3 in docs/CNT_ARCHITECTURE.md, mirrors
    // StandardKineticModel.cc:195). When the phase has a CNT block
    // configured AND currently holds zero mass, refuse to compute a
    // rate here so nothing sneaks in via the area = 1.0 fallback below.
    // CNT placement is then the sole entry point for this phase, and the
    // KineticController's DC-limit lock upstream still holds.
    if (nucleation_.has_value() && scaledMass <= 0.0) {
      massDissolved = 0.0;
      return;
    }

    double area = lattice_->getSurfaceArea(microPhaseId_);
    if (area <= 1.0e-9)
      area = 1.0;
    area *= surfaceAreaMultiplier_;

    double S = chemSys_->getMicroPhaseSI(microPhaseId_);

    if (!doTweak)
      std::clog << "      SaturatingRateModel::calculateKineticStep - "
                   "microPhaseId_/mPhName/SI : "
                << std::setw(3) << std::right << microPhaseId_ << " / "
                << std::setw(15) << std::left << name_ << " / " << S << endl;

    // Sign convention (matches StandardKineticModel):
    //   dissolution  -> positive dissrate  (solid mass decreases)
    //   precipitation-> negative dissrate  (solid mass increases)
    if (S < 1.0) {
      const double r = sat::dissolutionRate(dissolution_->rateConstant,
                                            dissolution_->B, dissolution_->n,
                                            S);
      dissrate = r * rhFactor_ * area;
    } else if (S > 1.0) {
      // Fall back to dissolution parameters if the input JSON did not
      // supply a precipitation block. Log once per phase so the working
      // symmetric-default assumption is auditable.
      const SaturatingRateParameters &pp =
          precipitation_.has_value() ? *precipitation_ : *dissolution_;
      if (!precipitation_.has_value() && !precipitationFallbackLogged_) {
        std::clog << "      SaturatingRateModel::calculateKineticStep - "
                     "precipitation block absent for "
                  << name_
                  << "; using dissolution parameters (microscopic reversibility"
                     " default)."
                  << endl;
        precipitationFallbackLogged_ = true;
      }
      const double r = sat::precipitationRate(pp.rateConstant, pp.B, pp.n, S);
      dissrate = -r * rhFactor_ * area;
    } else {
      dissrate = 0.0;
    }

    dissrate *= arrhenius_;

    // dissrate: mol of phase per 100 g of total solid per hour
    // timestep: hours
    // getDCMolarMass: g/mol
    // -> massDissolved: g of phase per 100 g of total solid
    massDissolved = dissrate * timestep * chemSys_->getDCMolarMass(DCId_);

    if (chemSys_->getKeepDCLowerLimit(DCId_) > 0) {
      if (massDissolved < 0) {
        chemSys_->setKeepDCLowerLimitZero(DCId_);
      } else {
        doNotModif = true;
        return;
      }
    }

    scaledMass = scaledMass_ - massDissolved;

    if (verbose_) {
      std::clog << endl
                << "    SaturatingRateModel::calculateKineticStep "
                   "dissrate/massDissolved/scaledMass_/scaledMass : "
                << dissrate << " / " << massDissolved << " / " << scaledMass_
                << " / " << scaledMass << endl;
    }

    if (scaledMass < 0.0) {
      massDissolved = scaledMass_;
      scaledMass = 0.0;
    }
    scaledMass_ = scaledMass;

    if (verbose_) {
      std::clog << endl
                << "  ****************** SRM_hT = " << timestep
                << "    cyc = " << cyc
                << "    microPhaseId_ = " << microPhaseId_
                << "    microPhase = " << name_
                << "    GEMPhaseIndex = " << GEMPhaseId_
                << " ******************" << endl;
      std::clog << "   SRM_hT   " << "rhFactor_: " << rhFactor_
                << "\tarrhenius_: " << arrhenius_
                << "\tS: " << S
                << "\tarea: " << area << endl;
      std::clog << "   SRM_hT   " << "dissrate: " << dissrate << endl;
      std::clog << "   SRM_hT   " << "initScaledMass_: " << initScaledMass_
                << "\tscaledMass_: " << scaledMass_
                << "\tmassDissolved: " << massDissolved << endl;
      std::clog.flush();
    }

  } catch (EOBException eex) {
    eex.printException();
    exit(1);
  } catch (DataException dex) {
    dex.printException();
    exit(1);
  } catch (FloatException flex) {
    flex.printException();
    exit(1);
  } catch (std::out_of_range &oor) {
    EOBException ex("SaturatingRateModel", "calculateKineticStep", oor.what(),
                    0, 0);
    ex.printException();
    exit(1);
  }

  return;
}

double SaturatingRateModel::estimateInitialDissolutionRate() const {
  // Initial estimate for adaptive-timestep sizing at cycle 0.
  // Assumes S -> 0 (maximum driving force). In the S -> 0 limit,
  // sat::dissolutionRate saturates at k, so the driving-force term is
  // just k. Surface area is estimated as ~1 m^2/g of phase, mirroring
  // StandardKineticModel::estimateInitialDissolutionRate.
  if (initScaledMass_ <= 0.0 || !dissolution_.has_value() ||
      dissolution_->rateConstant <= 0.0) {
    return 0.0;
  }

  const double estimatedArea = initScaledMass_ * surfaceAreaMultiplier_;
  const double dissrate =
      dissolution_->rateConstant * rhFactor_ * estimatedArea * arrhenius_;

  const double molarMass = chemSys_->getDCMolarMass(DCId_);
  return (dissrate * molarMass) / initScaledMass_;
}

double SaturatingRateModel::estimateInitialDissolutionRate(
    double saturationIndex) const {
  // Same estimate as the no-arg version but with the actual S from
  // pre-equilibration. Used by AdaptiveTimeController to size the
  // first cycle when SI is not exactly zero.
  if (initScaledMass_ <= 0.0 || !dissolution_.has_value() ||
      dissolution_->rateConstant <= 0.0) {
    return 0.0;
  }

  double drivingForce = dissolution_->rateConstant; // S <= 0 saturation
  if (saturationIndex > 0.0 && saturationIndex < 1.0) {
    drivingForce = sat::dissolutionRate(dissolution_->rateConstant,
                                        dissolution_->B, dissolution_->n,
                                        saturationIndex);
  } else if (saturationIndex >= 1.0) {
    const SaturatingRateParameters &pp =
        precipitation_.has_value() ? *precipitation_ : *dissolution_;
    drivingForce = sat::precipitationRate(pp.rateConstant, pp.B, pp.n,
                                          saturationIndex);
  }

  const double estimatedArea = initScaledMass_ * surfaceAreaMultiplier_;
  const double dissrate =
      drivingForce * rhFactor_ * estimatedArea * arrhenius_;

  const double molarMass = chemSys_->getDCMolarMass(DCId_);
  return (dissrate * molarMass) / initScaledMass_;
}

double SaturatingRateModel::getCurrentMolarRate(double scaledMass) const {
  // Current signed molar rate (mol/100g/h):
  //   > 0 dissolution, < 0 precipitation.
  // Does NOT touch state; safe from anywhere.
  if (scaledMass <= 0.0 || !dissolution_.has_value() ||
      dissolution_->rateConstant <= 0.0) {
    return 0.0;
  }

  double area = lattice_->getSurfaceArea(microPhaseId_);
  if (area <= 1.0e-9)
    area = 1.0;
  area *= surfaceAreaMultiplier_;

  const double S = chemSys_->getMicroPhaseSI(microPhaseId_);
  double dissrate = 0.0;

  if (S < 1.0) {
    const double r = sat::dissolutionRate(dissolution_->rateConstant,
                                          dissolution_->B, dissolution_->n, S);
    dissrate = r * rhFactor_ * area;
  } else if (S > 1.0) {
    const SaturatingRateParameters &pp =
        precipitation_.has_value() ? *precipitation_ : *dissolution_;
    const double r = sat::precipitationRate(pp.rateConstant, pp.B, pp.n, S);
    dissrate = -r * rhFactor_ * area;
  }

  dissrate *= arrhenius_;
  return dissrate;
}

//
// CNT hooks — symmetric with StandardKineticModel / PozzolanicModel.
//
double SaturatingRateModel::computeNucleationVoxels(double dt_hours) const {
  if (!nucleation_.has_value()) return 0.0;

  double S = chemSys_->getMicroPhaseSI(microPhaseId_);
  if (S <= 1.0) return 0.0;

  double v_m = chemSys_->getDCMolarVolume(DCId_);
  double T_K = temperature_;
  double V_voxel = lattice_->getVolumePerVoxel();
  double V_electrolyte =
      lattice_->getVolumeFraction(ELECTROLYTEID) *
      static_cast<double>(lattice_->getNumSites()) * V_voxel;
  double dt_seconds = dt_hours * 3600.0;

  return cnt::voxelsPerCycle(*nucleation_, v_m, T_K, S,
                             dt_seconds, V_electrolyte, V_voxel);
}

void SaturatingRateModel::accumulateNucleation(double dN) {
  if (dN > 0.0) nucleationAccumulator_ += dN;
}

int SaturatingRateModel::drainNucleationInteger() {
  if (nucleationAccumulator_ < 1.0) return 0;
  // Clamp double->int cast to INT_MAX-1 to avoid undefined behavior when
  // J is extreme (e.g., θ≈1° heterogeneous CNT at high S). Lattice can't
  // hold more voxels than getNumSites() anyway; anything above that is
  // capped again in nucleatePhaseRnd.
  constexpr double kIntMaxSafe =
      static_cast<double>(std::numeric_limits<int>::max() - 1);
  const double drained = std::min(std::floor(nucleationAccumulator_),
                                  kIntMaxSafe);
  int n = static_cast<int>(drained);
  nucleationAccumulator_ -= static_cast<double>(n);
  return n;
}

double SaturatingRateModel::getNucleationRate(double S) const {
  if (!nucleation_.has_value()) return 0.0;
  const double v_m = chemSys_->getDCMolarVolume(DCId_);
  return cnt::nucleationRate(*nucleation_, v_m, temperature_, S);
}

double SaturatingRateModel::getGrowthVelocity(double S) const {
  // Only precipitation-side (S > 1) is a growth event; dissolution
  // consumes mass rather than growing the phase, so returns 0.0
  // for S <= 1.
  if (S <= 1.0) return 0.0;
  if (!dissolution_.has_value()) return 0.0;

  // Fall back to dissolution params for precipitation if precipitation
  // block is absent — same microscopic-reversibility default as
  // calculateKineticStep.
  const SaturatingRateParameters &pp =
      precipitation_.has_value() ? *precipitation_ : *dissolution_;

  // Surface-normal molar rate: mol / m^2 / s. sat::precipitationRate
  // returns the magnitude (positive) at Omega > 1.
  const double r = sat::precipitationRate(pp.rateConstant, pp.B, pp.n, S);

  // Convert to linear velocity: m/s = (mol/m^2/s) * (m^3/mol).
  // Uses the DC's physical molar volume from GEMS (m^3/mol).
  const double v_m = chemSys_->getDCMolarVolume(DCId_);
  return r * v_m;
}
