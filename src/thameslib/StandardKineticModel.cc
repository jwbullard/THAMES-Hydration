/**
@file  StandardKineticModel.cc
@brief Method definitions for the StandardKineticModel class.

*/
#include "StandardKineticModel.h"

#include <cmath>

#include "NucleationRate.h"

using std::cout;
using std::endl;

StandardKineticModel::StandardKineticModel() {

  ///
  /// Default temperature in the PK model is 20 C (or 293 K)
  ///

  temperature_ = 293.15; // default temperature (K)
  refT_ = 293.15;        // default temperature (K)

  ///
  /// Default values for the rate constants
  ///

  dissolutionRateConst_ = 0.0;

  ///
  /// Default value for the exponents in the rate equation
  ///

  siexp_ = 1.0;
  dfexp_ = 1.0;
  lossOnIgnition_ = 0.0;

  name_ = "";
  microPhaseId_ = 2;
  DCId_ = 2;
  GEMPhaseId_ = 2;
  activationEnergy_ = 0.0;
  scaledMass_ = 0.0;
  initScaledMass_ = 0.0;

  temperature_ = lattice_->getTemperature();
  double critporediam = lattice_->getLargestSaturatedPore(); // in nm
  critporediam *= 1.0e-9;                                    // in m
  rh_ = exp(-6.23527e-7 / critporediam / temperature_);
  rh_ = rh_ > 0.55 ? rh_ : 0.551;
  rhFactor_ = rh_;

  arrhenius_ = exp((activationEnergy_ / GASCONSTANT) *
                   ((1.0 / refT_) - (1.0 / temperature_)));

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion hours or 114,000 years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;

  return;
}

StandardKineticModel::StandardKineticModel(ChemicalSystem *cs, Lattice *lattice,
                                           struct KineticData &kineticData,
                                           const bool verbose,
                                           const bool warning) {

  // Set the verbose and warning flags

  verbose_ = verbose;
  warning_ = warning;
#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  chemSys_ = cs;
  lattice_ = lattice;

  setDissolutionRateConst(kineticData.dissolutionRateConst);
  setSurfaceAreaMultiplier(kineticData.surfaceAreaMultiplier);
  setDissolvedUnits(kineticData.dissolvedUnits);
  setSiexp(kineticData.siexp);
  setDfexp(kineticData.dfexp);
  lossOnIgnition_ = kineticData.loi;

  ///
  /// Default initial solid mass is 100 g
  ///

  initSolidMass_ = 100.0;

  temperature_ = kineticData.temperature;
  refT_ = kineticData.reftemperature;

  modelName_ = "StandardKineticModel";
  name_ = kineticData.name;
  microPhaseId_ = kineticData.microPhaseId;
  DCId_ = kineticData.DCId;
  GEMPhaseId_ = kineticData.GEMPhaseId;
  activationEnergy_ = kineticData.activationEnergy;
  scaledMass_ = kineticData.scaledMass;
  initScaledMass_ = kineticData.scaledMass;

  // Copy in CNT parameters if present in the input JSON (empty otherwise).
  // std::optional propagates the empty state naturally when CNT is not
  // configured for this phase.
  nucleation_ = kineticData.nucleation;

  double critporediam = lattice_->getLargestSaturatedPore(); // in nm
  critporediam *= 1.0e-9;                                    // in m
  rh_ = exp(-6.23527e-7 / critporediam / temperature_);
  rh_ = rh_ > 0.55 ? rh_ : 0.551;
  rhFactor_ = rh_;

  arrhenius_ = exp((activationEnergy_ / GASCONSTANT) *
                   ((1.0 / refT_) - (1.0 / temperature_)));

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion hours or 114,000 years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;

  return;
}

void StandardKineticModel::calculateKineticStep(const double timestep,
                                                double &scaledMass,
                                                double &massDissolved, int cyc,
                                                double totalDOR, bool doTweak,
                                                bool &doNotModif) {
  ///
  /// Initialize local variables
  ///

  double dissrate = 1.0e9; // Nucleation and growth rate

  ///
  /// Determine if this is a normal step or a necessary
  /// tweak from a failed GEM_run call
  ///

  try {

    // First step each iteration is to equilibrate gas phase
    // with the electrolyte, while forbidding anything new
    // from precipitating.

    // RH factor is the same for all clinker phases
    // double vfvoid = lattice_->getVolumeFraction(VOIDID);
    // double vfh2o = lattice_->getVolumeFraction(ELECTROLYTEID);

    /// This is a big kluge for internal relative humidity
    /// @note Using new gel and interhydrate pore size distribution model
    ///       which is currently contained in the Lattice object.
    ///
    /// Surface tension of water is gamma = 0.072 J/m2
    /// Molar volume of water is Vm = 1.8e-5 m3/mole
    /// The Kelvin equation is
    ///    p/p0 = exp (-4 gamma Vm / d R T) = exp (-6.23527e-7 / (d T))
    ///
    ///    where d is the pore diameter in meters and T is absolute
    ///    temperature

    /// Assume a zero contact angle for now.
    /// @todo revisit the contact angle issue

    scaledMass_ = scaledMass;

    // ---------- CNT bypass (guard G3 in docs/CNT_ARCHITECTURE.md) ----------
    // If this phase has a nucleation block AND currently has zero mass,
    // skip the entire rate calculation. Two consequences:
    //   * The classical "area = 1.0" fallback at line ~199 below (which
    //     lets Standard produce a nonzero rate from zero surface area)
    //     cannot spontaneously precipitate the phase.
    //   * DCMoles_ is not touched here, so the pre-loop CNT-lock in
    //     KineticController::calculateKineticStep (which set DC upper and
    //     lower limits to 0 for this phase) still holds. CNT is then the
    //     ONLY code path that can bring the phase into existence.
    // Once CNT places nuclei this cycle, scaledMass becomes > 0 in the
    // next cycle and the bypass no longer fires — normal Standard growth
    // kinetics apply from then on.
    if (nucleation_.has_value() && scaledMass <= 0.0) {
      massDissolved = 0.0;
      return;
    }

    // JWB BEWARE: The new definition of area is truly a geometric calculation
    // made on the microstructure. It does not catch BET or internal surface
    // area if that ends up being important. The units of area are m2 per 100
    // g of intial total solid

    double area = lattice_->getSurfaceArea(microPhaseId_);
    if (area <= 1.0e-9)
      area = 1.0;

    // surfaceAreaMultiplier_ is a way to account for the influence
    // of subvoxel porosity that increases the total surface area in a voxel.

    area *= surfaceAreaMultiplier_;

    // Saturation index , but be sure that there is only one GEM Phase
    /// @note Assumes there is only one phase in this microstructure
    /// component
    /// @todo Generalize to multiple phases in a component (how?)

    double saturationIndex = chemSys_->getMicroPhaseSI(microPhaseId_);

    if (!doTweak)
      std::clog << "      StandardKineticModel::calculateKineticStep - "
                   "microPhaseId_/mPhName/SI : "
                << std::setw(3) << std::right << microPhaseId_ << " / "
                << std::setw(15) << std::left << name_ << " / "
                << chemSys_->getMicroPhaseSI(microPhaseId_) << endl;

    // dissolutionRateConst_ has units of mol/m2/h
    // area has units of m2 of phase per 100 g of total solid
    // Therefore dissrate has units of mol of phase per 100 g of all solid
    // per h

    // This equation basically implements the Dove and Crerar rate
    // equation for quartz.  Needs to be calibrated for silica fume, but
    // hopefully the BET area and LOI will help do that.

    if (saturationIndex < 1.0) {
      dissrate = dissolutionRateConst_ * rhFactor_ * area *
                 pow((1.0 - pow(saturationIndex, siexp_)), dfexp_);
    } else {
      dissrate = -dissolutionRateConst_ * rhFactor_ * area *
                 pow((pow(saturationIndex, siexp_) - 1.0), dfexp_);
    }

    dissrate *= arrhenius_;

    // dissrate has units of mol of phase per 100 g of solid per hour
    // timestep has units of hours
    // molar mass has units of grams per mole
    // Therefore, massDissolved has units of grams of phase per 100 g of all
    // solid
    massDissolved = dissrate * timestep * chemSys_->getDCMolarMass(DCId_);

    if (chemSys_->getKeepDCLowerLimit(DCId_) > 0) {
      if (massDissolved < 0) {
        // std::clog << "      StandardKineticModel::calculateKineticStep - 1 -
        // cyc = "
        //      << cyc << "   DCId_ = " << DCId_ << "   keepDCLowerLimit = "
        //      << chemSys_->getKeepDCLowerLimit(DCId_)
        //      << "   massDissolved(negative)/scaledMass_ : "
        //      << massDissolved << " / " << scaledMass_
        //      << "   doNotModif = " << doNotModif
        //      << endl;
        chemSys_->setKeepDCLowerLimitZero(DCId_);
      } else {
        // std::clog << "      StandardKineticModel::calculateKineticStep - 2 -
        // cyc = "
        //      << cyc << "   DCId_ = " << DCId_ << "   keepDCLowerLimit = "
        //      << chemSys_->getKeepDCLowerLimit(DCId_)
        //      << "   massDissolved(positive)/scaledMass_ : "
        //      << massDissolved << " / " << scaledMass_
        //      << "   doNotModif = " << doNotModif
        //      << endl;
        doNotModif = true;
        return;
      }
      // } else {
      //   std::clog << "      StandardKineticModel::calculateKineticStep - 3 -
      //   cyc = "
      //        << cyc << "   DCId_ = " << DCId_ << "   keepDCLowerLimit = "
      //        << chemSys_->getKeepDCLowerLimit(DCId_)
      //        << "   massDissolved(n-p)/scaledMass_ : "
      //        << massDissolved << " / " << scaledMass_
      //        << "   doNotModif = " << doNotModif
      //        << endl;
    }

    scaledMass = scaledMass_ - massDissolved;

    if (verbose_) {
      std::clog << endl
                << "    StandardKineticModel::calculateKineticStep "
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
                << "  ****************** SKM_hT = " << timestep
                << "    cyc = " << cyc
                << "    microPhaseId_ = " << microPhaseId_
                << "    microPhase = " << name_
                << "    GEMPhaseIndex = " << GEMPhaseId_
                << " ******************" << endl;
      std::clog << "   SKM_hT   " << "rhFactor_: " << rhFactor_
                << "\tarrhenius_: " << arrhenius_
                << "\tsaturationIndex: " << saturationIndex
                << "\tarea: " << area << endl;
      std::clog << "   SKM_hT   " << "dissrate: " << dissrate << endl;
      std::clog << "   SKM_hT   " << "initScaledMass_: " << initScaledMass_
                << "\tscaledMass_: " << scaledMass_
                << "\tmassDissolved: " << massDissolved << endl;
      std::clog << "   cyc = " << cyc << "    microPhaseId_ = " << microPhaseId_
                << "    microPhaseName = " << name_
                << "    saturationIndex = " << saturationIndex
                << "   Dc_a = " << chemSys_->getNode()->DC_a(DCId_) << endl;
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
  } catch (out_of_range &oor) {
    EOBException ex("StandardKineticModel", "calculateKineticStep", oor.what(),
                    0, 0);
    ex.printException();
    exit(1);
  }

  return;
}

double StandardKineticModel::estimateInitialDissolutionRate() const {
  //
  // Estimate the initial dissolution rate for adaptive timestep sizing.
  //
  // The Standard model rate equation is:
  //   dissrate [mol/100g/h] = k * rh * area * f(SI) * arrhenius
  //   massDissolved [g/100g] = dissrate * dt * molarMass
  //
  // For initial estimation, we assume:
  //   - SI ≈ 0 (far from equilibrium), so f(SI) ≈ 1
  //   - Initial surface area estimated from scaled mass and typical specific area
  //

  if (initScaledMass_ <= 0.0 || dissolutionRateConst_ <= 0.0) {
    return 0.0;
  }

  // Estimate initial surface area: assume ~1 m2 per gram of phase
  // (this is a rough estimate; actual value depends on particle size)
  double estimatedArea = initScaledMass_ * surfaceAreaMultiplier_;

  // dissrate in mol/100g/h (assuming SI=0, so driving force = 1)
  double dissrate = dissolutionRateConst_ * rhFactor_ * estimatedArea * arrhenius_;

  // Convert to mass fraction rate [1/h]
  // dissrate [mol/100g/h] * molarMass [g/mol] / initScaledMass [g/100g] = [1/h]
  double molarMass = chemSys_->getDCMolarMass(DCId_);
  double rate = (dissrate * molarMass) / initScaledMass_;

  return rate;
}

double StandardKineticModel::estimateInitialDissolutionRate(
    double saturationIndex) const {
  //
  // Estimate the initial dissolution rate using actual saturation index.
  //
  // This method uses the SI from pre-equilibration to compute the driving
  // force term rather than assuming SI ≈ 0 (maximum driving force).
  //
  // The driving force is: f(SI) = (1 - SI^siexp)^dfexp for dissolution (SI < 1)
  //                       f(SI) = (SI^siexp - 1)^dfexp for precipitation (SI > 1)
  //

  if (initScaledMass_ <= 0.0 || dissolutionRateConst_ <= 0.0) {
    return 0.0;
  }

  // Calculate driving force from actual SI
  double drivingForce = 1.0; // Default: maximum driving force

  if (saturationIndex > 0.0 && saturationIndex < 1.0) {
    // Undersaturated: dissolution occurs
    // f(SI) = (1 - SI^siexp)^dfexp
    drivingForce = pow((1.0 - pow(saturationIndex, siexp_)), dfexp_);
  } else if (saturationIndex >= 1.0) {
    // Supersaturated: precipitation occurs
    // Use absolute value of driving force for rate magnitude
    // f(SI) = (SI^siexp - 1)^dfexp
    drivingForce = pow((pow(saturationIndex, siexp_) - 1.0), dfexp_);
  }
  // If SI ≈ 0 or invalid, drivingForce remains 1.0 (maximum)

  // Estimate initial surface area
  double estimatedArea = initScaledMass_ * surfaceAreaMultiplier_;

  // dissrate in mol/100g/h with actual driving force
  double dissrate = dissolutionRateConst_ * rhFactor_ * estimatedArea *
                    drivingForce * arrhenius_;

  // Convert to mass fraction rate [1/h]
  double molarMass = chemSys_->getDCMolarMass(DCId_);
  double rate = (dissrate * molarMass) / initScaledMass_;

  return rate;
}

double StandardKineticModel::getCurrentMolarRate(double scaledMass) const {
  //
  // Calculate the current molar rate of dissolution/precipitation.
  //
  // The Standard model rate equation is:
  //   dissrate [mol/100g/h] = k * rh * area * f(SI) * arrhenius
  //
  // Returns positive for dissolution, negative for precipitation.
  //

  if (scaledMass <= 0.0 || dissolutionRateConst_ <= 0.0) {
    return 0.0;
  }

  // Get current surface area from lattice
  double area = lattice_->getSurfaceArea(microPhaseId_);
  if (area <= 1.0e-9) {
    area = 1.0;
  }
  area *= surfaceAreaMultiplier_;

  // Get current saturation index
  double saturationIndex = chemSys_->getMicroPhaseSI(microPhaseId_);

  // Calculate dissolution rate
  double dissrate = 0.0;

  if (saturationIndex < 1.0) {
    // Undersaturated: dissolution (positive rate)
    dissrate = dissolutionRateConst_ * rhFactor_ * area *
               pow((1.0 - pow(saturationIndex, siexp_)), dfexp_);
  } else {
    // Supersaturated: precipitation (negative rate)
    dissrate = -dissolutionRateConst_ * rhFactor_ * area *
               pow((pow(saturationIndex, siexp_) - 1.0), dfexp_);
  }

  dissrate *= arrhenius_;

  return dissrate;
}

double StandardKineticModel::computeNucleationVoxels(double dt_hours) const {
  if (!nucleation_.has_value()) return 0.0;

  // ChemicalSystem stores S already as the linear ratio IAP/Ksp; the log10
  // conversion from GEMS happens inside ChemicalSystem::calculateSI
  // (ChemicalSystem.cc:3244 via pow(10, Ph_SatInd(...))).
  double S = chemSys_->getMicroPhaseSI(microPhaseId_);
  if (S <= 1.0) return 0.0;

  double v_m = chemSys_->getDCMolarVolume(DCId_);           // [m^3/mol]
  double T_K = temperature_;                                 // [K]
  double V_voxel = lattice_->getVolumePerVoxel();            // [m^3]
  double V_electrolyte =
      lattice_->getVolumeFraction(ELECTROLYTEID) *
      static_cast<double>(lattice_->getNumSites()) * V_voxel; // [m^3]
  double dt_seconds = dt_hours * 3600.0;

  return cnt::voxelsPerCycle(*nucleation_, v_m, T_K, S,
                             dt_seconds, V_electrolyte, V_voxel);
}

void StandardKineticModel::accumulateNucleation(double dN) {
  if (dN > 0.0) nucleationAccumulator_ += dN;
}

int StandardKineticModel::drainNucleationInteger() {
  if (nucleationAccumulator_ < 1.0) return 0;
  int n = static_cast<int>(std::floor(nucleationAccumulator_));
  nucleationAccumulator_ -= static_cast<double>(n);
  return n;
}
