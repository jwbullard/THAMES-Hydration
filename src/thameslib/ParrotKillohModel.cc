/**
@file  ParrotKillohModel.cc
@brief Method definitions for the ParrotKillohModel class.

*/
#include "ParrotKillohModel.h"

using std::cout; using std::endl;

ParrotKillohModel::ParrotKillohModel() {
  ///
  /// Default value for w/c ratio in PK model is 0.45
  ///

  wsRatio_ = 0.45;

  ///
  /// Default value for specific surface area in PK model is 385
  /// m<sup>2</sup>/kg
  ///

  refSpecificSurfaceArea_ = RefSpecificSurfaceArea; // reference SSA (m2/kg)

  ///
  /// Default temperature in the PK model is 20 C (or 293 K)
  ///

  temperature_ = 293.15; // default temperature (K)
  refT_ = 293.15;        // default temperature (K)

  ///
  /// Clear out the vectors so they can be populated with values from the
  /// JSON input file
  ///

  name_ = "";
  microPhaseId_ = 2;
  DCId_ = 2;
  GEMPhaseId_ = 2;
  k1_ = 1.0;
  k2_ = 1.0;
  k3_ = 1.0;
  pfk_ = 1.0;
  n1_ = 1.0;
  n3_ = 1.0;
  activationEnergy_ = 0.0;
  scaledMass_ = 0.0;
  initScaledMass_ = 0.0;
  critDOR_ = 100.0;
  degreeOfReaction_ = 0.0;

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

ParrotKillohModel::ParrotKillohModel(ChemicalSystem *cs, Lattice *lattice,
                                     struct KineticData &kineticData,
                                     const bool verbose, const bool warning) {
  // Set the verbose and warning flags

  verbose_ = verbose;
  warning_ = warning;
#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
  cout << "ParrotKillohModel::ParrotKillohModel Constructor" << endl;
  cout.flush();
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  chemSys_ = cs;
  lattice_ = lattice;

  ///
  /// Default value for specific surface area in PK model is 385
  /// m<sup>2</sup>/kg
  ///

  specificSurfaceArea_ = lattice_->getCementSpecificSurfaceArea();

  refSpecificSurfaceArea_ = RefSpecificSurfaceArea;
  ssaFactor_ = specificSurfaceArea_ / refSpecificSurfaceArea_;

  ///
  /// Default value for w/c ratio in PK model is 0.45
  ///

  wsRatio_ = lattice_->getWsRatio();
  wcRatio_ = lattice_->getWcRatio();

  ///
  /// Default initial solid mass is 100 g
  ///

  initSolidMass_ = 100.0;

  ///
  /// Default temperature in the PK model is 20 C (or 293 K)
  ///

  lossOnIgnition_ = kineticData.loi;
  temperature_ = kineticData.temperature;
  refT_ = kineticData.reftemperature;

  modelName_ = "ParrotKillohModel";
  microPhaseId_ = kineticData.microPhaseId;
  name_ = kineticData.name;
  DCId_ = kineticData.DCId;
  GEMPhaseId_ = kineticData.GEMPhaseId;

  k1_ = kineticData.k1;
  k2_ = kineticData.k2;
  k3_ = kineticData.k3;

  /// Default to no pozzolanic influence on clinker phases

  pfk_ = 1.0;

  n1_ = kineticData.n1;
  n3_ = kineticData.n3;

  dorHcoeff_ = kineticData.dorHcoeff;

  activationEnergy_ = kineticData.activationEnergy;
  scaledMass_ = kineticData.scaledMass;
  initScaledMass_ = kineticData.scaledMass;

  critDOR_ = dorHcoeff_ * wcRatio_;
  degreeOfReaction_ = 0.0;

  double critporediam = lattice_->getLargestSaturatedPore(); // in nm
  critporediam *= 1.0e-9;                                    // in m
  rh_ = exp(-6.23527e-7 / critporediam / temperature_);
  rh_ = rh_ > 0.55 ? rh_ : 0.551;
  rhFactor_ = pow(((rh_ - 0.55) / 0.45), 4.0);

  arrhenius_ = exp((activationEnergy_ / GASCONSTANT) *
                   ((1.0 / refT_) - (1.0 / temperature_)));

  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion hours or 114,000 years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;

  return;
}

void ParrotKillohModel::calculateKineticStep(const double timestep,
                                             double &scaledMass, double &massDissolved,
                                             int cyc, double totalDOR, bool doTweak,
                                             bool &doNotModif) {

  ///
  /// Initialize local variables
  ///

  double DOR, newDOR;

  double wcFactor = 1.0;
  double ngrate = 1.0e-10;   // Nucleation and growth rate
  double hsrate = 1.0e-10;   // Hydration shell rate
  double diffrate = 1.0e-10; // Diffusion rate
  double rate = 1.0e-10;     // Selected rate

  ///
  /// Determine if this is a normal step or a necessary
  /// tweak from a failed GEM_run call
  ///

  try {

    // @todo BULLARD PLACEHOLDER
    // Still need to implement constant gas phase composition
    // Will involve equilibrating gas with aqueous solution
    //
    // First step each iteration is to equilibrate gas phase
    // with the electrolyte, while forbidding anything new
    // from precipitating.

    if (verbose_) {
      cout << "ParrotKillohModel::calculateKineticStep for " << name_ << endl;
      cout.flush();
    }

    if (!doTweak)
      cout << "      ParrotKillohModel::calculateKineticStep    - "
              "microPhaseId_/mPhName/SI : "       
           << setw(3) << right << microPhaseId_ << " / "
           << setw(15) << left << name_ << " / "
           << chemSys_->getMicroPhaseSI(microPhaseId_)
           << " (SI not used in PK model)" << endl;

    // RH factor is the same for all clinker phases

    /// This is a big kluge for internal relative humidity
    /// @note Using new gel and interhydrate pore size distribution model
    ///       which is currently contained in the Lattice object.
    ///
    /// Surface tension of water is gamma = 0.072 J/m2
    /// Molar volume of water is Vm = 1.8e-5 m3/mole
    /// The Kelvin equation is
    ///    p/p0 = exp (-4 gamma Vm / d R T) = exp (-6.23527e-7 / (d T))
    ///
    ///    where d is the pore diameter in meters and T is absolute temperature

    /// Assume a zero contact angle for now.
    /// @todo revisit the contact angle issue

    scaledMass_ = scaledMass;

    if (initScaledMass_ > 0.0) { // DOR @ t-1
      DOR = (initScaledMass_ - scaledMass_) / initScaledMass_;

      if (verbose_) {
        cout << "~~~~>DOR for " << name_ << " = " << DOR
             << "   initScaledMass_/scaledMass_/ssaFactor_ : " << initScaledMass_ << " / "
             << scaledMass_ << " / " << ssaFactor_ << endl;
        cout.flush();
      }
    } else {
      throw FloatException("ParrotKillohModel", "calculateKineticStep",
                           "initScaledMass_ = 0.0");
    }

    if (DOR < 1.0) {
      // Normal Parrott and Killoh implementation here

      if (fabs(n1_) > 0.0) {
        ngrate =
            (k1_ / n1_) * (1.0 - DOR) * pow((-log(1.0 - DOR)), (1.0 - n1_));
        ngrate *= (ssaFactor_); // only used for the N+G rate

        if (ngrate < 1.0e-10)
          ngrate = 1.0e-10;
      } else {
        throw FloatException("ParrotKillohModel", "calculateKineticStep",
                             "n1_ = 0.0");
      }

      hsrate = k3_ * pow((1.0 - DOR), n3_);
      if (hsrate < 1.0e-10)
        hsrate = 1.0e-10;

      if (DOR > 0.0) {
        diffrate = (k2_ * pow((1.0 - DOR), (2.0 / 3.0))) /
                   (1.0 - pow((1.0 - DOR), (1.0 / 3.0)));
        if (diffrate < 1.0e-10)
          diffrate = 1.0e-10;
      } else {
        diffrate = 1.0e9;
      }

      rate = (ngrate < hsrate) ? ngrate : hsrate;

      if (diffrate < rate)
        rate = diffrate;

      // double rate_ini = rate;

      rate *= (pfk_ * rhFactor_ * arrhenius_ / H_PER_DAY); // rate is R @ t-1

      double prod = rate * timestep;

      newDOR = DOR + prod;

      wcFactor = 1.0;

      if (totalDOR > critDOR_) {
        wcFactor = pow((1 + 3.333 * (critDOR_ - newDOR)), 4);
        prod = timestep * rate * wcFactor;
        newDOR = DOR + prod; // Could exceed 1.0 but we address that below
      }

      if (newDOR >= 1.0) {

        massDissolved = scaledMass_;

        scaledMass_ = 0.0;

      } else {

        scaledMass_ = initScaledMass_ * (1.0 - newDOR);

        // massDissolved = (newDOR - DOR) * initScaledMass_;
        massDissolved = initScaledMass_ * prod;
      }

      scaledMass = scaledMass_;

      if (verbose_) {
        cout << endl << "    ParrotKillohModel::calculateKineticStep "
                "rate/wcFactor/massDissolved : "
             << rate << " / " << wcFactor << " / " << massDissolved << endl;
        cout << "  ****************** PKM_hT = " << timestep
             << "    cyc = " << cyc << "    microPhaseId_ = " << microPhaseId_
             << "    microPhase = " << name_
             << "    GEMPhaseIndex = " << GEMPhaseId_ << " ******************"
             << endl;
        cout << "   PKM_hT   " << "pfk_: " << pfk_ << "    k1 = " << k1_
             << "    n1 = " << n1_ << "    k2 = " << k2_ << "    k3 = " << k3_
             << "    n3 = " << n3_ << endl;
        cout << "   PKM_hT   " << "dorHcoeff_: " << dorHcoeff_
             << "    Ea = " << activationEnergy_ << endl;
        cout << "   PKM_hT   "
             << "specificSurfaceArea_ = " << specificSurfaceArea_
             << "    refSpecificSurfaceArea_ = " << refSpecificSurfaceArea_
             << "    ssaFactor_ = " << ssaFactor_ << endl;
        cout << "   PKM_hT   " << "wcRatio_: " << wcRatio_
             << "\twsRatio_: " << wsRatio_ << endl;
        cout << "   PKM_hT   " << "ngrate: " << ngrate << "\thsrate: " << hsrate
             << "\tdiffrate: " << diffrate << "\trate: " << rate << endl;
        cout << "   PKM_hT   " << "wcFactor: " << wcFactor
             << "\trhFactor_: " << rhFactor_ << "\tarrhenius_: " << arrhenius_
             << endl;
        cout << "   PKM_hT   " << "DOR: " << DOR << "\tnewDOR: " << newDOR
             << "\ttotalDOR: " << totalDOR
             << "\tinitScaledMass_: " << initScaledMass_
             << "\tscaledMass_: " << scaledMass_
             << "\tmassDissolved: " << massDissolved << endl;
        cout.flush();
      }

    } else {
      throw DataException("ParrotKillohModel", "calculateKineticStep",
                          "DOR >= 1.0");
    }

  } // End of try block

  catch (EOBException eex) {
    eex.printException();
    exit(1);
  } catch (DataException dex) {
    dex.printException();
    exit(1);
  } catch (FloatException flex) {
    flex.printException();
    exit(1);
  } catch (out_of_range &oor) {
    EOBException eex("ParrotKillohModel", "calculateKineticStep", oor.what(), 0,
                     0);
    eex.printException();
    exit(1);
  }

  return;
}
