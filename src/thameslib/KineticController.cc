/**
@file  KineticController.cc
@brief Method definitions for the KineticController class.

*/
#include "KineticController.h"

using std::cout; using std::endl;
using std::string; using std::vector;

KineticController::KineticController() {
  temperature_ = 293.15;

  // default temperature (K)
  refT_ = 293.15;

  ///
  /// Clear out the vectors so they can be populated with values from the
  ///

  numPhases_ = 0;
  chemSys_ = NULL;
  lattice_ = NULL;
  phaseKineticModel_.clear();
  name_.clear();
  initScaledMass_.clear();
  scaledMass_.clear();
  specificSurfaceArea_.clear();
  refSpecificSurfaceArea_.clear();
  isKinetic_.clear();
  // ICName_.clear();
  DCNum_ = 0;
  // DCName_.clear();
  GEMPhaseNum_ = 0;

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion hours or 114K years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;
  beginAttackTime_ = 1.0e10;

  verbose_ = warning_ = false;

  return;
}

KineticController::KineticController(ChemicalSystem *cs, Lattice *lattice,
                                     const string &jsonFileName,
                                     const bool verbose, const bool warning)
    : chemSys_(cs), lattice_(lattice) {
  ///
  /// Clear out the vectors so they can be populated with values from the
  ///

  numPhases_ = 0;
  phaseKineticModel_.clear();
  name_.clear();
  isKinetic_.clear();

  // Set the verbose and warning flags

#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
  cout << "KineticController::KineticController Constructor" << endl;
  cout.flush();
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  ///
  /// Default temperature in the PK model is 20 C (or 293 K)
  ///

  temperature_ = 293.15;
  refT_ = 293.15;

  ///
  /// Clear out the vectors so they can be populated with values from the
  /// JSON input file
  ///

  name_.clear();
  microPhaseId_.clear();

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion hours or 114K years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;
  beginAttackTime_ = 1.0e10;

  ///
  /// Open the input JSON file for kinetic data and parse it
  ///

  string jsonext = ".json";
  size_t foundjson;
  foundjson = jsonFileName.find(jsonext);
  try {
    if (foundjson != string::npos) {
      if (verbose_) {
        cout << "KineticModel data file is a JSON file" << endl;
      }
      parseDoc(jsonFileName);
    } else {
      throw FileException("KineticModel", "KineticModel", jsonFileName,
                          "NOT in JSON format");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  } catch (DataException dex) {
    dex.printException();
    exit(1);
  }

  int microPhaseId;

  if (verbose_) {
    cout << "KineticController::KineticController Finished reading "
            "simparams.json "
         << endl;
    int size = microPhaseId_.size();
    for (int i = 0; i < size; ++i) {
      microPhaseId = microPhaseId_[i];
      if (isKinetic_[i]) {
        cout << "KineticController::KineticController kinetic phase "
             << microPhaseId << endl;
        cout << "KineticController::KineticController     name = "
             << chemSys_->getMicroPhaseName(microPhaseId) << endl;
      }
    }
    cout.flush();
  }

  // Assign the DC index for water

  DCNum_ = chemSys_->getNumDCs();
  // ICName_ = chemSys_->getICName();
  // DCName_ = chemSys_->getDCName();
  GEMPhaseNum_ = chemSys_->getNumGEMPhases();

  // ICMolesTot_.resize(ICNum_, 0.0);
  DCMoles_.resize(DCNum_, 0.0);
  DCMolesIni_.resize(DCNum_, 0.0);

  calcPhaseMasses();

  waterDCId_ = chemSys_->getDCId("H2O@");

  pKMsize_ = phaseKineticModel_.size();
  impurity_K2O_.resize(pKMsize_, 0);
  impurity_Na2O_.resize(pKMsize_, 0);
  impurity_Per_.resize(pKMsize_, 0);
  impurity_SO3_.resize(pKMsize_, 0);

  impurityDCID_.clear();
  impurityDCID_.push_back(chemSys_->getDCId("K2O"));
  impurityDCID_.push_back(chemSys_->getDCId("Na2O"));
  impurityDCID_.push_back(chemSys_->getDCId("Per"));
  impurityDCID_.push_back(chemSys_->getDCId("SO3"));

  // initScaledMass_, scaledMass_ & scaledMassIni_ are
  // initialized in KineticController::parseMicroPhase :
  // initScaledMass_.push_back(0.0);
  // scaledMass_.push_back(0.0);
  // scaledMassIni_.push_back(0.0);

  string modelName;
  int phID;
  initScaledCementMass_ = 0;
  cout << endl << "KineticController::KineticController(...) :" << endl;
  cout << "  - only these phases (controlled by the ParrotKillohModel or C3SKineticModel) "
          "contribute to initScaledCementMass_ & scaledCementMass_ :"
       << endl;

  for (int i = 0; i < pKMsize_; i++) {
    modelName = phaseKineticModel_[i]->getModelName();
    // cout << endl << "    modelName = " << modelName << endl;
    if (modelName == "ParrotKillohModel" || modelName == "C3SKineticModel") {
      phID = phaseKineticModel_[i]->getMicroPhaseId();
      initScaledCementMass_ += chemSys_->getMicroPhaseMass(phID);
      cout << "      microPhaseID/microPhaseName/microPhaseMass : " << setw(3)
           << right << phID << " / " << setw(15) << left
           << phaseKineticModel_[i]->getName() << " / "
           << chemSys_->getMicroPhaseMass(phID) << " g" << endl;
      chemSys_->setIsParrotKilloh(phID);
    }
  }
  cout << endl
       << "      initScaledCementMass_ = " << initScaledCementMass_
       << " g (same value for scaledCementMass_)" << endl;
  chemSys_->setInitScaledCementMass(initScaledCementMass_);

  hydTimeIni_ = 0;

  return;
}

KineticController::~KineticController() {
  for (int midx = 0; midx < pKMsize_; ++midx) {
    delete phaseKineticModel_[midx];
  }
}

void KineticController::parseDoc(const string &docName) {
  int numEntry = -1; // Tracks number of solid phases

  ///
  /// The kineticData structure is used to temporarily hold parsed data
  /// for a given phase before the data are loaded permanently into class
  /// members.
  ///

  struct KineticData kineticData;

  /// Test for JSON existence

  ifstream f(docName.c_str());

  /// Parse the JSON file all at once
  json data = json::parse(f);
  f.close();

  try {

    /// Get an iterator to the root node of the JSON file
    /// @todo Add a better JSON validity check.

    json::iterator it = data.find("environment");
    json::iterator cdi = it.value().begin();

    // Test for non-emptiness
    if (cdi == it.value().end() || it == data.end()) {
      throw FileException("Controller", "parseDoc", docName, "Empty JSON file");
    }

    cdi = it.value().find("temperature");
    temperature_ = cdi.value();
    cdi = it.value().find("reftemperature");
    refT_ = cdi.value();

    it = data.find("microstructure");

    /// Each phase is a more complicated grouping of data that
    /// has a separate method for parsing.

    try {
      parseMicroPhases(it, numEntry, kineticData);
    } catch (DataException dex) {
      throw dex;
    }

    /// Push a copy of the isKinetic vector to the ChemicalSystem

    chemSys_->setIsKinetic(isKinetic_);

  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  /// All kinetic components have been parsed now.  Next, this block tries
  /// to handle pozzolanic effects (loi, SiO2 content, etc.) on any other
  /// kinetic phases

  setPozzEffectOnPK();

  return;
}

void KineticController::parseMicroPhases(const json::iterator it, int &numEntry,
                                         struct KineticData &kineticData) {
  string testname;
  bool kineticfound = false;
  bool ispozz = false;
  bool isParrotKilloh = false;

  json::iterator cdi = it.value().find("phases");
  json::iterator p = cdi.value().begin();

  for (int i = 0; i < static_cast<int>(cdi.value().size()); ++i) {
    initKineticData(kineticData);
    isKinetic_.push_back(false);
    p = cdi.value()[i].find("thamesname");
    testname = p.value();
    kineticData.name = testname;
    kineticData.microPhaseId = chemSys_->getMicroPhaseId(testname);
    kineticfound = ispozz = isParrotKilloh = false;

    p = cdi.value()[i].find("kinetic_data");
    if (p != cdi.value()[i].end()) {
      numEntry += 1;
      kineticfound = true;
      isKinetic_[isKinetic_.size() - 1] = true;
      kineticData.GEMPhaseId =
          chemSys_->getMicroPhaseToGEMPhase(kineticData.microPhaseId, 0);
      kineticData.DCId =
          chemSys_->getMicroPhaseDCMembers(kineticData.microPhaseId, 0);

      ///
      /// Kinetic data are grouped together,
      /// so there is a method written just for parsing that grouping
      ///

      try {
        parseKineticData(p, kineticData);
      } catch (DataException dex) {
        throw dex;
      }
    }
    if (kineticfound) {
      kineticData.scaledMass =
          chemSys_->getMicroPhaseMass(kineticData.microPhaseId);
      kineticData.temperature = temperature_;
      kineticData.reftemperature = refT_;
      makeModel(kineticData);
    }

    /// Some items should be added to vectors whether kinetically controlled or
    /// not

    name_.push_back(kineticData.name);
    microPhaseId_.push_back(kineticData.microPhaseId);
    initScaledMass_.push_back(0.0);
    scaledMass_.push_back(0.0);
    scaledMassIni_.push_back(0.0);
  }

  return;
}

void KineticController::parseKineticData(const json::iterator p,
                                         struct KineticData &kineticData) {
  bool typefound = false;

  try {
    json::iterator pp = p.value().find("type");
    kineticData.type = pp.value();
    if (kineticData.type == ParrotKillohType) {
      typefound = true;
      try {
        parseKineticDataForParrotKilloh(p, kineticData);
      } catch (DataException dex) {
        throw dex;
      }
    } else if (kineticData.type == StandardType) {
      typefound = true;
      try {
        parseKineticDataForStandard(p, kineticData);
      } catch (DataException dex) {
        throw dex;
      }
    } else if (kineticData.type == PozzolanicType) {
      typefound = true;
      try {
        parseKineticDataForPozzolanic(p, kineticData);
      } catch (DataException dex) {
        throw dex;
      }
    } else if (kineticData.type == C3SType) {
      typefound = true;
      try {
        parseKineticDataForC3S(p, kineticData);
      } catch (DataException dex) {
        throw dex;
      }
    } else {
      throw HandleException("KineticController", "parseKineticData", "type",
                            "Model type not found");
    }

    if (!typefound) {
      throw HandleException("KineticController", "parseKineticData", "type",
                            "Model type not specified");
    }
  } catch (HandleException hex) {
    hex.printException();
  }

  return;
}

void KineticController::parseKineticDataForParrotKilloh(
    const json::iterator p, struct KineticData &kineticData) {

  if (verbose_) {
    cout << "--->Parsing PK data for " << kineticData.name << endl;
    cout.flush();
  }

  // Parrot-Killoh k1 parameter is like a rate const
  // Conventionally given in units of g per day
  // Immediately convert to units of g per h within model
  json::iterator pp = p.value().find("k1");
  kineticData.k1 = pp.value();

  // Parrot-Killoh k2 parameter
  // Conventionally given in units of g per day
  // Immediately convert to units of g per h within model
  pp = p.value().find("k2");
  kineticData.k2 = pp.value();

  // Parrot-Killoh k3 parameter
  // Conventionally given in units of g per day
  // Immediately convert to units of g per h within model
  pp = p.value().find("k3");
  kineticData.k3 = pp.value();

  // Parrot-Killoh n1 parameter
  pp = p.value().find("n1");
  kineticData.n1 = pp.value();

  // Parrot-Killoh n3 parameter
  pp = p.value().find("n3");
  kineticData.n3 = pp.value();

  // Parrot-Killoh DOR_Hcoeff parameter
  pp = p.value().find("dorHcoeff");
  kineticData.dorHcoeff = pp.value();

  // Activation energy
  pp = p.value().find("activationEnergy");
  kineticData.activationEnergy = pp.value();

  return;
}

void KineticController::parseKineticDataForC3S(
    const json::iterator p, struct KineticData &kineticData) {

  if (verbose_) {
    cout << "--->Parsing C3S kinetic data for " << kineticData.name
         << endl;
    cout.flush();
  }

  // How much to multiply the microstructure phase's surface
  // area to account for unresolved structure, roughness, etc.
  // This is an optional input, default value is 1.0
  json::iterator pp = p.value().find("surfaceAreaMultiplier");
  if (pp != p.value().end()) {
    kineticData.surfaceAreaMultiplier = pp.value();
  } else {
    kineticData.surfaceAreaMultiplier = 1.0;
  }

  // C3S constants
  // Dissolution rate constant input with units (mol/m2/s)
  // But immediately convert it to mol/m2/h within model
  pp = p.value().find("kC3S");
  if (pp != p.value().end()) {
    kineticData.kC3S = pp.value();
    kineticData.kC3S *= S_PER_H;
  } else {
    throw DataException("KineticController", "parseKineticDataForC3S",
                        "kC3S not found");
  }

  pp = p.value().find("C1");
  if (pp != p.value().end()) {
    kineticData.C1 = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForC3S",
                        "C1 not found");
  }

  pp = p.value().find("r");
  if (pp != p.value().end()) {
    kineticData.r = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForC3S",
                        "r not found");
  }

  // Activation energy
  pp = p.value().find("activationEnergy");
  if (pp != p.value().end()) {
    kineticData.activationEnergy  = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForC3S",
                        "activationEnergy not found");
  }

  return;
}

void KineticController::parseKineticDataForStandard(
    const json::iterator p, struct KineticData &kineticData) {

  if (verbose_) {
    cout << "--->Parsing standard kinetic data for " << kineticData.name
         << endl;
    cout.flush();
  }

  // How much to multiply the microstructure phase's surface
  // area to account for unresolved structure, roughness, etc.
  // This is an optional input, default value is 1.0
  json::iterator pp = p.value().find("surfaceAreaMultiplier");
  if (pp != p.value().end()) {
    kineticData.surfaceAreaMultiplier = pp.value();
  } else {
    kineticData.surfaceAreaMultiplier = 1.0;
  }

  // Dissolution rate constant input with units (mol/m2/s)
  // But immediately convert it to mol/m2/h within model
  pp = p.value().find("dissolutionRateConst");
  if (pp != p.value().end()) {
    kineticData.dissolutionRateConst = pp.value();
    kineticData.dissolutionRateConst *= S_PER_H;
  } else {
    throw DataException("KineticController", "parseKineticDataForStandard",
                        "dissolutionRateConst not found");
  }

  // Rate constant for early-age diffusion input with units (mol/m2/s)
  // But immediately convert it to mol/m2/h within model
  pp = p.value().find("diffusionRateConstEarly");
  if (pp != p.value().end()) {
    kineticData.diffusionRateConstEarly = pp.value();
    kineticData.diffusionRateConstEarly *= S_PER_H;
  } else {
    throw DataException("KineticController", "parseKineticDataForStandard",
                        "diffusionRateConstEarly not found");
  }

  // Dissolution rate constant input with units (mol/m2/s)
  // But immediately convert it to mol/m2/h within model
  pp = p.value().find("diffusionRateConstLate");
  if (pp != p.value().end()) {
    kineticData.diffusionRateConstLate = pp.value();
    kineticData.diffusionRateConstLate *= S_PER_H;
  } else {
    kineticData.diffusionRateConstLate = kineticData.diffusionRateConstEarly;
    cout << "WARNING: For " << kineticData.name
         << " diffusionRateConstLate not found; setting it to "
            "diffusionRateConstEarly"
         << endl;
  }

  // Number of DC units produced in dissociation reaction
  pp = p.value().find("dissolvedUnits");
  if (pp != p.value().end()) {
    kineticData.dissolvedUnits = pp.value();
  } else {
    kineticData.dissolvedUnits = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " dissolvedUnits not found; setting it to 1" << endl;
  }

  // Exponent on  the saturation index in the rate equation
  pp = p.value().find("siexp");
  if (pp != p.value().end()) {
    kineticData.siexp = pp.value();
  } else {
    kineticData.sio2 = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " sio2 not found; setting it to 1" << endl;
  }

  // Exponent on  the driving force term in the rate equation
  pp = p.value().find("dfexp");
  if (pp != p.value().end()) {
    kineticData.dfexp = pp.value();
  } else {
    kineticData.dfexp = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " dfexp not found; setting it to 1" << endl;
  }

  // Loss on ignition of the material
  pp = p.value().find("loi");
  if (pp != p.value().end()) {
    kineticData.loi = pp.value();
  } else {
    kineticData.loi = 1.0e-6;
    cout << "WARNING: For " << kineticData.name
         << " loi not found; setting it to 1.0e-6" << endl;
  }

  // Activation energy for dissolution
  pp = p.value().find("activationEnergy");
  if (pp != p.value().end()) {
    kineticData.activationEnergy = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForStandard",
                        "activationEnergy not found");
  }

  return;
}

void KineticController::parseKineticDataForPozzolanic(
    const json::iterator p, struct KineticData &kineticData) {

  if (verbose_) {
    cout << "--->Parsing pozzolanic data for " << kineticData.name << endl;
    cout.flush();
  }

  // How much to multiply the microstructure phase's surface
  // area to account for unresolved structure, roughness, etc.
  // This is an optional input, default value is 1.0
  json::iterator pp = p.value().find("surfaceAreaMultiplier");
  if (pp != p.value().end()) {
    kineticData.surfaceAreaMultiplier = pp.value();
  } else {
    kineticData.surfaceAreaMultiplier = 1.0;
  }

  // Dissolution rate constant input with units (mol/m2/s)
  // But immediately convert it to mol/m2/h within model
  pp = p.value().find("dissolutionRateConst");
  if (pp != p.value().end()) {
    kineticData.dissolutionRateConst = pp.value();
    kineticData.dissolutionRateConst *= S_PER_H;
  } else {
    throw DataException("KineticController", "parseKineticDataForPozzolanic",
                        "dissolutionRateConst not found");
  }

  // Early-age diffusion rate constant input with units (mol/m2/s)
  // But immediately convert it to mol/m2/h within model
  pp = p.value().find("diffusionRateConstEarly");
  if (pp != p.value().end()) {
    kineticData.diffusionRateConstEarly = pp.value();
    kineticData.diffusionRateConstEarly *= S_PER_H;
  } else {
    throw DataException("KineticController", "parseKineticDataForPozzolanic",
                        "diffusionRateConstEarly not found");
  }

  // Later-age diffusion rate constant input with units (mol/m2/s)
  // But immediately convert it to mol/m2/h within model
  pp = p.value().find("diffusionRateConstLate");
  if (pp != p.value().end()) {
    kineticData.diffusionRateConstLate = pp.value();
    kineticData.diffusionRateConstLate *= S_PER_H;
  } else {
    kineticData.diffusionRateConstLate = kineticData.diffusionRateConstEarly;
    cout << "WARNING: For " << kineticData.name
         << " diffusionRateConstLate not found; setting it to "
            "diffusionRateConstEarly"
         << endl;
  }

  // Number of DC units produced in dissociation reaction
  pp = p.value().find("dissolvedUnits");
  if (pp != p.value().end()) {
    kineticData.dissolvedUnits = pp.value();
  } else {
    kineticData.dissolvedUnits = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " dissolvedUnits not found; setting it to 1" << endl;
  }

  // Exponent on the saturation index in the rate equation
  pp = p.value().find("siexp");
  if (pp != p.value().end()) {
    kineticData.siexp = pp.value();
  } else {
    kineticData.siexp = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " siexp not found; setting it to 1" << endl;
  }

  // Exponent on the driving force term in the rate equation
  pp = p.value().find("dfexp");
  if (pp != p.value().end()) {
    kineticData.dfexp = pp.value();
  } else {
    kineticData.dfexp = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " dfexp not found; setting it to 1" << endl;
  }

  // Exponent on the degree of reaction term in the diffusion rate equation
  pp = p.value().find("dorexp");
  if (pp != p.value().end()) {
    kineticData.dorexp = pp.value();
  } else {
    kineticData.dorexp = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " dorexp not found; setting it to 1" << endl;
  }

  // Exponent on the hydroxy ion activity in the rate equation
  pp = p.value().find("ohexp");
  if (pp != p.value().end()) {
    kineticData.ohexp = pp.value();
  } else {
    kineticData.ohexp = 1.0;
    cout << "WARNING: For " << kineticData.name
         << " ohexp not found; setting it to 1" << endl;
  }

  // SiO2 mass fraction in the material
  pp = p.value().find("sio2");
  if (pp != p.value().end()) {
    kineticData.sio2 = pp.value();
    // } else {
    //   throw DataException("KineticController",
    //   "parseKineticDataForPozzolanic",
    //                       "sio2 not found");
  }

  // Al2O3 mass fraction in the material
  pp = p.value().find("al2o3");
  if (pp != p.value().end()) {
    kineticData.al2o3 = pp.value();
    // } else {
    //   throw DataException("KineticController",
    //   "parseKineticDataForPozzolanic",
    //                       "al2o3 not found");
  }

  // CaO mass fraction in the material
  pp = p.value().find("cao");
  if (pp != p.value().end()) {
    kineticData.cao = pp.value();
    // } else {
    //   throw DataException("KineticController",
    //   "parseKineticDataForPozzolanic",
    //                       "cao not found");
  }

  // Loss on ignition of the material
  pp = p.value().find("loi");
  if (pp != p.value().end()) {
    kineticData.loi = pp.value();
    // } else {
    //   throw DataException("KineticController",
    //   "parseKineticDataForPozzolanic",
    //                       "loi not found");
  }

  pp = p.value().find("activationEnergy");
  if (pp != p.value().end()) {
    kineticData.activationEnergy = pp.value();
  } else {
    throw DataException("KineticController", "parseKineticDataForPozzolanic",
                        "activationEnergy not found");
  }

  return;
}

void KineticController::calcPhaseMasses(void) {
  int microPhaseId;
  double pscaledMass = 0.0;

  int size = microPhaseId_.size();

  for (int i = 0; i < size; i++) {
    microPhaseId = microPhaseId_[i];
    if (microPhaseId != VOIDID && microPhaseId != ELECTROLYTEID) {
      pscaledMass = chemSys_->getMicroPhaseMass(microPhaseId);
      scaledMass_[i] = pscaledMass;
      initScaledMass_[i] = pscaledMass;
      scaledMassIni_[i] = pscaledMass;

      // Setting the phase mass will also automatically calculate the phase
      // volume

      if (verbose_) {
        cout
            << "KineticController::getPhaseMasses reads solid micphase mass of "
            << chemSys_->getMicroPhaseName(microPhaseId) << " as "
            << initScaledMass_[i] << endl;
        cout.flush();
      }
    }
  }

  return;
}

double KineticController::getSolidMass(void) {
  int microPhaseId;
  double totmass = 0.0;
  int size = microPhaseId_.size();
  for (int i = 0; i < size; i++) {
    microPhaseId = microPhaseId_[i];
    if (microPhaseId != VOIDID && microPhaseId != ELECTROLYTEID) {
      totmass += chemSys_->getMicroPhaseMass(microPhaseId);
    }
  }

  return (totmass);
}

void KineticController::makeModel(struct KineticData &kineticData) {
  KineticModel *km = NULL;

  if (kineticData.type == ParrotKillohType) {
    // Read remaining Parrot and Killoh model parameters
    km = new ParrotKillohModel(chemSys_, lattice_, kineticData, verbose_,
                               warning_);
  } else if (kineticData.type == StandardType) {
    // Read remaining Standard model parameters
    km = new StandardKineticModel(chemSys_, lattice_, kineticData, verbose_,
                                  warning_);
  } else if (kineticData.type == PozzolanicType) {
    // Read remaining Pozzolanic model parameters
    km = new PozzolanicModel(chemSys_, lattice_, kineticData, verbose_,
                             warning_);
  } else if (kineticData.type == C3SType) {
    // Read remaining C3S model parameters
    km = new C3SKineticModel(chemSys_, lattice_, kineticData, verbose_,
                                  warning_);
  }

  phaseKineticModel_.push_back(km);

  return;
}

void KineticController::setPozzEffectOnPK(void) {

  /// @todo This is the block where the influence of some components on the
  /// kinetic parameters of other components can be set.

  double refloi = 0.8;
  double loi = refloi;
  double maxloi = refloi;
  double sio2val = 0.94;
  double refsio2val = 0.94;
  double phaseSA = 29.0;
  double refPhaseSA = 29.0;
  double minpozzeffect = 1.0;
  double pozzeffect = 1.0;

  int size = phaseKineticModel_.size();

  for (int midx = 0; midx < size; ++midx) {
    loi = phaseKineticModel_[midx]->getLossOnIgnition();
    if (loi > maxloi)
      maxloi = loi;
    if (phaseKineticModel_[midx]->getType() == PozzolanicType) {
      sio2val = phaseKineticModel_[midx]->getSio2();
      phaseSA =
          lattice_->getSurfaceArea(phaseKineticModel_[midx]->getMicroPhaseId());
      pozzeffect = pow((sio2val / refsio2val), 2.0) * (phaseSA / refPhaseSA);
      if (pozzeffect < minpozzeffect)
        minpozzeffect = pozzeffect > 0.8 ? pozzeffect : 0.8;
      cout << endl
           << "KineticController::setPozzEffectOnPK for midx = " << midx
           << " (microPhaseId =  "
           << phaseKineticModel_[midx]->getMicroPhaseId()
           << ", microPhaseName = " << phaseKineticModel_[midx]->getName()
           << endl;

      cout << "  Ref LOI = " << refloi << endl;
      cout << "  LOI     = " << loi << endl;
      cout << "  Max LOI = " << maxloi << endl;
      cout << "  SiO2     = " << sio2val << endl;
      cout << "  Ref SiO2 = " << refsio2val << endl;
      cout << "  BET      = " << phaseSA << endl;
      cout << "  Ref BET  = " << refPhaseSA << endl;
      cout << "  Pozz Effect     = " << pozzeffect << endl;
      cout << "  Min Pozz Effect = " << minpozzeffect << endl;
      cout.flush();
    }
  }

  minpozzeffect *= (refloi / maxloi);

  /// The way this is set up, 0.0 <= refloi / maxloi <= 1.0
  for (int midx = 0; midx < size; ++midx) {
    if (phaseKineticModel_[midx]->getType() == ParrotKillohType) {
      phaseKineticModel_[midx]->setPfk(minpozzeffect);
    }
  }

  return;
}

void KineticController::calculateKineticStep(double time, const double timestep,
                                             int cyc) {
  ///
  /// Initialize local variables
  ///
  ///

  int i;

  // double massDissolved = 0.0;
  cout << scientific << setprecision(15);
  ///
  /// Determine if this is a normal step or a necessary
  /// tweak from a failed GEM_run call
  ///

  // vector<int> impurityDCID;
  // impurityDCID.clear();
  // impurityDCID.push_back(chemSys_->getDCId("K2O"));
  // impurityDCID.push_back(chemSys_->getDCId("Na2O"));
  // impurityDCID.push_back(chemSys_->getDCId("Per")); // 170
  // impurityDCID.push_back(chemSys_->getDCId("SO3"));

  // cout << endl << "impurityDCID : " << endl;
  // for(i = 0; i < chemSys_->getNumMicroImpurities(); i++){
  //     cout << i << "\t" << impurityDCID[i] << endl; cout.flush();
  // }
  // cout << endl ;

  double totMassImpurity, massImpurity;

  int DCId;
  // int pKMsize = phaseKineticModel_.size();
  // static vector<double> scaledMassIni;
  double keepNumDCMoles;
  vector<int> phaseDissolvedId;
  phaseDissolvedId.resize(pKMsize_, 0);
  double numDCMolesDissolved, scaledMass, massDissolved;

  double hyd_time = hydTimeIni_ + timestep;

  chemSys_->initDCLowerLimit(0); // check !

  bool doTweak = (chemSys_->getTimesGEMFailed() > 0) ? true : false;

  if (doTweak) {
    // hyd_time = hydTimeIni_ + timestep;
    if (verbose_)
      cout << endl
           << "    KineticController::calculateKineticStep - tweak cyc = " << cyc
           << " :  hyd_time = " << hyd_time
           << "   hydTimeIni_ = " << hydTimeIni_ << "   timestep = " << timestep
           << endl;

    for (int midx = 0; midx < pKMsize_; ++midx) {
      phaseDissolvedId[midx] = phaseKineticModel_[midx]->getMicroPhaseId();
      chemSys_->setMicroPhaseMass(phaseDissolvedId[midx], scaledMassIni_[midx]);
      if (verbose_) {
        cout << "       midx = " << midx
             << "     scaledMassIni[midx] = " << scaledMassIni_[midx]
             << "     microPhaseName = " << phaseKineticModel_[midx]->getName()
             << endl;
      }
    }

    for (i = 0; i < DCNum_; i++) {
      DCMoles_[i] = DCMolesIni_[i];
    }

    lattice_->resetSurfaceArea(surfaceAreaIni_);

  } else {

    // hyd_time = hydTimeIni_ + timestep;
    cout << "    KineticController::calculateKineticStep - cyc = " << cyc
         << " :  hyd_time = " << hyd_time << "   hydTimeIni_ = " << hydTimeIni_
         << "   timestep = " << timestep << endl;

    for (int midx = 0; midx < pKMsize_; ++midx) {
      phaseDissolvedId[midx] = phaseKineticModel_[midx]->getMicroPhaseId();
      scaledMassIni_[midx] =
          chemSys_->getMicroPhaseMass(phaseDissolvedId[midx]);
      if (verbose_) {
        cout << "      midx = " << midx
             << "     scaledMassIni[midx] = " << scaledMassIni_[midx]
             << "     microPhaseName = " << phaseKineticModel_[midx]->getName()
             << endl;
      }
    }

    for (i = 0; i < DCNum_; i++) {
      DCMoles_[i] = chemSys_->getDCMoles(i);
      DCMolesIni_[i] = DCMoles_[i];
    }
    surfaceAreaIni_ = lattice_->getSurfaceArea();
  }

  if (hyd_time <= beginAttackTime_) {

    try {
      // cout << "  KineticController::calculateKineticStep     hyd_time = "
      //      << hyd_time << "\tcyc = " << cyc << endl;

      // if (!doTweak) {
      //  @todo BULLARD PLACEHOLDER
      //  Still need to implement constant gas phase composition
      //  Will involve equilibrating gas with aqueous solution
      //
      //  First step each iteration is to equilibrate gas phase
      //  with the electrolyte, while forbidding anything new
      //  from precipitating.

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

      /// Loop over all kinetic models

      //*******
      double totalDOR = 0;

      /// @note The totalDOR is defined only as the combined degree of hydration
      /// of "cement" components, which the user defines. This was intended to
      /// be only portland cement clinker components; now is depending on user decision
      /// (arcanite, thenardite, gypsum, bassanite and hemihydrate can belong to this
      /// category)

      if (initScaledCementMass_ > 0) {
        totalDOR = (initScaledCementMass_ - chemSys_->getScaledCementMass()) /
                   initScaledCementMass_;
      } else {
        int numPKMorC3SKMphases = 0;
        for (int midx = 0; midx < pKMsize_; ++midx) {
          if (phaseKineticModel_[midx]->getModelName() == "ParrotKillohModel" ||
              phaseKineticModel_[midx]->getModelName() == "C3SKineticModel")
            numPKMorC3SKMphases++;
        }

        // This next block only if there are ParrotKilloh or C3S kinetic model phases
        if (numPKMorC3SKMphases > 0) {

          cout << endl
               << "   KineticController::calculateKineticStep error - "
                  "initScaledCementMass_ = 0  while numPKMorC3SKMphases = "
               << numPKMorC3SKMphases << " :" << endl;
          for (int midx = 0; midx < pKMsize_; ++midx) {
            phaseDissolvedId[midx] =
                phaseKineticModel_[midx]->getMicroPhaseId();
            scaledMassIni_[midx] =
                chemSys_->getMicroPhaseMass(phaseDissolvedId[midx]);
            cout << "     midx = " << midx
                 << "     scaledMassIni[midx] = " << scaledMassIni_[midx]
                 << "     microPhaseName = "
                 << phaseKineticModel_[midx]->getName() << endl;
          }
          cout << endl
               << "        cyc/doTweak/timesGEMFailed : " << cyc << " / "
               << doTweak << " / " << chemSys_->getTimesGEMFailed() << endl;
          throw FloatException("KineticController", "calculateKineticStep",
                               "initScaledCementMass_ = 0");
        }
      }

      if (totalDOR < 0) {
        cout << endl
             << "     KineticController::calculateKineticStep error : totalDOR "
                "< 0" << endl;
        cout << endl
             << "        cyc/doTweak/timesGEMFailed : " << cyc << " / "
             << doTweak << " / " << chemSys_->getTimesGEMFailed() << endl;
        cout << endl
             << "        initScaledCementMass_/scaledCementMass/totalDOR : "
             << initScaledCementMass_ << " / "
             << chemSys_->getScaledCementMass() << " / " << totalDOR << endl;
        throw DataException("KineticController", "calculateKineticStep",
                            "totalDOR < 0");
      }

      if (!doTweak) {
        cout << "    KineticController::calculateKineticStep - cyc = " << cyc
             << " :  scaledCementMass = " << chemSys_->getScaledCementMass()
             << "   totalDOR = " << totalDOR << endl;
      }

      //*******

      double dcmoles;
      bool runKM = true;

      for (int midx = 0; midx < pKMsize_; ++midx) {
        scaledMass = scaledMassIni_[midx];
        runKM = true;
        if (scaledMass == 0.0 &&
            ((phaseKineticModel_[midx]->getModelName() == "ParrotKillohModel") ||
             (phaseKineticModel_[midx]->getModelName() == "C3SKineticModel"))) {
          runKM = false;
        }

        if (runKM) {
          bool doNotModif = false;

          DCId = phaseKineticModel_[midx]->getDCId();

          // ONLY HERE IS WHERE WE CALCULATE THE SURFACE AREA EACH CYCLE
          // THIS WAY WE ONLY CALCULATE THE ONES WE NEED
          /// @todo Maybe we can do this more efficiently, like just updating
          /// the
          ///         local changes instead of recalculating each time.

          // We calculate surface area even for PK model phases although
          // they currently don't use this information but instead rely
          // on the user-prescribed Blaine fineness of the starting powder only.
          // We do this in anticipation of eventually discarding the PK model

          lattice_->calcSurfaceArea(phaseDissolvedId[midx]);
          massDissolved = 0.0;

          /// @note The totalDOR that is passed to the kinetic model is
          /// based SOLELY on the combined degree of reaction of
          /// "cement" components. It is intended for the Parrot-Killoh
          /// model usage.

          // cout << endl << "    KineticController::calculateKineticStep - 0 - cyc = " << cyc
          //      << " :  doNotModif = false <-> doNotModif = " << doNotModif
          //      << "   midx = " << midx << "   phName = " << phaseKineticModel_[midx]->getName()
          //      << "   DCId = " << DCId << "   scaledMass = " << scaledMass
          //      << "   massDissolved = " << massDissolved
          //      << "   keepDCLowerLimit = " << chemSys_->getKeepDCLowerLimit(DCId)
          //      << "   DCMoles_ = " << DCMoles_[DCId]
          //      << "   DCMoles_cs = " << chemSys_->getDCMoles(DCId)
          //     << endl;

          phaseKineticModel_[midx]->calculateKineticStep(
              timestep, scaledMass, massDissolved, cyc, totalDOR, doTweak, doNotModif);

          /// @note may want to change the condition of next block
          /// because it is possible for the scaled mass of a kinetic phase to
          /// be zero but not negative. We could allow a threshold number
          /// like -1.0e-9 or something like that to handle "nearly zero"

          if (scaledMass < 0.0) {
            cout << endl
                 << "KineticController::calculateKineticStep error for cyc = "
                 << cyc << " - scaledMass = " << scaledMass
                 << "   massDissolved = " << massDissolved << endl;
            cout << "   midx/phName/scaledMassIni_[midx] : " << midx << " / "
                 << phaseKineticModel_[midx]->getName() << " / "
                 << scaledMassIni_[midx] << endl;
            cout << endl << "end program" << endl;
            exit(0);
          }

          if (doNotModif) {
            // cout << "    KineticController::calculateKineticStep - 1 - cyc = " << cyc
            //      << " :  doNotModif = true <-> doNotModif = " << doNotModif
            //      << "   midx = " << midx << "   phName = " << phaseKineticModel_[midx]->getName()
            //      << "   DCId = " << DCId << "   scaledMass = " << scaledMass
            //      << "   massDissolved = " << massDissolved
            //      << "   keepDCLowerLimit = " << chemSys_->getKeepDCLowerLimit(DCId)
            //      << "   DCMoles_ = " << DCMoles_[DCId]
            //      << "   DCMoles_cs = " << chemSys_->getDCMoles(DCId)
            //      << endl;
            continue;
          // } else {
          //   cout << "    KineticController::calculateKineticStep - 2 - cyc = " << cyc
          //        << " :  doNotModif = false <-> doNotModif = " << doNotModif
          //        << "   midx = " << midx << "   phName = " << phaseKineticModel_[midx]->getName()
          //        << "   DCId = " << DCId << "   scaledMass = " << scaledMass
          //        << "   massDissolved = " << massDissolved
          //        << "   keepDCLowerLimit = " << chemSys_->getKeepDCLowerLimit(DCId)
          //        << "   DCMoles_ = " << DCMoles_[DCId]
          //        << "   DCMoles_cs = " << chemSys_->getDCMoles(DCId)
          //        << endl;
          }

          chemSys_->updateMicroPhaseMasses(phaseDissolvedId[midx], scaledMass, 0);

          if (verbose_) {
            cout << "New scaled mass = "
                 << chemSys_->getMicroPhaseMass(phaseDissolvedId[midx])
                 << " and new volume = "
                 << chemSys_->getMicroPhaseVolume(phaseDissolvedId[midx])
                 << endl;
            cout.flush();
          }

          /// @todo Allow any other component to be an impurity, not just the
          /// ones prescribed here
          totMassImpurity = 0;
          keepNumDCMoles = 0;
          numDCMolesDissolved = 0;

          massImpurity =
              massDissolved * chemSys_->getK2o(phaseDissolvedId[midx]);
          totMassImpurity += massImpurity;
          dcmoles = massImpurity / chemSys_->getDCMolarMass("K2O");
          DCMoles_[impurityDCID_[0]] += dcmoles;
          impurity_K2O_[midx] = dcmoles;

          massImpurity =
              massDissolved * chemSys_->getNa2o(phaseDissolvedId[midx]);
          totMassImpurity += massImpurity;
          dcmoles = massImpurity / chemSys_->getDCMolarMass("Na2O");
          DCMoles_[impurityDCID_[1]] += dcmoles;
          impurity_Na2O_[midx] = dcmoles;

          massImpurity =
              massDissolved * chemSys_->getMgo(phaseDissolvedId[midx]);
          totMassImpurity += massImpurity;
          dcmoles = massImpurity / chemSys_->getDCMolarMass("Per");
          DCMoles_[impurityDCID_[2]] += dcmoles;
          impurity_Per_[midx] = dcmoles;

          massImpurity =
              massDissolved * chemSys_->getSo3(phaseDissolvedId[midx]);
          totMassImpurity += massImpurity;
          dcmoles = massImpurity / chemSys_->getDCMolarMass("SO3");
          DCMoles_[impurityDCID_[3]] += dcmoles;
          impurity_SO3_[midx] = dcmoles;

          if (scaledMass > 0) {
            numDCMolesDissolved = (massDissolved - totMassImpurity) /
                                   chemSys_->getDCMolarMass(DCId);
            keepNumDCMoles = DCMoles_[DCId] - numDCMolesDissolved;

            if (keepNumDCMoles < 0) {
              cout << endl
                   << "KineticController::calculateKineticStep error for cyc = "
                   << cyc << " : keepNumDCMoles < 0  !!!" << endl;
              cout << "midx/DCId/DCMoles_/numDCMolesDissolved/keepNumDCMoles : "
                   << midx << " / " << DCId << " / " << DCMoles_[DCId] << " / "
                   << numDCMolesDissolved << " / " << keepNumDCMoles << endl;
              cout << "scaledMass/massDissolved/totMassImpurity/"
                      "massDissolved - totMassImpurity : "
                   << scaledMass << " / " << massDissolved << " / "
                   << totMassImpurity << " / " << massDissolved - totMassImpurity
                   << endl;
              cout << endl << "end program" << endl;
              exit(0);
            }
          } else if (scaledMass < 0) {
            cout << endl
                 << "KineticController::calculateKineticStep error for cyc = "
                 << cyc << " : scaledMass < 0  i.e. scaledMass = " << scaledMass
                 << "   massDissolved = " << massDissolved << endl;
            cout << "   midx/phName/scaledMassIni_[midx] : " << midx << " / "
                 << phaseKineticModel_[midx]->getName() << " / "
                 << scaledMassIni_[midx] << endl;
            cout << endl << "end program" << endl;
            exit(0);
          } else {
            keepNumDCMoles = 0.0;
          }

          // cout << "   :   DCMoles_/numDCMolesDissolved/keepNumDCMoles : "
          //      <<  DCMoles_[DCId] << " / " << numDCMolesDissolved << " / "
          //      << keepNumDCMoles << endl;

          chemSys_->setDCLowerLimit(DCId, keepNumDCMoles);
          chemSys_->setDCUpperLimit(DCId, keepNumDCMoles);

          if (verbose_) {
            cout << "    calculateKineticStep - "
                    "midx/DCId/DCMoles_/numDCMolesDissolved/keepNumDCMoles : "
                 << midx << " / " << DCId << " / " << DCMoles_[DCId] << " / "
                 << numDCMolesDissolved << " / " << keepNumDCMoles << endl;
            cout << "    calculateKineticStep - scaledMass/massDissolved/"
                    "totMassImpurity/massDissolved - totMassImpurity : "
                 << scaledMass << " / " << massDissolved << " / "
                 << totMassImpurity << " / " << massDissolved - totMassImpurity
                 << endl;
          }
        }
      }

      if (verbose_ && doTweak) {
        cout << endl
             << "  KineticController::calculateKineticStep "
                "- tweak after for cyc = "
             << cyc << endl;
      }

    } catch (EOBException eex) {
      eex.printException();
      exit(1);
    } catch (DataException dex) {
      dex.printException();
      exit(1);
    } catch (FloatException fex) {
      fex.printException();
      exit(1);
    } catch (out_of_range &oor) {
      EOBException ex("KineticController", "calculateKineticStep", oor.what(),
                      0, 0);
      ex.printException();
      exit(1);
    }
  } else {
    cout << endl
         << "     KineticController::calculateKineticStep : time >= "
            "beginAttackTime_ -> "
         << time << " >= " << beginAttackTime_ << " (hyd_time = " << hyd_time
         << ")" << endl;
    cout << "     KineticController::calculateKineticStep 0 : count_[VOIDID] = "
         << lattice_->getCount()[VOIDID] << "   &   count_[ELECTROLYTEID] = "
         << lattice_->getCount()[ELECTROLYTEID]
         << "  =>  waterMoles = " << DCMoles_[waterDCId_] << endl;

    double waterAddMoles = lattice_->fillAllPorosity(cyc);
    DCMoles_[waterDCId_] += waterAddMoles;

    cout << "     KineticController::calculateKineticStep 1 : count_[VOIDID] = "
         << lattice_->getCount()[VOIDID] << "   &   count_[ELECTROLYTEID] = "
         << lattice_->getCount()[ELECTROLYTEID]
         << "  =>  waterMoles = " << DCMoles_[waterDCId_] << endl;
  }

  for (i = 0; i < DCNum_; i++) {
    // cout << " " << i << "\t" << DCName_[i] << ": " << DCMoles_[i] << " mol"
    // << endl;
    chemSys_->setDCMoles(i, DCMoles_[i]);
    // cout << "          " << DCName_[i] << ": " << chemSys_->getDCMoles(i) <<
    // " mol" << endl;
  }

  return;
}

void KineticController::updateKineticStep(int cyc, int pId, double scaledMass,
                                          double timestep) {
  string modelName;
  double totMassImpurity, massImpurity;
  double keepNumDCMoles;
  // int phaseDissolvedId;
  double numDCMolesDissolved, massDissolved;

  double hyd_time = hydTimeIni_ + timestep;

  int midx;
  int DCId = -1;
  for (midx = 0; midx < pKMsize_; ++midx) {
    // phaseDissolvedId = phaseKineticModel_[midx]->getMicroPhaseId();
    // if (pId == phaseDissolvedId) {
    if (pId == phaseKineticModel_[midx]->getMicroPhaseId()) {
      DCId = phaseKineticModel_[midx]->getDCId();
      break;
    }
  }
  if (DCId == -1) {
    cout << endl
         << "  KineticController::updateKineticStep - error for cyc = " << cyc
         << " & pId = " << pId << "  =>  DCId = " << DCId << " !!!" << endl;
    cout << "    scaledMass = " << scaledMass << endl;
    cout << endl << "  >>> program stop <<<" << endl;
    exit(0);

  } else {
    chemSys_->setMicroPhaseMass(pId, scaledMassIni_[midx]);
    modelName = phaseKineticModel_[midx]
                    ->getModelName(); // updateKineticStep(scaledMass
                                      // , massDissolved, timestep);
    cout << "  KineticController::updateKineticStep - for cyc = " << cyc
         << " & phaseId = " << pId << " ["
         << phaseKineticModel_[midx]->getName()
         << " / DCId:" << chemSys_->getMicroPhaseDCMembers(pId, 0) << "]"
         << endl;
    cout << "    midx = " << midx << "   modelName : " << modelName
         << "   scaledMassIni[midx] = " << scaledMassIni_[midx]
         << "   scaledMass = " << scaledMass << endl;

    // DCMoles_[DCId] = DCMolesIni_[DCId];
    DCMoles_[impurityDCID_[0]] -= impurity_K2O_[midx];
    DCMoles_[impurityDCID_[1]] -= impurity_Na2O_[midx];
    DCMoles_[impurityDCID_[2]] -= impurity_Per_[midx];
    DCMoles_[impurityDCID_[3]] -= impurity_SO3_[midx];

    /// for this kinetic model

    massDissolved = scaledMassIni_[midx] - scaledMass;

    // chemSys_->setMicroPhaseMass(phaseDissolvedId, scaledMass);
    // chemSys_->setMicroPhaseMassDissolved(phaseDissolvedId, massDissolved);
    chemSys_->updateMicroPhaseMasses(pId, scaledMass, 1);

    totMassImpurity = 0;
    keepNumDCMoles = 0;
    numDCMolesDissolved = 0;

    double dcmoles;
    massImpurity = massDissolved * chemSys_->getK2o(pId);
    totMassImpurity += massImpurity;
    dcmoles = massImpurity / chemSys_->getDCMolarMass("K2O");
    DCMoles_[impurityDCID_[0]] += dcmoles;
    impurity_K2O_[midx] = dcmoles;

    massImpurity = massDissolved * chemSys_->getNa2o(pId);
    totMassImpurity += massImpurity;
    dcmoles = massImpurity / chemSys_->getDCMolarMass("Na2O");
    DCMoles_[impurityDCID_[1]] += dcmoles;
    impurity_Na2O_[midx] = dcmoles;

    massImpurity = massDissolved * chemSys_->getMgo(pId);
    totMassImpurity += massImpurity;
    dcmoles = massImpurity / chemSys_->getDCMolarMass("Per");
    DCMoles_[impurityDCID_[2]] += dcmoles;
    impurity_Per_[midx] = dcmoles;

    massImpurity = massDissolved * chemSys_->getSo3(pId);
    totMassImpurity += massImpurity;
    dcmoles = massImpurity / chemSys_->getDCMolarMass("SO3");
    DCMoles_[impurityDCID_[3]] += dcmoles;
    impurity_SO3_[midx] = dcmoles;

    numDCMolesDissolved =
        (massDissolved - totMassImpurity) / chemSys_->getDCMolarMass(DCId);
    keepNumDCMoles = DCMoles_[DCId] - numDCMolesDissolved;

    chemSys_->setDCLowerLimit(DCId, keepNumDCMoles);
    chemSys_->setDCUpperLimit(DCId, keepNumDCMoles);

    cout << "      massDissolved/totMassImpurity/massDissolved - "
            "totMassImpurity "
            ": "
         << massDissolved << " / " << totMassImpurity << " / "
         << massDissolved - totMassImpurity << endl;
    cout << "      DCMoles_/numDCMolesDissolved/keepNumDCMoles : "
         << DCMoles_[DCId] << " / " << numDCMolesDissolved << " / "
         << keepNumDCMoles << endl;
    // *****************

    if (hyd_time > beginAttackTime_) {
      cout << endl
           << "     KineticController::updateKineticStep : time >= "
              "beginAttackTime_ -> "
           << hyd_time << " >= " << beginAttackTime_ << " (hyd_time)" << endl;
      cout << "     KineticController::updateKineticStep 0 : count_[VOIDID] "
              "= "
           << lattice_->getCount()[VOIDID] << "   &   count_[ELECTROLYTEID] = "
           << lattice_->getCount()[ELECTROLYTEID]
           << "  =>  waterMoles = " << DCMoles_[waterDCId_] << endl;

      double waterAddMoles = lattice_->fillAllPorosity(cyc);
      DCMoles_[waterDCId_] += waterAddMoles;

      // if (waterAddMoles > 0)
      //   cout << "     KineticController::updateKineticStep : check if OK!"
      //        << endl;

      cout << "     KineticController::updateKineticStep 1 : count_[VOIDID] "
              "= "
           << lattice_->getCount()[VOIDID] << "   &   count_[ELECTROLYTEID] = "
           << lattice_->getCount()[ELECTROLYTEID]
           << "  =>  waterMoles = " << DCMoles_[waterDCId_] << endl;
    }

    for (int i = 0; i < DCNum_; i++) {
      // cout << " " << i << "\t" << DCName_[i] << ": " << DCMoles_[i] << " mol"
      // << endl;
      chemSys_->setDCMoles(i, DCMoles_[i]);
      // cout << "          " << DCName_[i] << ": " << chemSys_->getDCMoles(i)
      // << " mol" << endl;
    }
  }
}
