/**
@file  KineticController.cc
@brief Method definitions for the KineticController class.

*/
#include "KineticController.h"

#include <algorithm>

#include "SaturatingRateModel.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

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

  // JMAK per-phase state — parallel to name_, scaledMass_, etc.
  // Populated alongside them in parseKineticData.
  jmakEnabled_.clear();
  jmakParams_.clear();
  jmakGlobals_.clear();
  jmakGenerations_.clear();
  jmakSeedAccumulator_.clear();
  jmakVoxelsInLattice_.clear();
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

  // Initialize SI data tracking
  initialMicroPhaseSI_.clear();
  hasSIData_ = false;

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
  std::clog << "KineticController::KineticController Constructor" << endl;
  std::clog.flush();
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
        std::clog << "KineticModel data file is a JSON file" << endl;
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
    std::clog << "KineticController::KineticController Finished reading "
                 "simparams.json "
              << endl;
    int size = microPhaseId_.size();
    for (int i = 0; i < size; ++i) {
      microPhaseId = microPhaseId_[i];
      if (isKinetic_[i]) {
        std::clog << "KineticController::KineticController kinetic phase "
                  << microPhaseId << endl;
        std::clog << "KineticController::KineticController     name = "
                  << chemSys_->getMicroPhaseName(microPhaseId) << endl;
      }
    }
    std::clog.flush();
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
  std::clog << endl << "KineticController::KineticController(...) :" << endl;
  std::clog << "  - only these phases (controlled by the P-K model) "
               "contribute to initScaledCementMass_ & scaledCementMass_ :"
            << endl;

  for (int i = 0; i < pKMsize_; i++) {
    modelName = phaseKineticModel_[i]->getModelName();
    // std::clog << endl << "    modelName = " << modelName << endl;
    if (modelName == "ParrotKillohModel") {
      phID = phaseKineticModel_[i]->getMicroPhaseId();
      initScaledCementMass_ += chemSys_->getMicroPhaseMass(phID);
      std::clog << "      microPhaseID/microPhaseName/microPhaseMass : "
                << setw(3) << right << phID << " / " << setw(15) << left
                << phaseKineticModel_[i]->getName() << " / "
                << chemSys_->getMicroPhaseMass(phID) << " g" << endl;
      chemSys_->setIsParrotKilloh(phID);
    }
  }
  std::clog << endl
            << "      initScaledCementMass_ = " << initScaledCementMass_
            << " g (same value for scaledCementMass_)" << endl;
  chemSys_->setInitScaledCementMass(initScaledCementMass_);

  hydTimeIni_ = 0;

  // Initialize SI data tracking (will be set by setInitialSaturationIndices
  // if pre-equilibration is performed)
  initialMicroPhaseSI_.clear();
  hasSIData_ = false;

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

    // JMAK per-phase state — allocate parallel entries for every parsed
    // phase. A phase is jmak-enabled iff kineticData.jmak was set AND
    // kineticData.nucleation was set (JMAK is a growth extension of CNT
    // — makes no sense without a nucleation source). See
    // KineticController.h "JMAK per-voxel growth state" comment for the
    // per-cycle semantics.
    const bool jmakOn = kineticData.jmak.has_value() &&
                        kineticData.nucleation.has_value();
    jmakEnabled_.push_back(jmakOn);
    jmakParams_.push_back(jmakOn ? *kineticData.jmak
                                 : JMAKParameters{4.0, 4.0 * Pi / 3.0});
    jmakGlobals_.push_back(jmak::GlobalMoments{0.0, 0.0, 0.0, 0.0, 0.0});
    jmakGenerations_.push_back(std::vector<JMAKGeneration>{});
    jmakSeedAccumulator_.push_back(0.0);
    jmakVoxelsInLattice_.push_back(0);
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
    } else if (kineticData.type == SaturatingRateType) {
      typefound = true;
      try {
        parseKineticDataForSaturatingRate(p, kineticData);
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
    std::clog << "--->Parsing PK data for " << kineticData.name << endl;
    std::clog.flush();
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

void KineticController::parseKineticDataForStandard(
    const json::iterator p, struct KineticData &kineticData) {

  if (verbose_) {
    std::clog << "--->Parsing standard kinetic data for " << kineticData.name
              << endl;
    std::clog.flush();
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
    std::clog << "WARNING: For " << kineticData.name
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
    std::clog << "WARNING: For " << kineticData.name
              << " dissolvedUnits not found; setting it to 1" << endl;
  }

  // Exponent on  the saturation index in the rate equation
  pp = p.value().find("siexp");
  if (pp != p.value().end()) {
    kineticData.siexp = pp.value();
  } else {
    kineticData.sio2 = 1.0;
    std::clog << "WARNING: For " << kineticData.name
              << " sio2 not found; setting it to 1" << endl;
  }

  // Exponent on  the driving force term in the rate equation
  pp = p.value().find("dfexp");
  if (pp != p.value().end()) {
    kineticData.dfexp = pp.value();
  } else {
    kineticData.dfexp = 1.0;
    std::clog << "WARNING: For " << kineticData.name
              << " dfexp not found; setting it to 1" << endl;
  }

  // Loss on ignition of the material
  pp = p.value().find("loi");
  if (pp != p.value().end()) {
    kineticData.loi = pp.value();
  } else {
    kineticData.loi = 1.0e-6;
    std::clog << "WARNING: For " << kineticData.name
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

  // Optional CNT nucleation block; leaves kineticData.nucleation empty if absent.
  parseNucleationBlock(p, kineticData);

  return;
}

void KineticController::parseKineticDataForPozzolanic(
    const json::iterator p, struct KineticData &kineticData) {

  if (verbose_) {
    std::clog << "--->Parsing pozzolanic data for " << kineticData.name << endl;
    std::clog.flush();
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
    std::clog << "WARNING: For " << kineticData.name
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
    std::clog << "WARNING: For " << kineticData.name
              << " dissolvedUnits not found; setting it to 1" << endl;
  }

  // Exponent on the saturation index in the rate equation
  pp = p.value().find("siexp");
  if (pp != p.value().end()) {
    kineticData.siexp = pp.value();
  } else {
    kineticData.siexp = 1.0;
    std::clog << "WARNING: For " << kineticData.name
              << " siexp not found; setting it to 1" << endl;
  }

  // Exponent on the driving force term in the rate equation
  pp = p.value().find("dfexp");
  if (pp != p.value().end()) {
    kineticData.dfexp = pp.value();
  } else {
    kineticData.dfexp = 1.0;
    std::clog << "WARNING: For " << kineticData.name
              << " dfexp not found; setting it to 1" << endl;
  }

  // Exponent on the degree of reaction term in the diffusion rate equation
  pp = p.value().find("dorexp");
  if (pp != p.value().end()) {
    kineticData.dorexp = pp.value();
  } else {
    kineticData.dorexp = 1.0;
    std::clog << "WARNING: For " << kineticData.name
              << " dorexp not found; setting it to 1" << endl;
  }

  // Exponent on the hydroxy ion activity in the rate equation
  pp = p.value().find("ohexp");
  if (pp != p.value().end()) {
    kineticData.ohexp = pp.value();
  } else {
    kineticData.ohexp = 1.0;
    std::clog << "WARNING: For " << kineticData.name
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

  // Optional CNT nucleation block; leaves kineticData.nucleation empty if absent.
  parseNucleationBlock(p, kineticData);

  return;
}

// Parses an optional "nucleation" sub-block from a phase's "kinetic_data"
// object. Absent block = CNT disabled for this phase.
//
// Expected JSON (all three fields required when the block is present):
//   "nucleation": {
//     "gamma": {"value": 0.044,  "range": [0.030, 0.070], "provenance": "..."},
//     "theta": {"value": 180,    "range": [1, 180],       "provenance": "..."},
//     "A0":    {"value": 1.0e30, "range": [1.0e28, 1.0e32], "provenance": "..."}
//   }
//
// Only the "value" field is read into runtime memory. "range" and
// "provenance" are curation metadata for humans and sweep tooling; they
// stay in the JSON file as documentation but never enter the simulation.
//
// theta is validated to [1, 180]. theta = 180 designates homogeneous
// placement (uniform-random over electrolyte voxels via
// Lattice::nucleatePhaseRnd). theta in [1, 179] is reserved for future
// heterogeneous placement on substrate voxels.
void KineticController::parseNucleationBlock(const json::iterator p,
                                             struct KineticData &kineticData) {
  json::iterator nucIt = p.value().find("nucleation");
  if (nucIt == p.value().end()) {
    return;  // Absent block = CNT disabled; kineticData.nucleation stays empty
  }

  // Each parameter is a {value, range, provenance} object; runtime reads only
  // the value field. range and provenance are user-facing metadata (documented
  // in the JSON for what-if / sweep tooling) and are ignored here.
  auto readValue = [&](const std::string &key) -> json::iterator {
    auto it = nucIt.value().find(key);
    if (it == nucIt.value().end()) {
      throw DataException("KineticController", "parseNucleationBlock",
                          key + " not found in nucleation block");
    }
    auto valIt = it.value().find("value");
    if (valIt == it.value().end()) {
      throw DataException("KineticController", "parseNucleationBlock",
                          key + ".value not found");
    }
    return valIt;
  };

  NucleationParameters np;
  np.gamma = readValue("gamma").value().get<double>();
  np.theta_deg = readValue("theta").value().get<int>();
  np.A0 = readValue("A0").value().get<double>();

  if (np.theta_deg < 1 || np.theta_deg > 180) {
    throw DataException("KineticController", "parseNucleationBlock",
                        "theta must be integer in [1, 180]");
  }

  kineticData.nucleation = np;

  if (verbose_) {
    std::clog << "--->Parsed nucleation block for " << kineticData.name
              << ": gamma=" << np.gamma << " J/m^2, theta=" << np.theta_deg
              << " deg, A0=" << np.A0 << " /(m^3 s)" << endl;
    std::clog.flush();
  }

  // ------ Optional JMAK sub-block ------
  // If present, the phase uses per-voxel JMAK growth instead of the
  // classical "1 voxel = 1 nucleation event" placement path. See
  // docs/POST_ALPHA_TODOS.md "CNT growth model needs JMAK-per-voxel"
  // for the physical rationale.
  //
  // Expected JSON (both fields optional; defaults apply if absent):
  //   "jmak": {
  //     "n":     {"value": 4.0, "range": [2.5, 4.0], "provenance": "..."},
  //     "alpha": {"value": 4.18879, "range": [1.0, 10.0], "provenance": "..."}
  //   }
  //
  // n = Avrami exponent (4 = 3D continuous nucleation; range [2.5, 4]
  //     defensible for CNT-nucleated phases in a 1 um^3 voxel with
  //     non-ideal habitus).
  // alpha = morphology coefficient (4*pi/3 for 3D isotropic spheres).
  {
    json::iterator jmakIt = nucIt.value().find("jmak");
    if (jmakIt != nucIt.value().end()) {
      JMAKParameters jp;
      // Defaults if a field is absent.
      jp.n = 4.0;
      jp.alpha = 4.0 * Pi / 3.0;

      auto nField = jmakIt.value().find("n");
      if (nField != jmakIt.value().end()) {
        auto v = nField.value().find("value");
        if (v == nField.value().end()) {
          throw DataException("KineticController", "parseNucleationBlock",
                              "jmak.n.value not found");
        }
        jp.n = v.value().get<double>();
      }
      auto alphaField = jmakIt.value().find("alpha");
      if (alphaField != jmakIt.value().end()) {
        auto v = alphaField.value().find("value");
        if (v == alphaField.value().end()) {
          throw DataException("KineticController", "parseNucleationBlock",
                              "jmak.alpha.value not found");
        }
        jp.alpha = v.value().get<double>();
      }

      if (jp.n < 2.5 || jp.n > 4.0) {
        throw DataException("KineticController", "parseNucleationBlock",
                            "jmak.n must be in [2.5, 4.0]");
      }
      if (jp.alpha <= 0.0) {
        throw DataException("KineticController", "parseNucleationBlock",
                            "jmak.alpha must be positive");
      }

      kineticData.jmak = jp;

      if (verbose_) {
        std::clog << "--->Parsed jmak sub-block for " << kineticData.name
                  << ": n=" << jp.n << ", alpha=" << jp.alpha << endl;
        std::clog.flush();
      }
    }
  }
}

void KineticController::updateJMAKPhase(int midx, double timestep, int cyc) {
  // Precondition: caller has verified useNucleationKinetics_ AND
  // jmakEnabled_[midx] AND phaseKineticModel_[midx]->hasNucleation().
  //
  // Per JMAK-enabled phase, per cycle:
  //   1. Get S, J, G from the rate-law subclass.
  //   2. jmak::advanceMoments updates M_0..M_3 and G_acc.
  //   3. Seed new generations: expected new-seed voxels this cycle
  //      = J * V_unseeded_electrolyte * dt_seconds. Accumulate
  //      fractionally; when the accumulator crosses 1.0 create a new
  //      Generation with floor(acc) voxels and a moment snapshot.
  //   4. Iterate active generations, compute Y_g and X_g.
  //   5. Total transformed = sum_g N_g * X_g. Delta against
  //      jmakVoxelsInLattice_[midx] tells us how many new voxels to
  //      place in the lattice this cycle.
  //   6. Place via Lattice::nucleatePhaseRnd, update DCMoles_,
  //      KineticModel::scaledMass_, ChemicalSystem::microPhaseMass_
  //      using the same per-100g scaling as the classical CNT
  //      placement block.
  //   7. Prune generations that have reached X_g >= 0.999 to keep
  //      the vector length bounded.

  KineticModel *km = phaseKineticModel_[midx];
  const int microPhId = km->getMicroPhaseId();
  const int dcId = km->getDCId();

  // ---- 1. Query rate law for S, J, G ----
  const double S = chemSys_->getMicroPhaseSI(microPhId);
  const double J = km->getNucleationRate(S);          // [1/(m^3 s)]
  const double G = km->getGrowthVelocity(S);          // [m/s]
  const double dt_seconds = timestep * 3600.0;

  // ---- 2. Advance moments (idempotent when J = 0 or G = 0) ----
  jmak::advanceMoments(jmakGlobals_[midx], J, G, dt_seconds);

  // If nothing to do (both J = 0 AND no existing generations), skip the
  // rest. J may be > 0 with no generations yet (build-up phase) or
  // J = 0 with generations still growing (Ostwald tail).
  if (J <= 0.0 && jmakGenerations_[midx].empty()) return;

  const double V_voxel = lattice_->getVolumePerVoxel();          // [m^3]

  // ---- 3. Seed new generation? ----
  //
  // "Unseeded" volume = current electrolyte volume, minus the volume
  // already assigned to existing generations. In the JMAK model each
  // seeded voxel is a whole voxel worth of "reserved" electrolyte
  // (contains fractional-to-full Portlandite over time), so unseeded
  // voxel count = (electrolyte voxels in lattice) - (voxels reserved
  // by generations that are not yet fully transformed).
  int reservedInFlight = 0;
  for (const auto &g : jmakGenerations_[midx]) {
    reservedInFlight += g.Nvoxels;
  }
  // Subtract voxels already placed in the lattice — those are no
  // longer electrolyte, they're solid Portlandite. Reserving them
  // twice would inflate the unseeded count.
  reservedInFlight -= jmakVoxelsInLattice_[midx];
  if (reservedInFlight < 0) reservedInFlight = 0;

  const int nElectrolyte = lattice_->getCount(ELECTROLYTEID);
  int nUnseededVoxels = nElectrolyte - reservedInFlight;
  if (nUnseededVoxels < 0) nUnseededVoxels = 0;

  const double V_unseeded = static_cast<double>(nUnseededVoxels) * V_voxel;
  // Expected new-seed voxels this cycle (small-lambda Poisson approx).
  const double expectedNewSeeds = J * V_unseeded * dt_seconds;
  jmakSeedAccumulator_[midx] += expectedNewSeeds;

  if (jmakSeedAccumulator_[midx] >= 1.0) {
    const int nSeedInt =
        static_cast<int>(std::floor(jmakSeedAccumulator_[midx]));
    JMAKGeneration g;
    g.Nvoxels = nSeedInt;
    g.seed = jmak::snapshotSeed(jmakGlobals_[midx]);
    jmakGenerations_[midx].push_back(g);
    jmakSeedAccumulator_[midx] -= static_cast<double>(nSeedInt);
  }

  // ---- 4-5. Evaluate all generations; compute total transformed ----
  double transformedVoxelsCumulative = 0.0;
  for (auto &g : jmakGenerations_[midx]) {
    const double Y_over_V = jmak::extendedVolumePerVoxel(
        g.seed, jmakGlobals_[midx], jmakParams_[midx]);
    const double X = jmak::transformedFraction(Y_over_V);
    transformedVoxelsCumulative += static_cast<double>(g.Nvoxels) * X;
  }
  const int transformedIntTotal =
      static_cast<int>(std::floor(transformedVoxelsCumulative));
  const int deltaVoxels = transformedIntTotal - jmakVoxelsInLattice_[midx];

  // ---- 6. Place delta voxels in the lattice ----
  if (deltaVoxels > 0) {
    const int nPlaced = lattice_->nucleatePhaseRnd(microPhId, deltaVoxels);
    if (nPlaced > 0) {
      // Same per-100g scaling as the classical CNT placement block.
      // See KineticController.cc block comment above the classical
      // CNT placement for the derivation.
      const double molarMass = chemSys_->getDCMolarMass(dcId);
      const double vMolar = chemSys_->getDCMolarVolume(dcId);
      const double vfracPlaced =
          static_cast<double>(nPlaced) /
          static_cast<double>(lattice_->getNumSites());
      const double microPhaseMassPlacedPerCC =
          vfracPlaced * molarMass / vMolar / 1.0e6;          // g/cm^3
      const double placedMass =
          microPhaseMassPlacedPerCC * 100.0 /
          lattice_->getInitSolidMass();                       // g/100g solid
      const double moles = placedMass / molarMass;            // mol/100g solid

      DCMoles_[dcId] += moles;
      chemSys_->setDCLowerLimit(dcId, DCMoles_[dcId]);
      chemSys_->setDCUpperLimit(dcId, DCMoles_[dcId]);

      const double newScaledMass = km->getScaledMass() + placedMass;
      km->setScaledMass(newScaledMass);
      chemSys_->updateMicroPhaseMasses(microPhId, newScaledMass, 1);
      scaledMassIni_[midx] = newScaledMass;

      jmakVoxelsInLattice_[midx] += nPlaced;

      if (verbose_) {
        std::clog << "  JMAK: " << km->getName()
                  << " cyc=" << cyc
                  << " transformedTotal=" << transformedVoxelsCumulative
                  << " voxelsPlacedThisCycle=" << nPlaced
                  << " voxelsInLatticeCumul=" << jmakVoxelsInLattice_[midx]
                  << " (+" << placedMass << " g/100g)" << endl;
      }
    }
  }

  // ---- 7. Prune completed generations ----
  // Generations with X ~ 1 have contributed everything they can; drop
  // them so the vector stays bounded (avoids O(cycles) growth over
  // long runs). Threshold 0.999 is generous — we still keep some
  // late-life generations to avoid a discontinuity from pruning.
  auto &gens = jmakGenerations_[midx];
  gens.erase(std::remove_if(gens.begin(), gens.end(),
                            [&](const JMAKGeneration &g) {
                              const double Y = jmak::extendedVolumePerVoxel(
                                  g.seed, jmakGlobals_[midx],
                                  jmakParams_[midx]);
                              return jmak::transformedFraction(Y) >= 0.999;
                            }),
             gens.end());
}

void KineticController::parseSaturatingRateSubBlock(
    const json::iterator p, const std::string &blockName,
    std::optional<SaturatingRateParameters> &target) {
  json::iterator blockIt = p.value().find(blockName);
  if (blockIt == p.value().end()) {
    return;  // Absent block leaves target unchanged.
  }

  auto readValue = [&](const std::string &key) -> double {
    auto it = blockIt.value().find(key);
    if (it == blockIt.value().end()) {
      throw DataException("KineticController", "parseSaturatingRateSubBlock",
                          blockName + "." + key + " not found");
    }
    return it.value().get<double>();
  };

  SaturatingRateParameters params;
  params.rateConstant = readValue("rateConstant");
  params.B = readValue("B");
  params.n = readValue("n");
  target = params;
}

// Parser for a SaturatingRate-type phase (Bullard 2015 / Han 2025 Eq. 7).
//
// Expected JSON (dissolution required; precipitation optional):
//   "kinetic_data": {
//     "type": "SaturatingRate",
//     "surfaceAreaMultiplier": <double>,      (optional, default 1.0)
//     "dissolvedUnits":       <double>,        (required)
//     "activationEnergy":     <double>,        (required, J/mol)
//     "dissolution": {"rateConstant": ..., "B": ..., "n": ...},
//     "precipitation": {"rateConstant": ..., "B": ..., "n": ...},
//     "nucleation": {...optional CNT block...}
//   }
//
// Missing precipitation block => model uses dissolution parameters as
// symmetric fallback with a one-time log line noting the assumption.
void KineticController::parseKineticDataForSaturatingRate(
    const json::iterator p, struct KineticData &kineticData) {

  if (verbose_) {
    std::clog << "--->Parsing SaturatingRate kinetic data for "
              << kineticData.name << endl;
    std::clog.flush();
  }

  // surfaceAreaMultiplier (optional; default 1.0)
  json::iterator pp = p.value().find("surfaceAreaMultiplier");
  if (pp != p.value().end()) {
    kineticData.surfaceAreaMultiplier = pp.value();
  } else {
    kineticData.surfaceAreaMultiplier = 1.0;
  }

  // dissolvedUnits (required)
  pp = p.value().find("dissolvedUnits");
  if (pp != p.value().end()) {
    kineticData.dissolvedUnits = pp.value();
  } else {
    throw DataException("KineticController",
                        "parseKineticDataForSaturatingRate",
                        "dissolvedUnits not found");
  }

  // activationEnergy (required)
  pp = p.value().find("activationEnergy");
  if (pp != p.value().end()) {
    kineticData.activationEnergy = pp.value();
  } else {
    throw DataException("KineticController",
                        "parseKineticDataForSaturatingRate",
                        "activationEnergy not found");
  }

  // dissolution block (required)
  parseSaturatingRateSubBlock(p, "dissolution",
                              kineticData.saturatingDissolution);
  if (!kineticData.saturatingDissolution.has_value()) {
    throw DataException("KineticController",
                        "parseKineticDataForSaturatingRate",
                        "dissolution block not found; required for "
                        "SaturatingRate-type phase");
  }

  // precipitation block (optional; empty triggers the symmetric fallback
  // inside the model's rate calculation)
  parseSaturatingRateSubBlock(p, "precipitation",
                              kineticData.saturatingPrecipitation);

  // Optional CNT nucleation block; leaves kineticData.nucleation empty if absent.
  parseNucleationBlock(p, kineticData);
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
        std::clog
            << "KineticController::getPhaseMasses reads solid micphase mass of "
            << chemSys_->getMicroPhaseName(microPhaseId) << " as "
            << initScaledMass_[i] << endl;
        std::clog.flush();
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
    // Read remaining pozzolanic model parameters
    km = new StandardKineticModel(chemSys_, lattice_, kineticData, verbose_,
                                  warning_);
  } else if (kineticData.type == PozzolanicType) {
    // Read remaining pozzolanic model parameters
    km = new PozzolanicModel(chemSys_, lattice_, kineticData, verbose_,
                             warning_);
  } else if (kineticData.type == SaturatingRateType) {
    // Bullard 2015 / Han 2025 Eq. 7 saturating rate model
    km = new SaturatingRateModel(chemSys_, lattice_, kineticData, verbose_,
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
      std::clog << endl
                << "KineticController::setPozzEffectOnPK for midx = " << midx
                << " (microPhaseId =  "
                << phaseKineticModel_[midx]->getMicroPhaseId()
                << ", microPhaseName = " << phaseKineticModel_[midx]->getName()
                << endl;

      std::clog << "  Ref LOI = " << refloi << endl;
      std::clog << "  LOI     = " << loi << endl;
      std::clog << "  Max LOI = " << maxloi << endl;
      std::clog << "  SiO2     = " << sio2val << endl;
      std::clog << "  Ref SiO2 = " << refsio2val << endl;
      std::clog << "  BET      = " << phaseSA << endl;
      std::clog << "  Ref BET  = " << refPhaseSA << endl;
      std::clog << "  Pozz Effect     = " << pozzeffect << endl;
      std::clog << "  Min Pozz Effect = " << minpozzeffect << endl;
      std::clog.flush();
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
  std::clog << scientific << setprecision(15);
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

  // std::clog << endl << "impurityDCID : " << endl;
  // for(i = 0; i < chemSys_->getNumMicroImpurities(); i++){
  //     std::clog << i << "\t" << impurityDCID[i] << endl; std::clog.flush();
  // }
  // std::clog << endl ;

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
      std::clog << endl
                << "    KineticController::calculateKineticStep - tweak cyc = "
                << cyc << " :  hyd_time = " << hyd_time
                << "   hydTimeIni_ = " << hydTimeIni_
                << "   timestep = " << timestep << endl;

    for (int midx = 0; midx < pKMsize_; ++midx) {
      phaseDissolvedId[midx] = phaseKineticModel_[midx]->getMicroPhaseId();
      chemSys_->setMicroPhaseMass(phaseDissolvedId[midx], scaledMassIni_[midx]);
      if (verbose_) {
        std::clog << "       midx = " << midx
                  << "     scaledMassIni[midx] = " << scaledMassIni_[midx]
                  << "     microPhaseName = "
                  << phaseKineticModel_[midx]->getName() << endl;
      }
    }

    for (i = 0; i < DCNum_; i++) {
      DCMoles_[i] = DCMolesIni_[i];
    }

    lattice_->resetSurfaceArea(surfaceAreaIni_);

  } else {

    // hyd_time = hydTimeIni_ + timestep;
    std::clog << "    KineticController::calculateKineticStep - cyc = " << cyc
              << " :  hyd_time = " << hyd_time
              << "   hydTimeIni_ = " << hydTimeIni_
              << "   timestep = " << timestep << endl;

    for (int midx = 0; midx < pKMsize_; ++midx) {
      phaseDissolvedId[midx] = phaseKineticModel_[midx]->getMicroPhaseId();
      scaledMassIni_[midx] =
          chemSys_->getMicroPhaseMass(phaseDissolvedId[midx]);
      if (verbose_) {
        std::clog << "      midx = " << midx
                  << "     scaledMassIni[midx] = " << scaledMassIni_[midx]
                  << "     microPhaseName = "
                  << phaseKineticModel_[midx]->getName() << endl;
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
      // std::clog << "  KineticController::calculateKineticStep     hyd_time =
      // "
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
      /// be only portland cement clinker components; now is depending on user
      /// decision (arcanite, thenardite, gypsum, bassanite and hemihydrate can
      /// belong to this category)

      if (initScaledCementMass_ > 0) {
        totalDOR = (initScaledCementMass_ - chemSys_->getScaledCementMass()) /
                   initScaledCementMass_;
      } else {
        int numPKMphases = 0;
        for (int midx = 0; midx < pKMsize_; ++midx) {
          if (phaseKineticModel_[midx]->getModelName() == "ParrotKillohModel")
            numPKMphases++;
        }

        // This next block only if there are ParrotKilloh model phases
        if (numPKMphases > 0) {

          std::clog << endl
                    << "   KineticController::calculateKineticStep error - "
                       "initScaledCementMass_ = 0  while numPKMphases = "
                    << numPKMphases << " :" << endl;
          for (int midx = 0; midx < pKMsize_; ++midx) {
            phaseDissolvedId[midx] =
                phaseKineticModel_[midx]->getMicroPhaseId();
            scaledMassIni_[midx] =
                chemSys_->getMicroPhaseMass(phaseDissolvedId[midx]);
            std::clog << "     midx = " << midx
                      << "     scaledMassIni[midx] = " << scaledMassIni_[midx]
                      << "     microPhaseName = "
                      << phaseKineticModel_[midx]->getName() << endl;
          }
          std::clog << endl
                    << "        cyc/doTweak/timesGEMFailed : " << cyc << " / "
                    << doTweak << " / " << chemSys_->getTimesGEMFailed()
                    << endl;
          throw FloatException("KineticController", "calculateKineticStep",
                               "initScaledCementMass_ = 0");
        }
      }

      if (totalDOR < 0) {
        std::clog
            << endl
            << "     KineticController::calculateKineticStep error : totalDOR "
               "< 0"
            << endl;
        std::clog << endl
                  << "        cyc/doTweak/timesGEMFailed : " << cyc << " / "
                  << doTweak << " / " << chemSys_->getTimesGEMFailed() << endl;
        std::clog
            << endl
            << "        initScaledCementMass_/scaledCementMass/totalDOR : "
            << initScaledCementMass_ << " / " << chemSys_->getScaledCementMass()
            << " / " << totalDOR << endl;
        throw DataException("KineticController", "calculateKineticStep",
                            "totalDOR < 0");
      }
      if (!doTweak) {
        std::clog << "    KineticController::calculateKineticStep - cyc = "
                  << cyc << " :  scaledCementMass = "
                  << chemSys_->getScaledCementMass()
                  << "   totalDOR = " << totalDOR << endl;
      }

      //*******

      // ---------- CNT pre-loop lock (guard G4 in docs/CNT_ARCHITECTURE.md) ----------
      // For every phase that has (a) CNT configured, (b) currently zero mass,
      // set both DC limits to 0. This coordinates with three other pieces:
      //   1. Standard/Pozzolanic::calculateKineticStep has an early-return
      //      bypass (guard G3) that skips their rate calculation entirely
      //      for CNT-controlled zero-mass phases, so DCMoles_ is not
      //      incremented by the kinetic model.
      //   2. The CNT placement block at the tail of this same phase loop
      //      raises both DC limits to DCMoles_[dcId] AFTER placing nuclei
      //      (guard G5) — GEMS then respects the new floor and does not
      //      dissolve the placed nuclei back.
      //   3. GEMS at end-of-cycle sees limits = 0 for phases not yet
      //      bootstrapped by CNT and cannot precipitate them of its own
      //      accord. Only CNT can create such a phase from mass = 0.
      // Without this pre-loop lock, GEMS is free to precipitate portlandite
      // (or any CNT-configured phase) at S > 1 the moment it appears, and
      // CNT's role as the sole nucleation mechanism is defeated.
      if (useNucleationKinetics_) {
        for (int i = 0; i < pKMsize_; ++i) {
          if (phaseKineticModel_[i] != nullptr &&
              phaseKineticModel_[i]->hasNucleation() &&
              phaseKineticModel_[i]->getScaledMass() <= 0.0) {
            int dcId = phaseKineticModel_[i]->getDCId();
            chemSys_->setDCLowerLimit(dcId, 0.0);
            chemSys_->setDCUpperLimit(dcId, 0.0);
            if (verbose_) {
              std::clog << "  CNT-lock: "
                        << phaseKineticModel_[i]->getName()
                        << " has nucleation block and zero mass; DC["
                        << dcId << "] locked to 0 until CNT places nuclei"
                        << endl;
            }
          }
        }
      }

      double dcmoles;
      bool runKM = true;

      for (int midx = 0; midx < pKMsize_; ++midx) {
        scaledMass = scaledMassIni_[midx];
        runKM = true;
        if (scaledMass == 0.0 &&
            phaseKineticModel_[midx]->getModelName() == "ParrotKillohModel") {
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

          // std::clog << endl << "    KineticController::calculateKineticStep -
          // 0 - cyc = " << cyc
          //      << " :  doNotModif = false <-> doNotModif = " << doNotModif
          //      << "   midx = " << midx << "   phName = " <<
          //      phaseKineticModel_[midx]->getName()
          //      << "   DCId = " << DCId << "   scaledMass = " << scaledMass
          //      << "   massDissolved = " << massDissolved
          //      << "   keepDCLowerLimit = " <<
          //      chemSys_->getKeepDCLowerLimit(DCId)
          //      << "   DCMoles_ = " << DCMoles_[DCId]
          //      << "   DCMoles_cs = " << chemSys_->getDCMoles(DCId)
          //     << endl;

          phaseKineticModel_[midx]->calculateKineticStep(
              timestep, scaledMass, massDissolved, cyc, totalDOR, doTweak,
              doNotModif);

          /// @note may want to change the condition of next block
          /// because it is possible for the scaled mass of a kinetic phase to
          /// be zero but not negative. We could allow a threshold number
          /// like -1.0e-9 or something like that to handle "nearly zero"

          if (scaledMass < 0.0) {
            std::clog
                << endl
                << "KineticController::calculateKineticStep error for cyc = "
                << cyc << " - scaledMass = " << scaledMass
                << "   massDissolved = " << massDissolved << endl;
            std::clog << "   midx/phName/scaledMassIni_[midx] : " << midx
                      << " / " << phaseKineticModel_[midx]->getName() << " / "
                      << scaledMassIni_[midx] << endl;
            std::clog << endl << "end program" << endl;
            exit(0);
          }

          if (doNotModif) {
            // std::clog << "    KineticController::calculateKineticStep - 1 -
            // cyc = " << cyc
            //      << " :  doNotModif = true <-> doNotModif = " << doNotModif
            //      << "   midx = " << midx << "   phName = " <<
            //      phaseKineticModel_[midx]->getName()
            //      << "   DCId = " << DCId << "   scaledMass = " << scaledMass
            //      << "   massDissolved = " << massDissolved
            //      << "   keepDCLowerLimit = " <<
            //      chemSys_->getKeepDCLowerLimit(DCId)
            //      << "   DCMoles_ = " << DCMoles_[DCId]
            //      << "   DCMoles_cs = " << chemSys_->getDCMoles(DCId)
            //      << endl;
            continue;
            // } else {
            //   std::clog << "    KineticController::calculateKineticStep - 2 -
            //   cyc = " << cyc
            //        << " :  doNotModif = false <-> doNotModif = " <<
            //        doNotModif
            //        << "   midx = " << midx << "   phName = " <<
            //        phaseKineticModel_[midx]->getName()
            //        << "   DCId = " << DCId << "   scaledMass = " <<
            //        scaledMass
            //        << "   massDissolved = " << massDissolved
            //        << "   keepDCLowerLimit = " <<
            //        chemSys_->getKeepDCLowerLimit(DCId)
            //        << "   DCMoles_ = " << DCMoles_[DCId]
            //        << "   DCMoles_cs = " << chemSys_->getDCMoles(DCId)
            //        << endl;
          }

          chemSys_->updateMicroPhaseMasses(phaseDissolvedId[midx], scaledMass,
                                           0);

          if (verbose_) {
            std::clog << "New scaled mass = "
                      << chemSys_->getMicroPhaseMass(phaseDissolvedId[midx])
                      << " and new volume = "
                      << chemSys_->getMicroPhaseVolume(phaseDissolvedId[midx])
                      << endl;
            std::clog.flush();
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
              std::clog
                  << endl
                  << "KineticController::calculateKineticStep error for cyc = "
                  << cyc << " : keepNumDCMoles < 0  !!!" << endl;
              std::clog
                  << "midx/DCId/DCMoles_/numDCMolesDissolved/keepNumDCMoles : "
                  << midx << " / " << DCId << " / " << DCMoles_[DCId] << " / "
                  << numDCMolesDissolved << " / " << keepNumDCMoles << endl;
              std::clog << "scaledMass/massDissolved/totMassImpurity/"
                           "massDissolved - totMassImpurity : "
                        << scaledMass << " / " << massDissolved << " / "
                        << totMassImpurity << " / "
                        << massDissolved - totMassImpurity << endl;
              std::clog << endl << "end program" << endl;
              exit(0);
            }
          } else if (scaledMass < 0) {
            std::clog
                << endl
                << "KineticController::calculateKineticStep error for cyc = "
                << cyc << " : scaledMass < 0  i.e. scaledMass = " << scaledMass
                << "   massDissolved = " << massDissolved << endl;
            std::clog << "   midx/phName/scaledMassIni_[midx] : " << midx
                      << " / " << phaseKineticModel_[midx]->getName() << " / "
                      << scaledMassIni_[midx] << endl;
            std::clog << endl << "end program" << endl;
            exit(0);
          } else {
            keepNumDCMoles = 0.0;
          }

          // std::clog << "   :   DCMoles_/numDCMolesDissolved/keepNumDCMoles :
          // "
          //      <<  DCMoles_[DCId] << " / " << numDCMolesDissolved << " / "
          //      << keepNumDCMoles << endl;

          chemSys_->setDCLowerLimit(DCId, keepNumDCMoles);
          chemSys_->setDCUpperLimit(DCId, keepNumDCMoles);

          if (verbose_) {
            std::clog
                << "    calculateKineticStep - "
                   "midx/DCId/DCMoles_/numDCMolesDissolved/keepNumDCMoles : "
                << midx << " / " << DCId << " / " << DCMoles_[DCId] << " / "
                << numDCMolesDissolved << " / " << keepNumDCMoles << endl;
            std::clog << "    calculateKineticStep - scaledMass/massDissolved/"
                         "totMassImpurity/massDissolved - totMassImpurity : "
                      << scaledMass << " / " << massDissolved << " / "
                      << totMassImpurity << " / "
                      << massDissolved - totMassImpurity << endl;
          }
        }

        // ---------- CNT placement block (see docs/CNT_ARCHITECTURE.md §3) ----------
        // Runs at the natural END of each per-phase iteration so that any
        // early-continue path (doNotModif = true from GEMS "don't modify"
        // decisions, negative-mass exit, etc.) naturally skips CNT for
        // that phase this cycle.
        //
        // Sequence per phase:
        //   1. dN = model->computeNucleationVoxels(dt)  -- fractional voxels
        //   2. model->accumulateNucleation(dN)          -- per-phase accumulator
        //   3. nWant = floor(accumulator), fractional remainder carries forward
        //   4. Lattice::nucleatePhaseRnd places nWant voxels at uniform-random
        //      electrolyte sites; returns nPlaced (may be < nWant if sites exhausted)
        //   5. Convert placed voxel count to normalized-per-100g mass, moles,
        //      and volume (see block comment below for the scale factors).
        //   6. Update DCMoles_, KineticModel::scaledMass_, and
        //      ChemicalSystem::microPhaseMass_ so downstream code (GEMS
        //      constraint bounds, kinetic bookkeeping, and
        //      Lattice::changeMicrostructure via getMicroPhaseVolume) all see
        //      a consistent post-placement state.
        //   7. Raise BOTH DC lower and upper limits to the new DCMoles_[dcId] —
        //      lower to prevent GEMS from dissolving the nuclei back to the
        //      stale kinetic-computed floor, upper because the pre-loop CNT-lock
        //      may have zeroed it and DCMoles > 0 would violate UpperLimit = 0.
        if (useNucleationKinetics_) {
          // Two paths, selected per-phase by whether the config supplied
          // a jmak sub-block alongside the nucleation block:
          //
          //   JMAK PATH   (jmak block present): updateJMAKPhase applies the
          //     per-voxel Avrami-Poisson model — J·V·dt events distribute
          //     into generations, each generation's transformed fraction
          //     X_g evolves via the moment decomposition, and the lattice
          //     is synced to floor(sum_g N_g * X_g). This is the correct
          //     physical model at Portland-paste supersaturation regimes.
          //     See docs/POST_ALPHA_TODOS.md.
          //
          //   CLASSICAL PATH (no jmak block, or jmak disabled): the
          //     pre-2026-07-28 "1 voxel = 1 nucleation event" placement.
          //     Preserved for backward compatibility with existing
          //     configs and for phases where the JMAK per-voxel model
          //     is not needed.
          if (jmakEnabled_[midx]) {
            updateJMAKPhase(midx, timestep, cyc);
          } else {
          double dN =
              phaseKineticModel_[midx]->computeNucleationVoxels(timestep);
          if (dN > 0.0) {
            phaseKineticModel_[midx]->accumulateNucleation(dN);
            int nWant = phaseKineticModel_[midx]->drainNucleationInteger();
            if (nWant > 0) {
              int microPhId = phaseKineticModel_[midx]->getMicroPhaseId();
              int dcId = phaseKineticModel_[midx]->getDCId();
              int nPlaced = lattice_->nucleatePhaseRnd(microPhId, nWant);
              if (nPlaced > 0) {
                // Convert nPlaced voxels to normalized-per-100g state.
                //
                // THAMES stores DCMoles_[], microPhaseMass_[],
                // microPhaseVolume_[], and KineticModel::scaledMass_ in a
                // "per 100 g of total solid" reference frame, not physical
                // units. See Lattice::normalizePhaseMasses (Lattice.cc:815,
                // 834, 857) which establishes the convention. Any code that
                // combines physical getVolumePerVoxel()/getDCMolarVolume()
                // with these state variables needs the scale factor
                //     100 / (numSites * 1e6 * initSolidMass_)
                // to bridge the two frames (see reference conversion below).
                //
                // The formula mirrors Controller.cc:1296-1300, which does the
                // analogous conversion in the recall (numSitesNotAvailable)
                // path. Correct output is on the order of 10^-2 mol per 100g
                // for a 92k-voxel Portlandite placement in a 200³ Portland
                // paste; the pre-2026-07-25 formula (physical vVoxel/vMolar)
                // produced 10^-9 — six orders of magnitude too small — and
                // caused GEMS to hold Portlandite at picogram levels while
                // Lattice::changeMicrostructure dissolved the placement.
                double molarMass = chemSys_->getDCMolarMass(dcId);
                double vMolar = chemSys_->getDCMolarVolume(dcId);
                double vfracPlaced =
                    static_cast<double>(nPlaced) /
                    static_cast<double>(lattice_->getNumSites());
                double microPhaseMassPlacedPerCC =
                    vfracPlaced * molarMass / vMolar / 1.0e6; // g/cm³
                double placedMass =
                    microPhaseMassPlacedPerCC * 100.0 /
                    lattice_->getInitSolidMass();             // g per 100g solid
                double moles = placedMass / molarMass;        // mol per 100g solid

                DCMoles_[dcId] += moles;
                // Lock the newly-placed nuclei against GEMS redissolving
                // them back this cycle. Without this, GEMS uses the
                // kinetic-computed DCLowerLimit set at line ~1271 (pre-CNT
                // state) as its floor, and immediately "dissolves" all the
                // placed voxels back to that stale floor.
                // Also raise the UpperLimit — otherwise, if the pre-loop
                // "CNT-lock" block zeroed both limits (zero-mass CNT phase),
                // DCMoles > 0 would violate UpperLimit=0 and GEMS would
                // reject the placement.
                chemSys_->setDCLowerLimit(dcId, DCMoles_[dcId]);
                chemSys_->setDCUpperLimit(dcId, DCMoles_[dcId]);

                // Sync KineticModel::scaledMass_ AND
                // ChemicalSystem::microPhaseMass_/microPhaseVolume_ with the
                // placed mass. Reason: for kinetic phases (Standard/
                // Pozzolanic/SaturatingRate with a nucleation block),
                // ChemicalSystem::calculateState skips the
                // GEMPhaseVolume->microPhaseVolume fill (guarded by
                // `if (!isKinetic_[i])` at ChemicalSystem.cc:2763 and :2822).
                // Kinetic-phase volumes are maintained by KineticController
                // via updateMicroPhaseMasses. Without this update,
                // microPhaseVolume_[phase] stays at 0 all cycle,
                // Lattice::changeMicrostructure reads vfrac_next=0, and
                // dissolves every CNT-placed voxel the same cycle it was
                // placed.
                //
                // scaledMassIni_[midx] refresh matters if a lattice recall
                // path fires later in this cycle: Controller.cc:1302 calls
                // updateKineticStep which uses scaledMassIni_[midx] as the
                // reset baseline; a stale zero there would undo the CNT
                // placement in the recall path.
                double newScaledMass =
                    phaseKineticModel_[midx]->getScaledMass() + placedMass;
                phaseKineticModel_[midx]->setScaledMass(newScaledMass);
                chemSys_->updateMicroPhaseMasses(microPhId, newScaledMass, 1);
                scaledMassIni_[midx] = newScaledMass;

                if (verbose_) {
                  std::clog << "  CNT: "
                            << phaseKineticModel_[midx]->getName()
                            << " requested " << nWant << " voxels, placed "
                            << nPlaced << " (+" << moles << " mol DC["
                            << dcId << "], +" << placedMass
                            << " g/100g solid)" << endl;
                }
              }
              // If nPlaced < nWant, the residual is absorbed by exhausted
              // sites reality; do not re-accumulate. The fix for exhausted
              // electrolyte is a smaller dt next cycle, not a queue.
            }
          }
          }  // end classical-path else branch (JMAK is in the if branch above)
        }
      }

      if (verbose_ && doTweak) {
        std::clog << endl
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
    std::clog << endl
              << "     KineticController::calculateKineticStep : time >= "
                 "beginAttackTime_ -> "
              << time << " >= " << beginAttackTime_
              << " (hyd_time = " << hyd_time << ")" << endl;
    std::clog
        << "     KineticController::calculateKineticStep 0 : count_[VOIDID] = "
        << lattice_->getCount()[VOIDID] << "   &   count_[ELECTROLYTEID] = "
        << lattice_->getCount()[ELECTROLYTEID]
        << "  =>  waterMoles = " << DCMoles_[waterDCId_] << endl;

    double waterAddMoles = lattice_->fillAllPorosity(cyc);
    DCMoles_[waterDCId_] += waterAddMoles;

    std::clog
        << "     KineticController::calculateKineticStep 1 : count_[VOIDID] = "
        << lattice_->getCount()[VOIDID] << "   &   count_[ELECTROLYTEID] = "
        << lattice_->getCount()[ELECTROLYTEID]
        << "  =>  waterMoles = " << DCMoles_[waterDCId_] << endl;
  }

  for (i = 0; i < DCNum_; i++) {
    // std::clog << " " << i << "\t" << DCName_[i] << ": " << DCMoles_[i] << "
    // mol"
    // << endl;
    chemSys_->setDCMoles(i, DCMoles_[i]);
    // std::clog << "          " << DCName_[i] << ": " <<
    // chemSys_->getDCMoles(i) << " mol" << endl;
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
    std::clog << endl
              << "  KineticController::updateKineticStep - error for cyc = "
              << cyc << " & pId = " << pId << "  =>  DCId = " << DCId << " !!!"
              << endl;
    std::clog << "    scaledMass = " << scaledMass << endl;
    std::clog << endl << "  >>> program stop <<<" << endl;
    exit(0);

  } else {
    chemSys_->setMicroPhaseMass(pId, scaledMassIni_[midx]);
    modelName = phaseKineticModel_[midx]
                    ->getModelName(); // updateKineticStep(scaledMass
                                      // , massDissolved, timestep);
    std::clog << "  KineticController::updateKineticStep - for cyc = " << cyc
              << " & phaseId = " << pId << " ["
              << phaseKineticModel_[midx]->getName()
              << " / DCId:" << chemSys_->getMicroPhaseDCMembers(pId, 0) << "]"
              << endl;
    std::clog << "    midx = " << midx << "   modelName : " << modelName
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

    std::clog << "      massDissolved/totMassImpurity/massDissolved - "
                 "totMassImpurity "
                 ": "
              << massDissolved << " / " << totMassImpurity << " / "
              << massDissolved - totMassImpurity << endl;
    std::clog << "      DCMoles_/numDCMolesDissolved/keepNumDCMoles : "
              << DCMoles_[DCId] << " / " << numDCMolesDissolved << " / "
              << keepNumDCMoles << endl;
    // *****************

    if (hyd_time > beginAttackTime_) {
      std::clog << endl
                << "     KineticController::updateKineticStep : time >= "
                   "beginAttackTime_ -> "
                << hyd_time << " >= " << beginAttackTime_ << " (hyd_time)"
                << endl;
      std::clog
          << "     KineticController::updateKineticStep 0 : count_[VOIDID] "
             "= "
          << lattice_->getCount()[VOIDID] << "   &   count_[ELECTROLYTEID] = "
          << lattice_->getCount()[ELECTROLYTEID]
          << "  =>  waterMoles = " << DCMoles_[waterDCId_] << endl;

      double waterAddMoles = lattice_->fillAllPorosity(cyc);
      DCMoles_[waterDCId_] += waterAddMoles;

      // if (waterAddMoles > 0)
      //   std::clog << "     KineticController::updateKineticStep : check if
      //   OK!"
      //        << endl;

      std::clog
          << "     KineticController::updateKineticStep 1 : count_[VOIDID] "
             "= "
          << lattice_->getCount()[VOIDID] << "   &   count_[ELECTROLYTEID] = "
          << lattice_->getCount()[ELECTROLYTEID]
          << "  =>  waterMoles = " << DCMoles_[waterDCId_] << endl;
    }

    for (int i = 0; i < DCNum_; i++) {
      // std::clog << " " << i << "\t" << DCName_[i] << ": " << DCMoles_[i] << "
      // mol"
      // << endl;
      chemSys_->setDCMoles(i, DCMoles_[i]);
      // std::clog << "          " << DCName_[i] << ": " <<
      // chemSys_->getDCMoles(i)
      // << " mol" << endl;
    }
  }
}

void KineticController::setInitialSaturationIndices(
    const std::vector<double> &microPhaseSI) {
  //
  // Store the saturation indices from pre-equilibration GEMS calculation.
  // These are used by getMaxInitialDissolutionRate() to provide more
  // accurate rate estimates based on actual thermodynamic driving forces.
  //

  initialMicroPhaseSI_ = microPhaseSI;
  hasSIData_ = !microPhaseSI.empty();

  if (verbose_ && hasSIData_) {
    std::clog << "KineticController::setInitialSaturationIndices: "
              << "Stored SI data for " << microPhaseSI.size()
              << " microPhases" << std::endl;
  }
}

double KineticController::getMaxInitialDissolutionRate() const {
  //
  // Find the maximum initial dissolution rate across all kinetic models.
  // This determines the fastest-reacting phase, which constrains the
  // initial timestep for adaptive time stepping.
  //
  // If SI data from pre-equilibration is available, use the SI-aware
  // rate estimation for more accurate results.
  //

  double maxRate = 0.0;

  for (int i = 0; i < pKMsize_; ++i) {
    if (phaseKineticModel_[i] != nullptr) {
      double rate = 0.0;
      int microPhaseId = phaseKineticModel_[i]->getMicroPhaseId();

      // Use SI-aware estimation if SI data is available for this phase
      if (hasSIData_ && microPhaseId >= 0 &&
          microPhaseId < static_cast<int>(initialMicroPhaseSI_.size())) {
        double si = initialMicroPhaseSI_[microPhaseId];
        rate = phaseKineticModel_[i]->estimateInitialDissolutionRate(si);

        if (verbose_) {
          std::clog << "  KineticController::getMaxInitialDissolutionRate: "
                    << phaseKineticModel_[i]->getName()
                    << " SI=" << si << " rate=" << rate << " [1/h]"
                    << std::endl;
        }
      } else {
        // Fall back to no-SI estimation (assumes maximum driving force)
        rate = phaseKineticModel_[i]->estimateInitialDissolutionRate();
      }

      if (rate > maxRate) {
        maxRate = rate;
      }
    }
  }

  return maxRate;
}

double KineticController::computeKineticsBasedMaxTimestep(
    double maxRelativeChange) const {
  //
  // Compute the maximum timestep that limits DC mole changes to maxRelativeChange.
  //
  // For each kinetically-controlled DC:
  //   rate [mol/100g/h] = getCurrentMolarRate(scaledMass)
  //   moles = current DC moles in system
  //   dt_max = maxRelativeChange * moles / |rate|
  //
  // Return the minimum dt_max across all DCs.
  //

  const double largeTimestep = 1.0e6; // Return this if no constraint applies
  double minTimestep = largeTimestep;

  if (verbose_) {
    std::clog << "KineticController::computeKineticsBasedMaxTimestep: "
              << "maxRelativeChange=" << maxRelativeChange << std::endl;
  }

  for (int i = 0; i < pKMsize_; ++i) {
    if (phaseKineticModel_[i] != nullptr) {
      // Get the DC ID for this kinetic model
      int dcId = phaseKineticModel_[i]->getDCId();

      // Get the current scaled mass for this phase
      double scaledMass = phaseKineticModel_[i]->getScaledMass();

      // Skip phases with no mass remaining
      if (scaledMass <= 0.0) {
        continue;
      }

      // Get the current molar rate [mol/100g/h]
      double molarRate = phaseKineticModel_[i]->getCurrentMolarRate(scaledMass);

      // Skip phases with negligible rate
      if (fabs(molarRate) < 1.0e-20) {
        continue;
      }

      // Get the current DC moles from the chemical system
      double dcMoles = chemSys_->getDCMoles(dcId);

      // Skip DCs with negligible moles (avoid division issues)
      if (fabs(dcMoles) < 1.0e-15) {
        continue;
      }

      // Compute the timestep that would cause maxRelativeChange
      // |delta_moles| = |rate| * dt
      // |delta_moles| / |moles| = maxRelativeChange
      // dt = maxRelativeChange * |moles| / |rate|
      double dtConstraint = maxRelativeChange * fabs(dcMoles) / fabs(molarRate);

      if (verbose_) {
        std::clog << "  Phase " << phaseKineticModel_[i]->getName()
                  << ": dcMoles=" << dcMoles
                  << " molarRate=" << molarRate
                  << " dt_constraint=" << dtConstraint << " h" << std::endl;
      }

      if (dtConstraint < minTimestep) {
        minTimestep = dtConstraint;
      }
    }
  }

  if (verbose_) {
    if (minTimestep < largeTimestep) {
      std::clog << "KineticController::computeKineticsBasedMaxTimestep: "
                << "minTimestep=" << minTimestep << " h" << std::endl;
    } else {
      std::clog << "KineticController::computeKineticsBasedMaxTimestep: "
                << "No kinetics constraint applies" << std::endl;
    }
  }

  return minTimestep;
}

double KineticController::computeNucleationBasedMaxTimestep(
    double dtProposedHours, double capFraction) const {
  // The CNT rate J is fixed within a cycle (uses S from the previous
  // GEMS equilibration). Voxel count for a given dt scales linearly with dt.
  // If dt is too large, the fixed-J-in-cycle approximation breaks down
  // because the placement itself would depress solute concentration enough
  // to change J materially during the cycle.
  //
  // Two caps apply per phase; the tighter one wins:
  //
  //   (1) Electrolyte-fraction cap: capFraction * count_[ELECTROLYTEID].
  //       Prevents CNT from placing more voxels than a small fraction of
  //       the currently-saturated pore volume can host in one cycle.
  //
  //   (2) Mass-availability cap: for each IC that the phase's DC has
  //       non-zero stoichiometry, compute how many voxels the currently-
  //       aqueous IC moles could support. Take the min across ICs. This
  //       prevents CNT from placing voxels that GEMS then rolls back in
  //       Lattice::changeMicrostructure because Ca (or another IC) in
  //       solution can't feed them. This was added 2026-07-24 after S4
  //       validation of SaturatingRateModel showed CNT requesting ~92k
  //       Portlandite voxels/cycle at Portlandite SI ~ 10 while GEMS
  //       could only support ~125; see docs/POST_ALPHA_TODOS.md entry
  //       "CNT vs. Lattice::changeMicrostructure mass-balance mismatch"
  //       and docs/SATURATING_RATE.md §6.
  //
  // The per-phase effective cap = min(nCap_electrolyte, nCap_mass); if
  // computeNucleationVoxels(dtProposed) exceeds that, shrink dt so the
  // projected N would exactly land at the effective cap. Return the
  // minimum shrunk dt across all phases.
  //
  // Composes with the kinetics cap via sequential min at the two Controller
  // dt-selection sites (success path and post-failure path). Final dt is
  // min(adaptive-controller-proposed, kinetics cap, this CNT cap).
  if (!useNucleationKinetics_) return dtProposedHours;

  const double nElec =
      static_cast<double>(lattice_->getCount(ELECTROLYTEID));
  const double nCapElec = capFraction * nElec;
  double dtMax = dtProposedHours;

  // Snapshot the aqueous IC moles once per call. getSolution() iterates
  // over all DCs of aqueous class ('S','T','W') and accumulates
  // DCMoles * stoich, so O(numDCs * numICs). One snapshot is fine.
  std::vector<double> aqICMoles = chemSys_->getSolution();
  const double vVoxel = lattice_->getVolumePerVoxel();

  for (int i = 0; i < pKMsize_; ++i) {
    if (phaseKineticModel_[i] == nullptr) continue;
    if (!phaseKineticModel_[i]->hasNucleation()) continue;

    double nPhase =
        phaseKineticModel_[i]->computeNucleationVoxels(dtProposedHours);
    if (nPhase <= 0.0) continue;

    // Mass-availability cap for this phase.
    // nCapMass = min over ICs of (aqICMoles[ic] / molesPerVoxel_of_ic).
    // molesPerVoxel_of_ic = (vVoxel / vMolar_DC) * DCStoich[DCId][ic].
    int dcId = phaseKineticModel_[i]->getDCId();
    double vMolarDC = chemSys_->getDCMolarVolume(dcId);
    double nCapMass = std::numeric_limits<double>::infinity();
    if (vMolarDC > 0.0) {
      const double molesDCPerVoxel = vVoxel / vMolarDC;
      const int numICs = static_cast<int>(aqICMoles.size());
      // Last IC is charge (skip; DCStoich for charge would misbehave here).
      for (int ic = 0; ic < numICs - 1; ++ic) {
        const double stoich = chemSys_->getDCStoich(dcId, ic);
        if (stoich <= 0.0) continue;
        const double icPerVoxel = molesDCPerVoxel * stoich;
        if (icPerVoxel <= 0.0) continue;
        const double nFromThisIC = aqICMoles[ic] / icPerVoxel;
        if (nFromThisIC < nCapMass) nCapMass = nFromThisIC;
      }
    }

    const double effCap = std::min(nCapElec, nCapMass);
    if (nPhase > effCap && effCap >= 0.0) {
      // Guard against effCap == 0 which would send dtShrunk to 0.
      // In that case skip; the accumulator will hold the residual.
      if (effCap <= 0.0) continue;
      double dtShrunk = dtProposedHours * (effCap / nPhase);
      if (dtShrunk < dtMax) {
        dtMax = dtShrunk;
        if (verbose_) {
          const char *which =
              (nCapMass < nCapElec) ? "mass-availability" : "electrolyte-fraction";
          std::clog << "  CNT cap (" << which << "): phase "
                    << phaseKineticModel_[i]->getName()
                    << " would produce " << nPhase << " voxels at dt="
                    << dtProposedHours << " h (elec_cap=" << nCapElec
                    << ", mass_cap=" << nCapMass
                    << "); shrinking dt to " << dtShrunk << " h" << endl;
        }
      }
    }
  }
  return dtMax;
}
