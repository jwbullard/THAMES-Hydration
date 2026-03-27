/**
@file Controller.cc
@brief Definition of Controller class methods
*/

#include "Controller.h"

using std::cout;
using std::endl;
using std::map;
using std::string;
using std::vector;

Controller::Controller(Lattice *msh, KineticController *kc, ChemicalSystem *cs,
                       ThermalStrain *thmstr, const int simtype,
                       const string &jsonFileName, const string &jobname,
                       const bool verbose, const bool warning, const bool xyz) {

  stepTimeTHR_ = 1.e-5;  // Minimum timestep: 0.036 seconds (lowered from 1e-3)
  maxRelativeChange_ = 0.05; // Default: 5% max DC mole change per timestep

  // Initialize adaptive time stepping with unified parameters.
  // The kinetics-based timestep constraint (computeKineticsBasedMaxTimestep)
  // provides dynamic, physics-based adaptation that responds to actual
  // dissolution/precipitation rates regardless of model type.
  // These defaults may be overridden by simparams.json in parseDoc().
  useAdaptiveTimeStepping_ = true;
  AdaptiveTimeConfig adaptiveConfig;
  adaptiveConfig.dt_min = stepTimeTHR_;
  adaptiveConfig.dt_initial = 0.001;       // ~3.6 seconds (middle ground)
  adaptiveConfig.dt_max = 4.0;             // 4 hours max
  adaptiveConfig.growth_factor = 1.5;      // 50% growth per success
  adaptiveConfig.successes_for_growth = 2; // Grow after 2 consecutive successes
  adaptiveConfig.verbose = verbose;

  std::clog << "Controller: Using unified adaptive time stepping parameters"
            << " (kinetics constraint provides dynamic adaptation)" << endl;

  adaptiveTimeController_ = std::make_unique<AdaptiveTimeController>(adaptiveConfig);

  // Pre-equilibration: Run GEMS equilibration with initial electrolyte
  // composition to get realistic saturation indices for kinetic phases.
  // This allows the kinetic models to estimate dissolution rates based on
  // actual thermodynamic driving forces rather than assuming SI ≈ 0.
  std::clog << "Controller: Running pre-equilibration GEMS calculation..." << endl;
  try {
    cs->calculateSI(0, 0.0);  // cyc=0, time=0.0 for pre-equilibration

    // Pass the SI values to the KineticController
    std::vector<double> initialSI = cs->getMicroPhaseSI();
    kc->setInitialSaturationIndices(initialSI);

    if (verbose) {
      std::clog << "Controller: Pre-equilibration complete. "
                << "SI values for " << initialSI.size() << " phases available."
                << endl;
    }
  } catch (GEMException &gex) {
    std::clog << "WARNING: Pre-equilibration GEMS calculation failed. "
              << "Initial timestep will use conservative SI=0 assumption."
              << endl;
    // Continue without SI data - getMaxInitialDissolutionRate() will use
    // the no-argument version that assumes maximum driving force
  }

  // Set initial timestep based on kinetic rates
  // This provides a physics-based initial timestep rather than hard-coded default
  double maxKineticRate = kc->getMaxInitialDissolutionRate();
  if (maxKineticRate > 0.0) {
    adaptiveTimeController_->setInitialTimestepFromKinetics(maxKineticRate, maxRelativeChange_);
    std::clog << "Controller: Initial timestep set from kinetics, maxRate="
              << maxKineticRate << " 1/h, maxRelChange=" << (maxRelativeChange_ * 100)
              << "%, dt_initial="
              << adaptiveTimeController_->getCurrentTimestep() << " h" << endl;
  }

  // xyz_ = xyz;
  xyzMovie_ = xyz;
  xyzFiles_ = false;

  simType_ = simtype;
  chemSys_ = cs;
  kineticController_ = kc;
  lattice_ = msh;
  thermalstr_ = thmstr;
  jobRoot_ = jobname;

#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  ///
  /// Set default values for all parameters prior to any customization
  ///
  /// Setting the default times for outputting images, and for initiating the
  /// simulations for leaching or external sulfate attack.
  ///
  /// All times are given in hours, and the leaching and sulfate attack times
  /// are set to very high values so that they usually won't happen
  ///

  outputImageTimeInterval_ = -1.0;
  outputImageTime_.clear();

  leachTime_ = 1.0e10;
  sulfateAttackTime_ = 1.0e10;

  beginAttackTime_ = -1.0;
  endAttackTime_ = -1.0;
  attackTimeInterval_ = -1.0;

  oldDamageCount_ = 0;
  allDamageCount_ = 0;

  isParrotKilloh_ = chemSys_->getIsParrotKilloh();
  sizePK_ = isParrotKilloh_.size();

  numDCs_ = chemSys_->getNumDCs();
  numICs_ = chemSys_->getNumICs();
  numMicroPhases_ = chemSys_->getNumMicroPhases();
  numGEMPhases_ = chemSys_->getNumGEMPhases();

  ///
  /// Load up the pointers to the `ChemicalSystem` object and `Lattice` object
  ///

  // chemSys_ = lattice_->getChemSys();
  lattice_->setJobRoot(jobRoot_);

  ///
  /// Output the class codes for the solution and for DC components.
  /// Output the header for the microstructure phase stats file
  /// Output header for the file tracking pH
  /// Output header for the file tracking the C-S-H composition and Ca/Si ratios
  /// Output header for the file tracking the IC moles in the system
  ///

  ofstream outfs;
  string outfilename;

  try {
    outfilename = jobRoot_ + "_Solution.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    char cc;
    outfs << "Time(h)";
    for (int i = 0; i < numDCs_; i++) {
      cc = chemSys_->getDCClassCode(i);
      if (cc == 'S' || cc == 'T' || cc == 'W') {
        outfs << "," << chemSys_->getDCName(i);
      }
    }
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_DCVolumes.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }

    outfs << "Time(h)";
    for (int i = 0; i < numDCs_; i++) {
      cc = chemSys_->getDCClassCode(i);
      if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M' || cc == 'W') {
        outfs << "," << chemSys_->getDCName(i) << "(m3/100g)";
      }
    }
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_PhaseVolumes.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }

    outfs << "Time(h)";
    for (int i = 0; i < numGEMPhases_; i++) {
      outfs << "," << chemSys_->getGEMPhaseName(i) << "(m3/100g)";
    }
    outfs << ",Total Volume (m3/100g),Chemical Shrinkage Strain";
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_SurfaceAreas.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }

    outfs << "Time(h)";
    // JWB: Start with micro phase id 1 to avoid the Void phase
    for (int i = 1; i < numMicroPhases_; i++) {
      int dcid = chemSys_->getMicroPhaseDCMembers(i, 0);
      cc = chemSys_->getDCClassCode(dcid);
      if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M') {
        outfs << ",A_" << chemSys_->getMicroPhaseName(i) << "(m2/100g)";
      }
    }
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_SI.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }

    outfs << "Time(h)";
    // JWB: Start with micro phase id 1 to avoid the Void phase
    for (int i = 1; i < chemSys_->getNumMicroPhases(); i++) {
      int dcid = chemSys_->getMicroPhaseDCMembers(i, 0);
      cc = chemSys_->getDCClassCode(dcid);
      if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M') {
        outfs << ",SI_" << chemSys_->getMicroPhaseName(i);
      }
    }
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_CSH.csv";
    outfs.open(outfilename.c_str(), std::ios::app);
    if (!outfs) {
      throw FileException("Controller", "calculateState", outfilename,
                          "Could not append");
    }
    outfs << "Time(h)";
    for (int i = 0; i < chemSys_->getNumICs(); i++) {
      outfs << "," << chemSys_->getICName(i);
    }
    outfs << ",Ca/Si" << endl;
    outfs.close();

    outfilename = jobRoot_ + "_CSratio_solid.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "calculateState", outfilename,
                          "Could not append");
    }
    outfs << "Time(h),Ca/Si Ratio" << endl;
    outfs.close();

    outfilename = jobRoot_ + "_Microstructure.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }

    outfs << "Time(h)";
    for (int i = 0; i < chemSys_->getNumMicroPhases(); i++) {
      outfs << "," << chemSys_->getMicroPhaseName(i);
    }
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_pH.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    outfs << "Time(h),pH" << endl;
    outfs.close();

    outfilename = jobRoot_ + "_Enthalpy.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    outfs << "Time(h),Enthalpy(J/100g)" << endl;
    outfs.close();

    outfilename = jobRoot_ + "_DOR.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    outfs << "Time(h),DOR";
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_GelSpaceRatio.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    outfs << "Time(h),DOR,GelVol,VoidElectrolyteVol,GelSpaceRatioSim,"
             "GelSpaceRatioCalc";
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_reactProdRatio.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    outfs << "Time(h)";
    for (int i = 0; i < numMicroPhases_; i++) {
      if ((i >= FIRST_SOLID) && (!chemSys_->isCementComponent(i))) {
        outfs << "," << chemSys_->getMicroPhaseName(i);
      }
    }
    outfs << endl;
    outfs.close();

  } catch (FileException fex) {
    throw fex;
  }

  temperature_ = chemSys_->getTemperature();
  lattice_->setTemperature(temperature_);
  waterDCId_ = chemSys_->getDCId("H2O@");
  waterMolarMass_ = chemSys_->getDCMolarMass("H2O@");
  numSites_ = lattice_->getNumSites();
  initMicroVolume_ = chemSys_->getInitMicroVolume();

  std::clog << endl << "Controller::Controller(...) :" << endl;

  std::clog << endl
            << "   numGEMPhases_   = " << std::setw(3) << std::right
            << numGEMPhases_ << endl;
  std::clog << "   numMicroPhases_ = " << std::setw(3) << std::right
            << numMicroPhases_ << endl;
  std::clog << "   numDCs_         = " << std::setw(3) << std::right << numDCs_
            << endl;
  std::clog << "   numICs_         = " << std::setw(3) << std::right << numICs_
            << endl;
  std::clog << "   waterDCId_      = " << std::setw(3) << std::right
            << waterDCId_ << " (waterDCName = \""
            << chemSys_->getDCName(waterDCId_) << "\")" << endl;

  ///
  /// Output a file that directly links the microstructure ids to their
  /// rgb color.  This is only for easier image processing after the simulation
  /// is finished so we don't have to read the xml file
  ///

  lattice_->writeMicroColors();

  ///
  /// Write the initial microstructure image and its png image
  ///

  std::clog << endl
            << "Controller::Controller(...) - write initial "
               "microstructure files (writeLattice(0.0), etc)"
            << endl;

  string curTimeString = getTimeString(0.0);
  lattice_->writeLattice(curTimeString);
  lattice_->writeLatticePNG(curTimeString);

  if (xyzFiles_)
    lattice_->writeLatticeXYZ(0.0, curTimeString);

  if (xyzMovie_)
    lattice_->appendXYZ(0.0);

  ///
  /// Open and read the Controller parameter file
  ///

  string jsonext = ".json";
  size_t foundjson;

  try {
    foundjson = jsonFileName.find(jsonext);

    if (foundjson != string::npos) {
      parseDoc(jsonFileName);
    } else {
      std::clog << "Parameter file must be JSON" << endl;
      throw FileException("Controller", "Controller", jsonFileName,
                          "NOT JSON FORMAT");
    }
  } catch (FileException fex) {
    throw fex;
  }

  for (int i = 0; i < static_cast<int>(time_.size() - 1); i++) {
    if (abs(time_[i] - time_[i + 1]) <= 1.0e-6) {
      time_.erase(time_.begin() + i);
    }
  }

  // Re-compute kinetics-based initial timestep if maxRelativeChange_ was
  // updated by parseDoc() (from simparams.json adaptive_stepping config)
  if (useAdaptiveTimeStepping_ && adaptiveTimeController_) {
    double maxKinRate = kc->getMaxInitialDissolutionRate();
    if (maxKinRate > 0.0) {
      adaptiveTimeController_->setInitialTimestepFromKinetics(maxKinRate, maxRelativeChange_);
      std::clog << "Controller: Re-computed initial timestep after config load, "
                << "dt_initial=" << adaptiveTimeController_->getCurrentTimestep()
                << " h" << endl;
    }
  }

  int time_Size = time_.size();
  int outputImageTimeSize = outputImageTime_.size();

  if (time_Size == 0) {
    std::clog << endl
              << endl
              << "Controller::Controller error : a final time and "
              << "at least one output time value must be present in "
              << "the simulation parameters file!" << endl;
    std::clog
        << endl
        << "check and modify the simulation parameters file and run thames "
           "again"
        << endl;
    std::clog << endl << "end program" << endl;
    // exit(0);
    throw FileException(
        "Controller", "Controller", "simparams.json",
        "a final time and at least one "
        "outtime value must be present into simparams.json file!");
  }

  outfilename = jobRoot_ + "-times_used.json";
  outfs.open(outfilename.c_str());
  if (!outfs) {
    throw FileException("Controller", "Controller", outfilename,
                        "Could not append");
  }

  outfs << "{" << endl;
  outfs << "  \"time_parameters\": {" << endl;
  outfs << "    \"calctimes\": [" << endl;

  int j = 0;
  for (int i = 0; i < time_Size; i++) {
    j++;
    if (i < (time_Size - 1)) {
      if (j == 1) {
        outfs << "        " << std::fixed << time_[i] << ", ";
      } else if (j < 7) {
        outfs << time_[i] << ", ";
      } else {
        j = 0;
        outfs << time_[i] << "," << endl;
      }
    } else {
      if (j == 1) {
        outfs << "        " << time_[i] << endl;
      } else {
        outfs << time_[i] << endl;
      }
    }
  }
  outfs << "    ]," << endl;
  outfs << "    \"outtimes\": [" << endl;

  j = 0;
  for (int i = 0; i < outputImageTimeSize; i++) {
    j++;
    if (i < (outputImageTimeSize - 1)) {
      if (j == 1) {
        outfs << "        " << outputImageTime_[i] << ", ";
      } else if (j < 7) {
        outfs << outputImageTime_[i] << ", ";
      } else {
        j = 0;
        outfs << outputImageTime_[i] << "," << endl;
      }
    } else {
      if (j == 1) {
        outfs << "        " << outputImageTime_[i] << endl;
      } else {
        outfs << outputImageTime_[i] << endl;
      }
    }
  }

  if (attack_) {
    outfs << "    ]," << endl;
    outfs << "    \"beginattacktime\": " << beginAttackTime_ << "," << endl;
    outfs << "    \"endattacktime\": " << endAttackTime_ << "," << endl;
    outfs << "    \"attacktimeinterval\": " << attackTimeInterval_ << endl;
    outfs << "  }" << endl;
    outfs << "}" << endl;
  } else {
    outfs << "    ]," << endl;
    outfs << "    \"beginattacktime\": -1," << endl;
    outfs << "    \"endattacktime\": -1," << endl;
    outfs << "    \"attacktimeinterval\": -1" << endl;
    outfs << "  }" << endl;
    outfs << "}" << endl;
  }
  outfs.close();

  std::clog
      << endl
      << "   => time values (calctime & outtime in hours) have been used and "
         "writen as :"
      << endl;
  std::clog << "         " << outfilename << endl;

  leachTime_ = beginAttackTime_; // not default leachTime_ = 1.0e10
  sulfateAttackTime_ =
      beginAttackTime_; // not default sulfateAttackTime_ = 1.0e10

  chemSys_->setIniAttackTime(sulfateAttackTime_);
  lattice_->setSulfateAttackTime(sulfateAttackTime_);
  lattice_->setLeachTime(leachTime_);

  kineticController_->setIniAttackTime(sulfateAttackTime_);

  if (simType_ == SULFATE_ATTACK) {
    /*
    std::clog << endl << "   => attack = " << attack_ << endl;
    std::clog << "   parameters in hours:" << endl;
    std::clog << "     -> beginattacktime = " << std::setw(7) << std::right
         << static_cast<int>(beginAttackTime_) << endl;
    std::clog << "     -> endattacktime = " << std::setw(7) << std::right
         << static_cast<int>(endAttackTime_) << endl;
    std::clog << "     -> attacktimeinterval = " << std::setw(7) << std::right
         << static_cast<int>(attackTimeInterval_) << endl;
    */
    lattice_->createGrowingVectSA();
  }
}

void Controller::doCycle(double elemTimeInterval) {

  RestoreSystem iniLattice;
  RestoreSite site_l;
  RestoreInterface interface_l;

  ///
  /// This block arbitrarily sets the leaching initiation time to 100 days if
  /// the leaching module is to be run, or sets the sulfate attack initiation
  /// time to 100 days if the sulfate attack module is to be run
  ///
  /// @todo Think about generalizing this more, or allowing combinations of more
  /// than one (LEACHING/SULFATE_ATTACK/...)
  ///

  // Initialize the list of all interfaces in the lattice

  std::clog << endl
            << "Controller::doCycle(...) Entering Lattice::findInterfaces()"
            << endl;

  lattice_->findInterfaces();

  // lattice_->checkSite(8);
  // std::clog << endl << " exit controller" << endl;// exit(0);

  std::clog << endl
            << "Controller::doCycle(...) Entering Main time loop" << endl;

  static double timestep = 0.0;

  int timeIndexIMG = 0;

  int timesGEMFailed_loc = 0;

  int DCId;
  for (int i = FIRST_SOLID; i < numMicroPhases_; i++) {
    if (chemSys_->isKinetic(i)) {
      DCId = chemSys_->getMicroPhaseDCMembers(i, 0);
      chemSys_->setIsDCKinetic(DCId, true);
    }
  }
  std::clog << endl
            << "   numGEMPhases_  = " << std::setw(3) << std::right
            << chemSys_->getNumGEMPhases() << endl;
  std::clog << "   numDCs_        = " << std::setw(3) << std::right
            << chemSys_->getNumDCs() << endl;
  std::clog << "   numICs_        = " << std::setw(3) << std::right
            << chemSys_->getNumICs() << endl;

  std::clog << endl
            << "   ICTHRESH = " << std::setprecision(1) << ICTHRESH << " mol"
            << std::setprecision(15) << endl;

  chemSys_->setInitialElectrolyteComposition();
  chemSys_->setInitialGasComposition();
  chemSys_->calculateSI(0, 0.0);

  bool writeICsDCs = true;
  if (writeICsDCs)
    writeTxtOutputFiles_onlyICsDCs(0); // to check the total ICs

  writeTxtOutputFiles(0);

  // variables used in DCLowerLimit computation
  double volMolDiff, molarMassDiff, vfracDiff;
  double microPhaseMassDiff, scaledMassDiff, numMolesDiff;
  int timesGEMFailed_recall;

  std::clog << endl << endl << "     ===== START SIMULATION =====" << endl;

  // JWB 2024-03-18: Florin had this set to 1.e-7
  // JWB 2024-03-18: I increased it due to going from days to hours for time
  // units
  deltaTime_ = elemTimeInterval; // 1.e-6; // 1.e-5;

  double old_timestep = 0.0;

  double thrTimeToWriteLattice = 0.0167; // threshold ~ 1 minute

  /// used to check if all DCMoles_ have right values according to corresponding
  ///   DCUpperLimit_ & DCUpperLimit_
  // vector<int> updateDCId, updatePHId;
  // updateDCId.clear();
  // updatePHId.clear();

  double dissTimeInterval = 0.0;

  int i = 0;
  int cyc = 0;
  int timeSize = time_.size();
  double initialLastTime = time_[timeSize - 1];

  int lastGoodI = 0;
  lastGoodTime_ = 0.0;
  double currTime = 0.0;
  bool nucStopPrg = false;
  std::ofstream progressFile;

  // Main computation cycle
  // Use initialLastTime instead of time_[timeSize-1] because the time_ vector
  // gets modified during the loop (entries overwritten with actual cycle times).
  while ((currTime <= initialLastTime) &&
         (chemSys_->getGEMPhaseVolume(ElectrolyteGEMName) > 0.0)) {

    // Update progress file if necessary

    if (cyc % 10 == 0) {
      progressFile.open("progress.json", std::ios::out);
      if (!progressFile) {
        throw FileException("Controller", "calculateState", "progress.json",
                            "Could not write");
      }
      progressFile << "{";
      progressFile << "\"cycle\": " << cyc << ", ";
      progressFile << "\"time_hours\": " << currTime << ", ";
      progressFile << "\"target_time_hours\": " << initialLastTime
                   << ", ";
      progressFile << "\"timestamp\": \"" << lgf::time_stamp() << "\"";
      progressFile << "}";
      progressFile.close();
    }

    ///
    /// Do not advance the time step if GEM_run failed the last time
    ///

    if (timesGEMFailed_loc == 0) {
      lastGoodTime_ = currTime;
      std::clog << endl
                << "Controller::doCycle - vector time_ before next cycle - "
                   "i/cyc/lastGoodI = "
                << i << " / " << cyc << " / " << lastGoodI
                << " : lastGoodTime_ = " << lastGoodTime_;
      if (i > 0) {
        std::clog << "  => modif: time_[" << lastGoodI
                  << "] = " << time_[lastGoodI] << "  --->  time_[" << lastGoodI
                  << "] = " << currTime << endl;
        time_[lastGoodI] = currTime;
        lastGoodI++;
      } else {
        std::clog << "  => time_[" << lastGoodI << "] = " << time_[lastGoodI]
                  << " / timeSize = " << timeSize
                  << " / time_[timeSize - 1] = " << time_[timeSize - 1]
                  << " / currTime = " << currTime << endl;
      }
    }

    i++;

    bool isFirst = (cyc == 0) ? true : false;

    cyc++;

    if (currTime < beginAttackTime_) {
      if (timesGEMFailed_loc > 0) {
        // GEMS failed - adjust timestep using adaptive controller or legacy halving
        i--;

        old_timestep = timestep;

        if (useAdaptiveTimeStepping_ && adaptiveTimeController_) {
          // Use adaptive time controller to get new timestep
          timestep = adaptiveTimeController_->recordFailure(timesGEMFailed_loc);

          // Check if we've exceeded max failures
          if (adaptiveTimeController_->hasExceededMaxFailures()) {
            std::clog << endl
                      << "##### Controller::doCycle - ADAPTIVE TIME STEPPING: "
                         "MAX FAILURES EXCEEDED #####"
                      << endl;
            std::clog << adaptiveTimeController_->getStatusString() << endl;

            // Write exit status for UI
            std::ostringstream diagStream;
            long int totalCycles = adaptiveTimeController_->getTotalSuccesses() +
                                   adaptiveTimeController_->getTotalFailures();
            diagStream << "Simulation stopped at " << lastGoodTime_
                       << " hours after " << totalCycles
                       << " cycles. " << adaptiveTimeController_->getStatusString();
            writeExitStatus(1, "GEMS solver exceeded maximum consecutive failures",
                            diagStream.str());
            break;
          }

          // Calculate new current time based on adaptive timestep
          currTime = lastGoodTime_ + timestep;

        } else {
          // Legacy behavior: halve the timestep
          timestep = timestep / 2;

          if (timestep > stepTimeTHR_) { // stepTimeTHR_ = 1.e-3
            currTime -= timestep;
          } else {
            if (i < timeSize) {
              currTime = time_[i];
              time_.erase(time_.begin() + (i - 1));
              timeSize = time_.size();
              timestep = currTime - lastGoodTime_;
            } else {
              std::clog << endl
                        << endl
                        << endl
                        << "##### Controller::doCycle - START NEW CYCLE - "
                           "CHANGE TIME NOT POSSIBLE #####"
                        << endl;
              std::clog
                  << endl
                  << "##### i/timeSize/lastGoodTime_/initialLastTime/currTime : "
                  << i << " / " << timeSize << " / " << lastGoodTime_ << " / "
                  << initialLastTime << " / " << currTime
                  << " (time in hours) #####" << endl;
              break;
            }
          }
        }

        dissTimeInterval = currTime - lastGoodTime_;

        std::clog << endl
                  << "Controller::doCycle GEM failed (" << timesGEMFailed_loc
                  << " times) - i/cyc = " << i << " / " << cyc
                  << "  =>  START NEW CYCLE - change time_[i] - "
                     "lastGoodI/lastGoodTime_/"
                     "time_[i]/currTime         : "
                  << lastGoodI << " / " << lastGoodTime_ << " / " << time_[i]
                  << " / " << currTime << endl;
        std::clog << "Controller::doCycle GEM failed (" << timesGEMFailed_loc
                  << " times) - i/cyc = " << i << " / " << cyc
                  << "  =>  START NEW CYCLE - change time_[i] - "
                     "lastGoodI/lastGoodTime_/"
                     "old_timestep/new_timestep : "
                  << lastGoodI << " / " << lastGoodTime_ << " / "
                  << old_timestep << " / " << timestep << endl;

      } else {
        // Previous GEMS call succeeded - compute next timestep

        if (useAdaptiveTimeStepping_ && adaptiveTimeController_) {
          // Use adaptive time controller to get next timestep
          timestep = adaptiveTimeController_->getNextTimestep();

          // Apply kinetics-based constraint to limit DC mole changes per step
          double kineticsMaxTimestep = kineticController_->computeKineticsBasedMaxTimestep(maxRelativeChange_);
          if (kineticsMaxTimestep < timestep) {
            std::clog << "    Kinetics constraint: reducing timestep from "
                      << timestep << " to " << kineticsMaxTimestep << " h" << std::endl;
            timestep = kineticsMaxTimestep;
          }

          // Ensure we don't overshoot final simulation time.
          // Use initialLastTime (saved before the time_ vector gets modified
          // by the cycle loop) rather than time_[timeSize-1] which gets
          // overwritten with lastGoodTime_ each cycle.
          if (lastGoodTime_ + timestep > initialLastTime) {
            timestep = initialLastTime - lastGoodTime_;
          }

          // Check if we've reached the final time (timestep would be zero or tiny)
          if (timestep < 1.0e-12) {
            std::clog << std::endl
                      << std::endl
                      << "##### Controller::doCycle - ADAPTIVE TIME STEPPING: "
                         "FINAL TIME REACHED #####"
                      << std::endl;
            std::clog << "##### lastGoodTime_ = " << lastGoodTime_
                      << " h, finalTime = " << initialLastTime << " h #####" << std::endl;

            // Write exit status for UI - normal completion
            std::ostringstream diagStream;
            long int totalCycles = adaptiveTimeController_->getTotalSuccesses() +
                                   adaptiveTimeController_->getTotalFailures();
            diagStream << "Simulation completed " << lastGoodTime_ << " hours in "
                       << totalCycles << " cycles. "
                       << adaptiveTimeController_->getStatusString();
            writeExitStatus(0, "Simulation completed successfully", diagStream.str());
            break;
          }

          // Ensure we hit output times exactly
          for (size_t ot = 0; ot < outputImageTime_.size(); ot++) {
            double outTime = outputImageTime_[ot];
            if (outTime > lastGoodTime_ && outTime < lastGoodTime_ + timestep) {
              // There's an output time before our proposed next time
              timestep = outTime - lastGoodTime_;
              break;
            }
          }

          currTime = lastGoodTime_ + timestep;
          dissTimeInterval = timestep;

          chemSys_->setAllKeepDCLowerLimitZero();

          std::clog << endl
                    << endl
                    << endl
                    << "##### Controller::doCycle - START NEW CYCLE (ADAPTIVE) "
                       "i/cyc/currTime/lastGoodTime_/timestep : "
                    << i << " / " << cyc << " / " << currTime << " / "
                    << lastGoodTime_ << " / " << timestep
                    << " (time in hours) #####" << endl;

        } else {
          // Legacy mode: use pre-generated time array
          if (i - 1 < timeSize) {
            currTime = time_[i - 1];

            timestep = currTime - lastGoodTime_;

            dissTimeInterval = timestep;

            chemSys_->setAllKeepDCLowerLimitZero();

            std::clog << endl
                      << endl
                      << endl
                      << "##### Controller::doCycle - START NEW CYCLE   "
                         "i/cyc/time_[i]/currTime/lastGoodTime_/timestep : "
                      << i << " / " << cyc << " / " << time_[i] << " / "
                      << currTime << " / " << lastGoodTime_ << " / " << timestep
                      << " (time in hours) #####" << endl;
          } else {
            std::clog << endl
                      << endl
                      << endl
                      << "##### Controller::doCycle - START NEW CYCLE : "
                         "ALL TIME VALUES HAVE BEEN USED !  #####"
                      << endl;
            std::clog
                << endl
                << "##### i/timeSize/lastGoodTime_/initialLastTime/currTime : "
                << i << " / " << timeSize << " / " << lastGoodTime_ << " / "
                << initialLastTime << " / " << currTime
                << " (time in hours) #####" << endl;
            break;
          }
        }
      }

    } else { // for SA
      if (timesGEMFailed_loc > 0) {
        std::clog << endl
                  << ">>>>> Controller::doCycle - SA - timesGEMFailed_loc > 0 "
                     "=> i/cyc/time_[i]/currTime/lastGoodTime_/timestep : "
                  << i << " / " << cyc << " / " << time_[i] << " / " << currTime
                  << " / " << lastGoodTime_ << " / " << timestep << endl;
        std::clog << endl
                  << ">>>>> \"normal exit\" (for now!) under sulfate attack "
                     "conditions"
                  << endl;
        // exit(0);
        bool is_Error = false;
        throw MicrostructureException(
            "Controller", "doCycle", "SA - first GEM_run failed (0)", is_Error);
      } else {
        if (i - 1 < timeSize) {
          currTime = time_[i - 1];

          timestep = currTime - lastGoodTime_;

          dissTimeInterval = timestep;
          std::clog << endl
                    << endl
                    << endl
                    << "##### Controller::doCycle - START NEW CYCLE - SA  "
                       "i/cyc/time_[i]/currTime/lastGoodTime_/timestep : "
                    << i << " / " << cyc << " / " << time_[i] << " / "
                    << currTime << " / " << lastGoodTime_ << " / " << timestep
                    << " (time in hours) #####" << endl;
        } else {
          std::clog << endl
                    << endl
                    << endl
                    << "##### Controller::doCycle - START NEW CYCLE - SA : "
                       "ALL TIME VALUES HAVE BEEN USED !  #####"
                    << endl;
          std::clog
              << endl
              << "##### i/timeSize/lastGoodTime_/initialLastTime/currTime : "
              << i << " / " << timeSize << " / " << lastGoodTime_ << " / "
              << initialLastTime << " / " << currTime
              << " (time in hours) #####" << endl;
          break;
        }
      }
    }

    /// Assume that only voxel-scale pore water is chemically reactive,
    /// while water in nanopores is chemically inert.
    ///
    ///
    /// This is the main step of the cycle; the calculateState method
    /// runs all the major steps of a computational cycle

    try {

      // std::clog << endl
      //      << "test_calculateState_ini - i/cyc/currTime/timesGEMFailed_loc: "
      //      << i << " / " << cyc << " / " << currTime << " / " <<
      //      timesGEMFailed_loc
      //      << endl;

      timesGEMFailed_loc =
          calculateState(currTime, dissTimeInterval, isFirst, cyc);

      if (timesGEMFailed_loc > 0) {
        if (currTime < beginAttackTime_) {
          std::clog << endl
                    << "Controller::doCycle GEM_run failed "
                       "i/cyc/time[i]/currTime/getTimesGEMFailed_loc : "
                    << i << " / " << cyc << " / " << time_[i] << " / "
                    << currTime << " / " << timesGEMFailed_loc << endl;
          continue;
        } else {
          std::clog << endl
                    << "Controller::doCycle => SA - first GEM_run failed "
                       "i/cyc/time[i]/currTime/getTimesGEMFailed_loc : "
                    << i << " / " << cyc << " / " << time_[i] << " / "
                    << currTime << " / " << timesGEMFailed_loc << endl;

          std::clog << endl
                    << ">>>>> \"normal exit\" (for now!) under sulfate attack "
                       "conditions"
                    << endl;
          // exit(0);
          bool is_Error = false;
          throw MicrostructureException("Controller", "doCycle",
                                        "SA - first GEM_run failed (1)",
                                        is_Error);
        }
      } else {

        std::clog << endl
                  << "Controller::doCycle GEM_run OK "
                     "i/cyc/time[i]/currTime/getTimesGEMFailed_loc: "
                  << i << " / " << cyc << " / " << time_[i] << " / " << currTime
                  << " / " << timesGEMFailed_loc << endl;

        // Record successful GEMS calculation for adaptive time stepping
        if (useAdaptiveTimeStepping_ && adaptiveTimeController_) {
          long int iterDone = chemSys_->getIterDone();
          double pci = chemSys_->getPCI();
          double dxm = chemSys_->getDXM();
          adaptiveTimeController_->recordSuccess(iterDone, pci, dxm);

          if (verbose_) {
            std::clog << "AdaptiveTime status: "
                      << adaptiveTimeController_->getStatusString() << endl;
          }
        }
      }

    } catch (GEMException gex) {

      string curTimeString = getTimeString(currTime);
      lattice_->writeLattice(curTimeString);
      lattice_->writeLatticePNG(curTimeString);

      if (xyzFiles_)
        lattice_->writeLatticeXYZ(currTime, curTimeString);

      if (xyzMovie_)
        lattice_->appendXYZ(currTime);

      throw gex;
    }

    ///
    /// Once the change in state is determined, propagate the consequences
    /// to the 3D microstructure only if the GEM_run calculation succeeded.
    /// Otherwise we adjust the time step and try again
    ///

    // if (time_[i] < beginAttackTime_) {
    //   if (timesGEMFailed_loc > 0) {
    //     std::clog << endl
    //          << "  Controller::doCycle first GEM_run failed "
    //             "i/cyc/time[i]/getTimesGEMFailed_loc : "
    //          << i << " / " << cyc << " / " << time_[i] << " / "
    //          << timesGEMFailed_loc << endl;
    //
    //    //**********************
    //    ...
    //
    //   this old algorithm to search for a good time value using ten small
    //   intervals has been replaced by another one and moved into
    //   Controller::calculateState(...)
    //
    //    ...
    // }

    if (verbose_) {
      std::clog << "Controller::doCycle Entering Lattice::changeMicrostructure"
                << endl;
      std::clog.flush();
    }

    ///
    /// Next function can encounter EOB exceptions within but they
    /// are caught there and the program will then exit from within
    /// this function rather than throwing an exception itself
    ///

    try {

      // set iniLattice i.e. a copy of initial lattice/system configuration
      // (including RNG state and all DCs values)
      // from ChemicalSystem:
      iniLattice.DCMoles = kineticController_->getDCMoles();

      // from Lattice:
      iniLattice.count = lattice_->getCount();

      iniLattice.growthInterfaceSize = lattice_->getGrowthInterfaceSize();

      iniLattice.dissolutionInterfaceSize =
          lattice_->getDissolutionInterfaceSize();

      iniLattice.site.clear();
      iniLattice.site.shrink_to_fit();

      // RestoreSite site_l; // only one declaration
      for (int ij = 0; ij < numSites_; ij++) {
        site_l.microPhaseId = (lattice_->getSite(ij))->getMicroPhaseId();
        site_l.growth = (lattice_->getSite(ij))->getGrowthPhases();
        site_l.wmc = (lattice_->getSite(ij))->getWmc();
        site_l.wmc0 = (lattice_->getSite(ij))->getWmc0();
        site_l.visit = 0;
        site_l.inGrowInterfacePos =
            (lattice_->getSite(ij))->getInGrowInterfacePosVector();
        site_l.inDissInterfacePos =
            (lattice_->getSite(ij))->getInDissInterfacePos();
        iniLattice.site.push_back(site_l);
      }

      iniLattice.interface.clear();
      iniLattice.interface.shrink_to_fit();

      // RestoreInterface interface_l; // only one declaration
      int dimLatticeInterface = lattice_->getInterfaceSize();
      for (int ij = 0; ij < dimLatticeInterface; ij++) {
        interface_l.microPhaseId = lattice_->getInterface(ij).getMicroPhaseId();
        interface_l.growthSites = lattice_->getInterface(ij).getGrowthSites();
        interface_l.dissolutionSites =
            lattice_->getInterface(ij).getDissolutionSites();
        iniLattice.interface.push_back(interface_l);
      }

      iniLattice.numRNGcall_0 = lattice_->getNumRNGcall_0();
      iniLattice.numRNGcallLONGMAX = lattice_->getNumRNGcallLONGMAX();
      iniLattice.lastRNG = lattice_->getLastRNG();

      int phId = 0;
      int changeLattice = -100;
      int whileCount = 0;

      vector<int> numSitesNotAvailable;
      numSitesNotAvailable.clear();
      vector<int> vectPhIdDiff;
      vectPhIdDiff.clear();
      vector<string> vectPhNameDiff;
      vectPhNameDiff.clear();
      nucStopPrg = false;

      changeLattice = lattice_->changeMicrostructure(
          currTime, simType_, numSitesNotAvailable, vectPhIdDiff,
          vectPhNameDiff, whileCount, nucStopPrg, cyc);

      // if not all the voxels requested by KM/GEM for a certain microphase
      //  phDiff (DCId)can be dissolved because of the system configuration
      //  (lattice):
      //   - comeback to the initial system configuration : iniLattice
      //   - re-run GEM with restrictions impossed by the system configuration
      //       (DC distribution on the lattice sites) i.e. the primal solution
      //        must contain a number of moles corresponding to
      //        numSitesNotAvailable ("numDiff" in
      //        lattice_->changeMicrostructure" - sites that cannot be
      //        dissolved) lattice sites for the microphase phDiff

      if (changeLattice == 0) {

        timesGEMFailed_recall = -1;
        bool testDiff = true;
        int numSitesNotAvailableSize;
        while (changeLattice == 0) { // - for many phases!
          whileCount++;
          numSitesNotAvailableSize = numSitesNotAvailable.size();
          std::clog << endl
                    << "Controller::doCycle - cyc = " << cyc
                    << " :  changeLattice = " << changeLattice
                    << "  =>  whileCount = " << whileCount << endl;

          while (timesGEMFailed_recall != 0) {

            // reset for ChemicalSystem:
            for (int ij = 0; ij < numDCs_; ij++) {
              chemSys_->setDCMoles(ij, iniLattice.DCMoles[ij]);
            }

            // reset for Lattice:
            lattice_->setCount(iniLattice.count);
            for (int ij = 0; ij < numSites_; ij++) {
              (lattice_->getSite(ij))
                  ->setMicroPhaseId(iniLattice.site[ij].microPhaseId);
              (lattice_->getSite(ij))
                  ->setGrowthPhases(iniLattice.site[ij].growth);
              (lattice_->getSite(ij))->setWmc(iniLattice.site[ij].wmc);
              (lattice_->getSite(ij))->setWmc0(iniLattice.site[ij].wmc0);
              (lattice_->getSite(ij))
                  ->setVisit(iniLattice.site[ij].visit); // or 0!
              (lattice_->getSite(ij))
                  ->setInGrowInterfacePosVector(
                      iniLattice.site[ij].inGrowInterfacePos);
              (lattice_->getSite(ij))
                  ->setInDissInterfacePos(
                      iniLattice.site[ij].inDissInterfacePos);
            }
            for (int ij = 0; ij < dimLatticeInterface; ij++) {
              lattice_->setGrowthSites(ij,
                                       iniLattice.interface[ij].growthSites);
              lattice_->setDissolutionSites(
                  ij, iniLattice.interface[ij].dissolutionSites);
            }
            lattice_->setGrowthInterfaceSize(iniLattice.growthInterfaceSize);
            lattice_->setDissolutionInterfaceSize(
                iniLattice.dissolutionInterfaceSize);
            lattice_->resetRNG(iniLattice.numRNGcall_0,
                               iniLattice.numRNGcallLONGMAX,
                               iniLattice.lastRNG);

            std::clog << "Controller::doCycle - cyc = " << cyc
                      << " :  reset system OK & GEM_run recall for "
                         "i/whileCount/numSitesNotAvailable.size() = "
                      << i << " / " << whileCount << " / "
                      << numSitesNotAvailableSize << endl;
            std::clog << "Controller::doCycle - cyc = " << cyc
                      << " :  reset DCLowerLimits :" << endl;

            // updateDCId.clear();
            // updatePHId.clear();
            for (int ij = 0; ij < numSitesNotAvailableSize; ij++) {

              phId = vectPhIdDiff[ij]; // microPhase Id

              if (chemSys_->isKinetic(phId) && (time_[i] < beginAttackTime_)) {
                // each KC microPhase must correspond to a single GEM phase and
                // more important: to a single DC !!! attention to bassanite!!!

                vfracDiff = (static_cast<double>(numSitesNotAvailable[ij])) /
                            (static_cast<double>(numSites_));

                DCId = chemSys_->getMicroPhaseDCMembers(phId, 0);
                // updateDCId.push_back(DCId);
                // updatePHId.push_back(phId);

                volMolDiff = chemSys_->getDCMolarVolume(DCId);  // m3/mol
                molarMassDiff = chemSys_->getDCMolarMass(DCId); // g/mol

                microPhaseMassDiff =
                    vfracDiff * molarMassDiff / volMolDiff / 1.0e6; // g/cm3

                scaledMassDiff =
                    microPhaseMassDiff * 100.0 / lattice_->getInitSolidMass();

                kineticController_->updateKineticStep(cyc, phId, scaledMassDiff,
                                                      dissTimeInterval);
              } else {
                // to a nKC microPhase can correspond one or more GEM phases so,
                // one or more DCs!

                std::clog
                    << endl
                    << "    Controller::doCycle - not a KM phase - for cyc = "
                    << cyc << " & phaseId = " << phId << " ["
                    << chemSys_->getMicroPhaseName(phId) << " / DCId(phId,0):"
                    << chemSys_->getMicroPhaseDCMembers(phId, 0) << "]" << endl;

                double numTotSites_phId = lattice_->getCount()[phId];
                double numTotMoles = 0;
                int numTotCompNotZero = 0;
                vector<int> compDC = chemSys_->getMicroPhaseDCMembers(phId);
                int size = static_cast<int>(compDC.size());

                std::clog << "      DCs components:" << endl;
                for (int k = 0; k < size; k++) {
                  numMolesDiff = numSitesNotAvailable[ij] *
                                 chemSys_->getDCMoles(compDC[k]) /
                                 numTotSites_phId;

                  chemSys_->setDCLowerLimit(compDC[k], numMolesDiff);
                  chemSys_->setDCUpperLimit(compDC[k], 1.0e6); // numMolesDiff);
                  // updateDCId.push_back(compDC[k]);
                  // updatePHId.push_back(phId);

                  std::clog
                      << "          DCId = " << std::setw(3) << std::right
                      << compDC[k] << "   DCName = " << std::setw(15)
                      << std::left << chemSys_->getDCName(compDC[k])
                      << "   DCMoles = " << chemSys_->getDCMoles(compDC[k])
                      << "   DCLowerLimit = "
                      << chemSys_->getDCLowerLimit(compDC[k]) << endl;

                  numTotMoles += chemSys_->getDCMoles(compDC[k]);
                  if (chemSys_->getDCMoles(compDC[k]) > 0) {
                    numTotCompNotZero++;
                  }
                }

                std::clog << "        numTotCompNotZero = " << numTotCompNotZero
                          << "    numTotMoles = " << numTotMoles << endl;
              }
            }

            std::clog << endl;

            // timesGEMFailed_recall =
            //     chemSys_->calculateState(currTime, updateDCId, updatePHId,
            //                              isFirst, cyc);
            timesGEMFailed_recall =
                chemSys_->calculateState(currTime, isFirst, cyc);

            std::clog << endl
                      << "Controller::doCycle - cyc = " << cyc
                      << " :  i/time[i]/currTime/getTimesGEMFailed_recall = "
                      << i << " / " << time_[i] << " / " << currTime << " / "
                      << timesGEMFailed_recall << endl;

            if (timesGEMFailed_recall > 0) {
              std::clog
                  << "  Controller::doCycle - GEM_run failed for whileCount = "
                  << whileCount << endl;
              // std::clog.flush();
              timesGEMFailed_loc = timesGEMFailed_recall;
              testDiff = false;
              for (int iii = 0; iii < numSitesNotAvailableSize; iii++) {
                phId = vectPhIdDiff[iii];
                if (lattice_->getCount(phId) > numSitesNotAvailable[iii]) {
                  std::clog
                      << "  Controller::doCycle - for i/cyc/phId/iii = " << i
                      << " / " << cyc << " / " << phId << " / " << iii
                      << " : ini numSitesNotAvailable[" << iii
                      << "] = " << numSitesNotAvailable[iii];
                  numSitesNotAvailable[iii]++;
                  std::clog << "   =>   fin numSitesNotAvailable[" << iii
                            << "] = " << numSitesNotAvailable[iii] << endl;
                  // std::clog.flush();
                  testDiff = true;
                  break;
                } else {
                  vector<int> compDC = chemSys_->getMicroPhaseDCMembers(phId);
                  int size = static_cast<int>(compDC.size());

                  std::clog
                      << "  Controller::doCycle - for i/cyc/phId/iii = " << i
                      << " / " << cyc << " / " << phId << " / " << iii
                      << " - keepDCLowerLimit :" << endl;

                  for (int k = 0; k < size; k++) {
                    chemSys_->setKeepDCLowerLimit(compDC[k]);

                    std::clog
                        << "          DCId = " << std::setw(3) << std::right
                        << compDC[k] << "   DCName = " << std::setw(15)
                        << std::left << chemSys_->getDCName(compDC[k])
                        << "   DCMoles = " << chemSys_->getDCMoles(compDC[k])
                        << "   DCLowerLimit = "
                        << chemSys_->getDCLowerLimit(compDC[k])
                        << "   keepDCLowerLimit = "
                        << chemSys_->getKeepDCLowerLimit(compDC[k]) << endl;
                  }
                }
              }

              std::clog.flush();
              if (!testDiff) {
                std::clog
                    << "Controller::doCycle - do not update the microstructure "
                    << endl;
                std::clog.flush();
                break;
              }

            } else {
              timesGEMFailed_loc = timesGEMFailed_recall;
              std::clog << "Controller::doCycle - cyc = " << cyc
                        << " :  GEM_run OK for whileCount = " << whileCount
                        << endl;
              testDiff = true;
            }
          }

          if (testDiff) {

            // numSitesNotAvailable.clear(); // clear in
            // changeMicrostructure(...) vectPhIdDiff.clear();         // clear
            // in changeMicrostructure(...) vectPhNameDiff.clear();       //
            // clear in changeMicrostructure(...)
            nucStopPrg = false;
            changeLattice = lattice_->changeMicrostructure(
                currTime, simType_, numSitesNotAvailable, vectPhIdDiff,
                vectPhNameDiff, whileCount, nucStopPrg, cyc);
            std::clog << endl
                      << "Controller::doCycle - cyc = " << cyc
                      << "  &  whileCount = " << whileCount
                      << "  :  timesGEMFailed_recall = "
                      << timesGEMFailed_recall
                      << "  &  changeLattice = " << changeLattice << endl;
            if (timesGEMFailed_recall == 0 && changeLattice == 0) {
              std::clog << "     => redo adjustment for cyc/whileCount = "
                        << cyc << " / " << whileCount << endl;
              timesGEMFailed_recall = -1;
              // testDiff = true;
            }
          } else {
            break;
          }
        } // while (changeLattice == 0) { // - for many phases!

        if (timesGEMFailed_loc > 0) {
          continue;
        } else {
          kineticController_->setHydTimeIni(currTime);
          std::clog
              << endl
              << "Controller::doCycle => normal end after reset system for "
                 "cyc = "
              << cyc << " (i = " << i << ")" << endl;
        }

      } else {

        kineticController_->setHydTimeIni(currTime);

        std::clog << endl
                  << "Controller::doCycle - hydration & lattice update - cyc = "
                  << cyc << " (i = " << i
                  << "  &  time_.size()/time_[size - 1] = " << time_.size()
                  << " / " << time_[time_.size() - 1] << ")   =>   normal end"
                  << endl;
      }

    } catch (DataException dex) {

      string curTimeString = getTimeString(currTime);
      dex.printException();
      lattice_->writeLattice(curTimeString);
      lattice_->writeLatticePNG(curTimeString);

      if (xyzFiles_)
        lattice_->writeLatticeXYZ(currTime, curTimeString);

      if (xyzMovie_)
        lattice_->appendXYZ(currTime);

      throw dex;
    } catch (EOBException ex) {

      string curTimeString = getTimeString(currTime);
      ex.printException();
      lattice_->writeLattice(curTimeString);
      lattice_->writeLatticePNG(curTimeString);

      if (xyzFiles_)
        lattice_->writeLatticeXYZ(currTime, curTimeString);

      if (xyzMovie_)
        lattice_->appendXYZ(currTime);

      throw ex;
    } catch (MicrostructureException mex) {

      string curTimeString = getTimeString(currTime);
      std::clog
          << endl
          << "Controller::doCycle MicroEx from Lattice::changeMicrostructure "
             "- cyc = "
          << cyc << endl;
      mex.printException();
      lattice_->writeLattice(curTimeString);
      lattice_->writeLatticePNG(curTimeString);

      if (xyzFiles_)
        lattice_->writeLatticeXYZ(currTime, curTimeString);

      if (xyzMovie_)
        lattice_->appendXYZ(currTime);

      // write output .txt files
      writeTxtOutputFiles(currTime);
      if (writeICsDCs)
        writeTxtOutputFiles_onlyICsDCs(currTime);

      throw mex;
    }

    /// After one cycle we have removed artificial phase information
    /// introduced by the GEMS data files, so we can now calculate
    /// the initial system volume according to GEMS
    if (isFirst)
      chemSys_->setInitGEMVolume();

    // write output .txt files
    writeTxtOutputFiles(currTime);

    if (writeICsDCs)
      writeTxtOutputFiles_onlyICsDCs(currTime);

    ///
    /// Calculate the pore size distribution and saturation
    ///

    // lattice_->calculatePoreSizeDistribution();

    // thrTimeToWriteLattice threshold ~ 1 minute i.e 0.0167 hours
    if ((timeIndexIMG < static_cast<int>(outputImageTime_.size())) &&
        ((currTime >= outputImageTime_[timeIndexIMG]) ||
         (abs(currTime - outputImageTime_[timeIndexIMG]) <
          thrTimeToWriteLattice))) {

      double writeTime = currTime;
      if (abs(currTime - outputImageTime_[timeIndexIMG]) <
          thrTimeToWriteLattice)
        writeTime = outputImageTime_[timeIndexIMG];

      string curTimeString = getTimeString(writeTime);

      // if (verbose_)
      std::clog << endl
                << "Controller::doCycle - write microstructure files at i = "
                << i << " / currTime = " << currTime << " / outputImageTime_["
                << timeIndexIMG << "] = " << outputImageTime_[timeIndexIMG]
                << ", writeTime = " << writeTime << endl;

      lattice_->writeLattice(curTimeString);
      lattice_->writeLatticePNG(curTimeString);

      if (xyzFiles_)
        lattice_->writeLatticeXYZ(writeTime, curTimeString);

      if (xyzMovie_)
        lattice_->appendXYZ(writeTime);

      lattice_->calculatePoreSizeDistribution();
      lattice_->writePoreSizeDistribution(writeTime, curTimeString);

      timeIndexIMG++;
    }

    ///
    /// Check if there is any voxel-scale pore water remaining.  If not then
    /// we ASSUME hydration has stopped.
    ///
    /// @todo Generalize this idea to allow nanopore water to react by taking
    /// into account its lower chemical potential.

    // double watervolume = chemSys_->getMicroPhaseVolume(ELECTROLYTEID);

    // if (watervolume < 2.0e-18) { // Units in m3, so this is about two voxels,
    //   // we will stop hydration
    //   if (warning_) {
    //     cout << "Controller::doCycle WARNING: System is out of voxel-scale
    //     pore
    //     "
    //             "water."
    //          << endl;
    //     cout << "Controller::doCycle          This version of code assumes "
    //             "that only voxel-scale"
    //          << endl;
    //     cout << "Controller::doCycle          water is chemically reactive,
    //     so "
    //             "the system is"
    //          << endl;
    //     cout << "Controller::doCycle          is assumed to be incapable of "
    //             "further hydration."
    //          << endl;
    //     cout.flush();
    //   }
    // }

    ///
    /// The following block executes only for sulfate attack simulations
    ///

    if (currTime >
        sulfateAttackTime_) { // check! change time[i] with currTime!!!

      std::clog
          << endl
          << endl
          << "Controller::doCycle - STRAIN module (sulfate attack) : cyc = "
          << cyc << "   i = " << i
          << "   sulfateAttackTime_ = " << sulfateAttackTime_
          << "   currTime = " << currTime
          << "   time_.size() = " << time_.size()
          << "   time_[size - 1] = " << time_[time_.size() - 1] << endl;

      // std::clog << endl
      //      << " Controller::doCycle - for sulfate attack, check conditions
      //      for "
      //         "addDissolutionSites & coordination sphere "
      //      << endl;

      if (verbose_) {
        std::clog << "Controller::doCycle STRAIN module (sulfate attack)"
                  << endl;
        std::clog.flush();
      }
      map<int, vector<double>> expansion;
      expansion = lattice_->getExpansion();

      ifstream instopexp("stopexp.dat");
      if (!instopexp) {
        // if (verbose_)
        // std::clog << endl
        //      << "Controller::doCycle - STRAIN module (sulfate attack) : cyc =
        //      " << cyc
        //      << "   =>   expansion.size() = " << expansion.size() << endl;
        double maxExpVal = -1.0e6;
        double minExpVal = 1.0e6;
        int expindex;
        int maxExpInd = -1, minExpInd = -1;
        vector<double> expanval;
        for (map<int, vector<double>>::iterator it = expansion.begin();
             it != expansion.end(); it++) {
          expindex = it->first;
          expanval = it->second;
          if (expanval[0] > maxExpVal) {
            maxExpVal = expanval[0];
            maxExpInd = expindex;
          }
          if (expanval[0] < minExpVal) {
            minExpVal = expanval[0];
            minExpInd = expindex;
          }
        }
        std::clog
            << endl
            << "Controller::doCycle - STRAIN module (sulfate attack) : cyc = "
            << cyc << "   =>   expansion.size() = " << expansion.size()
            << "   minExpInd/minExpVal = " << minExpInd << " / " << minExpVal
            << "   maxExpInd/maxExpVal = " << maxExpInd << " / " << maxExpVal
            << endl;

      } else {
        expansion.clear();
        std::clog << endl
                  << "Controller::doCycle - STRAIN module (sulfate attack) = "
                  << cyc
                  << "   =>   expansion has been stopped due to the "
                     "percolation of damage"
                  << endl;
      }
      std::clog.flush();

      ///
      /// Stop FM temporarily
      /// @todo What is this?  The following if block will never be run if
      /// uncommented!
      ///

      // expansion.clear();

      newDamageCount_ = 0;

      std::clog << "  cyc = " << cyc
                << " -> damaged sites before set damage :  "
                   "oldDamageCount_/newDamageCount_/allDamageCount_ = "
                << oldDamageCount_ << " / " << newDamageCount_ << " / "
                << allDamageCount_ << endl;

      // /* // SA!!!
      if (expansion.size() > 0) {
        // double strxx, stryy, strzz;
        vector<double> locEleStress;
        double locTstrength;

        // aliId => aliId_ etc
        // int aliId = chemSys_->getMicroPhaseId_SA("Alite");
        // int belId = chemSys_->getMicroPhaseId_SA("Belite");
        // int aluId = chemSys_->getMicroPhaseId_SA("Aluminate");
        // int ferId = chemSys_->getMicroPhaseId_SA("Ferrite");
        // vector<int> isParrotKilloh_ = chemSys_->getIsParrotKilloh();
        // int sizePK_ = isParrotKilloh_.size();
        bool notPKPhase = true;

        int oldDamageCount = 0;
        // int newDamageCount = 0;
        // double poreintroduce = 0.5;

        if (verbose_) {
          std::clog
              << "Controller::doCycle STRAIN module (sulfate attack) writing "
              << endl;
          std::clog << "Controller::doCycle lattice at time_[" << i
                    << "] = " << time_[i] << ", currTime = " << currTime
                    << endl;
          std::clog << "controller::doCycle outputImageTime_[" << timeIndexIMG
                    << "] = " << outputImageTime_[timeIndexIMG] << endl;
          std::clog.flush();
        }

        // string ofileName(jobRoot_);

        // ostringstream ostrT;
        // ostrT << setprecision(3) << temperature_;
        // string tempstr(ostrT.str());

        // ostringstream ostrY, ostrD, ostrH, ostrM;
        // ostrY << setfill('0') << std::setw(3) << formattedTime.years;
        // string timestrY(ostrY.str());
        // ostrD << setfill('0') << std::setw(3) << formattedTime.days;
        // string timestrD(ostrD.str());
        // ostrH << setfill('0') << std::setw(2) << formattedTime.hours;
        // string timestrH(ostrH.str());
        // ostrM << setfill('0') << std::setw(2) << formattedTime.minutes;
        // string timestrM(ostrM.str());

        // string timeString = timestrY + "y" + timestrD + "d" +
        //                     timestrH + "h" + timestrM + "m";

        string curTimeString = getTimeString(currTime);
        // ofileName = ofileName + "." + curTimeString + "." + tempstr +
        // "K.img";

        ///
        /// In the sulfate attack algorithm, calculate the stress and strain
        /// distributions
        ///

        thermalstr_->setEigen();

        int expindex;
        vector<double> expanval;
        vector<int> expcoordin;

        for (map<int, vector<double>>::iterator it = expansion.begin();
             it != expansion.end(); it++) {

          expindex = it->first;
          expanval = it->second;
          // vector<int> expcoordin = lattice_->getExpansionCoordin(expindex);
          expcoordin = lattice_->getSite(expindex)->getXYZ();
          thermalstr_->setEigen(expindex, expanval[0], expanval[1], expanval[2],
                                0.0, 0.0, 0.0);
          thermalstr_->setExpansionCoord(expindex, expcoordin);
        }

        ///
        /// Calculate the stress-free strain (thermal strain) state in the
        /// microstructure, and then write the displacement field
        ///

        // thermalstr_->Calc(time_[i], ofileName, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        vector<int> sitePhIdVect = lattice_->getAllSitesPhId();
        vector<int> *p_sitePhIdVect = &sitePhIdVect;

        thermalstr_->calc(cyc, currTime, p_sitePhIdVect, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0);

        // thermalstr_ -> writeStress(jobRoot_,time_[i],0); //write strxx
        // thermalstr_ -> writeStrainEngy(jobRoot_,time_[i]);

        // thermalstr_->writeDisp(jobRoot_, time_[i]);
        thermalstr_->writeDisp(jobRoot_, curTimeString);

        ///
        /// Get the true volume of each voxel after FEM calculation
        ///

        //// truevolume not used!!!
        // double otruevolume = 0.0;
        // double truevolume = 0.0;
        // for (int ii = 0; ii < numSites_; ii++) {
        //   Site *ste;
        //   ste = lattice_->getSite(ii);
        //   otruevolume = ste->getTrueVolume();
        //   for (int j = 0; j < 3; j++) {
        //     truevolume += thermalstr_->getEleStrain(ii, j);
        //   }
        //   truevolume = otruevolume * (1 + truevolume);
        //   ste->setTrueVolume(truevolume);
        // }

        // double dwmcval = poreintroduce;
        double poreincrease = 0.2;
        double damageexp_0 = 1.0 / 3.0 * poreincrease;
        vector<double> damageexp(3, damageexp_0);
        vector<double> damageexpo;
        Site *ste; // , *stenb;
        int pid;

        for (int index = 0; index < numSites_; index++) {

          ste = lattice_->getSite(index);
          pid = ste->getMicroPhaseId();

          if ((ste->IsDamage())) {

            oldDamageCount++;
            // std::clog << endl << "SA-test: oldDamageCount_ = " <<
            // oldDamageCount_
            //      << "  index = " << index << "  pid = " << pid << endl;

            // double strxx, stryy, strzz;
            // strxx = stryy = strzz = 0.0;
            // strxx = thermalstr_->getEleStress(index, 0);
            // stryy = thermalstr_->getEleStress(index, 1);
            // strzz = thermalstr_->getEleStress(index, 2);
            locEleStress = thermalstr_->getEleStress(index);
            // strxx = locEleStress[0];
            // stryy = locEleStress[1];
            // strzz = locEleStress[2];
            // if ((strxx >= 1.0) || (stryy >= 1.0) || (strzz >= 1.0)) {
            if ((locEleStress[0] >= 1.0) || (locEleStress[1] >= 1.0) ||
                (locEleStress[2] >= 1.0)) {
              damageexpo = lattice_->getExpansion(index);
              damageexpo[0] += damageexp_0;
              damageexpo[1] += damageexp_0;
              damageexpo[2] += damageexp_0;
              lattice_->setExpansion(index, damageexpo);

              ///
              /// Set expansion site to be damaged if there is one, as
              /// determined by the setEigen function returning every site above
              /// damage stress threshold
              ///

              // lattice_->dWaterChange(poreintroduce);

              // lattice_->setWmc0(index,dwmcval); // *
              // lattice_->dWmc(index, dwmcval);
              // for (int j = 0; j < NN_NNN; j++) { //ste->nbSize(2)
              //   stenb = ste->nb(j);
              //   stenb->dWmc(dwmcval);
              // }
              // need interface update???

              // lattice_->dWaterChange(poreincrease);
              /// JWB: This next line must be from some earlier version
              /// The called method does not exist any longer
              ///
              /// @todo: Determine whether it is necessary to add this back
              /// in for crystallization pressure calculations
              //
              // ste->setVolume(VOIDID,(ste->getVolume(VOIDID) + poreincrease));
              //
            }
          } else {

            ///
            /// The next block gets the stress in each voxel that does NOT
            /// contain a clinker phase (C3S, C2S, C3A, or C4AF), then determine
            /// if the voxel should be damaged as a result

            /// Prefer to make this independent of whether or not there is C4AF
            /// in the phase definitions.  What if this is a white cement or
            /// something?
            ///
            /// @note Associating the last clinker phase with id 5 is a kluge
            /// @todo Give each phase a calcstress property or something like
            /// that
            ///       that can be checked instead of hardwiring phase ids

            notPKPhase = true;
            for (int i = 0; i < sizePK_; i++) {
              if (pid == isParrotKilloh_[i]) {
                notPKPhase = false;
                break;
              }
            }

            if ((pid > ELECTROLYTEID) && notPKPhase &&
                (chemSys_->isPorous(pid) || chemSys_->isWeak(pid))) {
              // double strxx, stryy, strzz;
              // strxx = stryy = strzz = 0.0;
              // strxx = thermalstr_->getEleStress(index, 0);
              // stryy = thermalstr_->getEleStress(index, 1);
              // strzz = thermalstr_->getEleStress(index, 2);
              locEleStress = thermalstr_->getEleStress(index);
              // strxx = locEleStress[0];
              // stryy = locEleStress[1];
              // strzz = locEleStress[2];
              // if ((strxx >= thermalstr_->getTstrength(index)) ||
              //     (stryy >= thermalstr_->getTstrength(index)) ||
              //     (strzz >= thermalstr_->getTstrength(index))) {
              locTstrength = thermalstr_->getTstrength(index);
              // if ((strxx >= locTstrength) ||
              //     (stryy >= locTstrength) ||
              //     (strzz >= locTstrength)) {
              if ((locEleStress[0] >= locTstrength) ||
                  (locEleStress[1] >= locTstrength) ||
                  (locEleStress[2] >= locTstrength)) {

                newDamageCount_++;
                ste->setDamage();
                lattice_->setExpansion(index, damageexp);
                // lattice_->dWaterChange(poreintroduce);

                // lattice_->setWmc0(index,dwmcval); // *
                // lattice_->dWmc(index, dwmcval);
                // for (int j = 0; j < NN_NNN; j++) {
                //   stenb = ste->nb(j);
                //   stenb->dWmc(dwmcval);
                //   if ((stenb->getWmc() > 0.0) &&
                //       (stenb->getMicroPhaseId() > ELECTROLYTEID)) {
                //     lattice_->addDissolutionSite(stenb,
                //                                  stenb->getMicroPhaseId());
                //   }
                // }

                // vector<double> damageexp;
                // damageexp.clear();
                // double poreindamage = 0.6;
                // damageexp.resize(3,(1.0 / 3.0 * poreindamage));
                // lattice_->setExpansion(index,damageexp);
                // vector<int> coordin;
                // coordin.clear();
                // coordin.resize(3,0);
                // coordin[0] = ste->getX();
                // coordin[1] = ste->getY();
                // coordin[2] = ste->getZ();
                // lattice_->setExpansionCoordin(index,coordin);
                // lattice_->dWaterChange(poreindamage);
                // ste->setVolume(VOIDID,poreindamage);
                // }
              }
            }
          }
        } // End of loop over all voxels
        if (oldDamageCount_ != oldDamageCount) {
          std::clog << endl
                    << "Controller::doCycle STRAIN module (sulfate attack) "
                       "-  error : oldDamageCount_ != "
                       "oldDamageCount"
                       " <-> cyc/oldDamageCount_/oldDamageCount = "
                    << cyc << " / " << oldDamageCount_ << " / "
                    << oldDamageCount << endl;
          std::clog << endl
                    << "         newDamageCount_/allDamageCount_ = "
                    << newDamageCount_ << " / " << allDamageCount_ << endl;
          std::clog << endl << "exit" << endl;
          exit(0);
        }

        if (newDamageCount_ > 0) { // or allDamageCount_ > 0 ? // check!
          // to see whether new damage is generated
          lattice_->writeDamageLattice(curTimeString);
          lattice_->writeDamageLatticePNG(curTimeString);
        }
      }
      // */ // SA!!!

      std::clog
          << endl
          << "Controller::doCycle - STRAIN module (sulfate attack) - cyc = "
          << cyc << " (i = " << i << ")   =>   normal end" << endl;

      allDamageCount_ = newDamageCount_ + oldDamageCount_;
      std::clog << "  cyc = " << cyc
                << " -> damaged sites after set damage  :  "
                   "oldDamageCount_/newDamageCount_/allDamageCount_ = "
                << oldDamageCount_ << " / " << newDamageCount_ << " / "
                << allDamageCount_
                << " (if newDamageCount_ = 0 => no damage file written!)"
                << endl;
      oldDamageCount_ = allDamageCount_;
    }

    if (nucStopPrg)
      break;
  }

  ///
  /// Write the final lattice state to an ASCII file and to a PNG file for
  /// visualization
  ///
  if (currTime > outputImageTime_[outputImageTime_.size() - 1] || nucStopPrg) {

    string curTimeString = getTimeString(currTime);

    lattice_->writeLattice(curTimeString);
    lattice_->writeLatticePNG(curTimeString);

    if (xyzFiles_)
      lattice_->writeLatticeXYZ(currTime, curTimeString);

    if (xyzMovie_)
      lattice_->appendXYZ(currTime);

    // if (currTime > sulfateAttackTime_ && allDamageCount_ > 0) {
    //   lattice_->writeDamageLattice(curTimeString);
    //   lattice_->writeDamageLatticePNG(curTimeString);
    // }

    lattice_->calculatePoreSizeDistribution();
    lattice_->writePoreSizeDistribution(currTime, curTimeString);

    std::clog << endl
              << "Controller::doCycle - program stops (currTime = " << currTime
              << " hours) : " << endl;
    if (nucStopPrg) {
      std::clog << endl
                << "                      number of nucleation events < nuclei "
                   "requested number!"
                << endl;
    } else {
      std::clog << endl
                << "                      outputImageTime_.size() = "
                << outputImageTime_.size() << endl;
      std::clog << "                      outputImageTime_["
                << outputImageTime_.size() - 1
                << "] = " << outputImageTime_[outputImageTime_.size() - 1]
                << endl;
      std::clog << endl
                << "                      currTime > "
                   "outputImageTime_[outputImageTime_.size() - 1]"
                   "  =>  normal end"
                << endl;
    }
  }

  return;
}

int Controller::calculateState(double &currTime, double dt, bool isFirst,
                               int cyc) {

  int timesGEMFailed = 0;

  try {

    std::clog << endl
              << "  Controller::calculateState - enter  : currTime = "
              << currTime << "   currTimeLoc = " << currTime
              << "   timesGEMFailed = " << timesGEMFailed << endl;

    ///
    /// We must pass some vectors to the `calculateKineticStep` method that
    /// will hold the amounts of impurity elements released from the clinker
    /// phases.  These values do not go in the `ChemicalSystem` object, but will
    /// still need to be processed afterward.
    ///

    // vector<double> impurityrelease;
    // impurityrelease.clear();
    // impurityrelease.resize(chemSys_->getNumMicroImpurities(), 0.0);

    ///
    /// Get the number of moles of each IC dissolved from kinetically controlled
    /// phases
    ///

    /// The thermodynamic calculation returns the saturation index of phases,
    /// which is needed for calculations of driving force for dissolution
    /// or growth.
    ///
    /// 2024-05-29:  At the moment, if a microstructure phase is defined
    /// to be one or more GEM phases, the SI of the microstructure phase
    /// is calculated as the mole-weighted average of the SIs of the
    /// constituent GEM CSD phases.  Can't think of a better way to do this
    /// except to prohibit users from defining mixtures of CSD phases
    /// as microstructure phases.

    chemSys_->initDCLowerLimit(0.0);
    chemSys_->initDCUpperLimit(1.0e6);
    kineticController_->calculateKineticStep(currTime, dt, cyc);

    ///
    /// Now that the method is done determining the change in moles of each IC,
    /// launch a thermodynamic calculation to determine new equilibrium state
    ///
    /// The `ChemicalSystem` object provides an interface for these calculations
    ///

    /// used to check if all DCMoles_ have right values according to
    /// corresponding
    ///   DCUpperLimit_ & DCUpperLimit_
    // vector<int> updateDCId, updatePHId;
    // updateDCId.clear();
    // updatePHId.clear();

    try {

      // timesGEMFailed = chemSys_->calculateState(currTime, updateDCId,
      // updatePHId,
      //                                           isFirst, cyc);
      timesGEMFailed = chemSys_->calculateState(currTime, isFirst, cyc);
      if (verbose_) {
        std::clog << "*Returned from ChemicalSystem::calculateState" << endl;
        std::clog << "*called by function Controller::calculateState" << endl;
        std::clog << "*timesGEMFailed = " << timesGEMFailed << endl;
        std::clog.flush();
      }
    } catch (GEMException gex) {
      gex.printException();
      exit(1);
    }

    ///
    /// The thermodynamic calculation returns the saturation index of phases,
    /// which is needed for calculations of driving force for dissolution
    /// or growth.  Assign this to the lattice in case crystallization pressures
    /// should be calculated.
    ///

    // Log result - retry logic is now handled by doCycle() with adaptive time stepping
    if (timesGEMFailed > 0) {
      std::clog << "  Controller::calculateState - GEMS failed : currTime = "
                << currTime << "   timesGEMFailed = " << timesGEMFailed
                << " (retry handled by adaptive time stepping in doCycle)"
                << endl;
    } else {
      std::clog << "  Controller::calculateState - OK : currTime = " << currTime
                << "   timesGEMFailed = " << timesGEMFailed << endl;
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  } catch (FloatException flex) {
    flex.printException();
    exit(1);
  }
  return timesGEMFailed;
}

void Controller::writeTxtOutputFiles(double time) {
  ///
  /// Set the kinetic DC moles.  This adds the clinker components to the DC
  /// moles.
  ///
  /// @todo Find out what this is and why it needs to be done
  ///
  ///
  // int numICs = chemSys_->getNumICs();
  // int numDCs = chemSys_->getNumDCs();
  int i;

  // Output to files the solution composition data, phase data, DC data,
  // microstructure data, pH, and C-S-H composition and Ca/Si ratio

  string outfilename = jobRoot_ + "_Solution.csv";
  ofstream outfs(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }

  outfs << setprecision(5) << time;
  char cc;
  for (i = 0; i < numDCs_; i++) {
    cc = chemSys_->getDCClassCode(i);
    if (cc == 'S' || cc == 'T' || cc == 'W') {
      outfs << "," << (chemSys_->getNode())->Get_cDC((long int)i); // molality
    }
  }
  outfs << endl;
  outfs.close();

  outfilename = jobRoot_ + "_DCVolumes.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }

  outfs << setprecision(5) << time;
  for (i = 0; i < numDCs_; i++) {
    if (chemSys_->getDCMolarMass(i) > 0.0) {
      cc = chemSys_->getDCClassCode(i);
      if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M' || cc == 'W') {
        string dcname = chemSys_->getDCName(i);
        double V0 =
            chemSys_->getDCMoles(dcname) * chemSys_->getDCMolarVolume(dcname);
        outfs << "," << V0;
      }
    } else {
      string msg = "Divide by zero error for DC " + chemSys_->getDCName(i);
      outfs.close();
      throw FloatException("Controller", "calculateState", msg);
    }
  }
  outfs << endl;
  outfs.close();

  outfilename = jobRoot_ + "_PhaseVolumes.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }

  outfs << setprecision(5) << time;
  for (i = 0; i < numGEMPhases_; i++) {
    outfs << "," << chemSys_->getGEMPhaseVolume(i);
  }
  if (time >= time_[1]) {
    double sysvol = chemSys_->getGEMVolume();
    double initsysvol = chemSys_->getInitGEMVolume();
    outfs << "," << sysvol << "," << ((initsysvol - sysvol) / initsysvol);
  } else {
    outfs << ", ,0.0";
  }
  outfs << endl;
  outfs.close();

  outfilename = jobRoot_ + "_SurfaceAreas.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "Controller", outfilename,
                        "Could not append");
  }
  outfilename = jobRoot_ + "_SI.csv";
  ofstream outfs01;
  outfs01.open(outfilename.c_str(), ios::app);
  if (!outfs01) {
    throw FileException("Controller", "Controller", outfilename,
                        "Could not append");
  }

  outfs << setprecision(5) << time;
  outfs01 << setprecision(5) << time;
  // JWB: Start with micro phase id 1 to avoid the Void phase
  for (int i = 1; i < chemSys_->getNumMicroPhases(); i++) {
    int dcid = chemSys_->getMicroPhaseDCMembers(i, 0);
    cc = chemSys_->getDCClassCode(dcid);
    if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M') {
      outfs << "," << lattice_->getSurfaceArea(i);
      ;
      outfs01 << "," << chemSys_->getMicroPhaseSI(i);
    }
  }
  outfs << endl;
  outfs01 << endl;
  outfs.close();
  outfs01.close();

  outfilename = jobRoot_ + "_Microstructure.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }

  outfs << setprecision(5) << time;
  for (i = 0; i < numMicroPhases_; i++) {
    outfs << "," << (lattice_->getVolumeFraction(i));
  }
  outfs << endl;
  outfs.close();

  outfilename = jobRoot_ + "_pH.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }
  outfs << setprecision(5) << time;
  outfs << "," << (chemSys_->getPH()) << endl;
  outfs.close();

  chemSys_->setGEMPhaseStoich();

  double *CSHcomp;
  try {
    CSHcomp = chemSys_->getPGEMPhaseStoich(chemSys_->getGEMPhaseId(CSHGEMName));
  } catch (EOBException eex) {
    eex.printException();
    exit(1);
  }
  if (verbose_) {
    std::clog << "Done!" << endl;
    std::clog.flush();
  }

  outfilename = jobRoot_ + "_CSH.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }
  outfs << setprecision(5) << time;
  for (i = 0; i < numICs_; i++) {
    outfs << "," << CSHcomp[i];
  }
  int id_Ca = chemSys_->getICId("Ca");
  int id_Si = chemSys_->getICId("Si");
  double CaMoles = CSHcomp[id_Ca];
  double SiMoles = CSHcomp[id_Si];
  if (CaMoles < 1.0e-16)
    CaMoles = 1.0e-16;
  if (SiMoles < 1.0e-16)
    SiMoles = 1.0e-16;
  double CaSiRatio = CaMoles / SiMoles;
  outfs << "," << CaSiRatio << endl;
  outfs.close();

  double *phaseRecord;
  CaMoles = SiMoles = 0.0;
  outfilename = jobRoot_ + "_CSratio_solid.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }
  outfs << setprecision(5) << time;
  for (i = 0; i < numGEMPhases_; i++) {
    cc = chemSys_->getGEMPhaseClassCode(i);
    if (cc == 's') {
      phaseRecord = chemSys_->getPGEMPhaseStoich(i);
      // ICIndex = chemSys_->getICId("Ca");
      // CaMoles += phaseRecord[ICIndex];
      // ICIndex = chemSys_->getICId("Si");
      // SiMoles += phaseRecord[ICIndex];
      CaMoles += phaseRecord[id_Ca];
      SiMoles += phaseRecord[id_Si];
    }
  }
  if (SiMoles != 0) {
    CaSiRatio = CaMoles / SiMoles;
    outfs << "," << CaSiRatio << endl;
  } else {
    outfs << ",Si_moles is ZERO" << endl;
  }
  outfs.close();

  outfilename = jobRoot_ + "_Enthalpy.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }

  double enth = 0.0;
  for (i = 0; i < numDCs_; i++) {
    enth += (chemSys_->getDCEnthalpy(i));
  }

  outfs << setprecision(5) << time;
  outfs << "," << enth << endl;
  outfs.close();

  outfilename = jobRoot_ + "_DOR.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "Controller", outfilename,
                        "Could not append");
  }
  double dor = chemSys_->getTotalDOR();
  outfs << setprecision(5) << time;
  outfs << "," << dor << endl;
  outfs.close();

  outfilename = jobRoot_ + "_GelSpaceRatio.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "Controller", outfilename,
                        "Could not append");
  }
  double gelVolume = 0, spaceVolume = 0;
  // outfs << "Time(h),DOR,GelVol,VoidElectrolyteVol,GelSpaceRatio,GSRcalc";
  // outfs << setprecision(5) << time << ",";
  outfs << setprecision(5) << time << "," << dor << ",";
  for (i = 0; i < numMicroPhases_; i++) {
    if (i < FIRST_SOLID) {
      spaceVolume += lattice_->getVolumeFraction(i);
    } else {
      if (!chemSys_->isCementComponent(i)) {
        gelVolume += lattice_->getVolumeFraction(i);
      }
    }
  }

  double gsrSim = -1.0;
  if ((spaceVolume + gelVolume) > 0.0)
    gsrSim = gelVolume / (spaceVolume + gelVolume);
  double gsrCalc = 0.68 * dor / (0.32 * dor + 0.45);
  outfs << gelVolume << "," << spaceVolume << "," << gsrSim << "," << gsrCalc
        << endl;
  outfs.close();

  outfilename = jobRoot_ + "_reactProdRatio.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "Controller", outfilename,
                        "Could not append");
  }
  outfs << setprecision(5) << time;
  for (i = 0; i < numMicroPhases_; i++) {
    if ((i >= FIRST_SOLID) && (!chemSys_->isCementComponent(i))) {
      if (gelVolume == 0) {
        outfs << ",0";
      } else {
        outfs << "," << lattice_->getVolumeFraction(i) / gelVolume;
      }
    }
  }
  outfs << endl;
  outfs.close();
}

void Controller::writeTxtOutputFiles_onlyICsDCs(double time) {

  int i, j;
  vector<double> ICMoles;
  ICMoles.resize(numICs_, 0.0);
  vector<double> DCMoles;
  DCMoles.resize(numDCs_, 0.0);
  for (i = 0; i < numDCs_; i++) {
    DCMoles[i] = chemSys_->getDCMoles(i);
  }

  string outfilenameIC = jobRoot_ + "_icmoles.csv";
  string outfilenameDC = jobRoot_ + "_dcmoles.csv";
  ofstream outfs;
  if (time < 1.e-10) {
    outfs.open(outfilenameIC.c_str());
    if (!outfs) {
      throw FileException("Controller", "writeTxtOutputFiles_onlyICsDCs",
                          outfilenameIC, "Could not append");
    }
    outfs << "Time(h)";
    for (i = 0; i < numICs_; i++) {
      outfs << "," << chemSys_->getICName(i);
    }
    outfs << endl;
    outfs.close();

    outfs.open(outfilenameDC.c_str());
    if (!outfs) {
      throw FileException("Controller", "writeTxtOutputFiles_onlyICsDCs",
                          outfilenameDC, "Could not append");
    }
    outfs << "Time(h)";
    for (i = 0; i < numDCs_; i++) {
      outfs << "," << chemSys_->getDCName(i);
    }
    outfs << endl;
    outfs.close();
  }

  vector<int> impurityDCID;
  impurityDCID.clear();
  impurityDCID.push_back(chemSys_->getDCId("K2O"));
  impurityDCID.push_back(chemSys_->getDCId("Na2O"));
  impurityDCID.push_back(chemSys_->getDCId("Per"));
  impurityDCID.push_back(chemSys_->getDCId("SO3"));

  double scMass;
  int mPhId;
  double massImpurity, totMassImpurity;

  // std::clog << endl << "getIsDCKinetic: " << endl;
  for (j = 0; j < numDCs_; j++) {
    if (chemSys_->getIsDCKinetic(j)) {
      // molMass = chemSys_->getDCMolarMass(j);
      mPhId = chemSys_->getDC_to_MPhID(j);
      scMass = chemSys_->getMicroPhaseMass(mPhId);

      totMassImpurity = 0;

      massImpurity = scMass * chemSys_->getK2o(mPhId);
      totMassImpurity += massImpurity;
      DCMoles[impurityDCID[0]] +=
          massImpurity / chemSys_->getDCMolarMass("K2O");

      massImpurity = scMass * chemSys_->getNa2o(mPhId);
      totMassImpurity += massImpurity;
      DCMoles[impurityDCID[1]] +=
          massImpurity / chemSys_->getDCMolarMass("Na2O");

      massImpurity = scMass * chemSys_->getMgo(mPhId);
      totMassImpurity += massImpurity;
      DCMoles[impurityDCID[2]] +=
          massImpurity / chemSys_->getDCMolarMass("Per"); // MgO

      massImpurity = scMass * chemSys_->getSo3(mPhId);
      totMassImpurity += massImpurity;
      DCMoles[impurityDCID[3]] +=
          massImpurity / chemSys_->getDCMolarMass("SO3");

      DCMoles[j] = (scMass - totMassImpurity) / chemSys_->getDCMolarMass(j);
    }
  }

  for (j = 0; j < numDCs_; j++) {
    for (i = 0; i < numICs_; i++) {
      ICMoles[i] += DCMoles[j] * chemSys_->getDCStoich(j, i);
    }
  }

  outfs.open(outfilenameIC.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "writeTxtOutputFiles_onlyICsDCs",
                        outfilenameIC, "Could not append");
  }

  outfs << setprecision(5) << time;
  for (i = 0; i < numICs_; i++) {
    outfs << "," << ICMoles[i];
  }
  outfs << endl;
  outfs.close();

  outfs.open(outfilenameDC.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "writeTxtOutputFiles_onlyICsDCs",
                        outfilenameDC, "Could not append");
  }

  outfs << setprecision(5) << time;
  for (i = 0; i < numDCs_; i++) {
    outfs << "," << DCMoles[i];
  }
  outfs << endl;
  outfs.close();
}

void Controller::parseDoc(const string &docName) {

  /// check if the JSON file exists

  ifstream f(docName.c_str());
  if (!f.is_open()) {
    std::clog << endl
              << "Controller::parseDoc: JSON " << docName << " file not found"
              << endl;
    throw FileException("Controller", "parseDoc", docName, "File not found");
  } else {
    std::clog << endl
              << "Controller::parseDoc: JSON " << docName
              << " file found => start reading" << endl;
  }

  /// Parse the JSON file all at once

  json data = json::parse(f);
  f.close();

  try {

    /// Get an iterator to the root node of the JSON file
    /// @todo Add a better JSON validity check.

    json::iterator it = data.find("time_parameters");
    json::iterator cdi = it.value().find("finaltime");

    if (cdi == it.value().end() || it == data.end()) {
      throw FileException("Controller", "parseDoc", docName, "Empty JSON file");
    }

    // Input times are conventionally in days
    // Immediately convert to hours within model
    double finalTime = cdi.value();
    finalTime *= (H_PER_DAY);

    /// Users may specify either individual output times
    /// or a single output frequency for simplicity.
    /// If both are specified the output times are ignored

    // Output times are conventionally in days
    // Immediately convert to hours within model

    outputImageTime_.clear();
    outputImageTimeInterval_ = 0.0;

    cdi = it.value().find("outfreq");
    double testTime = 0.0;
    if (cdi != it.value().end()) {
      outputImageTimeInterval_ = cdi.value();
      outputImageTimeInterval_ *= (H_PER_DAY);

      // Knowing the time interval, construct the output times
      if (outputImageTimeInterval_ > 1.0e-9) {
        while (testTime < finalTime) {
          testTime += outputImageTimeInterval_;
          if (testTime < finalTime) {
            outputImageTime_.push_back(testTime);
          } else {
            outputImageTime_.push_back(finalTime);
          }
        }
      }
    } else {
      cdi = it.value().find("outtimes");
      // double testTime = 0.0;
      int outtimenum = cdi.value().size();
      for (int i = 0; i < outtimenum; ++i) {
        testTime = cdi.value()[i];
        testTime *= (H_PER_DAY);
        outputImageTime_.push_back(testTime);
      }
    }

    //
    // Now populate the calculation times on a natural log scale
    // and fold in the output times in order
    time_.clear();
    testTime = 0.0;
    // time_.push_back(testTime);
    while (testTime <= finalTime) {
      // GODZILLA
      // GODZILLA We need better adaptive time stepping based on kinetic models
      // GODZILLA and rate constants
      // GODZILLA Next line is exponential
      // testTime += (0.05 * (testTime + 0.0024));
      // GODZILLA Next line is linear for a VERY fast rate
      testTime += 0.0006;
      // testTime += 0.0024;
      if (testTime >= finalTime) {
        time_.push_back(finalTime);
        break;
      }
      time_.push_back(testTime);
    }

    int outputImageTimeSize = static_cast<int>(outputImageTime_.size());
    int i_j = 0;
    for (int i = 0; i < outputImageTimeSize; i++) {
      for (int j = i_j; j < static_cast<int>(time_.size()); j++) {
        if (time_[j] >= outputImageTime_[i]) {
          if (abs(outputImageTime_[i] - time_[j]) >= 1.0e-3) {
            time_.insert(time_.begin() + j, outputImageTime_[i]);
          } else {
            time_[j] = outputImageTime_[i];
          }
          i_j = j;
          break;
        }
      }
    }

    // Read adaptive time stepping configuration if present in JSON.
    // If absent, the hard-coded defaults from the constructor are kept.
    auto asi = it.value().find("adaptive_stepping");
    if (asi != it.value().end()) {
      std::clog << "Controller::parseDoc: Found adaptive_stepping configuration"
                << endl;

      // Read the enabled flag
      auto enabledIt = asi.value().find("enabled");
      if (enabledIt != asi.value().end()) {
        useAdaptiveTimeStepping_ = enabledIt.value().get<bool>();
      }

      if (useAdaptiveTimeStepping_ && adaptiveTimeController_) {
        // Read parameters and apply to the existing adaptive controller
        AdaptiveTimeConfig config = adaptiveTimeController_->getConfig();

        auto valIt = asi.value().find("dt_initial");
        if (valIt != asi.value().end()) {
          config.dt_initial = valIt.value().get<double>();
        }

        valIt = asi.value().find("dt_max");
        if (valIt != asi.value().end()) {
          config.dt_max = valIt.value().get<double>();
        }

        valIt = asi.value().find("growth_factor");
        if (valIt != asi.value().end()) {
          config.growth_factor = valIt.value().get<double>();
        }

        valIt = asi.value().find("shrink_factor");
        if (valIt != asi.value().end()) {
          config.shrink_factor = valIt.value().get<double>();
        }

        valIt = asi.value().find("successes_for_growth");
        if (valIt != asi.value().end()) {
          config.successes_for_growth = valIt.value().get<int>();
        }

        valIt = asi.value().find("max_consecutive_failures");
        if (valIt != asi.value().end()) {
          config.max_consecutive_failures = valIt.value().get<int>();
        }

        adaptiveTimeController_->setConfig(config);
        adaptiveTimeController_->reset();

        std::clog << "  Adaptive params loaded: dt_initial=" << config.dt_initial
                  << ", dt_max=" << config.dt_max
                  << ", growth=" << config.growth_factor
                  << ", shrink=" << config.shrink_factor
                  << ", successes=" << config.successes_for_growth
                  << ", max_failures=" << config.max_consecutive_failures << endl;
      } else if (!useAdaptiveTimeStepping_) {
        std::clog << "  Adaptive time stepping DISABLED by configuration" << endl;
      }

      // max_relative_change applies regardless of adaptive stepping enabled/disabled
      auto mrcIt = asi.value().find("max_relative_change");
      if (mrcIt != asi.value().end()) {
        maxRelativeChange_ = mrcIt.value().get<double>();
        std::clog << "  maxRelativeChange set to " << (maxRelativeChange_ * 100)
                  << "%" << endl;
      }
    } else {
      std::clog << "Controller::parseDoc: No adaptive_stepping section found, "
                << "using defaults" << endl;
    }

    // Read suppressed phases list (top-level key in simparams.json).
    // For each suppressed GEMS phase, set DCUpperLimit_ to a very small
    // value to prevent GEMS from precipitating it.
    auto spIt = data.find("suppressed_phases");
    if (spIt != data.end() && spIt.value().is_array()) {
      int suppressedCount = 0;
      int skippedCount = 0;
      // Get the full DC members map once for safe lookup
      int numGEMPhases = chemSys_->getNumGEMPhases();

      // Build a complete phase-to-DC mapping from the GEMS nDCinPH array.
      // DCs are assigned sequentially: phase 0 gets DCs 0..nDCinPH[0]-1,
      // phase 1 gets the next nDCinPH[1] DCs, etc.
      auto *node = chemSys_->getNode();
      std::map<int, std::vector<int>> allPhaseDCs;
      int dcOffset = 0;
      for (int ph = 0; ph < numGEMPhases; ph++) {
        int nDC = node->pCSD()->nDCinPH[ph];
        std::vector<int> dcIds;
        for (int d = 0; d < nDC; d++) {
          dcIds.push_back(dcOffset + d);
        }
        allPhaseDCs[ph] = dcIds;
        dcOffset += nDC;
      }

      for (const auto &phaseName : spIt.value()) {
        std::string name = phaseName.get<std::string>();
        int gemPhaseIdx = chemSys_->getGEMPhaseId(name);
        if (gemPhaseIdx >= numGEMPhases) {
          skippedCount++;
          continue;
        }
        auto &dcIds = allPhaseDCs[gemPhaseIdx];
        for (int dcId : dcIds) {
          chemSys_->addSuppressedDC(dcId);
          chemSys_->setDCUpperLimit(dcId, 0.0);
        }
        suppressedCount++;
      }
      std::clog << endl
                << "Controller::parseDoc: Suppressed " << suppressedCount
                << " GEMS phases (DCUpperLimit set to 0)";
      if (skippedCount > 0)
        std::clog << " (" << skippedCount << " skipped)";
      std::clog << endl;
    }

    // There may be times associated with chemical attack
    // Next three blocks search for this
    cdi = it.value().find("beginattacktime");
    if (cdi != it.value().end()) {
      beginAttackTime_ = cdi.value();
      // beginAttackTime_ *= (H_PER_DAY);
    }

    // Input times are conventionally in days
    // Immediately convert to hours within model
    cdi = it.value().find("endattacktime");
    if (cdi != it.value().end()) {
      endAttackTime_ = cdi.value();
      // endAttackTime_ *= (H_PER_DAY);
    }

    // Input times are conventionally in days
    // Immediately convert to hours within model
    cdi = it.value().find("attacktimeinterval");
    if (cdi != it.value().end()) {
      attackTimeInterval_ = cdi.value();
      // attackTimeInterval_ *= (H_PER_DAY);
    }

    // Done searching for chemical attack times
    attack_ = false;
    bool errorAttackTime = false;
    bool errorAttackComp = false;
    if (simType_ == SULFATE_ATTACK || simType_ == LEACHING) {
      attack_ = true;
      if (beginAttackTime_ < 0 || endAttackTime_ < 0 ||
          attackTimeInterval_ < 0 || endAttackTime_ < beginAttackTime_) {
        errorAttackTime = true;
      } else {
        if (endAttackTime_ == beginAttackTime_) {
          if (attackTimeInterval_ > 0)
            errorAttackTime = true;
        } else {
          if (attackTimeInterval_ == 0)
            errorAttackTime = true;
        }
      }
      if (chemSys_->getAttackSolCompSize() <= 0)
        errorAttackComp = true;
    } else {
      if (chemSys_->getAttackSolCompSize() != 0)
        errorAttackComp = true;
    }

    std::clog << endl
              << "   simType_ = " << simType_
              << "  =>  you decided to simulate a ";
    if (simType_ == LEACHING) {
      std::clog << "LEACHING !" << endl;
    } else if (simType_ == SULFATE_ATTACK) {
      std::clog << "SULFATE ATTACK !" << endl;
    } else if (simType_ == HYDRATION) {
      std::clog << "HYDRATION !" << endl;
    }

    if (attack_) {

      if (errorAttackTime) {

        std::clog << endl
                  << endl
                  << "   ************ ERROR ************" << endl;
        std::clog << endl
                  << "   Controller::parseDoc - error : "
                  << "your time parameters are not set accordingly !!!" << endl;
        std::clog
            << endl
            << "   the specific controll parameters "
               "must fulfill some additional conditions (simparams.json) :"
            << endl;
        std::clog << "     -> beginattacktime >= 0" << endl;
        std::clog << "     -> endattacktime >= beginattacktime" << endl;
        std::clog << "     -> attacktimeinterval >= 0" << endl;
        std::clog << "          -> if endattacktime = beginattacktime => "
                     "attacktimeinterval = 0 !"
                  << endl;
        std::clog << "          -> if endattacktime > beginattacktime => "
                     "attacktimeinterval > 0 !"
                  << endl;

        // std::clog << endl
        //      << "   by default all these variables are set to -1.0" << endl;

        std::clog << endl
                  << "   for the current simulation their values are:" << endl;
        std::clog << "     -> beginattacktime    = " << beginAttackTime_
                  << endl;
        std::clog << "     -> endattacktime      = " << endAttackTime_ << endl;
        std::clog << "     -> attacktimeinterval = " << attackTimeInterval_
                  << endl;

        std::clog << endl
                  << "   => before to restart the program, please modify these "
                     "values "
                     "into simparams.json file!"
                  << endl;

        std::clog << endl
                  << "   ************ STOP! ************" << endl
                  << endl;

        throw DataException(
            "Controller", "Controller",
            "leaching or sulfate attack time parameters setting");
      }
      if (errorAttackComp) {

        std::clog << endl
                  << endl
                  << "   ************ ERROR ************" << endl;
        std::clog << endl
                  << "   Controller::parseDoc - error : "
                  << "you didn't decide the composition of attack solution !!!"
                  << endl;
        std::clog << endl
                  << "   => for LEACHING or SULFATE ATTACK, simparams.json file"
                     " must contain \"electrolyte_conditions\" information !!!"
                  << endl;
        std::clog << endl << "   THIS IS AN EXAMPLE :" << endl;

        std::clog << endl << "      \"environment\": {" << endl;
        std::clog << "        \"temperature\": 298.15," << endl;
        std::clog << "        \"reftemperature\": 298.15," << endl;
        std::clog << "        \"saturated\": 1," << endl;
        std::clog << "        \"electrolyte_conditions\": [" << endl;
        std::clog << "          { \"DCname\": \"Na+\"," << endl;
        std::clog << "            \"condition\": \"attack\"," << endl;
        std::clog << "            \"concentration\": 0.2" << endl;
        std::clog << "          }," << endl;
        std::clog << "          { \"DCname\": \"SO4-2\"," << endl;
        std::clog << "            \"condition\": \"attack\"," << endl;
        std::clog << "            \"concentration\": 0.1" << endl;
        std::clog << "          }" << endl;
        std::clog << "        ]" << endl;
        std::clog << "      }," << endl;
        std::clog << endl << "      ..." << endl;

        std::clog << endl
                  << "   => before to restart the program, please add"
                     " such information !!!"
                  << endl;

        std::clog << endl
                  << "   ************ STOP! ************" << endl
                  << endl;
        throw DataException("Controller", "Controller",
                            "sulfate attack missing information");
      }

      if (simType_ == LEACHING) {
        leachTime_ = beginAttackTime_; // not default leachTime_ = 1.0e10
        // std::clog << "leaching ";
      } else if (simType_ == SULFATE_ATTACK) {
        sulfateAttackTime_ =
            beginAttackTime_; // not default sulfateAttackTime_ = 1.0e10
        // std::clog << "sulfate attack ";
      }

      std::clog
          << endl
          << "   => you use these time parameters (only beginAttackTime_ has a"
             " time meaning!):"
          << endl;
      std::clog << "     -> beginattacktime    = " << std::setw(5) << std::right
                << static_cast<int>(beginAttackTime_) << "d / "
                << static_cast<int>(beginAttackTime_ * H_PER_DAY) << "h"
                << endl;
      std::clog << "     -> endattacktime      = " << std::setw(5) << std::right
                << static_cast<int>(endAttackTime_) << "d / "
                << static_cast<int>(endAttackTime_ * H_PER_DAY) << "h" << endl;
      std::clog << "     -> attacktimeinterval = " << std::setw(5) << std::right
                << static_cast<int>(attackTimeInterval_) << "d / "
                << static_cast<int>(attackTimeInterval_ * H_PER_DAY) << "h"
                << endl;
      std::clog << "     -> finalTime          = " << time_[time_.size() - 1]
                << " h" << endl;

      // Input times for beginAttackTime_/endAttackTime_/attackTimeInterval_
      // are conventionally in days => convert to hours within model
      beginAttackTime_ *= (H_PER_DAY);
      endAttackTime_ *= (H_PER_DAY);
      attackTimeInterval_ *= (H_PER_DAY);

      double tp = beginAttackTime_;
      int tempSize;

      tempSize = outputImageTime_.size();
      for (int i = 0; i < tempSize; i++) {
        if (outputImageTime_[i] > beginAttackTime_) {
          outputImageTime_.erase(outputImageTime_.begin() + i,
                                 outputImageTime_.begin() + tempSize);
          break;
        }
      }

      tempSize = time_.size();
      for (int i = 0; i < tempSize; i++) {
        if (time_[i] > beginAttackTime_) {
          time_.erase(time_.begin() + i, time_.begin() + tempSize);
          break;
        }
      }

      // tp = beginAttackTime_
      while (tp < endAttackTime_) {
        tp += attackTimeInterval_;
        if (tp > endAttackTime_)
          tp = endAttackTime_;
        outputImageTime_.push_back(tp);
        time_.push_back(tp);
      }

    } else if (simType_ == HYDRATION) {

      if (chemSys_->getAttackSolCompSize() != 0) {
        // attack_ is false & errorAttackComp is true
        std::clog << endl
                  << endl
                  << "   ************ ERROR ************" << endl;
        std::clog << endl
                  << "   Controller::parseDoc - error : "
                  << "your simparams.json file contains  <<\"condition\": "
                     "\"attack\">> !!!"
                  << endl;
        std::clog << endl
                  << "   => for HYDRATION, \"condition\" can have only two "
                     "values: \"initial\""
                     " and/or \"fixed\" !!!"
                  << endl;

        std::clog << endl << "   THIS IS AN EXAMPLE :" << endl;

        std::clog << endl << "      \"environment\": {" << endl;
        std::clog << "        \"temperature\": 298.15," << endl;
        std::clog << "        \"reftemperature\": 298.15," << endl;
        std::clog << "        \"saturated\": 1," << endl;
        std::clog << "        \"electrolyte_conditions\": [" << endl;
        std::clog << "          { \"DCname\": \"Ca(SO4)@\"," << endl;
        std::clog << "            \"condition\": \"initial\"," << endl;
        std::clog << "            \"concentration\": 1.0e-6" << endl;
        std::clog << "          }," << endl;
        std::clog << "          { \"DCname\": \"Na+\"," << endl;
        std::clog << "            \"condition\": \"fixed\"," << endl;
        std::clog << "            \"concentration\": 0.1" << endl;
        std::clog << "          }" << endl;
        std::clog << "          { \"DCname\": \"HCO3-\"," << endl;
        std::clog << "            \"condition\": \"fixed\"," << endl;
        std::clog << "            \"concentration\": 0.1" << endl;
        std::clog << "          }" << endl;
        std::clog << "        ]" << endl;
        std::clog << "      }," << endl;
        std::clog << "      ..." << endl;

        std::clog << endl << "   => before to restart the program: " << endl;
        std::clog << "         -> delete \"electrolyte_conditions\": [...] "
                  << endl;
        std::clog << "           or" << endl;
        std::clog << "         -> change <<\"condition\": \"attack\">>" << endl;

        std::clog << endl
                  << "   ************ STOP! ************" << endl
                  << endl;
        throw DataException("Controller", "Controller",
                            "only hydration with \"condition\": \"attack\"");
      }
      beginAttackTime_ = 1.e10;
      endAttackTime_ = 1.e10;
      attackTimeInterval_ = 1.e10;

      std::clog << endl << "   => you use these time parameters:" << endl;
      // std::clog << "             -> beginattacktime    = 1.e10 days" << endl;
      // std::clog << "             -> endattacktime      = 1.e10 days" << endl;
      // std::clog << "             -> attacktimeinterval = 1.e10 days" << endl;
      std::clog << "     -> finalTime          = " << time_[time_.size() - 1]
                << " hours" << endl;
    }

  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  return;
}

string Controller::getTimeString(const double curtime) {
  int s_per_h = static_cast<int>(S_PER_H);
  int s_per_year = static_cast<int>(S_PER_YEAR);
  int s_per_day = static_cast<int>(S_PER_DAY);
  int s_per_minute = static_cast<int>(S_PER_MINUTE);

  // Convert curtime (currently in h) into nearest second
  double curtime_in_s_dbl = curtime * S_PER_H;
  int curtime_s = static_cast<int>(curtime_in_s_dbl + 0.5);

  int years, days, hours, mins, secs;
  // How many years is this?
  years = curtime_s / s_per_year;
  curtime_s -= (years * s_per_year);
  // Convert remaining time into days
  days = curtime_s / s_per_day;
  curtime_s -= (days * s_per_day);
  // Convert remaining time into hours
  hours = curtime_s / s_per_h;
  curtime_s -= (hours * s_per_h);
  // Convert remaining time into minutes
  mins = curtime_s / s_per_minute;
  curtime_s -= (mins * s_per_minute);
  // Remaining time is seconds
  secs = curtime_s;

  ostringstream ostrY, ostrD, ostrH, ostrM, ostrS;
  ostrY << setfill('0') << std::setw(3) << years;
  string timestrY(ostrY.str());
  ostrD << setfill('0') << std::setw(3) << days;
  string timestrD(ostrD.str());
  ostrH << setfill('0') << std::setw(2) << hours;
  string timestrH(ostrH.str());
  ostrM << setfill('0') << std::setw(2) << mins;
  string timestrM(ostrM.str());
  ostrS << setfill('0') << std::setw(2) << secs;
  string timestrS(ostrS.str());

  // Include seconds in filename when non-zero for sub-minute precision
  // Format: 000y000d00h00m (backward compatible) or 000y000d00h00m00s
  string timeString =
      timestrY + "y" + timestrD + "d" + timestrH + "h" + timestrM + "m";
  if (secs > 0) {
    timeString += timestrS + "s";
  }

  return timeString;
}

void Controller::writeExitStatus(int exit_code, const std::string &exit_reason,
                                  const std::string &diagnostics) {
  //
  // Write exit status to JSON file for UI consumption.
  // The file is written to the Result directory alongside other output files.
  //
  std::string filename = "Result/exit_status.json";
  std::ofstream ofs(filename.c_str());

  if (!ofs) {
    std::clog << "WARNING: Could not write exit_status.json" << std::endl;
    return;
  }

  // Get current timestamp
  time_t now = time(nullptr);
  char timestamp[64];
  strftime(timestamp, sizeof(timestamp), "%c", localtime(&now));

  // Write JSON format
  ofs << "{" << std::endl;
  ofs << "  \"exit_code\": " << exit_code << "," << std::endl;
  ofs << "  \"exit_reason\": \"" << exit_reason << "\"," << std::endl;
  ofs << "  \"success\": " << (exit_code == 0 ? "true" : "false") << ","
      << std::endl;
  ofs << "  \"timestamp\": \"" << timestamp << "\"";

  if (!diagnostics.empty()) {
    ofs << "," << std::endl;
    ofs << "  \"diagnostics\": \"" << diagnostics << "\"" << std::endl;
  } else {
    ofs << std::endl;
  }

  ofs << "}" << std::endl;
  ofs.close();

  std::clog << "Exit status written to " << filename << ": "
            << (exit_code == 0 ? "SUCCESS" : "FAILURE") << " - " << exit_reason
            << std::endl;
}
