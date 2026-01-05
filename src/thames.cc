/**
@file thames.cc

*/

#include "thames.h"

using std::cin;
using std::cout;
using std::endl;
using std::string;

/**
@brief The main block for running THAMES.

@return 0 on successful completion, non-zero otherwise
*/
int main(int argc, char **argv) {

  // Set up logfile name
  const char logFileName[] = "thames.log";

  // Instantiate the object that manages the redirection of std::clog output to
  // file

  lgf::init_helper log_to_file(logFileName);

  // Check command line arguments

  string OutputFolder;
  if (checkArgs(argc, argv, OutputFolder)) {
    exit(0);
  }

  // Set up the strainenergy vector.  We allow no more than 156 phases,
  // but this can be changed below.

  // strainenergy.clear();
  // strainenergy.resize(156, 0.0);

  int simtype;
  string buff = "";
  ChemicalSystem *ChemSys = NULL;
  Lattice *Mic = NULL;
  ThermalStrain *ThermalStrainSolver = NULL;
  AppliedStrain *AppliedStrainSolver = NULL;
  KineticController *KController = NULL;
  Controller *Ctrl = NULL;
  RanGen *RNG = NULL;

  bool errorProgram = false;

  //
  // testSimparamsFile is used in ChemicalSystem to check (or not) if a DC
  //   identified by GEMS exists in simparams.json input file (true/false)
  //
  bool testSimparamsFile = false;

  //
  // Main menu where user decides what kind of simulation this will be.
  //

  std::clog << endl << "Enter simulation type: " << endl;
  std::clog << "  " << QUIT_PROGRAM << ") Exit program " << endl;
  std::clog << "  " << HYDRATION << ") Hydration " << endl;
  std::clog << "  " << LEACHING << ") Leaching " << endl;
  std::clog << "  " << SULFATE_ATTACK << ") Sulfate attack " << endl;
  std::clog << "  " << ELASTIC_CALC << ") Elastic calculation " << endl;
  cin >> simtype;

  if (simtype == SULFATE_ATTACK) {
    std::clog << endl
              << "simtype = " << simtype << " (Sulfate attack / SA)" << endl;
  } else if (simtype == ELASTIC_CALC) {
    std::clog << endl
              << "simtype = " << simtype << " (Elastic calculation)" << endl;
  } else {
    std::clog << endl << "simtype = " << simtype << endl;
  }

  // std::clog << "epsilon for double : \t" << numeric_limits<double>::epsilon()
  // << endl; std::clog << "epsilon for int : \t" <<
  // numeric_limits<int>::epsilon() << endl; std::clog << "epsilon for float :
  // \t" << numeric_limits<float>::epsilon()
  // << endl;

  if (simtype <= QUIT_PROGRAM || simtype > ELASTIC_CALC) {

    std::clog << "Exiting program now." << endl << endl;
    exit(1);
  }

  int seedRNG = -25943; // -142234;
  // cin >> seedRNG;
  std::clog << endl
            << "The RNG seed is                 : seedRNG = " << seedRNG
            << endl;

  double elemTimeInterval = (24.0 * 1.e-5); // 24 factor to convert from days to
                                            // hours
  // cin >> elemTimeInterval;
  std::clog << "The elementary time interval is : elemTimeInterval = "
            << setprecision(3) << elemTimeInterval
            << " hours (used in Parrot-Killoh model)" << endl;

  double corPorCSHQ = 1.0;
  std::clog << "correction CSHQ porrosity : corPorCSHQ = " << corPorCSHQ
            << endl;

  std::clog << scientific << setprecision(15) << endl;

  time_t lt = time(NULL);
  struct tm *inittime;
  inittime = localtime(&lt);
  std::clog << asctime(inittime);
  clock_t starttime = clock();

  //
  // User must provide the name of the GEM chemical system definition (CSD) file
  // for the aqueous solution
  //

  // Read the newline character.  Wish there was a better way!
  getline(cin, buff);

  //
  // User must provide the name of the GEM CSD for the whole system
  //

  std::clog << endl << "What is the name of the GEM input file? " << endl;
  getline(cin, buff);
  const string gemInputName(buff);
  std::clog << "   - gemInputName      :  " << gemInputName << endl;
  std::clog.flush();

  // Set up the strainenergy vector dimension according to
  //  GEMS input files (...-dch.dat) - strainenergy is used in THAMES
  //  only for SA but GEMS has to know about it in any case (any simtype)
  string dchName(buff);
  int pos = dchName.find("-dat.lst");
  // std::clog << endl << "pos = " << pos << endl;
  dchName.replace(pos, pos + 7, "-dch.dat");
  std::clog << endl << "     - dchName = " << dchName << endl;
  ifstream f(dchName.c_str());
  int nDC = -1;
  if (!f.is_open()) {
    std::clog << "     - file " << dchName << " not found!" << endl;
    std::clog << "exit" << endl;
    exit(0);
  } else {
    std::clog << "     - file " << dchName << " found => " << endl;
    string ch;
    while (nDC == -1) {
      f >> ch;
      if (ch == "<nDC>")
        f >> nDC;
    }
  }
  std::clog << "            number of dependent components: nDC = " << nDC
            << endl;

  strainenergy.clear();
  strainenergy.resize(nDC, 0.0);

  //
  // User must provide the name of the file specifying the microstructre
  // phase data
  //

  std::clog << endl
            << "What is the name of the simulation parameter file? " << endl;
  getline(cin, buff);
  const string simParamName(buff);
  std::clog << "   - simParamName        :  " << simParamName << endl;

  //
  // Create the ChemicalSystem object
  //

  try {
    ChemSys = new ChemicalSystem(gemInputName, simParamName, VERBOSE, WARNING,
                                 testSimparamsFile);
  } catch (bad_alloc &ba) {
    std::clog << "Bad memory allocation in ChemicalSystem constructor: "
              << ba.what() << endl;
    errorProgram = true;
  } catch (FileException fex) {
    fex.printException();
    errorProgram = true;
  } catch (GEMException gex) {
    gex.printException();
    errorProgram = true;
  } catch (DataException dex) {
    dex.printException();
    errorProgram = true;
  }
  if (errorProgram) {
    deleteDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                      AppliedStrainSolver, KController, Ctrl, starttime, lt,
                      errorProgram, OutputFolder);
  }

  ChemSys->setCorPorCSHQ(corPorCSHQ);

  //
  // Create the random number generator
  //

  RNG = new RanGen(seedRNG);

  //
  // User must specifiy the file containing the 3D microstructure itself
  //

  std::clog << endl << "What is the name of the MICROSTRUCTURE file? " << endl;
  getline(cin, buff);
  const string initMicName(buff);
  std::clog << "   - initMicName       :  " << initMicName << endl;

  //
  // Create the Lattice object to hold the microstructure
  //

  try {
    Mic = new Lattice(ChemSys, RNG, seedRNG, initMicName, VERBOSE, WARNING);
    std::clog << endl << "Lattice creation done... " << endl;
    std::clog << "X size of lattice is " << Mic->getXDim() << endl;
    std::clog << "Y size of lattice is " << Mic->getYDim() << endl;
    std::clog << "Z size of lattice is " << Mic->getZDim() << endl;
    std::clog << "Total number of sites is " << Mic->getNumSites() << endl;

  } catch (bad_alloc &ba) {
    std::clog << "Bad memory allocation in Lattice constructor: " << ba.what()
              << endl;
    std::cerr << "Bad memory allocation in Lattice constructor: " << ba.what()
              << endl;
    errorProgram = true;
  } catch (FileException fex) {
    fex.printException();
    errorProgram = true;
  } catch (GEMException gex) {
    gex.printException();
    errorProgram = true;
  } catch (EOBException eex) {
    eex.printException();
    errorProgram = true;
  } catch (FloatException flex) {
    flex.printException();
    errorProgram = true;
  } catch (DataException dex) {
    dex.printException();
    errorProgram = true;
  }

  if (errorProgram) {
    deleteDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                      AppliedStrainSolver, KController, Ctrl, starttime, lt,
                      errorProgram, OutputFolder);
  }

  if (simtype == SULFATE_ATTACK) {

    //
    // This block is executed only if simulating external sulfate attack,
    // in which case we need information about the elastic moduli of the
    // constituent phases, and will need to include a finite element solver
    //

    // std::clog << endl << "What is the name of the elastic modulus file?" <<
    // endl; buff = "";
    // // cin >> buff;  // C++ >> operator does not allow spaces
    // getline(cin, buff);
    // const string phasemod_fileName(buff);

    // std::clog << endl << "The name of the elastic modulus file is : " <<
    // endl; const string phasemod_fileName = "phasemod.txt"; std::clog << "   -
    // phasemod_fileName :  " << phasemod_fileName << endl;

    //
    // Create the ThermalStrain FE solver, which handles phase transformation
    // misfit
    //

    try {
      ThermalStrainSolver = new ThermalStrain(
          Mic->getXDim(), Mic->getYDim(), Mic->getZDim(),
          (Mic->getNumSites() + 2), ChemSys, 1, VERBOSE, WARNING);
      std::clog << "ThermalStrain object creation done... " << endl;
      // ThermalStrainSolver->setPhasemodfileName(phasemod_fileName);
    } catch (bad_alloc &ba) {
      std::clog << "Bad memory allocation in ThermalStrain constructor: "
                << ba.what() << endl;
      std::cerr << "Bad memory allocation in ThermalStrain constructor: "
                << ba.what() << endl;
      errorProgram = true;
    } catch (FileException fex) {
      fex.printException();
      errorProgram = true;
    } catch (GEMException gex) {
      gex.printException();
      errorProgram = true;
    }
    if (errorProgram) {
      deleteDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                        AppliedStrainSolver, KController, Ctrl, starttime, lt,
                        errorProgram, OutputFolder);
    }

    int nx, ny, nz;
    nx = ny = nz = 3;
    int ns = nx * ny * nz;

    //
    // Create the AppliedStrain FE solver, which handles applied external
    // strains
    //

    try {
      bool hasAggregateSlab = Mic->getHasAggregateSlab();
      AppliedStrainSolver = new AppliedStrain(
          nx, ny, nz, ns, ChemSys, hasAggregateSlab, VERBOSE, WARNING);
      std::clog << "AppliedStrain object creation done... " << endl;
      // AppliedStrainSolver->setPhasemodfileName(phasemod_fileName);
    } catch (bad_alloc &ba) {
      std::clog << "Bad memory allocation in AppliedStrain constructor: "
                << ba.what() << endl;
      errorProgram = true;
    } catch (FileException fex) {
      fex.printException();
      errorProgram = true;
    } catch (GEMException gex) {
      gex.printException();
      errorProgram = true;
    }
    if (errorProgram) {
      deleteDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                        AppliedStrainSolver, KController, Ctrl, starttime, lt,
                        errorProgram, OutputFolder);
    }

    Mic->setFEsolver(AppliedStrainSolver);
  }

  //
  // Handle standalone elastic calculation mode
  //
  if (simtype == ELASTIC_CALC) {
    std::clog << endl << "=== ELASTIC CALCULATION MODE ===" << endl;
    std::clog << "Microstructure dimensions: " << Mic->getXDim() << " x "
              << Mic->getYDim() << " x " << Mic->getZDim() << endl;

    // Create the AppliedStrain FE solver with actual microstructure dimensions
    int nx = Mic->getXDim();
    int ny = Mic->getYDim();
    int nz = Mic->getZDim();
    int ns = nx * ny * nz;

    try {
      AppliedStrainSolver =
          new AppliedStrain(nx, ny, nz, ns, ChemSys, 1, VERBOSE, WARNING);
      std::clog << "AppliedStrain object created for elastic calculation"
                << endl;
    } catch (bad_alloc &ba) {
      std::clog << "Bad memory allocation in AppliedStrain constructor: "
                << ba.what() << endl;
      errorProgram = true;
    } catch (FileException fex) {
      fex.printException();
      errorProgram = true;
    } catch (GEMException gex) {
      gex.printException();
      errorProgram = true;
    }

    if (errorProgram) {
      deleteDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                        AppliedStrainSolver, KController, Ctrl, starttime, lt,
                        errorProgram, OutputFolder);
    }

    // Get the phase vector from the microstructure
    std::vector<int> phaseIds = Mic->getAllSitesPhId();

    // Calculate effective bulk modulus using combined loading
    // getBulkModulus applies strains (0.1, 0.1, 0.1, 0.05, 0.05, 0.05)
    // which allows extraction of both bulk AND shear moduli from single solve
    std::clog << endl << "Calculating effective elastic moduli..." << endl;
    std::clog
        << "Applying combined strain state for bulk and shear calculation..."
        << endl;
    double bulkModulus = AppliedStrainSolver->getBulkModulus(
        &phaseIds, Mic->getResolution(), OutputFolder);
    std::clog << "Effective bulk modulus: " << bulkModulus << " GPa" << endl;

    // Extract shear modulus from the SAME solution (no second solve needed!)
    // The combined loading includes shear strains (exz=eyz=exy=0.05)
    // Average the three shear components as VCCTL does:
    // shear = (strxy/sxy + strxz/sxz + stryz/syz) / 3
    double strxy = AppliedStrainSolver->getStress(5); // xy stress
    double strxz = AppliedStrainSolver->getStress(3); // xz stress
    double stryz = AppliedStrainSolver->getStress(4); // yz stress
    double sxy = AppliedStrainSolver->getStrain(5);   // xy strain
    double sxz = AppliedStrainSolver->getStrain(3);   // xz strain
    double syz = AppliedStrainSolver->getStrain(4);   // yz strain
    double shearModulus = (strxy / sxy + strxz / sxz + stryz / syz) / 3.0;
    std::clog << "Effective shear modulus: " << shearModulus << " GPa" << endl;

    // Calculate derived quantities
    double youngsModulus = 0.0;
    double poissonsRatio = 0.0;
    if (bulkModulus > 0 || shearModulus > 0) {
      double denom = 3.0 * bulkModulus + shearModulus;
      if (denom > 0) {
        youngsModulus = 9.0 * bulkModulus * shearModulus / denom;
        poissonsRatio =
            (3.0 * bulkModulus - 2.0 * shearModulus) / (2.0 * denom);
      }
    }

    std::clog << endl << "=== ELASTIC CALCULATION RESULTS ===" << endl;
    std::clog << "Effective bulk modulus K:    " << bulkModulus << " GPa"
              << endl;
    std::clog << "Effective shear modulus G:   " << shearModulus << " GPa"
              << endl;
    std::clog << "Effective Young's modulus E: " << youngsModulus << " GPa"
              << endl;
    std::clog << "Effective Poisson's ratio:   " << poissonsRatio << endl;
    std::clog << "===================================" << endl;

    // Write results to output file
    string elasticResultsFile = OutputFolder + "/EffectiveModuli.csv";
    std::clog << endl
              << "Attempting to write results to: " << elasticResultsFile
              << endl;
    std::ofstream elasticOut(elasticResultsFile);
    if (elasticOut.is_open()) {
      elasticOut << "Property,Value,Units" << endl;
      elasticOut << "Microstructure," << initMicName << "," << endl;
      elasticOut << "X_Dimension," << nx << ",voxels" << endl;
      elasticOut << "Y_Dimension," << ny << ",voxels" << endl;
      elasticOut << "Z_Dimension," << nz << ",voxels" << endl;
      elasticOut << "Resolution, " << (1.0e6 * Mic->getResolution())
                 << ",um/voxel" << endl;
      elasticOut << "Bulk_modulus," << bulkModulus << ",GPa" << endl;
      elasticOut << "Shear_modulus," << shearModulus << ",GPa" << endl;
      elasticOut << "Youngs_modulus," << youngsModulus << ",GPa" << endl;
      elasticOut << "Poissons_ratio," << poissonsRatio << "," << endl;
      elasticOut.close();
      std::clog << "Results successfully written to: " << elasticResultsFile
                << endl;
    } else {
      std::clog << "ERROR: Could not open file for writing: "
                << elasticResultsFile << endl;
      std::clog << "Check that the directory '" << OutputFolder << "' exists."
                << endl;
    }

    // Write ITZ results to output file if they are available
    if (AppliedStrainSolver->getHasAggregateSlab()) {

      string elasticITZFile = OutputFolder + "/ITZModuli.csv";
      std::clog << endl
                << "Attempting to write results to: " << elasticITZFile << endl;
      elasticOut.open(elasticITZFile);
      if (!elasticOut.is_open()) {
        std::clog << "ERROR: Could not open file for writing: "
                  << elasticITZFile << endl;
        std::clog << "Check that the directory '" << OutputFolder << "' exists."
                  << endl;
      } else {
        double resInMicrometers = (1.0e6 * Mic->getResolution());
        elasticOut << "Property,Value,Units" << endl;
        elasticOut << "Microstructure," << initMicName << "," << endl;
        elasticOut << "X_Dimension," << nx << ",voxels" << endl;
        elasticOut << "Y_Dimension," << ny << ",voxels" << endl;
        elasticOut << "Z_Dimension," << nz << ",voxels" << endl;
        elasticOut << "Resolution," << resInMicrometers << ",um/voxel" << endl;

        ///
        /// Must know the surface position of the aggregate slab first
        ///

        int aggX = Mic->getAggregateSurfacePosition();
        std::clog << "Aggregate surface position = " << aggX << " voxels"
                  << endl;
        std::clog.flush();
        double distance = (-0.5) * resInMicrometers;
        double bulkmod_i, bulkmod_ni, shearmod_i, shearmod_ni;
        double bulkModulus_at_x, shearModulus_at_x, youngsModulus_at_x,
            poissonsRatio_at_x;
        int i_inverse;
        for (int i = aggX - 1; i > 0; i--) {
          i_inverse = i - aggX + 2;
          distance += resInMicrometers;
          try {
            bulkmod_i = AppliedStrainSolver->getLayerBulkModulus(i);
            bulkmod_ni = AppliedStrainSolver->getLayerBulkModulus(nx - i - 1);
            shearmod_i = AppliedStrainSolver->getLayerShearModulus(i);
            shearmod_ni = AppliedStrainSolver->getLayerShearModulus(nx - i - 1);
            std::clog << "K[" << i << "] = " << bulkmod_i << endl;
            std::clog.flush();
            std::clog << "K[" << nx - i - 1 << "] = " << bulkmod_ni << endl;
            std::clog.flush();
            std::clog << "G[" << i << "] = " << shearmod_i << endl;
            std::clog.flush();
            std::clog << "G[" << nx - i - 1 << "] = " << shearmod_ni << endl;
            std::clog.flush();
            bulkModulus_at_x = 0.50 * (bulkmod_i + bulkmod_ni);
            shearModulus_at_x = 0.50 * (shearmod_i + shearmod_ni);
            youngsModulus_at_x = 9.0 * bulkModulus_at_x * shearModulus_at_x /
                                 (3.0 * bulkModulus_at_x + shearModulus_at_x);
            poissonsRatio_at_x =
                (3.0 * bulkModulus_at_x - 2.0 * shearModulus_at_x) /
                (2.0 * (3.0 * bulkModulus_at_x + shearModulus_at_x));
            std::clog << distance << " " << bulkModulus_at_x << " "
                      << shearModulus_at_x << endl;
            std::clog.flush();

            elasticOut << "Layer_" << i_inverse << "_distance," << distance
                       << ",um" << endl;
            elasticOut << "Layer_" << i_inverse << "_Bulk_modulus,"
                       << bulkModulus_at_x << ",GPa" << endl;
            elasticOut << "Layer_" << i_inverse << "_Shear_modulus,"
                       << shearModulus_at_x << ",GPa" << endl;
            elasticOut << "Layer_" << i_inverse << "_Youngs_modulus,"
                       << youngsModulus_at_x << ",GPa" << endl;
            elasticOut << "Layer_" << i_inverse << "_Poissons_ratio,"
                       << poissonsRatio_at_x << ",";
            if (i > 1) {
              elasticOut << endl;
            }
          } catch (const std::out_of_range &e) {
            std::cerr << "ERROR: Out of bounds exception in layered moduli:"
                      << endl;
            std::cerr << "       " << e.what() << endl;
          }
        }
        elasticOut.close();
        std::clog << "Results successfully written to: " << elasticITZFile
                  << endl;
      }
    }

    // Clean up and exit
    std::clog << endl << "Elastic calculation complete." << endl;
    deleteDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                      AppliedStrainSolver, KController, Ctrl, starttime, lt,
                      errorProgram, OutputFolder);
    return 0; // Exit here for elastic-only mode
  }

  string jobRoot, statFileName;
  if (VERBOSE) {
    std::clog << endl << "About to enter KineticController constructor" << endl;
    std::clog << "exit" << endl;
    exit(0);
    std::clog.flush();
  }

  //
  // Create the KineticController object
  //

  try {
    KController =
        new KineticController(ChemSys, Mic, simParamName, VERBOSE, WARNING);
  } catch (bad_alloc &ba) {
    std::cerr << "Bad memory allocation in KineticController constructor: "
              << ba.what() << endl;
    errorProgram = true;
  } catch (FileException fex) {
    fex.printException();
    errorProgram = true;
  } catch (GEMException gex) {
    gex.printException();
    errorProgram = true;
  } catch (FloatException flex) {
    flex.printException();
    errorProgram = true;
  } catch (DataException dex) {
    dex.printException();
    errorProgram = true;
  }
  if (errorProgram) {
    deleteDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                      AppliedStrainSolver, KController, Ctrl, starttime, lt,
                      errorProgram, OutputFolder);
  }

  if (VERBOSE) {
    std::clog << "Finished constructing KineticController KController" << endl;
    std::clog.flush();
  }

  std::clog << endl
            << "What shall be the root name of all output files?" << endl;
  getline(cin, buff);
  jobRoot.assign(buff);
  std::clog << "   - files root name   :  " << jobRoot << endl;

  prepOutputFolder(OutputFolder, jobRoot, gemInputName, statFileName,
                   initMicName, simParamName);

  //
  // Create the Controller object to direct flow of the program
  //

  try {
    Ctrl =
        new Controller(Mic, KController, ChemSys, ThermalStrainSolver, simtype,
                       simParamName, jobRoot, VERBOSE, WARNING, XYZ);
  } catch (bad_alloc &ba) {
    std::cerr << "Bad memory allocation in Controller constructor: "
              << ba.what() << endl;
    errorProgram = true;
  } catch (FileException fex) {
    fex.printException();
    errorProgram = true;
  } catch (GEMException gex) {
    gex.printException();
    errorProgram = true;
  } catch (DataException dex) {
    dex.printException();
    errorProgram = true;
  }
  if (errorProgram) {
    deleteDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                      AppliedStrainSolver, KController, Ctrl, starttime, lt,
                      errorProgram, OutputFolder);
  }

  //
  // Write a formatted output of the simulation parameters for later reference
  //

  writeReport(jobRoot, inittime, initMicName, simParamName, gemInputName,
              ChemSys);

  //
  // Launch the main controller to run the simulation
  //

  if (VERBOSE) {
    std::clog << endl << "Going into Controller::doCycle" << endl;
    std::clog.flush();
  }

  try {

    Ctrl->doCycle(elemTimeInterval);

  } catch (GEMException gex) {
    gex.printException();
    errorProgram = true;
  } catch (DataException dex) {
    dex.printException();
    errorProgram = true;
  } catch (EOBException eex) {
    eex.printException();
    errorProgram = true;
  } catch (MicrostructureException mex) {
    mex.printException();
    errorProgram = true;
  }
  if (errorProgram) {
    deleteDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                      AppliedStrainSolver, KController, Ctrl, starttime, lt,
                      errorProgram, OutputFolder);
  }

  //
  // Simulation is finished.
  //

  deleteDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver, AppliedStrainSolver,
                    KController, Ctrl, starttime, lt, errorProgram,
                    OutputFolder);

  return 0;
}

void deleteDynAllocMem(ChemicalSystem *ChemSys, Lattice *Mic, RanGen *RNG,
                       ThermalStrain *ThermalStrainSolver,
                       AppliedStrain *AppliedStrainSolver,
                       KineticController *KController, Controller *Ctrl,
                       clock_t st_time, time_t lt, bool errorProgram,
                       const string &OutputFolder) {

  string buff = "";
  int resCallSystem;

  if (Ctrl) {
    delete Ctrl;
  }
  if (KController) {
    delete KController;
  }
  if (AppliedStrainSolver) {
    delete AppliedStrainSolver;
  }
  if (ThermalStrainSolver) {
    delete ThermalStrainSolver;
  }
  if (Mic) {
    writeLastCount(ChemSys, Mic);
    delete Mic;
  }
  if (ChemSys) {
    delete ChemSys;
  }
  if (RNG) {
    delete RNG;
  }

  //
  // Move remaining output files to output folder
  //

  string name = "ipmlog.txt";
  ifstream f(name.c_str());
  if (f.good()) {
    buff = "mv -f ipmlog.txt " + OutputFolder + "/.";
    resCallSystem = system(buff.c_str());
    if (resCallSystem == -1) {
      throw FileException("thames", "deleteDynAllocMem", buff, "FAILED");
    }
    // std::clog << buff << endl;
    f.close();
  }

  name = "IPM_dump.txt";
  f.open(name.c_str());
  if (f.good()) {
    buff = "mv -f IPM_dump.txt " + OutputFolder + "/.";
    resCallSystem = system(buff.c_str());
    if (resCallSystem == -1) {
      throw FileException("thames", "deleteDynAllocMem", buff, "FAILED");
    }
    // std::clog << buff << endl;
    f.close();
  }

  std::clog << endl
            << "=> IPM_dump.txt & ipmlog.txt are into " << OutputFolder
            << " folder" << endl;

  std::clog << endl << "STOP Program" << endl;
  timeCount(st_time, lt);

  if (errorProgram) {
    exit(1);
  } else {
    exit(0);
  }
}

void timeCount(clock_t time_, time_t lt_) {

  time_t lt1 = time(NULL);
  struct tm *inittime1;
  inittime1 = localtime(&lt1);
  std::clog << endl << asctime(inittime1);
  clock_t endtime = clock();

  double elapsedtime = static_cast<double>(endtime - time_) / CLOCKS_PER_SEC;
  double ltD = difftime(lt1, lt_);
  std::clog << endl << "Total time = " << ltD << " seconds" << endl;
  std::clog << endl
            << "Total time with clock = " << elapsedtime << " seconds" << endl;
}

void printHelp(void) {
  std::clog << endl;
  std::clog << "Usage: \"thames [--verbose|-v] [--suppress|-s] [--xyz|-x] "
               "[--outfolder|-o] folder [--help|-h]\""
            << endl;
  std::clog << "        --verbose [-v]      Produce verbose output" << endl;
  std::clog << "        --suppress [-s]     Suppress warning messages" << endl;
  std::clog << "        --xyz [-x]          Create 3D visualization movie"
            << endl;
  std::clog
      << "        --outfolder [-o]    Name of folder for output data (default "
         "is Result)"
      << endl;
  std::clog << endl;
  std::clog << "Note: thames --help [-h] will print this help message" << endl;

  return;
}

int checkArgs(int argc, char **argv, string &OutputFolder) {

  // Many of the variables here are defined in the getopts.h system header file
  // Can define more options here if we want

  const char *const short_opts = "vsxj:o:h";
  const option long_opts[] = {{"verbose", no_argument, nullptr, 'v'},
                              {"suppress", no_argument, nullptr, 's'},
                              {"xyz", no_argument, nullptr, 'x'},
                              {"outfolder", required_argument, nullptr, 'o'},
                              {"help", no_argument, nullptr, 'h'},
                              {nullptr, no_argument, nullptr, 0}};

  VERBOSE = false;
  WARNING = true;
  XYZ = false;

  // Default value of output folder unless user overrides it
  OutputFolder = "Result";

  while (true) {

    const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

    if (-1 == opt)
      break; // breaks out of this while loop

    switch (opt) {
    case 'v':
      VERBOSE = true; // Verbose defined in thameslib global.h
      std::clog << "**Will produce verbose output**" << endl;
      break;
    case 's':
      WARNING = false; // Verbose defined in thameslib global.h
      std::clog << "**Will suppress warning messages**" << endl;
      break;
    case 'x':
      XYZ = true; // Verbose defined in thameslib global.h
      std::clog << "**Will create 3D visualization file **" << endl;
      break;
    case 'o':
      OutputFolder = optarg; // Verbose defined in thameslib global.h
      if ((OutputFolder[0] == ' ') || (OutputFolder[0] == '-') ||
          (OutputFolder[0] == '\\')) {
        printHelp();
        return (1);
      }
      std::clog << "Output folder: " << OutputFolder << endl;
      break;
    case 'h': // -h or --help
    case '?': // Unrecognized option
    default:
      printHelp();
      return (1);
      break;
    }
  }
  return (0);
}

void prepOutputFolder(const string &OutputFolder, string &jobRoot,
                      const string &gemInputName, string &statFileName,
                      const string &initMicName, const string &simParamName) {

  int resCallSystem;

  string buff = "mkdir -p " + OutputFolder;
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }
  jobRoot = OutputFolder + "/" + jobRoot;
  std::clog << "   - jobRoot           :  " << jobRoot << endl;
  std::clog.flush();

  statFileName = jobRoot + ".stats";

  // Read the gem input master file to get file names
  // and copy them to the output folder

  ifstream in(gemInputName);
  string buff1;
  in >> buff1; // discard flag

  // DCH file
  buff1.clear();
  in >> buff1; // This is the dch file with quotes
  if (buff1[0] == '"' || buff1[0] == '\'') {
    buff1 = buff1.substr(1, buff1.size() - 2);
  }
  buff = "cp -f " + buff1 + " " + OutputFolder + "/.";
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }

  // IPM file
  buff1.clear();
  in >> buff1; // This is the ipm file with quotes
  if (buff1[0] == '"' || buff1[0] == '\'') {
    buff1 = buff1.substr(1, buff1.size() - 2);
  }
  buff = "cp -f " + buff1 + " " + OutputFolder + "/.";
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }

  // DBR file
  buff1.clear();
  in >> buff1; // This is the dbr file with quotes
  in.close();

  if (buff1[0] == '"' || buff1[0] == '\'') {
    buff1 = buff1.substr(1, buff1.size() - 2);
  }
  buff = "cp -f " + buff1 + " " + OutputFolder + "/.";
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }

  buff = "cp -f " + gemInputName + " " + OutputFolder + "/.";
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }

  buff = "cp -f " + initMicName + " " + OutputFolder + "/.";
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }

  buff = "cp -f " + simParamName + " " + OutputFolder + "/.";
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }

  std::clog << "     => All input files have been copied into " << OutputFolder
            << " folder" << endl;

  return;
}

void writeReport(const string &jobRoot, struct tm *itime,
                 const string &initMicName, const string &simParamName,
                 const string &csdName, ChemicalSystem *csys) {

  string statName = jobRoot + ".stats";
  string jFileName = jobRoot + ".report";
  ofstream out(jFileName.c_str());
  if (!out.is_open()) {
    if (WARNING)
      std::clog << "WARNING:  Could not open report file" << endl;
    std::cerr << "WARNING:  Could not open report file" << endl;
    return;
  }

  //
  // Write the time the job was executed
  //

  out << "THAMES simulation " << jobRoot;
  out << " initialized on " << asctime(itime) << endl;
  out << endl;
  out << "INPUT FILES USED:" << endl;
  out << "              Microstructure file name: " << initMicName << endl;
  out << "   Microstructure definition file name: " << simParamName << endl;
  out << "                   GEM input file name: " << csdName << endl;
  out << endl;
  out << "OUTPUT FILES GENERATED:" << endl;
  out << "                Global phase fractions: " << statName << endl;
  out << endl;
  out << "----------------------------------------------------------" << endl;
  out << endl;

  csys->writeChemSys(out);
  out << endl;

  out.close();

  return;
}

void writeLastCount(ChemicalSystem *chemsys, Lattice *mic) {
  std::vector<int> count;
  count = mic->getCount();
  int dim = count.size();
  int numPhases = chemsys->getNumMicroPhases();
  std::vector<string> phaseName = chemsys->getMicroPhaseName();
  long int sumCount = 0;

  if (dim == 0) {
    count.resize(numPhases, 0);
    dim = numPhases;
  }
  std::clog << endl
            << "   Last number of voxels for each microPhase <-> nVmPh(i) :"
            << endl;
  for (int i = 0; i < dim; i++) {
    sumCount += count[i];
    std::clog << "    " << std::setw(3) << std::right << i << " - "
              << std::setw(18) << std::left << phaseName[i] << " : "
              << std::setw(10) << std::right << count[i] << endl;
  }
  std::clog << endl << "   nVmPh_total = " << sumCount << endl;
  std::clog << "   system dimension is : dimX*dimY*dimZ = " << mic->getXDim()
            << " * " << mic->getYDim() << " * " << mic->getZDim() << " = "
            << mic->getXDim() * mic->getYDim() * mic->getZDim() << endl;
}
