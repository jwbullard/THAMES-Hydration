/**
@file thames.cc

*/

#include "thames.h"

using std::cout; using std::cin; using std::endl;
using std::string;

/**
@brief The main block for running THAMES.

@return 0 on successful completion, non-zero otherwise
*/
int main(int argc, char **argv) {

  // Check command line arguments

  string outputFolder;
  cout << scientific << setprecision(15);
  if (checkArgs(argc, argv, outputFolder)) {
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

  cout << endl << "Enter simulation type: " << endl;
  cout << "  " << QUIT_PROGRAM << ") Exit program " << endl;
  cout << "  " << HYDRATION << ") Hydration " << endl;
  cout << "  " << LEACHING << ") Leaching " << endl;
  cout << "  " << SULFATE_ATTACK << ") Sulfate attack " << endl;
  cin >> simtype;

  if (simtype > SULFATE_ATTACK) {
    cout << endl << "simtype = " << simtype << " (Sulfate attack / SA)" << endl;
  } else {
    cout << endl << "simtype = " << simtype << endl;
  }

  // cout << "epsilon for double : \t" << numeric_limits<double>::epsilon() <<
  // endl; cout << "epsilon for int : \t" << numeric_limits<int>::epsilon() <<
  // endl; cout << "epsilon for float : \t" << numeric_limits<float>::epsilon()
  // << endl;

  if (simtype <= QUIT_PROGRAM || simtype > SULFATE_ATTACK) {

    cout << "Exiting program now." << endl << endl;
    exit(1);
  }

  int seedRNG = -25943; // -142234;
  // cin >> seedRNG;
  cout << endl
       << "The RNG seed is                 : seedRNG = " << seedRNG << endl;

  double elemTimeInterval = (24.0 * 1.e-5); // 24 factor to convert from days to
                                            // hours
  // cin >> elemTimeInterval;
  cout << "The elementary time interval is : elemTimeInterval = "
       << setprecision(3) << elemTimeInterval
       << " hours (used in Parrot-Killoh model)" << endl;

  double corPorCSHQ = 1.0;
  cout << "correction CSHQ porrosity : corPorCSHQ = " << corPorCSHQ << endl;

  cout << scientific << setprecision(15) << endl;

  time_t lt = time(NULL);
  struct tm *inittime;
  inittime = localtime(&lt);
  cout << asctime(inittime);
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

  cout << endl << "What is the name of the GEM input file? " << endl;
  getline(cin, buff);
  const string gemInputName(buff);
  cout << "   - gemInputName      :  " << gemInputName << endl;
  cout.flush();

  // Set up the strainenergy vector dimension according to
  //  GEMS input files (...-dch.dat) - strainenergy is used in THAMES
  //  only for SA but GEMS has to know about it in any case (any simtype)
  string dchName(buff);
  int pos = dchName.find("-dat.lst");
  // cout << endl << "pos = " << pos << endl;
  dchName.replace(pos, pos + 7, "-dch.dat");
  cout << endl << "     - dchName = " << dchName << endl;
  ifstream f(dchName.c_str());
  int nDC = -1;
  if (!f.is_open()) {
    cout << "     - file " << dchName << " not found!" << endl;
    cout << "exit" << endl;
    exit(0);
  } else {
    cout << "     - file " << dchName << " found => " << endl;
    string ch;
    while (nDC == -1) {
      f >> ch;
      if (ch == "<nDC>")
        f >> nDC;
    }
  }
  cout << "            number of dependent components: nDC = " << nDC << endl;

  strainenergy.clear();
  strainenergy.resize(nDC, 0.0);

  //
  // User must provide the name of the file specifying the microstructre
  // phase data
  //

  cout << endl << "What is the name of the simulation parameter file? " << endl;
  getline(cin, buff);
  const string simParamName(buff);
  cout << "   - simParamName        :  " << simParamName << endl;

  //
  // Create the ChemicalSystem object
  //

  try {
    ChemSys = new ChemicalSystem(gemInputName, simParamName, VERBOSE, WARNING,
                                 testSimparamsFile);
  } catch (bad_alloc &ba) {
    cout << "Bad memory allocation in ChemicalSystem constructor: " << ba.what()
         << endl;
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
                      errorProgram, outputFolder);
  }

  ChemSys->setCorPorCSHQ(corPorCSHQ);

  //
  // Create the random number generator
  //

  RNG = new RanGen(seedRNG);

  //
  // User must specifiy the file containing the 3D microstructure itself
  //

  cout << endl << "What is the name of the MICROSTRUCTURE file? " << endl;
  getline(cin, buff);
  const string initMicName(buff);
  cout << "   - initMicName       :  " << initMicName << endl;

  //
  // Create the Lattice object to hold the microstructure
  //

  try {
    Mic = new Lattice(ChemSys, RNG, seedRNG, initMicName, VERBOSE, WARNING);
    cout << endl << "Lattice creation done... " << endl;
    cout << "X size of lattice is " << Mic->getXDim() << endl;
    cout << "Y size of lattice is " << Mic->getYDim() << endl;
    cout << "Z size of lattice is " << Mic->getZDim() << endl;
    cout << "Total number of sites is " << Mic->getNumSites() << endl;

  } catch (bad_alloc &ba) {
    cout << "Bad memory allocation in Lattice constructor: " << ba.what()
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
                      errorProgram, outputFolder);
  }

  if (simtype == SULFATE_ATTACK) {

    //
    // This block is executed only if simulating external sulfate attack,
    // in which case we need information about the elastic moduli of the
    // constituent phases, and will need to include a finite element solver
    //

    // cout << endl << "What is the name of the elastic modulus file?" << endl;
    // buff = "";
    // // cin >> buff;  // C++ >> operator does not allow spaces
    // getline(cin, buff);
    // const string phasemod_fileName(buff);

    // cout << endl << "The name of the elastic modulus file is : " << endl;
    // const string phasemod_fileName = "phasemod.txt";
    // cout << "   - phasemod_fileName :  " << phasemod_fileName << endl;

    //
    // Create the ThermalStrain FE solver, which handles phase transformation
    // misfit
    //

    try {
      ThermalStrainSolver = new ThermalStrain(
          Mic->getXDim(), Mic->getYDim(), Mic->getZDim(),
          (Mic->getNumSites() + 2), ChemSys, 1, VERBOSE, WARNING);
      cout << "ThermalStrain object creation done... " << endl;
      // ThermalStrainSolver->setPhasemodfileName(phasemod_fileName);
    } catch (bad_alloc &ba) {
      cout << "Bad memory allocation in ThermalStrain constructor: "
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
                        errorProgram, outputFolder);
    }

    int nx, ny, nz;
    nx = ny = nz = 3;
    int ns = nx * ny * nz;

    //
    // Create the AppliedStrain FE solver, which handles applied external
    // strains
    //

    try {
      AppliedStrainSolver =
          new AppliedStrain(nx, ny, nz, ns, ChemSys, 1, VERBOSE, WARNING);
      cout << "AppliedStrain object creation done... " << endl;
      // AppliedStrainSolver->setPhasemodfileName(phasemod_fileName);
    } catch (bad_alloc &ba) {
      cout << "Bad memory allocation in AppliedStrain constructor: "
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
                        errorProgram, outputFolder);
    }

    Mic->setFEsolver(AppliedStrainSolver);
  }

  string jobRoot, statFileName;
  if (VERBOSE) {
    cout << endl << "About to enter KineticController constructor" << endl;
    cout << "exit" << endl;
    exit(0);
    cout.flush();
  }

  //
  // Create the KineticController object
  //

  try {
    KController =
        new KineticController(ChemSys, Mic, simParamName, VERBOSE, WARNING);
  } catch (bad_alloc &ba) {
    cout << "Bad memory allocation in KineticController constructor: "
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
                      errorProgram, outputFolder);
  }

  if (VERBOSE) {
    cout << "Finished constructing KineticController KController" << endl;
    cout.flush();
  }

  cout << endl << "What shall be the root name of all output files?" << endl;
  getline(cin, buff);
  jobRoot.assign(buff);
  cout << "   - files root name   :  " << jobRoot << endl;

  prepOutputFolder(outputFolder, jobRoot, gemInputName, statFileName,
                   initMicName, simParamName);

  //
  // Create the Controller object to direct flow of the program
  //

  try {
    Ctrl =
        new Controller(Mic, KController, ChemSys, ThermalStrainSolver, simtype,
                       simParamName, jobRoot, VERBOSE, WARNING, XYZ);
  } catch (bad_alloc &ba) {
    cout << "Bad memory allocation in Controller constructor: " << ba.what()
         << endl;
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
                      errorProgram, outputFolder);
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
    cout << endl << "Going into Controller::doCycle" << endl;
    cout.flush();
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
                      errorProgram, outputFolder);
  }

  //
  // Simulation is finished.
  //

  deleteDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver, AppliedStrainSolver,
                    KController, Ctrl, starttime, lt, errorProgram,
                    outputFolder);

  return 0;
}

void deleteDynAllocMem(ChemicalSystem *ChemSys, Lattice *Mic, RanGen *RNG,
                       ThermalStrain *ThermalStrainSolver,
                       AppliedStrain *AppliedStrainSolver,
                       KineticController *KController, Controller *Ctrl,
                       clock_t st_time, time_t lt, bool errorProgram,
                       const string &outputFolder) {

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
    buff = "mv -f ipmlog.txt " + outputFolder + "/.";
    resCallSystem = system(buff.c_str());
    if (resCallSystem == -1) {
      throw FileException("thames", "deleteDynAllocMem", buff, "FAILED");
    }
    // cout << buff << endl;
    f.close();
  }

  name = "IPM_dump.txt";
  f.open(name.c_str());
  if (f.good()) {
    buff = "mv -f IPM_dump.txt " + outputFolder + "/.";
    resCallSystem = system(buff.c_str());
    if (resCallSystem == -1) {
      throw FileException("thames", "deleteDynAllocMem", buff, "FAILED");
    }
    // cout << buff << endl;
    f.close();
  }

  cout << endl
       << "=> IPM_dump.txt & ipmlog.txt are into " << outputFolder << " folder"
       << endl;

  cout << endl << "STOP Program" << endl;
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
  cout << endl << asctime(inittime1);
  clock_t endtime = clock();

  double elapsedtime = static_cast<double>(endtime - time_) / CLOCKS_PER_SEC;
  double ltD = difftime(lt1, lt_);
  cout << endl << "Total time = " << ltD << " seconds" << endl;
  cout << endl
       << "Total time with clock = " << elapsedtime << " seconds" << endl;
}

void printHelp(void) {
  cout << endl;
  cout << "Usage: \"thames [--verbose|-v] [--suppress|-s] [--xyz|-x] "
          "[--outfolder|-o] folder [--help|-h]\""
       << endl;
  cout << "        --verbose [-v]      Produce verbose output" << endl;
  cout << "        --suppress [-s]     Suppress warning messages" << endl;
  cout << "        --xyz [-x]          Create 3D visualization movie" << endl;
  cout << "        --outfolder [-o]    Name of folder for output data (default "
          "is Result)"
       << endl;
  cout << endl;
  cout << "Note: thames --help [-h] will print this help message" << endl;

  return;
}

int checkArgs(int argc, char **argv, string &outputFolder) {

  // Many of the variables here are defined in the getopts.h system header file
  // Can define more options here if we want

  const char *const short_opts = "vsxo:h";
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
  outputFolder = "Result";

  while (true) {

    const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

    if (-1 == opt)
      break; // breaks out of this while loop

    switch (opt) {
    case 'v':
      VERBOSE = true; // Verbose defined in thameslib global.h
      cout << "**Will produce verbose output**" << endl;
      break;
    case 's':
      WARNING = false; // Verbose defined in thameslib global.h
      cout << "**Will suppress warning messages**" << endl;
      break;
    case 'x':
      XYZ = true; // Verbose defined in thameslib global.h
      cout << "**Will create 3D visualization file **" << endl;
      break;
    case 'o':
      outputFolder = optarg; // Verbose defined in thameslib global.h
      if ((outputFolder[0] == ' ') || (outputFolder[0] == '-') ||
          (outputFolder[0] == '\\')) {
        printHelp();
        return (1);
      }
      cout << "Output folder: " << outputFolder << endl;
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

void prepOutputFolder(const string &outputFolder, string &jobRoot,
                      const string &gemInputName, string &statFileName,
                      const string &initMicName, const string &simParamName) {

  int resCallSystem;

  string buff = "mkdir -p " + outputFolder;
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }
  jobRoot = outputFolder + "/" + jobRoot;
  cout << "   - jobRoot           :  " << jobRoot << endl;
  cout.flush();

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
  buff = "cp -f " + buff1 + " " + outputFolder + "/.";
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
  buff = "cp -f " + buff1 + " " + outputFolder + "/.";
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
  buff = "cp -f " + buff1 + " " + outputFolder + "/.";
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }

  buff = "cp -f " + gemInputName + " " + outputFolder + "/.";
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }

  buff = "cp -f " + initMicName + " " + outputFolder + "/.";
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }

  buff = "cp -f " + simParamName + " " + outputFolder + "/.";
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    throw FileException("thames", "prepOutputFolder", buff, "FAILED");
  }

  cout << "     => All input files have been copied into " << outputFolder
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
      cout << "WARNING:  Could not open report file" << endl;
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
  cout << endl << "   Last number of voxels for each microPhase <-> nVmPh(i) :"
       << endl;
  for (int i = 0; i < dim; i++) {
    sumCount += count[i];
    cout << "    " << setw(3) << right << i
         << " - " << setw(18) << left << phaseName[i]
         << " : " << setw(10) << right << count[i]
         << endl;
  }
  cout << endl << "   nVmPh_total = " << sumCount << endl;
  cout << "   system dimension is : dimX*dimY*dimZ = "
       << mic->getXDim() << " * " << mic->getYDim() << " * " << mic->getZDim()
       << " = " << mic->getXDim() * mic->getYDim() * mic->getZDim()
       << endl;
}
