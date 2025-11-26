/**
@file Lattice.cc
@brief Definition of methods for the Lattice class.

*/
#include "Lattice.h"

using std::cout; using std::endl;
using std::string; using std::vector; using std::map;

Lattice::Lattice(ChemicalSystem *cs) : chemSys_(cs) {
  xdim_ = ydim_ = zdim_ = 0;
  time_ = 0.0;
  temperature_ = REFTEMP; // in Kelvin (see global.h)
  oldtemp_ = REFTEMP;     // in Kelvin (see global.h)
  numSites_ = 0;
  resolution_ = REFRES * 1.0e-6;                // in meters (see global.h)
  areaPerFace_ = resolution_ * resolution_;     // in m2 units
  volumePerVoxel_ = areaPerFace_ * resolution_; // in m3 units
  site_.clear();
  depthEffect_ = false;
  masterPoreSizeDist_.clear();
#ifdef DEBUG
  verbose_ = true;
#else
  verbose_ = false;
#endif
}

Lattice::Lattice(ChemicalSystem *cs, RanGen *rg, int seedRNG,
                 const string &fileName, const bool verbose,
                 const bool warning) {
  // : chemSys_(cs), rg_(rg) {

  bool seededCSH = false;
  bool seedMassCEM = false;
  bool seedMassC3S = true;
  bool seedMassC2S = false;
  double massFraction = 0.02;

  if (seededCSH) {
    cout << endl
         << "Lattice::Lattice : for this simulation the initial configuration contains CSH seeds!"
         << endl;
    cout << "  - the number of seeds is computed as a mass fraction ("
         << massFraction << ") of the total amount of ";
    if (seedMassCEM == true) {
      cout << "cement"<< endl;
    } else if (seedMassC3S == true && seedMassC2S == true) {
      cout << "alite and belite"<< endl;
    } else if (seedMassC3S == true) {
      cout << "alite"<< endl;
    } else if (seedMassC2S == true) {
      cout << "belite"<< endl;
    }
    cout << "  - the CSH seeds are placed randomly into the structure" << endl;
  } else {
    cout << endl
         << "Lattice::Lattice : for this simulation the initial configuration doesn't "
            "contain CSH seeds"
         << endl;
  }

  int i, j, k;
  int ii;
  string buff;
  int xn, yn, zn;
  int idn;
  int pid;
  string msg;

  chemSys_ = cs;
  rg_ = rg;

  convFactDbl2IntPor_ = chemSys_->getConvFactDbl2IntPor();
  microPhasePorosityInt_.clear();
  microPhasePorosityInt_ = chemSys_->getMicroPhasePorosityInt();

  affinityInt_.clear();
  affinityInt_ = chemSys_->getAffinityInt();

  numMicroPhases_ = chemSys_->getNumMicroPhases();
  waterDCId_ = chemSys_->getDCId("H2O@");
  waterMolarMass_ = chemSys_->getDCMolarMass(waterDCId_);
  waterMolarVol_ = chemSys_->getDCMolarVolume(waterDCId_);

  growthTemplate_.clear();
  bool tmpl;
  vector<bool> grVectTmpl;
  for (int i = 0; i < numMicroPhases_; i++) {
    grVectTmpl.clear();
    for (int j = 0; j < numMicroPhases_; j++) {
      tmpl =chemSys_->isGrowthTemplate(i, j);
      grVectTmpl.push_back(tmpl);
    }
    growthTemplate_.push_back(grVectTmpl);
  }

  ///
  /// Gather data from the periodic table for writing a .cfg file (alternative
  /// to .xyz)
  populateElementData();
  particRadius_ = 0.7;

  xdim_ = ydim_ = zdim_ = 0;
  time_ = 0.0;
  temperature_ = REFTEMP; // in Kelvin (see global.h)
  oldtemp_ = REFTEMP;     // in Kelvin (see global.h)
  numSites_ = 0;
  resolution_ = REFRES * 1.0e-6;                // in meters (see global.h)
  areaPerFace_ = resolution_ * resolution_;     // in m2 units
  volumePerVoxel_ = areaPerFace_ * resolution_; // in m3 units
  site_.clear();
  // depthEffect_ = true;
  depthEffect_ = false;
  masterPoreSizeDist_.clear();

#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
  cout << "Lattice::Lattice Constructor" << endl;
  cout.flush();
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  ///
  /// Open the microstructure input file and process it
  ///

  ifstream in(fileName.c_str());
  if (!in) {
    throw FileException("Lattice", "Lattice", fileName, "Could not open");
  }

  in >> buff;
  if (buff == VERSIONSTRING) {
    in >> version_;
    cout << endl
         << "Lattice::Lattice - " << fileName
         << " input file generated using "
            "THAMES Version : "
         << version_ << endl;
    in >> buff; // X size string identifier
    in >> xdim_;
    in >> buff; // Y size string identifier
    in >> ydim_;
    in >> buff; // Z size string identifier
    in >> zdim_;
    double testres;
    in >> buff; // Voxel resolution identifier
    in >> testres;
    setResolution(testres * 1.0e-6); // Convert to meters
  } else {

    ///
    /// Image file generated prior to VCCTL Version 3.0.
    /// Allow backward compatibility by defaulting system
    /// size to 100 and resolution to 1.0 micrometers
    ///
    cout << endl
         << "Lattice::Lattice - " << fileName
         << " input file generated using "
            "THAMES Version prior to THAMES Version 3.0.0"
         << endl;
    version_ = "2.0";
    double testres = 1.0e-6; // in meters
    setResolution(testres);
    xdim_ = ydim_ = zdim_ = 100;
  }

  if (xdim_ <= 0 || ydim_ <= 0 || zdim_ <= 0) {
    cout << endl
         << ">>>>> all xdim_, ydim_ & zdim_ must be grater than 0!" << endl;
    cout << endl << "      xdim_ = " << xdim_ << endl;
    cout << "      ydim_ = " << ydim_ << endl;
    cout << "      zdim_ = " << zdim_ << endl;
    cout << endl << ">>>>> stop <<<<<" << endl;
    exit(0);
  }

  ostringstream ostrMAJ, ostrMIN;
  ostrMAJ << VERSION_MAJOR;
  ostrMIN << VERSION_MINOR;
  string majVers(ostrMAJ.str());
  string minVers(ostrMIN.str());
  thamesVersion_ = majVers + "." + minVers + "." + VERSIONBUGFIX;
  cout << endl
       << "Lattice::Lattice - .img output files generated using "
          "THAMES Version: "
       << thamesVersion_ << endl;

  ///
  /// Print out the microstructure size and characteristics
  ///

  if (verbose_) {
    cout << "Lattice::Lattice Read microstructure file header..." << endl;
    cout << "Lattice::Lattice     THAMES Version = " << version_ << endl;
    cout << "Lattice::Lattice     xdim_ = " << xdim_ << endl;
    cout << "Lattice::Lattice     ydim_ = " << ydim_ << endl;
    cout << "Lattice::Lattice     zdim_ = " << zdim_ << endl;
    cout.flush();
  }
  numSites_ = xdim_ * ydim_ * zdim_;
  if (verbose_) {
    cout << "Lattice::Lattice    numSites_ = " << numSites_ << endl;
    cout << "Lattice::Lattice    resolution_ = " << resolution_ << endl;
    cout.flush();
  }

  ///
  /// Allocate a random number generator object and seed it
  ///

  latticeRNGseed_ = seedRNG; //-76241;
  // cout << endl << "Lattice::Lattice   latticeRNGseed_ = " << latticeRNGseed_
  // << endl; setRNGseed(latticeRNGseed_); try {
  //   rg_->setSeed(latticeRNGseed_);
  // } catch (bad_alloc &ba) {
  //   cout << "Lattice constructor failed when allocating rg_";
  //  cout.flush();
  //   exit(1);
  // }
  // //cout << endl << "Lattice::Lattice   rg_->getSeed() = " << rg_->getSeed()
  // << endl;
  numRNGcall_0_ = 0;
  numRNGcallLONGMAX_ = 0;
  lastRNG_ = 1.e-16;

  growthInterfaceSize_.clear();
  dissolutionInterfaceSize_.clear();
  growthInterfaceSize_.resize(numMicroPhases_, 0);
  dissolutionInterfaceSize_.resize(numMicroPhases_, 0);

  // vectors based on simparams.json file, calling ChemicalSystem::getVectorsSA
  growingVectSA_.clear();
  shrinking_.clear();
  volratios_.clear();

  count_.resize(numMicroPhases_, 0);

  expansion_.clear();
  // expansion_coordin_.clear();

  waterChange_ = 0.0;

  ///
  /// Add the Site objects to the `site_` array, defining the structure
  /// of the lattice.  Not yet assigning a phase to the site because
  /// the file has not yet been read
  ///

  for (k = 0; k < zdim_; k++) {
    for (j = 0; j < ydim_; j++) {
      for (i = 0; i < xdim_; i++) {
        addSite(i, j, k);
      }
    }
  }

  ///
  /// Set up the array of neighbors for each site
  ///

  for (ii = 0; ii < numSites_; ii++) {

    for (j = 0; j < NN_NNN; j++) {
      ///
      /// Store up to 2nd nearest neighbors
      ///

      switch (j) {
      case 0: // West neighbor
        xn = site_[ii].getX() - 1;
        if (xn < 0 && BC != 1)
          xn += xdim_;
        if (xn < 0 && BC == 1)
          xn = 0;
        yn = site_[ii].getY();
        zn = site_[ii].getZ();
        break;
      case 1: // East neighbor
        xn = site_[ii].getX() + 1;
        if (xn >= xdim_ && BC != 1)
          xn -= xdim_;
        if (xn >= xdim_ && BC == 1)
          xn -= 1;
        yn = site_[ii].getY();
        zn = site_[ii].getZ();
        break;
      case 2: // South neighbor
        yn = site_[ii].getY() - 1;
        if (yn < 0 && BC != 2)
          yn += ydim_;
        if (yn < 0 && BC == 2)
          yn = 0;
        xn = site_[ii].getX();
        zn = site_[ii].getZ();
        break;
      case 3: // North neighbor
        yn = site_[ii].getY() + 1;
        if (yn >= ydim_ && BC != 2)
          yn -= ydim_;
        if (yn >= ydim_ && BC == 2)
          yn -= 1;
        xn = site_[ii].getX();
        zn = site_[ii].getZ();
        break;
      case 4: // Down neighbor
        zn = site_[ii].getZ() - 1;
        if (zn < 0 && BC != 3)
          zn += zdim_;
        if (zn < 0 && BC == 3)
          zn = 0;
        xn = site_[ii].getX();
        yn = site_[ii].getY();
        break;
      case 5: // Up neighbor
        zn = site_[ii].getZ() + 1;
        if (zn >= zdim_ && BC != 3)
          zn -= zdim_;
        if (zn >= zdim_ && BC == 3)
          zn -= 1;
        xn = site_[ii].getX();
        yn = site_[ii].getY();
        break;
      case 6: // Southwest neighbor
        xn = site_[ii].getX() - 1;
        if (xn < 0 && BC != 1)
          xn += xdim_;
        if (xn < 0 && BC == 1)
          xn = 0;
        yn = site_[ii].getY() - 1;
        if (yn < 0 && BC != 2)
          yn += ydim_;
        if (yn < 0 && BC == 2)
          yn = 0;
        zn = site_[ii].getZ();
        break;
      case 7: // Northwest neighbor
        xn = site_[ii].getX() - 1;
        if (xn < 0 && BC != 1)
          xn += xdim_;
        if (xn < 0 && BC == 1)
          xn = 0;
        yn = site_[ii].getY() + 1;
        if (yn >= ydim_ && BC != 2)
          yn -= ydim_;
        if (yn >= ydim_ && BC == 2)
          yn -= 1;
        zn = site_[ii].getZ();
        break;
      case 8: // Northeast neighbor
        xn = site_[ii].getX() + 1;
        if (xn >= xdim_ && BC != 1)
          xn -= xdim_;
        if (xn >= xdim_ && BC == 1)
          xn -= 1;
        yn = site_[ii].getY() + 1;
        if (yn >= ydim_ && BC != 2)
          yn -= ydim_;
        if (yn >= ydim_ && BC == 2)
          yn -= 1;
        zn = site_[ii].getZ();
        break;
      case 9: // Southeast neighbor
        xn = site_[ii].getX() + 1;
        if (xn >= xdim_ && BC != 1)
          xn -= xdim_;
        if (xn >= xdim_ && BC == 1)
          xn -= 1;
        yn = site_[ii].getY() - 1;
        if (yn < 0 && BC != 2)
          yn += ydim_;
        if (yn < 0 && BC == 2)
          yn = 0;
        zn = site_[ii].getZ();
        break;
      case 10: // Downwest neighbor
        xn = site_[ii].getX() - 1;
        if (xn < 0 && BC != 1)
          xn += xdim_;
        if (xn < 0 && BC == 1)
          xn = 0;
        yn = site_[ii].getY();
        zn = site_[ii].getZ() - 1;
        if (zn < 0 && BC != 3)
          zn += zdim_;
        if (zn < 0 && BC == 3)
          zn = 0;
        break;
      case 11: // Downnorth neighbor
        xn = site_[ii].getX();
        yn = site_[ii].getY() + 1;
        if (yn >= ydim_ && BC != 2)
          yn -= ydim_;
        if (yn >= ydim_ && BC == 2)
          yn -= 1;
        zn = site_[ii].getZ() - 1;
        if (zn < 0 && BC != 3)
          zn += zdim_;
        if (zn < 0 && BC == 3)
          zn = 0;
        break;
      case 12: // Downeast neighbor
        xn = site_[ii].getX() + 1;
        if (xn >= xdim_ && BC != 1)
          xn -= xdim_;
        if (xn >= xdim_ && BC == 1)
          xn -= 1;
        yn = site_[ii].getY();
        zn = site_[ii].getZ() - 1;
        if (zn < 0 && BC != 3)
          zn += zdim_;
        if (zn < 0 && BC == 3)
          zn = 0;
        break;
      case 13: // Downsouth neighbor
        xn = site_[ii].getX();
        yn = site_[ii].getY() - 1;
        if (yn < 0 && BC != 2)
          yn += ydim_;
        if (yn < 0 && BC == 2)
          yn = 0;
        zn = site_[ii].getZ() - 1;
        if (zn < 0 && BC != 3)
          zn += zdim_;
        if (zn < 0 && BC == 3)
          zn = 0;
        break;
      case 14: // Upwest neighbor
        xn = site_[ii].getX() - 1;
        if (xn < 0 && BC != 1)
          xn += xdim_;
        if (xn < 0 && BC == 1)
          xn = 0;
        yn = site_[ii].getY();
        zn = site_[ii].getZ() + 1;
        if (zn >= zdim_ && BC != 3)
          zn -= zdim_;
        if (zn >= zdim_ && BC == 3)
          zn -= 1;
        break;
      case 15: // Upnorth neighbor
        xn = site_[ii].getX();
        yn = site_[ii].getY() + 1;
        if (yn >= ydim_ && BC != 2)
          yn -= ydim_;
        if (yn >= ydim_ && BC == 2)
          yn -= 1;
        zn = site_[ii].getZ() + 1;
        if (zn >= zdim_ && BC != 3)
          zn -= zdim_;
        if (zn >= zdim_ && BC == 3)
          zn -= 1;
        break;
      case 16: // Upeast neighbor
        xn = site_[ii].getX() + 1;
        if (xn >= xdim_ && BC != 1)
          xn -= xdim_;
        if (xn >= xdim_ && BC == 1)
          xn -= 1;
        yn = site_[ii].getY();
        zn = site_[ii].getZ() + 1;
        if (zn >= zdim_ && BC != 3)
          zn -= zdim_;
        if (zn >= zdim_ && BC == 3)
          zn -= 1;
        break;
      case 17: // Upsouth neighbor
        xn = site_[ii].getX();
        yn = site_[ii].getY() - 1;
        if (yn < 0 && BC != 2)
          yn += ydim_;
        if (yn < 0 && BC == 2)
          yn = 0;
        zn = site_[ii].getZ() + 1;
        if (zn >= zdim_ && BC != 3)
          zn -= zdim_;
        if (zn >= zdim_ && BC == 3)
          zn -= 1;
        break;
      }

      idn = xn + (xdim_ * yn) + (xdim_ * ydim_) * zn;
      site_[ii].setNb(j, &site_[idn]);
    }
  }

  setNeighborhoodSA();

  ///
  /// This loop reads the phase for each site and assigns it,
  /// also updating the count of the phase.
  ///

  for (i = 0; i < numSites_; i++) {
    in >> pid;
    site_[i].setMicroPhaseId(pid);
    count_[pid]++;
  }

  ///
  /// Done with the input microstructure file, so close it.
  ///

  in.close();

  ///
  /// With the phase counts known, calculate phase volume fractions
  /// We cannot calculate the mass fractions until we know the phase identities
  /// Phase identities are determined by the KineticModel
  ///

  vector<double> microPhaseMass(numMicroPhases_, 0.0);

  surfaceArea_.clear();
  surfaceArea_.resize(numMicroPhases_, 0.0);
  specificSurfaceArea_.clear();
  specificSurfaceArea_.resize(numMicroPhases_, 0.0);
  volumeFraction_.clear();
  volumeFraction_.resize(numMicroPhases_, 0.0);
  // initVolumeFraction_.clear();
  // initVolumeFraction_.resize(numMicroPhases_, 0.0);

  double vfrac, molarMass, molarVolume, density;
  double solidMass = 0.0;
  double cementMass = 0.0;
  int microPhaseId = 0;
  int DCId = 0;
  string myname;
  try {
    if (verbose_) {
      cout << "Lattice::Lattice Calculating volume fractions now ..." << endl;
      for (ii = 0; ii < numMicroPhases_; ii++) {
        cout << "Micro phase "
             << chemSys_->getMicroPhaseName(ii)
             << ", count = " << count_[ii] << " of " << site_.size() << endl;
      }
      cout.flush();
    }

    // DAMAGEID_ = chemSys_->getMicroPhaseId("DAMAGE");

    if (numSites_ > 0) {
      for (microPhaseId = 0; microPhaseId < numMicroPhases_; microPhaseId++) {
        myname = chemSys_->getMicroPhaseName(microPhaseId);
        vfrac = (static_cast<double>(count_[microPhaseId])) /
                (static_cast<double>(numSites_));

        if (vfrac > 0.0) {
          try {
            volumeFraction_[microPhaseId] = vfrac;
            // setVolumeFraction(microPhaseId, vfrac);
            // setInitVolumeFraction(microPhaseId, vfrac);
          } catch (EOBException eex) {
            throw eex;
          }
          if (verbose_) {
            cout << "Lattice::Lattice microPhaseId = " << microPhaseId
                 << ", microPhase = " << myname << ", count_[" << microPhaseId
                 << "] = " << count_[microPhaseId] << ", volume fraction = "
                 << vfrac << endl;
            cout.flush();
          }
          if (microPhaseId == ELECTROLYTEID) {
            DCId = waterDCId_;
          } else if (microPhaseId !=
                     VOIDID) { // && microPhaseId != DAMAGEID_) {
            DCId = chemSys_->getMicroPhaseDCMembers(microPhaseId, 0);
          }
          if (microPhaseId != VOIDID) { // && microPhaseId != DAMAGEID_) {
            molarMass = chemSys_->getDCMolarMass(DCId);     // g/mol
            molarVolume = chemSys_->getDCMolarVolume(DCId); // m3/mol
            density = 0.0;
            if (molarVolume > 1.0e-12)
              density = molarMass / molarVolume / 1.0e6; // g/cm3

            // microPhaseMass, solidMass, and cementMass here all have
            // units of g of phase per cm3 of whole microstructure
            microPhaseMass[microPhaseId] = vfrac * density;
            if (microPhaseId != ELECTROLYTEID) {
              solidMass += microPhaseMass[microPhaseId];
              if (chemSys_->isCementComponent(microPhaseId))
                cementMass += microPhaseMass[microPhaseId];
            }
            if (verbose_ && vfrac > 0.0) { // if (vfrac > 0.0) {
              cout << microPhaseId << "\tLattice::Lattice Phase "
                   << chemSys_->getMicroPhaseName(microPhaseId)
                   << "\tmicroPhaseId: " << microPhaseId
                   << ", DCName = " << chemSys_->getDCName(DCId) << endl;
              cout << "Lattice::Lattice     Molar mass = " << molarMass
                   << " g/mol" << endl;
              cout << "Lattice::Lattice     Molar volume = " << molarVolume
                   << " m3/mol" << endl;
              cout << "Lattice::Lattice     Density = " << density << " g/cm3"
                   << endl;
              cout << "Lattice::Lattice     Volume Fraction = " << vfrac
                   << endl;
              cout << "Lattice::Lattice     Mass density = "
                   << (vfrac * density) << " g/cm3 of system" << endl;
              cout.flush();
            }
          }
        }
      } // End of for loop over all microstructure phases
    } else {
      msg = "Divide by zero error:  site_.size() = 0";
      throw FloatException("Lattice", "Lattice", msg);
    }

    // for (ii = 0; ii < numMicroPhases_; ii++) {
    //   cout << "  ii : " << ii << "\tcementComponent : "
    //        << chemSys_->isCementComponent(ii) << endl;
    // }

    // Set the water-solids mass ratio based on the initial microstructure

    initSolidMass_ = solidMass;
    calcOneFaceAreaPerHundredGramSolid();
    chemSys_->setInitSolidMass(initSolidMass_);

    wsRatio_ = microPhaseMass[ELECTROLYTEID] / solidMass;
    wcRatio_ = microPhaseMass[ELECTROLYTEID] / cementMass;

    if (verbose_) {
      cout << "Lattice::Lattice Microstructure w/s = " << wsRatio_ << endl;
      cout << "Lattice::Lattice Mass of water = "
           << microPhaseMass[ELECTROLYTEID]
           << ", mass of solids = " << solidMass << endl;
      cout.flush();
    }

    // set to zero ALL ICs/DCs in the system (ChemicalSystem)
    int numICs = chemSys_->getNumICs();
    for (int i = 0; i < numICs; i++) {
      chemSys_->setICMoles(i, 0.0);
      // cout << "   " << i << "\t" << chemSys_->getICMoles(i)
      //     << "\t" << chemSys_->getICName(i) << endl;
    }

    int numDCs = chemSys_->getNumDCs();
    for (int i = 0; i < numDCs; i++) {
      chemSys_->setDCMoles(i, 0.0);
      chemSys_->setDCLowerLimit(i, 0.0);
      chemSys_->setDCUpperLimit(i, 1.0e6);
      // cout << "   " << i << "\t" << chemSys_->getDCMoles(i)
      //      << "\t" << chemSys_->getDCClassCode(i) << "\t" <<
      //      chemSys_->getDCName(i) << endl;
    }

    // Next we set the initial normalized phase masses, microstructure
    // phase volume, and subvoxel porosity.  This is all triggered when we set
    // the mass. All is normalized to 100 g of total solid in the initial
    // microstructure

    // calc porosity in
    // normalizePhaseMasses->setMicroPhaseMass->setMicroPhaseVolume->calcMicroPhasePorosity

    try {
      normalizePhaseMasses(microPhaseMass);
      // Next we also normalize the surface areas to the same 100 g of total
      // solid in the initial microstructure
      // Calculate the surface area per 100 g solid
      calcSurfaceAreas();
    } catch (FloatException flex) {
      throw flex;
    }

    // Set the initial total volume of the microstructure

    double totmicvol = 0.0;
    int totCount_ = 0;
    cout << endl << "Lattice::Lattice - volumes & count_ : " << endl;
    for (int i = 0; i < numMicroPhases_; i++) {
      microPhaseId = chemSys_->getMicroPhaseId(i);
      if (microPhaseId != VOIDID) {
        cout << "    mPhId = " << i << "\tvolume = "
             << chemSys_->getMicroPhaseVolume(microPhaseId)
             << "\t" << chemSys_->getMicroPhaseName(i) << " / count_ = "
             << count_[i] << endl ;
        totmicvol += chemSys_->getMicroPhaseVolume(microPhaseId);
        totCount_ += count_[i];
      }
    }
    cout << "Lattice::Lattice - totCount_/totmicvol: " << totCount_
         << " / " << totmicvol << endl;

    chemSys_->setInitMicroVolume(totmicvol);
    chemSys_->setInitGEMVolume();
    initialMicrostructureVolume_ = totmicvol;

    // Initially assume that all free water and void space
    // is in voxels

    // voxelPoreVolumeFraction_ =
    //    getVolumeFraction(ELECTROLYTEID) + getVolumeFraction(VOIDID);

    if (seededCSH) {
      addSeedCSHQ(seedMassCEM, seedMassC3S, seedMassC2S, massFraction);
      calcSurfaceAreas();
    }

  } catch (FloatException flex) {
    throw flex;
  } catch (EOBException eex) {
    throw eex;
  }

  // calc & set wmc
  electrolyteIntPorosity_ = microPhasePorosityInt_[ELECTROLYTEID];
  voidIntPorosity_ = microPhasePorosityInt_[VOIDID];
  int phId;
  // double rng;
  for (i = 0; i < numSites_; i++) {
    // stId = site_[i].getId();
    phId = site_[i].getMicroPhaseId();
    if (phId == ELECTROLYTEID) {
      site_[i].setWmc0(electrolyteIntPorosity_); // 1.e5
    } else if ((phId != VOIDID)) {
      // rng = callRNG();
      // if (rng >= thrPorosity) {
      //   site_[i].setWmc0(microPhasePorosityInt_[phId]);
      // } else {
      //   site_[i].setWmc0(0);
      // }
      site_[i].setWmc0(microPhasePorosityInt_[phId]);
    } else {
      site_[i].setWmc0(voidIntPorosity_); // VOID 1.e5
    }
  }

  int wmcLocInt;
  int nbIdLoc;
  for (i = 0; i < numSites_; i++) {
    site_[i].setVisit(0);
    wmcLocInt = site_[i].getWmc0();
    for (j = 0; j < NN_NNN; j++) {
      nbIdLoc = (site_[i].nb(j))->getId();
      wmcLocInt += site_[nbIdLoc].getWmc0();
    }
    site_[i].setWmc(wmcLocInt);
  }

  // voxels without contact with electrolyte i.e. voxels having low probability
  // to dissolve in this step
  // findIsolatedClusters();
}

Lattice::~Lattice() {
  interface_.clear();
  site_.clear();
  masterPoreSizeDist_.clear();
  // delete rg_;
}

void Lattice::addSite(int x, const int y, const int z) {
  // string msg;

  // try {
  //   if (x >= xdim_ || x < 0) {
  //     throw EOBException("Lattice", "addSite", "site_", xdim_, x);
  //   } else if (y >= ydim_ || y < 0) {
  //     throw EOBException("Lattice", "addSite", "site_", ydim_, y);
  //   } else if (z >= zdim_ || z < 0) {
  //     throw EOBException("Lattice", "addSite", "site_", zdim_, z);
  //   }
  // } catch (EOBException ex) {
  //   ex.printException();
  //   exit(1);
  // }

  // Assume that the site is WATER by default.  This will break if
  // there is ever a sytem without water.

  site_.push_back(Site(x, y, z, xdim_, ydim_, zdim_, NN_NNN, chemSys_));
}

void Lattice::normalizePhaseMasses(vector<double> microPhaseMass) {
  int microPhaseId, DCId;
  double pscaledMass = 0.0;
  double molarMass;
  double totalSolidMass = 0.0, totalCementMass = 0.0;

  for (int i = 0; i < numMicroPhases_; i++) {
    microPhaseId = chemSys_->getMicroPhaseId(i);
    if (microPhaseId == ELECTROLYTEID) {
      // int waterId = chemSys_->getDCId("H2O@");
      // double waterMolarMass = chemSys_->getDCMolarMass(waterId);
      // molarMass = chemSys_->getDCMolarMass(waterId);
      pscaledMass = wsRatio_ * 100.0; // Mass of solids scaled to 100 g now

      chemSys_->setDCMoles(waterDCId_, (pscaledMass / waterMolarMass_));
      if (verbose_) {
        cout << "Lattice::normalizePhaseMasses Setting initial micphase mass "
                "and volume of "
             << chemSys_->getMicroPhaseName(ELECTROLYTEID) << endl;
        cout.flush();
      }

      chemSys_->setMicroPhaseMass(ELECTROLYTEID, pscaledMass);
      chemSys_->setMicroPhaseMassDissolved(ELECTROLYTEID, 0.0);

    } else if (microPhaseId != VOIDID) { // && microPhaseId != DAMAGEID_) {
      // microPhaseMass has units of grams per cm3 of whole microstructure
      // initSolidMass_ is the sum of all microPhaseMass values
      // So pscaledMass (p stands for percent) is mass percent on a total solids
      // basis
      pscaledMass = microPhaseMass[microPhaseId] * 100.0 / initSolidMass_;
      DCId = chemSys_->getMicroPhaseDCMembers(microPhaseId, 0);
      molarMass = chemSys_->getDCMolarMass(DCId);
      if (verbose_) {
        cout << "Lattice::normalizePhaseMasses, 1 cm3 of microstructure has "
                "mass of "
             << i << "   " << chemSys_->getMicroPhaseName(microPhaseId) << " ("
             << microPhaseId << ") = " << microPhaseMass[microPhaseId]
             << " g out of " << initSolidMass_ << " g total" << endl;
        // Setting the phase mass will also automatically calculate the phase
        // volume
        cout.flush();
      }

      totalSolidMass += pscaledMass;
      if (chemSys_->isCementComponent(microPhaseId))
        totalCementMass += pscaledMass;

      // into the next method : setMicroPhaseVolume & calcMicroPhasePorosity
      chemSys_->setMicroPhaseMass(microPhaseId, pscaledMass);

      // ChemicalSystem DC moles is in units of moles per 100 g of total solid
      chemSys_->setDCMoles(DCId, (pscaledMass / molarMass));

      chemSys_->setMicroPhaseMassDissolved(microPhaseId, 0.0);
    }
  }

  // initScaledCementMass is percent of all solid mass that is cement
  // chemSys_->setInitScaledCementMass(cementMass * 100 / initSolidMass_);

  // Up to this point we could not really handle volume of void space, but
  // now we can do so in proportion to electrolyte volume

  double totalMass = 0;
  cout << endl << "  Lattice::normalizePhaseMasses - normalized masses:" << endl;
  for (int i = 1; i < numMicroPhases_; i++) {
    if (chemSys_->isCementComponent(i)) {
      cout << setw(6) << right << i << " : " << setw(15) << left
           << chemSys_->getMicroPhaseName(i) << chemSys_->getMicroPhaseMass(i)
           << " g (*)" << endl;
    } else {
      cout << setw(6) << right << i << " : " << setw(15) << left
           << chemSys_->getMicroPhaseName(i) << chemSys_->getMicroPhaseMass(i)
           << " g" << endl;
    }
    totalMass += chemSys_->getMicroPhaseMass(i);
  }
  cout << endl
       << "    totalMass       = " << totalMass
       << " g <all phases including water>" << endl;
  cout << "    totalSolidMass  = " << totalSolidMass << " g <all solid phases>"
       << endl;
  cout << "    totalCementMass = " << totalCementMass
       << " g <only cement phases (*)>" << endl;
  cout << "    wsRatio_        = " << wsRatio_
       << "   <waterMass/totalSolidMass>>" << endl;
  cout << "    wcRatio_        = " << wcRatio_
       << "   <waterMass/totalCementMass>" << endl;

  // chemSys_->setInitScaledCementMass(cementMass * 100 / solidMass);
  // cout << "normalizePhaseMasses totalSolidMass/totalCementMass = "
  //      << totalSolidMass << " / " << totalCementMass << endl;

  // Up to this point we could not really handle volume of void space, but
  // now we can do so in proportion to electrolyte volume

  double vfv = getVolumeFraction(VOIDID);
  double vfe = getVolumeFraction(ELECTROLYTEID);
  double ve = chemSys_->getMicroPhaseVolume(ELECTROLYTEID);

  chemSys_->setMicroPhaseVolume(VOIDID, (ve * vfv / vfe));

  // cout << endl << "normalizePhaseMasses DCs:" << endl;
  // for (int i = 0; i < chemSys_->getNumDCs(); i++){
  //   cout << "   " << i << "\t" << chemSys_->getDCMoles(i)
  //        << "\t" << chemSys_->getDCClassCode(i)
  //        << "\t" << chemSys_->getDCName(i) << endl;
  // }
  // cout << endl << "normalizePhaseMasses OK!" << endl << endl;

  return;
}

void Lattice::findInterfaces(void) {
  int i, kk;
  int k;
  vector<Site *> gsite, dsite;
  vector<Site *>::iterator beginLocation, endLocation;
  int stId;
  int gvsize = gsite.size();
  int dvsize = dsite.size();

  ///
  /// An interface must have at least one adjacent site that is water or void
  ///
  ///
  cout << endl << "   Lattice::findInterfaces() :" << endl;

  interface_.clear();
  for (i = 0; i < numMicroPhases_; i++) {
    if (i != ELECTROLYTEID && i != VOIDID) { // Solid phase of some kind
      gsite.clear();
      dsite.clear();
      for (k = 0; k < numSites_; k++) {
        if (site_[k].getWmc() > 0) { // Is there some water nearby?
          if ((site_[k].getMicroPhaseId() == i)) {
            dsite.push_back(&site_[k]);
            for (kk = 0; kk < NN_NNN; kk++) {
              if ((site_[k].nb(kk))->getMicroPhaseId() == ELECTROLYTEID) {
                gsite.push_back(site_[k].nb(kk));
                site_[k].nb(kk)->setGrowthSite(i);
              }
            }
            //}

            /// @note There is no reason to make the phase
            /// itself be a template for itself, nor water be
            /// a growth template for any phase, because we already
            /// tested for those above.

          } else if (chemSys_->isGrowthTemplate(i,
                                                site_[k].getMicroPhaseId())) {
            for (kk = 0; kk < NUM_NEAREST_NEIGHBORS;
                 kk++) { // site_[k].nbSize(1)
              if ((site_[k].nb(kk))->getMicroPhaseId() == ELECTROLYTEID) {
                gsite.push_back(site_[k].nb(kk));
                site_[k].nb(kk)->setGrowthSite(i);
              }
            }
          }
        }
      }

      // Eliminate duplicate values from the growth site vector
      gvsize = gsite.size();
      if (gvsize > 0) {
        sort(gsite.begin(), gsite.end());
        beginLocation = gsite.begin();
        endLocation = unique(gsite.begin(), gsite.end());
        gsite.erase(endLocation, gsite.end());
      }
      // Eliminate duplicate values from the dissolution site vector
      dvsize = dsite.size();
      if (dvsize > 0) {
        sort(dsite.begin(), dsite.end());
        beginLocation = dsite.begin();
        endLocation = unique(dsite.begin(), dsite.end());
        dsite.erase(endLocation, dsite.end());
      }

      gvsize = gsite.size();
      for (int j = 0; j < gvsize; j++) {
        stId = gsite[j]->getId();
        site_[stId].setInGrowInterfacePos(i, j);
      }
      dvsize = dsite.size();
      for (int j = 0; j < dvsize; j++) {
        stId = dsite[j]->getId();
        site_[stId].setInDissInterfacePos(j);
      }

      growthInterfaceSize_[i] = gvsize;
      dissolutionInterfaceSize_[i] = dvsize;

      interface_.push_back(Interface(chemSys_, gsite, dsite, i, verbose_));

    } else { // Not a solid phase (either electrolyte or void)

      interface_.push_back(Interface(verbose_));
    }
  }

  cout
      << endl
      << "                         ***** initial count_ & interface sizes *****"
      << endl;
  cout << "   numMicroPhases = " << numMicroPhases_ << endl;

  for (int i = 0; i < numMicroPhases_; i++) {
    cout << "  " << setw(3) << right << i << " : " << setw(15) << left
         << chemSys_->getMicroPhaseName(i) << setw(4) << right
         << " id:" << setw(3) << chemSys_->getMicroPhaseId(i)
         << "     count_ = " << setw(8) << count_[i]
         << "     dissolutionInterfaceSize_ =  " << setw(8)
         << dissolutionInterfaceSize_[i]
         << "     growthInterfaceSize_ =  " << setw(8)
         << growthInterfaceSize_[i]
         << "     porosity : " << chemSys_->getMicroPhasePorosity(i)
         << "     porosityInt : " << setw(8)
         << chemSys_->getMicroPhasePorosityInt(i)
         << "     templates : ";
    for (int j = 0; j < numMicroPhases_; j++) {
      if (chemSys_->isGrowthTemplate(i, j)) {
        cout << j << " ";
        //}else{
        //    cout <<"";
      }
    }
    cout << endl;
  }

  // test!
  // cout << endl << "  Lattice::findInterfaces - local growthTemplate_ :" << endl;
  // for (int i = 0; i < numMicroPhases_; i++ ) {
  //   cout << "     " << setw(3) << right << i << " : for " << setw(15) << left
  //        << chemSys_->getMicroPhaseName(i) << setw(4) << right
  //        << " id:" << setw(3) << right << chemSys_->getMicroPhaseId(i) << "   templates : ";
  //   for (int j = 0; j < numMicroPhases_; j++ ) {
  //     if (growthTemplate_[i][j])
  //       cout << j << " ";
  //   }
  //   cout << endl;
  // }

  /*
    {
      for (int i = 0; i < numMicroPhases_; i++) {
        int sizeIntLatt = dissolutionInterfaceSize_[i];
        if (sizeIntLatt > 0) {
          int sizeIntInte = interface_[i].getDissolutionSize();
          int intSiteId, intPhId, posSite;
          cout << endl
               << endl
               << "   Lattice::findInterfaces() PhId/sizeIntLatt/sizeIntInte : "
    << i << " / "
               << sizeIntLatt << " / " << sizeIntInte << endl
               << endl;

          for (int jj = 0; jj < sizeIntLatt; jj++) {
            intSiteId = interface_[i].getDissolutionSitesId(jj);
            intPhId = site_[intSiteId].getMicroPhaseId();
            posSite = site_[intSiteId].getInDissInterfacePos();
            if (jj == posSite) {
              cout << "      " << jj << "     intSiteId : " << intSiteId
                   << "     posSite : " << posSite << "     intPhId : " <<
    intPhId << endl; } else { cout << "      " << jj << "     intSiteId : " <<
    intSiteId
                   << "     posSite : " << posSite << "     intPhId : " <<
    intPhId << "   +++"
                   << endl;
            }
          }
        }
      }
    }
  */

  // cout << "***   findInterfaces exit" << endl; exit(0);

  /*
  // **************************** test for modifying affinities (contact angles)
  vector<int> nbC3S_tot, nbCSH_tot, nbC3S_woCSH, nbCSH_woC3S;
  nbC3S_tot.clear();
  nbCSH_tot.clear();
  nbC3S_woCSH.clear();
  nbCSH_woC3S.clear();
  nbC3S_tot.resize(19, 0);
  nbCSH_tot.resize(19, 0);
  nbC3S_woCSH.resize(19, 0);
  nbCSH_woC3S.resize(19, 0);
  int inbC3S_tot, inbCSH_tot; // , inbCSH_woC3S, inbC3S_woCSH;
  int siteID, nbPhId;
  Site *stenb;

  int cshqID = chemSys_->getMicroPhaseId("CSHQ");
  int aliteID = chemSys_->getMicroPhaseId("Alite");
  vector<Isite> isite;
  isite = interface_[cshqID].getGrowthSites();
  int dimIsite = isite.size();
  for (int jj = 0; jj < dimIsite; jj++) {
    siteID = isite[jj].getId();
    inbC3S_tot = 0;
    inbCSH_tot = 0;
    // inbC3S_woCSH = 0;
    // inbCSH_woC3S = 0;
    for (int i = 0; i < NN_NNN; i++) {
      stenb = site_[siteID].nb(i);
      nbPhId = stenb->getMicroPhaseId();
      if (nbPhId == cshqID){
        inbCSH_tot++;
      }
      if (nbPhId == aliteID){
        inbC3S_tot++;
      }
    }
    nbCSH_tot[inbCSH_tot]++;
    nbC3S_tot[inbC3S_tot]++;
    if (inbC3S_tot == 0) {
      nbCSH_woC3S[inbCSH_tot]++;
    }
    if (inbCSH_tot == 0) {
      nbC3S_woCSH[inbC3S_tot]++;
    }
  }

  ofstream outfs;
  string outfilename;

  outfilename = "nbCSH_tot.csv";
  outfs.open(outfilename.c_str());
  for (i = 0; i < 19; i++) {
    outfs << i << "," << nbCSH_tot[i] << endl;
  }
  outfs.close();

  outfilename = "nbC3S_tot.csv";
  outfs.open(outfilename.c_str());
  for (i = 0; i < 19; i++) {
    outfs << i << "," << nbC3S_tot[i] << endl;
  }
  outfs.close();

  outfilename = "nbCSH_woC3S.csv";
  outfs.open(outfilename.c_str());
  for (i = 0; i < 19; i++) {
    outfs << i << "," << nbCSH_woC3S[i] << endl;
  }
  outfs.close();

  outfilename = "nbC3S_woCSH.csv";
  outfs.open(outfilename.c_str());
  for (i = 0; i < 19; i++) {
    outfs << i << "," << nbC3S_woCSH[i] << endl;
  }
  outfs.close();

  // cout << endl << "exit Lattice::findInterfaces" << endl;
  // exit(0);
  // **************************** test for modifying affinities (contact angles)
  */

  return;
}

bool Lattice::hasPorousSolidNeighbor(const int siteID,
                                     const int neighborRange) {
  for (int i = 0; i < neighborRange; ++i) {
    if (site_[siteID].nb(i)->isPorousSolid())
      return (true);
  }
  return (false);
}

vector<int> Lattice::growPhase(vector<int> growPhaseIDVect,
                               vector<int> numSiteGrowVect,
                               vector<string> growPhNameVect, int &numadded_G,
                               bool &nucOK, int totalTRC) {

  //*** for controll
  int bcl = 0;
  int static trc_g;
  trc_g++;

  int i, j, jj;

  int growPhaseIDVectSize = growPhaseIDVect.size();
  vector<int> numChange(growPhaseIDVectSize, 0);
  vector<int> dim_isite(growPhaseIDVectSize, 0);
  vector<int> nucleated(growPhaseIDVectSize, 0);
  vector<bool> writeFirst(growPhaseIDVectSize, false);
  vector<int> numLeft = numSiteGrowVect; // numtotake
  vector<int> inGrowInterfacePos;
  vector<Isite> isite;
  vector<int> plist;
  int posGrPhId, plistSize;
  int phaseInGrowVect;

  bool needUpdate = false;
  bool needNucleation = true;

  int wmcIni, wmcEnd;
  int steWmc, stenbWmc, dwmcval;
  double rng, probRNG;

  int pid;
  int stId, pos, posGrowVect;
  int isitePos;
  int posVect;
  int siteID;
  int phaseID;
  int phaseIDn; // for nucleation

  Site *ste, *stenb;

  // int totInterfaceSize = 0;
  int numChangeTot = 0;
  int numLeftTot = 0;
  for (i = 0; i < growPhaseIDVectSize; i++) {
    numLeftTot += numSiteGrowVect[i];
    // totInterfaceSize += growthInterfaceSize_[growPhaseIDVect[i]];
  }
  if (numLeftTot == 0) {
    cout << endl << "     Lattice::growPhase error => numLeftTot = 0" << endl;
    cout << "     totalTRC/trc_g/bcl :  "
         << "   " << totalTRC << "/" << trc_g << "/" << bcl << endl;
    cout << endl << "     stop program" << endl;
    exit(0);
  }

  // growth probabilities based on affinities
  vector<structGrowVect> growthVector;
  structGrowVect growStruct;

  int posProbVect = 0;
  int aftyInt = 0;
  long int affSumInt = 0;
  for (i = 0; i < growPhaseIDVectSize; i++) {
    phaseID = growPhaseIDVect[i];
    isite = interface_[phaseID].getGrowthSites();
    dim_isite[i] = isite.size();
    for (jj = 0; jj < dim_isite[i]; jj++) {
      siteID = isite[jj].getId();
      aftyInt = isite[jj].getAffinityInt();
      affSumInt += aftyInt;
      growStruct.id = siteID;
      growStruct.posVect = i;
      growStruct.affinityInt = aftyInt;
      growthVector.push_back(growStruct);
      site_[siteID].setInGrowthVectorPos(phaseID, posProbVect);
      posProbVect++;
    }
  }
  double affSumDbl = static_cast<double>(affSumInt);
  int growthVectorSize = growthVector.size();

  // cout << endl
  //      << "    Lattice::growPhase GROW_INI totalTRC/trc_g/bcl/affSum " <<
  //      totalTRC
  //      << "/" << trc_g << "/" << bcl << "/" << affSum << endl;
  // cout << "      GROW_INI growPhaseIDVectSize = " << growPhaseIDVectSize
  //      << "   growthVectorSize = " << growthVectorSize
  //      << "   numLeftTot = " << numLeftTot
  //      << "   numChangeTot = " << numChangeTot << endl;
  cout << endl
       << "    Lattice::growPhase GROW_INI totalTRC/trc_g " << totalTRC << "/"
       << trc_g << " : growPhaseIDVectSize = " << growPhaseIDVectSize
       << "   growthVectorSize = " << growthVectorSize
       << "   numLeftTot = " << numLeftTot
       << "   numChangeTot = " << numChangeTot << endl;
  for (i = 0; i < growPhaseIDVectSize; i++) {
    phaseID = growPhaseIDVect[i];
    cout << "        GROW_INI for i = " << setw(3) << i
         << "  => phaseID phaseName count_ dim_isite numleft numchange  :  "
         << setw(3) << phaseID << "   " << setw(15) << left << growPhNameVect[i]
         << "   " << setw(8) << right << count_[phaseID] << "   " << setw(8)
         << dim_isite[i] << "   " << setw(8) << numLeft[i] << "   " << setw(8)
         << numChange[i] << endl;
  }
  cout << "        WAIT to grow " << numLeftTot << " voxels ..." << endl;
  cout.flush();

  bool noPosToGrow = false;

  if ((numLeftTot > 0) && (growthVectorSize == 0)) {
    int nucPhaseId;
    int numNuclei = 0;
    cout << endl << "    Lattice::growPhase growPhaseIDVectSize = 0" << endl;
    cout << "      => need an initial nucleation for at least one of the "
            "growing phases"
         << endl;
    cout << "      there is(are) " << growPhaseIDVectSize
         << " growing phase(s):" << endl;
    cout << "      ";
    for (int i = 0; i < growPhaseIDVectSize; i++) {
      cout << "  " << growPhaseIDVect[i];
    }
    while ((numLeftTot > 0) && (growthVectorSize == 0)) {
      nucPhaseId = -1;
      rng = callRNG();
      for (int i = 0; i < growPhaseIDVectSize; i++) {
        if (rng <= (i / static_cast<double>(growPhaseIDVectSize))) {
          if (numLeft[i] > 0) {
            nucPhaseId = i;
            break;
          }
        }
      }

      if (nucPhaseId > -1) {
        phaseID = growPhaseIDVect[nucPhaseId];
        cout << endl
             << "      one of them is chosen randomly for nucleation: "
             << phaseID << endl;
        dim_isite[nucPhaseId] = growthInterfaceSize_[phaseID];
        if (dim_isite[nucPhaseId] == 0) { // nucleation
          needNucleation = true;
          cout << endl
               << "    *** Lattice::growPhase - need nucleation for phaseID = "
               << phaseID << endl;
          cout << "      interface dimension dim_site[" << nucPhaseId
               << "] = " << dim_isite[nucPhaseId] << " while numLeft["
               << nucPhaseId << "] = " << numLeft[nucPhaseId] << endl;
          cout << "      => for this microPhase (" << growPhNameVect[nucPhaseId]
               << ") a number of " << numLeft[nucPhaseId]
               << " seed(s)/site(s) will be nucleated" << endl;

          // numNuclei = nucleatePhaseAff(phaseID, numLeft[i]);
          cout << "      WAIT for nucleation ..." << endl;
          cout.flush();

          numNuclei = nucleatePhaseRnd(phaseID, numLeft[nucPhaseId]);

          if (numNuclei < numLeft[nucPhaseId]) {
            noPosToGrow = true;
            nucOK = true;
          }

          nucleated[nucPhaseId] = numNuclei;

          numLeftTot = numLeftTot - numNuclei;
          numChangeTot = numChangeTot + numNuclei;

          numChange[nucPhaseId] = numChange[nucPhaseId] + numNuclei;
          numLeft[nucPhaseId] = numLeft[nucPhaseId] - numNuclei;

          writeFirst[i] = true;
          cout << endl
               << "    Lattice::growPhase GROW_END BY NUCLEATION for i = "
               << nucPhaseId << "   totalTRC/trc_g/bcl " << totalTRC << "/"
               << trc_g << "/" << bcl << endl;
          cout << "      GROW_END growPhaseIDVectSize = " << growPhaseIDVectSize
               << "   growthVectorSize = " << growthVectorSize
               << "   numLeftTot = " << numLeftTot
               << "   numChangeTot = " << numChangeTot << endl;
          cout << "        GROW_END phaseid count_ dim_isite numleft numchange "
                  " :  "
               << setw(3) << growPhaseIDVect[nucPhaseId] << "   " << setw(8)
               << right << count_[phaseID] << "   " << setw(8)
               << interface_[phaseID].getGrowthSites().size() << "   "
               << setw(8) << numLeft[nucPhaseId] << "   " << setw(8)
               << numChange[nucPhaseId] << endl;
          cout << "        WAIT to update the system ..." << endl;
          cout.flush();

          posProbVect = 0;
          affSumInt = 0;
          for (int j = 0; j < growPhaseIDVectSize; j++) {
            if (numLeft[j] > 0) {
              phaseID = growPhaseIDVect[j];
              isite = interface_[phaseID].getGrowthSites();
              dim_isite[j] = isite.size();

              for (int jj = 0; jj < dim_isite[j]; jj++) {
                aftyInt = isite[jj].getAffinityInt();
                affSumInt += aftyInt;
                growStruct.id = isite[jj].getId();
                growStruct.posVect = j;
                growStruct.affinityInt = aftyInt;
                growthVector.push_back(growStruct);
                site_[isite[jj].getId()].setInGrowthVectorPos(phaseID,
                                                              posProbVect);
                posProbVect++;
              }
            }
          }
          growthVectorSize = growthVector.size();
        }
      }
      if (noPosToGrow)
        break;
    }
    growPhaseIDVectSize = growPhaseIDVect.size();
    cout << endl
         << "    Lattice::growPhase growPhaseIDVectSize = " << growPhaseIDVectSize
         << endl;
    cout << "      => end initial nucleation" << endl;
    cout.flush();
  }

  if (!noPosToGrow) {
    while ((numLeftTot > 0) && (growthVectorSize >= 1)) {
      bcl++;

      // growth probabilities based on affinities
      // affAllPos - works with modified affinities in affinity_ table from
      // ChemicalSystem (all positives!)

      // calc probabilities & choose a site
      // affSum = 0.0;
      // for (j = 0; j < growthVectorSize; j++) {
      //   affSum += growthVector[j].affinity;
      // }

      rng = callRNG();
      if (affSumInt > 0) {
        affSumDbl = static_cast<double>(affSumInt);
        probRNG = growthVector[0].affinityInt / affSumDbl;
        if (rng <= probRNG) {
          isitePos = 0;
        } else {
          for (isitePos = 1; isitePos < growthVectorSize; isitePos++) {
            probRNG += (growthVector[isitePos].affinityInt / affSumDbl);
            if (rng <= probRNG)
              break;
          }
          if (isitePos == growthVectorSize)
            isitePos--;
        }
      } else {
        isitePos = static_cast<int>(rng * growthVectorSize);
      }

      ste = &site_[growthVector[isitePos].id];
      pid = ste->getMicroPhaseId(); // always ELECTROLYTEID !!
      posVect = growthVector[isitePos].posVect;
      phaseID = growPhaseIDVect[posVect]; // phase to grow

      if (pid != ELECTROLYTEID) {
        cout << endl
             << "Lattice::growPhase error: phaseid != ELECTROLYTEID" << endl;
        cout << "  phaseid totalTRC/trc_g/bcl count_ numLeftTot numChangeTot  "
                ":  "
             << phaseID << "   " << totalTRC << "/" << trc_g << "/" << bcl
             << "   " << count_[phaseID] << "   " << numLeftTot << "   "
             << numChangeTot << endl;
        cout << " posVect pid : " << posVect << "   " << pid << endl;
        cout << " ste.id numLeft numChange : " << ste->getId() << "   "
             << numLeft[posVect] << "   " << numChange[posVect] << endl;
        cout << "STOP" << endl;
        exit(0);
      }

      wmcIni = ste->getWmc0();

      plist = ste->getGrowthPhases(); // include phaseID i.e. faza care creste pe
                                      // acest voxel
      plistSize = plist.size();

      for (int j = 0; j < plistSize; j++) {
        for (int k = 0; k < growPhaseIDVectSize; k++) {
          if (growPhaseIDVect[k] == plist[j] && numLeft[k] > 0) {
            posGrPhId = ste->getInGrowthVectorPos(plist[j]);
            affSumInt -= growthVector[posGrPhId].affinityInt;
            if (posGrPhId != growthVectorSize - 1) {
              growthVector[posGrPhId] = growthVector[growthVectorSize - 1];
              pos = growthVector[posGrPhId].posVect;
              phaseInGrowVect = growPhaseIDVect[pos];
              site_[growthVector[posGrPhId].id].setInGrowthVectorPos(
                    phaseInGrowVect, posGrPhId);
            }
            growthVector.pop_back();
            growthVectorSize--;
            ste->setInGrowthVectorPos(plist[j], -1);
            break;
          }
        }
        removeGrowthSite_grow(ste, plist[j]);
      }
      if (affSumInt < 0) {
        cout << endl
             << "    Lattice::growPhase test:  affSumInt = " << affSumInt
             << "  - for growthVectorSize = " << growthVectorSize << endl;
        cout << "      phaseid totalTRC/trc_g/bcl count_ numLeftTot numChangeTot :  "
             << phaseID << "   " << totalTRC << "/" << trc_g << "/" << bcl
             << "   " << count_[phaseID] << "   " << numLeftTot << "   "
             << numChangeTot << endl;
        cout << "      posVect pid : " << posVect << "   " << pid << endl;
        cout << "      ste.id numLeft numChange : " << ste->getId() << "   "
             << numLeft[posVect] << "   " << numChange[posVect] << endl;
        cout << endl
             << "    Lattice::growPhase error: affSumInt < 0   =>   exit"
             << endl;
        exit(0);
      }

      ste->clearGrowth();
      setMicroPhaseId(ste, phaseID);

      ///
      /// Weighted mean curvature (wmc) is changed by the difference
      /// between the growing phase's porosity and the template's porosity.
      ///
      /// @todo Determine why the calculation works this way.
      ///

      // double por = chemSys_->getMicroPhasePorosity(phaseID);
      wmcEnd = microPhasePorosityInt_[phaseID];
      ste->setWmc0(wmcEnd);

      // dwmcval = chemSys_->getMicroPhasePorosity(phaseid) -
      //           chemSys_->getMicroPhasePorosity(pid);
      dwmcval = wmcEnd - wmcIni;
      ste->dWmc(dwmcval);

      ///
      /// Now that the site has been added, it is eligible for dissolution
      /// later on, so we add it to the list of dissolution sites.
      ///

      steWmc = ste->getWmc();

      if (steWmc > 0) {
        addDissolutionSite(ste, phaseID);
      }

      ///
      /// Update the wmc of the neighboring sites.  This can be done
      /// because the wmc is originally calculated within a box around each
      /// site, so any time the id of a site within that box changes, it
      /// will change the wmc of the site at the box's center.
      ///
      /// NN_NNN = NUM_NEAREST_NEIGHBORS +  NUM_SECONDNEAREST_NEIGHBORS;
      ///

      for (j = 0; j < NN_NNN; j++) {
        stenb = ste->nb(j);
        stenb->dWmc(dwmcval);
        stenbWmc = stenb->getWmc();
        if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
          inGrowInterfacePos = stenb->getInGrowInterfacePosVector();
          for (int k = FIRST_SOLID; k < numMicroPhases_; k++) {
            pos = inGrowInterfacePos[k]; // instenb->getInGrowInterfacePos(k);
            if (pos != -1) {
              aftyInt = chemSys_->getAffinityInt(k, phaseID);
              interface_[k].updateAffinityInt(pos, aftyInt);
              posGrowVect = stenb->getInGrowthVectorPos(k);
              if (posGrowVect != -1) {
                affSumInt += aftyInt;
                growthVector[posGrowVect].affinityInt += aftyInt;
              }

            } else {

              if (k == phaseID) {
                addGrowthSite(stenb, phaseID);

                if (numLeft[posVect] > 0) {
                  aftyInt = interface_[k].getAffinityInt(growthInterfaceSize_[k] - 1);
                  affSumInt += aftyInt;
                  growStruct.id = stenb->getId();
                  growStruct.posVect = posVect;
                  growStruct.affinityInt = aftyInt;
                  growthVector.push_back(growStruct);
                  stenb->setInGrowthVectorPos(phaseID, growthVectorSize);
                  growthVectorSize++;
                }
              }
            }
          }
        } else if ((stenbWmc == 0) &&
                   (stenb->getMicroPhaseId() > ELECTROLYTEID)) {
          removeDissolutionSite(stenb, stenb->getMicroPhaseId());
        }
      }

      for (j = 0; j < NUM_NEAREST_NEIGHBORS; j++) {
        stenb = ste->nb(j);
        if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
          for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
               phaseTmpl++) {
            // if (chemSys_->isGrowthTemplate(phaseTmpl, phaseID)) {
            if (growthTemplate_[phaseTmpl][phaseID]) {
              if (stenb->getInGrowInterfacePos(phaseTmpl) == -1) {
                addGrowthSite(stenb, phaseTmpl);
                for (int kk = 0; kk < growPhaseIDVectSize; kk++) {
                  if (growPhaseIDVect[kk] == phaseTmpl && numLeft[kk] > 0) {
                    aftyInt = interface_[phaseTmpl].getAffinityInt(
                          growthInterfaceSize_[phaseTmpl] - 1);
                    affSumInt += aftyInt;
                    growStruct.id = stenb->getId();
                    growStruct.posVect = kk;
                    growStruct.affinityInt = aftyInt;
                    growthVector.push_back(growStruct);
                    stenb->setInGrowthVectorPos(phaseTmpl, growthVectorSize);
                    growthVectorSize++;
                    break;
                  }
                }
              }
            }
          }
        }
      }

      numLeftTot--;
      numChangeTot++;

      numLeft[posVect]--;
      numChange[posVect]++;

      if (numLeft[posVect] < 0) {
        cout << endl
             << "error_numLeft totalTRC trc_g bcl numLeftTot numChangeTot  :  "
             << totalTRC << "   " << trc_g << "   " << bcl << "   " << numLeftTot
             << "   " << numChangeTot << endl;
        cout << "  phaseid totalTRC/trc_g/bcl count_ numLeftTot numChangeTot  "
                ":  "
             << phaseID << "   " << totalTRC << "/" << trc_g << "/" << bcl
             << "   " << count_[phaseID] << "   " << numLeftTot << "   "
             << numChangeTot << endl;
        cout << " posVect pid : " << posVect << "   " << pid << endl;
        cout << " ste.id numLeft numChange : " << ste->getId() << "   "
             << numLeft[posVect] << "   " << numChange[posVect] << endl;
        cout << "STOP" << endl;
        exit(0);
      }

      // growth probabilities based on affinities
      // test nucleatie - ini
      needUpdate = false;
      needNucleation = true;
      int numNuclei = 0;

      while (needNucleation) {
        needNucleation = false;
        for (i = 0; i < growPhaseIDVectSize; i++) {
          if (numLeft[i] > 0) {
            numNuclei = 0;
            phaseID = growPhaseIDVect[i];
            dim_isite[i] = growthInterfaceSize_[phaseID];
            if (dim_isite[i] == 0) { // nucleation
              needNucleation = true;
              cout
                  << endl
                  << "    *** Lattice::growPhase - need nucleation for phaseID = "
                  << phaseID << endl;
              cout << "      interface dimension dim_site[" << i
                   << "] = " << dim_isite[i] << " while numLeft[" << i
                   << "] = " << numLeft[i] << endl;
              cout << "      => for this microPhase (" << growPhNameVect[i]
                      << ") a number of " << numLeft[i]
                         << " seed(s)/site(s) will be nucleated" << endl;

              cout << endl
                   << "      *** Lattice::growPhase before nucleation -> "
                      "totalTRC/trc_g "
                   << totalTRC << "/" << trc_g
                   << " : growPhaseIDVectSize = " << growPhaseIDVectSize
                   << "   growthVectorSize = " << growthVectorSize
                   << "   numLeftTot = " << numLeftTot
                   << "   numChangeTot = " << numChangeTot << endl;
              for (int ij = 0; ij < growPhaseIDVectSize; ij++) {
                phaseIDn = growPhaseIDVect[ij];
                cout << "              GROW_INI for ij = " << setw(3) << ij
                     << "  => phaseID phaseName count_ dim_isite numleft "
                        "numchange  :  "
                     << setw(3) << phaseIDn << "   " << setw(15) << left
                     << growPhNameVect[ij] << "   " << setw(8) << right
                     << count_[phaseIDn] << "   " << setw(8) << dim_isite[ij]
                        << "   " << setw(8) << numLeft[ij] << "   " << setw(8)
                        << numChange[ij] << endl;
              }

              try {
                // numNuclei = nucleatePhaseAff(phaseID, numLeft[i]);
                numNuclei = nucleatePhaseRnd(phaseID, numLeft[i]);
              } catch (MicrostructureException mex) {
                cout << endl
                     << "    *** Lattice::growPhase - MicroEx from growPhase - "
                        "totalTRC/trc_g = "
                     << totalTRC << " / " << trc_g << endl;
                throw mex;
              }
              if (numNuclei < numLeft[i])
                noPosToGrow = true;

              nucleated[i] = numNuclei;

              numLeftTot = numLeftTot - numNuclei;
              numChangeTot = numChangeTot + numNuclei;

              numChange[i] = numChange[i] + numNuclei;
              numLeft[i] = numLeft[i] - numNuclei;

              // nucleated[i] = numLeft[i];
              //
              // numLeftTot = numLeftTot - numLeft[i];
              // numChangeTot = numChangeTot + numLeft[i];
              //
              // numChange[i] = numChange[i] + numLeft[i];
              // numLeft[i] = 0;
              writeFirst[i] = true;
              /*
            cout << endl
                 << "    Lattice::growPhase GROW_END BY NUCLEATION for i = "
                 << i << "   totalTRC/trc_g/bcl " << totalTRC << "/" << trc_g
                 << "/" << bcl << endl;
            cout << "      GROW_END growPhaseIDVectSize = "
                 << growPhaseIDVectSize
                 << "   growthVectorSize = " << growthVectorSize
                 << "   numLeftTot = " << numLeftTot
                 << "   numChangeTot = " << numChangeTot << endl;
            cout << "        GROW_END phaseid count_ dim_isite numleft "
                    "numchange  :  "
                 << setw(3) << growPhaseIDVect[i] << "   " << setw(8) << right
                 << count_[growPhaseIDVect[i]] << "   " << setw(8)
                 << dim_isite[i] << "   " << setw(8) << numLeft[i] << "   "
                 << setw(8) << numChange[i] << endl;
            cout.flush();
            */

              needUpdate = true;

              if (noPosToGrow) break;

            }
          }
        }
        if (noPosToGrow) break;
      }
      // test nucleatie - fin

      for (i = 0; i < growPhaseIDVectSize; i++) {
        if (numLeft[i] == 0) {
          if (writeFirst[i]) {
            continue;
          } else {
            writeFirst[i] = true;
            phaseID = growPhaseIDVect[i];
            dim_isite[i] = growthInterfaceSize_[phaseID];

            /*
          cout << "    Lattice::growPhase GROW_END for i = " << i << "
                  " totalTRC/trc_g/bcl "
               << totalTRC << "/" << trc_g << "/" << bcl << endl;
          cout << "      GROW_END growPhaseIDVectSize = " << growPhaseIDVectSize
               << "   growthVectorSize = " << growthVectorSize << "   numLeftTot
          = " << numLeftTot
               << "   numChangeTot = " << numChangeTot << endl;
          cout << "        GROW_END phaseid count_ dim_isite numleft numchange
          :  "
               << setw(3) << growPhaseIDVect[i] << "   "
               << setw(8) << right << count_[growPhaseIDVect[i]] << "   "
               << setw(8) << dim_isite[i] << "   "
               << setw(8) << numLeft[i] << "   "
               << setw(8) << numChange[i] << endl;
          cout.flush();
          */

            needUpdate = true;
          }
        }
      }

      if (needUpdate) {
        for (int i = 0; i < growthVectorSize; i++) {
          stId = growthVector[i].id;
          pid = growPhaseIDVect[growthVector[i].posVect];
          site_[stId].setInGrowthVectorPos(pid, -1);
          // site_[stId].resetInGrowthVectorPos();
        }
        growthVector.clear();
        posProbVect = 0;
        affSumInt = 0;
        for (int i = 0; i < growPhaseIDVectSize; i++) {
          if (numLeft[i] > 0) {
            phaseID = growPhaseIDVect[i];
            isite = interface_[phaseID].getGrowthSites();
            dim_isite[i] = isite.size();

            for (int jj = 0; jj < dim_isite[i]; jj++) {
              aftyInt = isite[jj].getAffinityInt();
              affSumInt += aftyInt;
              growStruct.id = isite[jj].getId();
              growStruct.posVect = i;
              growStruct.affinityInt = aftyInt;
              growthVector.push_back(growStruct);
              site_[isite[jj].getId()].setInGrowthVectorPos(phaseID, posProbVect);
              posProbVect++;
            }
          }
        }
        growthVectorSize = growthVector.size();
      }
      if (noPosToGrow) {
        nucOK = true;
        break;
      }
    }
  } // if (!noPosToGrow)

  numadded_G = numChangeTot;

  return (nucleated);
}

int Lattice::nucleatePhaseRnd(const int phaseID, const int numToNucleate) {

  int numRemaining = numToNucleate;

  vector<int> waterNucSites;
  vector<int> voidNucSites;
  int sizeWS = 0;
  int sizeVS = 0;
  double rng;
  int pos;
  int j, k;
  vector<int> seedID;

  int numThisPhase = 0;
  for (int i = 0; i < numSites_; i++) {
    if (site_[i].getMicroPhaseId() == phaseID)
      numThisPhase++;
  }
  int initNumThisPhase = numThisPhase;

  string namePhase = chemSys_->getMicroPhaseName(phaseID);
  cout << endl
       << "      Lattice::nucleatePhaseRnd INI for phaseID = " << phaseID
       << "   namePhase = " << namePhase << endl;
  cout << "        sites to add        : numToNucleate = " << numToNucleate
       << endl;
  cout << "        sites in the system : count_[" << phaseID
       << "] = " << count_[phaseID]
       << "   &   check numThisPhase = " << numThisPhase << endl;
  cout << "        pore sites          : count_[0] = " << count_[0] << endl;
  cout << "        electrolyte sites   : count_[1] = " << count_[1] << endl;
  cout << "        growthInterfaceSize_ = " << growthInterfaceSize_[phaseID]
       << endl;
  cout << "        dissolutionInterfaceSize_ = "
       << dissolutionInterfaceSize_[phaseID] << endl;

  seedID.clear();

  // Make a list of all saturated sites capable of hosting nucleation
  waterNucSites.clear();
  for (k = 0; k < numSites_; ++k) {
    if (site_[k].getMicroPhaseId() == ELECTROLYTEID) {
      waterNucSites.push_back(k);
    }
  }
  sizeWS = waterNucSites.size();
  int allPNS = sizeWS;

  if (numToNucleate > sizeWS) {
    // Now look for void sites that neighbor porous solids
    // We may need them if there are not enough saturated pore voxels

    voidNucSites.clear();
    for (k = 0; k < numSites_; ++k) {
      if ((site_[k].getMicroPhaseId() == VOIDID) &&
          (hasPorousSolidNeighbor(k, NUM_NEAREST_NEIGHBORS))) {
        voidNucSites.push_back(k);
      }
    }
    sizeVS = voidNucSites.size();
    allPNS += sizeVS;
  }

  cout << endl
       << "      Lattice::nucleatePhaseRnd numToNucleate = " << numToNucleate << endl;
  cout << "        eligible water sites : sizeWS = " << sizeWS << endl;
  cout << "        eligible void sites  : sizeVS = " << sizeVS << endl;
  cout << "        total eligible sites : allPNS = " << allPNS << endl;
  cout.flush();

  // Now randomly order the saturated pore voxels
  // Then append randomly orders void voxels that neighbor porous solids

  while ((sizeWS > 0) && (numRemaining > 0)) {
    rng = callRNG();
    pos = static_cast<int>(rng * sizeWS);
    if (pos == sizeWS)
      pos--;

    seedID.push_back(waterNucSites[pos]);

    numRemaining--;

    waterNucSites[pos] = waterNucSites[sizeWS - 1];
    waterNucSites.pop_back();

    sizeWS--;
    allPNS--;
  }

  while ((sizeVS > 0) && (numRemaining > 0)) {
    rng = callRNG();
    pos = static_cast<int>(rng * sizeVS);
    if (pos == sizeVS)
      pos--;

    seedID.push_back(voidNucSites[pos]);

    numRemaining--;

    voidNucSites[pos] = voidNucSites[sizeVS - 1];
    voidNucSites.pop_back();

    sizeVS--;
    allPNS--;
  }

  // We might still have run out of all possible nucleation sites
  /*// check!
  if (numRemaining > 0) {
    cout << endl
         << "     Lattice::nucleatePhaseRnd => requested nucleation "
            "numToNucleate is larger than all possible nucleation sites!"
         << endl;
    cout << "     numToNucleate > all nuc sites : " << numToNucleate << " > "
         << allPNS << endl;
    cout << endl
         << "     There is no room to nucleate => normal exit of the program"
         << endl;
    bool is_Error = false;
    throw MicrostructureException("Lattice", "nucleatePhaseRnd",
                                  "no room for nucleation", is_Error);
    // exit(0);
  }
  */

  int sizeSeedID = seedID.size();
  int sID, pid;
  Site *ste, *stenb;
  // double dwmcval, steWmc, stenbWmc;
  // double wmcIni, wmcEnd;
  int dwmcval, steWmc, stenbWmc;
  int wmcIni, wmcEnd;

  for (int i = 0; i < sizeSeedID; i++) {
    sID = seedID[i];
    ste = &site_[sID];
    pid = ste->getMicroPhaseId(); // always ELECTROLYTEID or possibly VOIDID !!

    wmcIni = ste->getWmc0();

    if (pid == ELECTROLYTEID) {
      removeGrowthSite_nucleation(ste);
    }
    setMicroPhaseId(ste, phaseID);

    ///
    /// Weighted mean curvature (wmc) is changed by the difference
    /// between the growing phase's porosity and the template's porosity.
    ///
    /// @todo Determine why the calculation works this way.
    ///

    // double por = chemSys_->getMicroPhasePorosity(phaseID);
    wmcEnd = microPhasePorosityInt_[phaseID];
    ste->setWmc0(wmcEnd);

    dwmcval = wmcEnd - wmcIni;
    ste->dWmc(dwmcval);

    ///
    /// Now that the site has been added, it is eligible for dissolution
    /// later on, so we add it to the list of dissolution sites.
    ///

    steWmc = ste->getWmc();

    if (steWmc > 0) {
      addDissolutionSite(ste, phaseID);
    }

    ///
    /// Update the wmc of the neighboring sites.  This can be done
    /// because the wmc is originally calculated within a box around each
    /// site, so any time the id of a site within that box changes, it
    /// will change the wmc of the site at the box's center.
    ///
    /// NN_NNN = NUM_NEAREST_NEIGHBORS +  NUM_SECONDNEAREST_NEIGHBORS;
    ///

    for (j = 0; j < NN_NNN; j++) {
      stenb = ste->nb(j);
      stenb->dWmc(dwmcval);
      stenbWmc = stenb->getWmc();

      if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
        // addGrowthSite(stenb, phaseID);
        if (stenb->getInGrowInterfacePos(phaseID) == -1) {
          // if (chemSys_->getAffinity(phaseID, phaseID) > 10) {//itself
          // affinity!
          addGrowthSite(stenb, phaseID);
          // }
        } else {
          int pos = stenb->getInGrowInterfacePos(phaseID);
          if (interface_[phaseID].getGrowthSitesId(pos) != stenb->getId()) {
            cout << endl
                 << "error Lattice::nucleatePhaseRnd interface_[phaseID]."
                    "getGrowthSitesId(pos) != stenb->getId() : "
                    "ste pid steWmc phaseID stenb"
                    " stenb_pid stenbWmc j : "
                 << ste->getId() << "   " << pid << "   " << ste->getWmc()
                 << "   " << phaseID << "   " << stenb->getId() << "   "
                 << stenb->getMicroPhaseId() << "   " << stenbWmc << "   " << j
                 << endl;
            cout << endl
                 << "pos/stenb->getInGrowInterfacePos(phaseID) : " << pos
                 << " / " << stenb->getInGrowInterfacePos(phaseID) << endl;

            cout << endl << "stop program" << endl;
            exit(1);
          }
        }

        if (j < NUM_NEAREST_NEIGHBORS) {
          for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
               phaseTmpl++) {
            // if (chemSys_->isGrowthTemplate(phaseTmpl, phaseID)) {
            if (growthTemplate_[phaseTmpl][phaseID]) {
              // addGrowthSite(stenb, nbgrowthtemp[jj]);
              // addGrowthSite(stenb, phaseTmpl);
              if (stenb->getInGrowInterfacePos(phaseTmpl) == -1) {
                addGrowthSite(stenb, phaseTmpl);
              } else {
                int pos = stenb->getInGrowInterfacePos(phaseTmpl);
                if (interface_[phaseTmpl].getGrowthSitesId(pos) !=
                    stenb->getId()) {
                  cout << endl
                       << "error Lattice:: growPhase interface_[phaseTmpl]."
                          "getGrowthSitesId(pos) != stenb->getId() : "
                          "ste pid steWmc phaseID stenb"
                          " stenb_pid stenbWmc j : "
                       << ste->getId() << "   " << pid << "   " << ste->getWmc()
                       << "   " << phaseID << "   "
                       << "   " << stenb->getId() << "   "
                       << stenb->getMicroPhaseId() << "   " << stenbWmc << "   "
                       << j << endl;
                  cout << endl
                       << "pos/phaseTmpl : " << pos << "   " << phaseID << endl;
                  cout << endl
                       << "stenb->getInGrowInterfacePos(phaseTmpl) : "
                       << stenb->getInGrowInterfacePos(phaseTmpl) << endl;

                  cout << endl << "stop program" << endl;
                  exit(1);
                }
              }
            }
          }
        }
      } else if ((stenbWmc == 0) &&
                 (stenb->getMicroPhaseId() > ELECTROLYTEID) &&
                 (stenb->getInDissInterfacePos() > -1)) {
        removeDissolutionSite(stenb, stenb->getMicroPhaseId());
      }
    }
  }
  /* // check!
  numThisPhase = 0;
  for (int i = 0; i < numSites_; i++) {
    if (site_[i].getMicroPhaseId() == phaseID)
      numThisPhase++;
  }

  if (numThisPhase != initNumThisPhase + numToNucleate) {
    cout << endl
         << "     error: numThisPhase != initNumThisPhase + numToNucleate -> "
         << numThisPhase << " != " << initNumThisPhase + numToNucleate << endl;
    cout << endl
         << "     error: initNumThisPhase = " << initNumThisPhase
         << "  &  numToNucleate = " << numToNucleate << endl;
    cout << endl
         << "     Lattice::nucleatePhaseRnd END for phaseID = " << phaseID
         << endl;

    cout << endl << "     exit" << endl;
    bool is_Error = false;
    throw MicrostructureException("Lattice", "nucleatePhaseRnd",
                                  "not all numToNucleate have been nucleated ",
                                  is_Error);
    // exit(1);
  }
  */

  if (numRemaining > 0) {
    cout << endl
         << "     Lattice::nucleatePhaseRnd : requested numToNucleate is"
            " larger than the number of all possible nucleation sites!"
         << endl;
    cout << "     numToNucleate > all possible nuc sites : "
         << numToNucleate << " > " << sizeSeedID << endl;
    cout << endl
         << "     There is no room for all nucleation events ("
         << numToNucleate << ")  =>  only " << sizeSeedID
         << " nucleation events have been realized  =>  normal exit"
         << endl;

    // bool is_Error = false;
    // throw MicrostructureException("Lattice", "nucleatePhaseRnd",
    //                               "no room for nucleation", is_Error);

  } else {

    numThisPhase = 0;
    for (int i = 0; i < numSites_; i++) {
      if (site_[i].getMicroPhaseId() == phaseID)
        numThisPhase++;
    }

    if (numThisPhase != initNumThisPhase + numToNucleate) {
      cout << endl
           << "     error: numThisPhase != initNumThisPhase + numToNucleate -> "
           << numThisPhase << " != " << initNumThisPhase + numToNucleate << endl;
      cout << endl
           << "     error: initNumThisPhase = " << initNumThisPhase
           << "  &  numToNucleate = " << numToNucleate << endl;
      cout << endl
           << "     Lattice::nucleatePhaseRnd END for phaseID = " << phaseID
           << endl;

      cout << endl << "     exit" << endl;
      bool is_Error = false;
      throw MicrostructureException("Lattice", "nucleatePhaseRnd",
                                    "not all numToNucleate have been nucleated ",
                                    is_Error);
      // exit(1);
    }

    cout << endl
         << "      Lattice::nucleatePhaseRnd END for phaseID = " << phaseID
         << endl;
    cout << "        numNewSites = " << numToNucleate << endl;
    cout << "        growthInterfaceSize_ = " << growthInterfaceSize_[phaseID]
         << endl;
    cout << "        dissolutionInterfaceSize_ = "
         << dissolutionInterfaceSize_[phaseID] << endl;
  }

  return(sizeSeedID);
}

int Lattice::nucleatePhaseAff(const int phaseID, const int numToNucleate) {

  int numRemaining = numToNucleate;
  struct localStruct {
    int id;
    int affInt;
    double prob;
  };
  localStruct un;
  vector<localStruct> waterNucSites, voidNucSites;
  int affInt;
  long int wAffSumInt, vAffSumInt;
  vector<Site *> localNb;
  int sizeWS = 0;
  int sizeVS = 0;
  double rng;
  int pos;
  int j, k;
  vector<int> seedID;

  int numSites = 0;
  for (int i = 0; i < numSites_; i++) {
    if (site_[i].getMicroPhaseId() == phaseID)
      numSites++;
  }
  int numSites0 = numSites;

  string namePhase = chemSys_->getMicroPhaseName(phaseID);
  cout << endl
       << "      Lattice::nucleatePhaseAff INI for phaseID = " << phaseID
       << "   namePhase = " << namePhase << endl;
  cout << "        sites to add        : numToNucleate = " << numToNucleate
       << endl;
  cout << "        sites in the system : count_[" << phaseID
       << "] = " << count_[phaseID] << "   &   check numSites = " << numSites
       << endl;
  cout << "        pore sites          : count_[0] = " << count_[0] << endl;
  cout << "        electrolyte sites   : count_[1] = " << count_[1] << endl;
  cout << "        growthInterfaceSize_ = " << growthInterfaceSize_[phaseID]
       << endl;
  cout << "        dissolutionInterfaceSize_ = "
       << dissolutionInterfaceSize_[phaseID] << endl;

  seedID.clear();
  /// JWB: Must be enough saturated pore voxels for all the
  /// nucleating phase to occupy.
  ///
  /// @todo Prefer saturated pore voxels, but then allow
  /// additional nucleation in void sites adjacent to saturated subvoxel
  /// porosity

  // Start filling in saturated pore voxels
  wAffSumInt = vAffSumInt = 0;
  waterNucSites.clear();
  for (k = 0; k < numSites_; ++k) {
    if (site_[k].getMicroPhaseId() == ELECTROLYTEID) {
      un.id = site_[k].getId();
      affInt = 0;
      localNb = site_[k].getNb();
      for (j = 0; j < NN_NNN; j++) {
        affInt += affinityInt_[phaseID][localNb[j]->getMicroPhaseId()];
      }
      un.affInt = affInt;
      waterNucSites.push_back(un);
      wAffSumInt += affInt;
    }
  }
  sizeWS = waterNucSites.size();
  int allPNS = sizeWS;

  cout << "      Lattice::nucleatePhaseAff -> sizeWS/wAffSumInt: " << sizeWS
       << " / " << wAffSumInt << endl;

  if (numToNucleate > sizeWS) {
    // We also will allow void sites that neighbor porous solids
    // because there are not enough saturated pore voxels

    voidNucSites.clear();
    for (k = 0; k < numSites_; ++k) {
      if ((site_[k].getMicroPhaseId() == VOIDID) &&
          (hasPorousSolidNeighbor(k, NUM_NEAREST_NEIGHBORS))) {
        un.id = site_[k].getId();
        affInt = 0;
        localNb = site_[k].getNb();
        for (j = 0; j < NUM_NEAREST_NEIGHBORS; j++) {
          affInt += affinityInt_[phaseID][localNb[j]->getMicroPhaseId()];
        }
        un.affInt = affInt;
        voidNucSites.push_back(un);
        vAffSumInt += affInt;
      }
    }
    sizeVS = voidNucSites.size();
    allPNS += sizeVS;
  }

  // Now randomly order the saturated pore voxels
  // Then append randomly orders void voxels that neighbor porous solids

  double wAffSumDbl;
  while (sizeWS > 0 && numRemaining > 0) {
    rng = callRNG();
    if (wAffSumInt > 0) {
      wAffSumDbl = static_cast<double>(wAffSumInt);
      waterNucSites[0].prob = (waterNucSites[0].affInt) / wAffSumDbl;
      for (k = 1; k < sizeWS; k++) {
        waterNucSites[k].prob =
            waterNucSites[k - 1].prob + (waterNucSites[k].affInt) / wAffSumDbl;
      }

      for (pos = 0; pos < sizeWS; pos++) {
        if (rng <= waterNucSites[pos].prob)
          break;
      }
    } else {

      pos = static_cast<int>(rng * sizeWS);
    }
    if (pos == sizeWS)
      pos--;

    seedID.push_back(waterNucSites[pos].id);

    numRemaining--;

    waterNucSites[pos] = waterNucSites[sizeWS - 1];
    waterNucSites.pop_back();

    sizeWS--;
    allPNS--;

    wAffSumInt = 0;
    for (k = 0; k < sizeWS; ++k) {
      wAffSumInt += waterNucSites[k].affInt;
    }
  }

  double vAffSumDbl;
  while (sizeVS > 0 && numRemaining > 0) {
    rng = callRNG();
    if (vAffSumInt > 0) {
      vAffSumDbl = static_cast<double>(vAffSumInt);
      voidNucSites[0].prob = (voidNucSites[0].affInt) / vAffSumDbl;
      for (k = 1; k < sizeVS; k++) {
        voidNucSites[k].prob =
            voidNucSites[k - 1].prob + (voidNucSites[k].affInt) / vAffSumDbl;
      }

      for (pos = 0; pos < sizeVS; pos++) {
        if (rng <= voidNucSites[pos].prob)
          break;
      }
    } else {

      pos = static_cast<int>(rng * sizeVS);
    }
    if (pos == sizeVS)
      pos--;

    seedID.push_back(voidNucSites[pos].id);

    numRemaining--;

    voidNucSites[pos] = voidNucSites[sizeVS - 1];
    voidNucSites.pop_back();

    sizeVS--;
    allPNS--;

    vAffSumInt = 0;
    for (k = 0; k < sizeVS; ++k) {
      vAffSumInt += voidNucSites[k].affInt;
    }
  }

  // We might still have run out of all possible nucleation sites
  /*// check!
  if (numRemaining > 0) {
    cout << endl
         << "      Lattice::nucleatePhaseAff => requested nucleation "
            "numToNucleate is larger than all possible nucleation sites!"
         << endl;
    cout << "        numToNucleate > all nuc sites : " << numToNucleate << " > "
         << allPNS << endl;
    cout << endl
         << "        There is no room to nucleate => normal exit of the program"
         << endl;
    bool is_Error = false;
    throw MicrostructureException("Lattice", "nucleatePhaseAff",
                                  "no room for nucleation", is_Error);
    // exit(0);
  }
  */

  int sizeSeedID = seedID.size();
  int sID, pid;
  Site *ste, *stenb;
  // double dwmcval;
  int dwmcval;
  // double steWmc, stenbWmc;
  int steWmc, stenbWmc;
  // double wmcIni, wmcEnd;
  int wmcIni, wmcEnd;

  for (int i = 0; i < sizeSeedID; i++) {
    sID = seedID[i];
    ste = &site_[sID];
    pid = ste->getMicroPhaseId(); // always ELECTROLYTEID or possibly VOIDID !!

    wmcIni = ste->getWmc0();

    if (pid == ELECTROLYTEID) {
      removeGrowthSite_nucleation(ste);
    }
    setMicroPhaseId(ste, phaseID);

    ///
    /// Weighted mean curvature (wmc) is changed by the difference
    /// between the growing phase's porosity and the template's porosity.
    ///
    /// @todo Determine why the calculation works this way.
    ///

    // double por = chemSys_->getMicroPhasePorosity(phaseID);
    wmcEnd = microPhasePorosityInt_[phaseID];
    ste->setWmc0(wmcEnd);

    dwmcval = wmcEnd - wmcIni;
    ste->dWmc(dwmcval);

    ///
    /// Now that the site has been added, it is eligible for dissolution
    /// later on, so we add it to the list of dissolution sites.
    ///

    steWmc = ste->getWmc();
    if (steWmc > 0) {
      addDissolutionSite(ste, phaseID);
    }

    ///
    /// Update the wmc of the neighboring sites.  This can be done
    /// because the wmc is originally calculated within a box around each
    /// site, so any time the id of a site within that box changes, it
    /// will change the wmc of the site at the box's center.
    ///
    /// NN_NNN = NUM_NEAREST_NEIGHBORS +  NUM_SECONDNEAREST_NEIGHBORS;
    ///

    for (j = 0; j < NN_NNN; j++) {
      stenb = ste->nb(j);
      stenb->dWmc(dwmcval);
      stenbWmc = stenb->getWmc();

      if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
        // addGrowthSite(stenb, phaseID);
        if (stenb->getInGrowInterfacePos(phaseID) == -1) {
          addGrowthSite(stenb, phaseID);
        } else {
          int pos = stenb->getInGrowInterfacePos(phaseID);
          if (interface_[phaseID].getGrowthSitesId(pos) != stenb->getId()) {
            cout << endl
                 << "error Lattice:: growPhase interface_[phaseID]."
                    "getGrowthSitesId(pos) != stenb->getId() : "
                    "ste pid steWmc phaseID stenb"
                    " stenb_pid stenbWmc j : "
                 << ste->getId() << "   " << pid << "   " << ste->getWmc()
                 << "   " << phaseID << "   " << stenb->getId() << "   "
                 << stenb->getMicroPhaseId() << "   " << stenbWmc << "   " << j
                 << endl;
            cout << endl
                 << "pos/stenb->getInGrowInterfacePos(phaseID) : " << pos
                 << " / " << stenb->getInGrowInterfacePos(phaseID) << endl;

            cout << endl << "stop program" << endl;
            exit(1);
          }
        }
        if (j < NUM_NEAREST_NEIGHBORS) {
          for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
               phaseTmpl++) {
            // if (chemSys_->isGrowthTemplate(phaseTmpl, phaseID)) {
            if (growthTemplate_[phaseTmpl][phaseID]) {
              // addGrowthSite(stenb, nbgrowthtemp[jj]);
              // addGrowthSite(stenb, phaseTmpl);
              if (stenb->getInGrowInterfacePos(phaseTmpl) == -1) {
                addGrowthSite(stenb, phaseTmpl);
              } else {
                int pos = stenb->getInGrowInterfacePos(phaseTmpl);
                if (interface_[phaseTmpl].getGrowthSitesId(pos) !=
                    stenb->getId()) {
                  cout << endl
                       << "error Lattice:: growPhase interface_[phaseTmpl]."
                          "getGrowthSitesId(pos) != stenb->getId() : "
                          "ste pid steWmc phaseID stenb"
                          " stenb_pid stenbWmc j : "
                       << ste->getId() << "   " << pid << "   " << ste->getWmc()
                       << "   " << phaseID << "   " << stenb->getId() << "   "
                       << stenb->getMicroPhaseId() << "   " << stenbWmc << "   "
                       << j << endl;
                  cout << endl
                       << "pos/phaseTmpl : " << pos << "   " << phaseID << endl;
                  cout << endl
                       << "stenb->getInGrowInterfacePos(phaseTmpl) : "
                       << stenb->getInGrowInterfacePos(phaseTmpl) << endl;

                  cout << endl << "stop program" << endl;
                  exit(1);
                }
              }
            }
          }
        }
      } else if ((stenbWmc == 0.0) &&
                 (stenb->getMicroPhaseId() > ELECTROLYTEID) &&
                 (stenb->getInDissInterfacePos() > -1)) {
        removeDissolutionSite(stenb, stenb->getMicroPhaseId());
      }
    }
  }

  /* // check!
  numSites = 0;
  for (int i = 0; i < numSites_; i++) {
    if (site_[i].getMicroPhaseId() == phaseID)
      numSites++;
  }
  if (numSites != numSites0 + numToNucleate) {
    cout << endl
         << "     error: numSites != numSites0 + numToNucleate -> " << numSites
         << " != " << numSites0 + numToNucleate << endl;
    cout << endl
         << "     error: numSites0 = " << numSites0
         << "  &  numToNucleate = " << numToNucleate << endl;
    cout << endl
         << "     Lattice::nucleatePhaseAff END for phaseID = " << phaseID
         << endl;

    cout << endl << "     exit" << endl;
    bool is_Error = false;
    throw MicrostructureException("Lattice", "nucleatePhaseAff",
                                  "not all numToNucleate have been nucleated ",
                                  is_Error);
    // exit(1);
  }
  */

  if (numRemaining > 0) {

    cout << endl
         << "     Lattice::nucleatePhaseAff : requested numToNucleate is"
            " larger than the number of all possible nucleation sites!"
         << endl;
    cout << "     numToNucleate > all possible nuc sites : "
         << numToNucleate << " > " << sizeSeedID << endl;
    cout << endl
         << "     There is no room for all nucleation events ("
         << numToNucleate << ")  =>  only " << sizeSeedID
         << " nucleation events have been realized  =>  normal exit"
         << endl;

    // bool is_Error = false;
    // throw MicrostructureException("Lattice", "nucleatePhaseAff",
    //                               "no room for nucleation", is_Error);
    // exit(0);

  } else {

    numSites = 0;
    for (int i = 0; i < numSites_; i++) {
      if (site_[i].getMicroPhaseId() == phaseID)
        numSites++;
    }

    if (numSites != numSites0 + numToNucleate) {
      cout << endl
           << "     error: numSites != numSites0 + numToNucleate -> "
           << numSites << " != " << numSites0 + numToNucleate << endl;
      cout << endl
           << "     error: numSites0 = " << numSites0
           << "  &  numToNucleate = " << numToNucleate << endl;
      cout << endl
           << "     Lattice::nucleatePhaseAff END for phaseID = " << phaseID
           << endl;

      cout << endl << "     exit" << endl;
      bool is_Error = false;
      throw MicrostructureException("Lattice", "nucleatePhaseAff",
                                    "not all numToNucleate have been nucleated ",
                                    is_Error);
      // exit(1);
    }

    cout << "      Lattice::nucleatePhaseAff END for phaseID = " << phaseID
         << endl;
    cout << "        numNewSites = " << numToNucleate << endl;
    cout << "        growthInterfaceSize_ = " << growthInterfaceSize_[phaseID]
         << endl;
    cout << "        dissolutionInterfaceSize_ = "
         << dissolutionInterfaceSize_[phaseID] << endl;
  }

  return(sizeSeedID);
}

vector<int> Lattice::dissolvePhase(vector<int> dissPhaseIDVect,
                                   vector<int> numSiteDissVect,
                                   vector<string> dissPhNameVect,
                                   int &numadded_D, int totalTRC) {

  //*** controll
  int bcl = 0;
  int static trc_d;
  trc_d++;

  int i, ii, jj;

  Site *ste, *stenb;
  int pid = 0;
  int isitePos, phaseID;
  int stId, nbid, nbpid;
  int posVect, posnb;
  double rng, probRNG;
  int wmcIni, wmcEnd, dwmcval;
  bool phaseid_exist;
  int stWmcInt;

  int dissPhaseIDVectSize = dissPhaseIDVect.size();
  vector<int> numChange(dissPhaseIDVectSize, 0);
  vector<int> dim_isite(dissPhaseIDVectSize, 0);
  vector<bool> writeFirst(dissPhaseIDVectSize, false);
  vector<int> numLeft = numSiteDissVect; // numtotake
  vector<Isite> isite;
  vector<int> growth_local;
  int grLocSize;

  int numChangeTot = 0;
  int numLeftTot = 0;
  for (i = 0; i < dissPhaseIDVectSize; i++) {
    numLeftTot += numSiteDissVect[i];
  }
  if (numLeftTot == 0) {
    cout << "Lattice::dissolvePhase error numLeftTot = 0" << endl;
    cout << "   totalTRC/trc_d/bcl :  "
         << "   " << totalTRC << "/" << trc_d << "/" << bcl << endl;
    cout << "stop program" << endl;
    exit(0);
  }

  // dissolution probabilities based on wmc
  vector<structDissVect> dissolutionVector;
  structDissVect dissStruct;
  int posProbVect = 0;
  long int sumWmcInt = 0;
  for (i = 0; i < dissPhaseIDVectSize; i++) {
    phaseID = dissPhaseIDVect[i];
    isite = interface_[phaseID].getDissolutionSites();
    dim_isite[i] = isite.size();
    for (jj = 0; jj < dim_isite[i]; jj++) {

      stId = isite[jj].getId();
      stWmcInt = site_[stId].getWmc();
      sumWmcInt += stWmcInt;
      dissStruct.id = stId;
      dissStruct.posVect = i;
      dissStruct.wmcInt = stWmcInt;

      dissolutionVector.push_back(dissStruct);
      // dissProbStruct(int id_i = 0, double instab_i = 0, int posVect_i = 0,
      // double prob_i = 0, double prob_01_i = 0) {}

      site_[stId].setInDissolutionVectorPos(posProbVect);
      posProbVect++;
    }
  }
  int dissolutionVectorSize = dissolutionVector.size();

  if (sumWmcInt < 0) {
    cout << endl
         << "    Lattice::dissolvePhase test_-1: sumWmcInt = " << sumWmcInt
         << "  - for dissolutionVectorSize = " << dissolutionVectorSize << endl;
    cout << "    Lattice::dissolvePhase test: totalTRC/trc_d/bcl  "
         << totalTRC << "/" << trc_d << "/" << bcl << endl;
    if (dissolutionVectorSize > 1) {
      cout << endl
           << "    Lattice::dissolvePhase error: dissolutionVectorSize > 1   =>   exit"
           << endl;
      exit(0);
    }
  }

  // cout << endl
  //      << "    Lattice::dissolvePhase DISS_INI totalTRC/trc_d/bcl/sumWmcInt " <<
  //      totalTRC
  //      << "/" << trc_d << "/" << bcl << "/" << sumWmcInt << endl;
  // cout << "      DISS_INI dissPhaseIDVectSize = " << dissPhaseIDVectSize
  //      << "   dissolutionVectorSize = " << dissolutionVectorSize
  //      << "   numLeftTot = " << numLeftTot
  //      << "   numChangeTot = " << numChangeTot << endl;
  cout << endl
       << "    Lattice::dissolvePhase DISS_INI totalTRC/trc_d " << totalTRC
       << "/" << trc_d << " : dissPhaseIDVectSize = " << dissPhaseIDVectSize
       << "   dissolutionVectorSize = " << dissolutionVectorSize
       << "   numLeftTot = " << numLeftTot
       << "   numChangeTot = " << numChangeTot << endl;
  for (i = 0; i < dissPhaseIDVectSize; i++) {
    phaseID = dissPhaseIDVect[i];
    // isite = interface_[phaseID].getDissolutionSites();
    // dim_isite = isite.size();
    cout << "        DISS_INI for i = " << setw(3) << i
         << "  => phaseID phaseName count_ dim_isite numleft numchange  :  "
         << setw(3) << phaseID << "   " << setw(15) << left << dissPhNameVect[i]
         << "   " << setw(8) << right << count_[phaseID] << "   " << setw(8)
         << dim_isite[i] << "   " << setw(8) << numLeft[i] << "   " << setw(8)
         << numChange[i] << endl;
  }
  cout << "        WAIT to dissolve " << numLeftTot << " voxels ..." << endl;
  cout.flush();

  // int isitePosError = 0;

  int callGEM = -1;

  double sumWmcDbl;
  while ((numLeftTot > 0) && (dissolutionVectorSize >= 1)) {
    try {
      bcl++;

      // dissolution probabilities based on wmc
      rng = callRNG();

      // sumWmc = 0.;
      // for (int i = 0; i < dissolutionVectorSize; i++) {
      //   sumWmc +=  dissolutionVector[i].wmc;
      // }
      if (sumWmcInt <= 0) {
        cout << endl
             << "    Lattice::dissolvePhase test_0: sumWmcInt = " << sumWmcInt
             << "   numLeftTot = " << numLeftTot
             << "   dissolutionVectorSize = " << dissolutionVectorSize
             << endl;
        cout << "    Lattice::dissolvePhase test: totalTRC/trc_d/bcl  "
             << totalTRC << "/" << trc_d << "/" << bcl << endl;
        cout << endl
             << "    Lattice::dissolvePhase error   =>   exit"
             << endl;
        exit(0);
      }

      sumWmcDbl = static_cast<double>(sumWmcInt);
      probRNG = dissolutionVector[0].wmcInt / sumWmcDbl;
      if (rng <= probRNG) {
        isitePos = 0;
      } else {
        for (isitePos = 1; isitePos < dissolutionVectorSize; isitePos++) {
          probRNG += (dissolutionVector[isitePos].wmcInt / sumWmcDbl);
          if (rng <= probRNG)
            break;
        }
        if (isitePos == dissolutionVectorSize) isitePos--;
      }

      ste = &site_[dissolutionVector[isitePos].id];
      pid = ste->getMicroPhaseId(); // intrebare pid diff phaseid ???
      posVect = dissolutionVector[isitePos].posVect;

      wmcIni = ste->getWmc0();

      sumWmcInt -= ste->getWmc();

      /*
      if (sumWmcInt < 0) {
        cout << endl
             << "    Lattice::dissolvePhase test_1: sumWmcInt = " << sumWmcInt
             << "  - for dissolutionVectorSize = " << dissolutionVectorSize << endl;
        cout << "    Lattice::dissolvePhase test: steId/pid/posVect/isitePos/wmcIni/ste->getWmc()  "
             << ste->getId() << "/" << pid << "/" << posVect << "/" << isitePos << "/" << wmcIni << "/"
             << ste->getWmc()
             << endl;
        cout << "    Lattice::dissolvePhase test: totalTRC/trc_d/bcl  "
             << totalTRC << "/" << trc_d << "/" << bcl << endl;
        if (dissolutionVectorSize > 1) {
          cout << endl
               << "    Lattice::dissolvePhase error: dissolutionVectorSize > 1   =>   exit"
               << endl;
          exit(0);
        }
        cout << "    =>  sumWmcInt = 0" << endl;
        sumWmcInt = 0;
      }

      if ((dissolutionVectorSize == 1) && (sumWmcInt != 0)) {
        cout << endl
             << "    Lattice::dissolvePhase test_2: sumWmcInt = " << sumWmcInt
             << "  - for dissolutionVectorSize = " << dissolutionVectorSize << endl;
        cout << "    Lattice::dissolvePhase test: steId/pid/posVect/isitePos/wmcIni  "
             << ste->getId() << "/" << pid << "/" << posVect << "/" << isitePos << "/" << wmcIni
             << endl;
        cout << "    Lattice::dissolvePhase test: totalTRC/trc_d/bcl  "
             << totalTRC << "/" << trc_d << "/" << bcl << endl;
        cout << endl
             << "    Lattice::dissolvePhase error: dissolutionVectorSize > 1   =>   exit"
             << endl;
        exit(0);
      }
      */

      if (ste->getInDissInterfacePos() == -1) {
        cout << endl
             << "    Lattice::dissolvePhase error: "
                "ste->getInDissInterfacePos() = -1 for dissolutionVectorSize = "
             << dissolutionVectorSize << endl;
        cout << "    Lattice::dissolvePhase error: steId/pid/posVect/isitePos  "
             << ste->getId() << "/" << pid << "/" << posVect << "/" << isitePos
             << endl;
        cout << "    Lattice::dissolvePhase error: totalTRC/trc_d/bcl  "
             << totalTRC << "/" << trc_d << "/" << bcl << endl;
        cout << "    Lattice::dissolvePhase error: exit" << endl;
        exit(0);
      }

      removeDissolutionSite(ste, pid);

      setMicroPhaseId(ste, ELECTROLYTEID);

      /// Weighted mean curvature (wmc) is changed by the difference
      /// between the growing phase's porosity and the template's porosity.
      ///
      /// @todo Determine why the calculation works this way.
      ///

      wmcEnd = electrolyteIntPorosity_; // ELECTROLYTEID;
      ste->setWmc0(wmcEnd);

      dwmcval = wmcEnd - wmcIni;
      ste->dWmc(dwmcval);

      for (i = 0; i < NN_NNN; i++) {
        stenb = ste->nb(i);
        stenb->dWmc(dwmcval);
        nbpid = stenb->getMicroPhaseId();

        if (nbpid > ELECTROLYTEID) {

          ///
          /// Now that the site has been dissolved, it is eligible for growth
          /// later on, so we add it to the list of growth sites for certain
          /// phases.
          ///

          if (stenb->getInDissInterfacePos() == -1) {
            addDissolutionSite(stenb, nbpid);

            if (stenb->getInDissolutionVectorPos() != -1) {
              // error
              nbid = stenb->getId();
              cout << endl
                   << "Lattice::dissolvePhase : "
                      "stenb->getInDissolutionVectorPos() != -1"
                   << endl;
              cout << endl
                   << "id/pid : " << ste->getId() << " / " << pid << endl;
              cout << endl
                   << "stenb->getInDissolutionVectorPos() = "
                   << stenb->getInDissolutionVectorPos() << endl;
              cout << endl
                   << "dissolutionVectorSize : " << dissolutionVectorSize
                   << endl;
              cout << endl
                   << "i/nbid/nbpid : " << i << " / " << nbid << " / " << nbpid
                   << endl;
              cout << endl
                   << "totalTRC trc_d bcl numLeftTot numChangeTot  :  "
                   << totalTRC << "   " << trc_d << "   " << bcl << "   "
                   << numLeftTot << "   " << numChangeTot << endl;
              cout << endl << "exit" << endl;
              exit(1);
            } else {

              for (int k = 0; k < dissPhaseIDVectSize; k++) {
                if (dissPhaseIDVect[k] == nbpid && numLeft[k] > 0) {
                  stWmcInt = stenb->getWmc();
                  sumWmcInt += stWmcInt;
                  stenb->setInDissolutionVectorPos(dissolutionVectorSize);
                  dissStruct.id = stenb->getId();
                  dissStruct.posVect = k;
                  dissStruct.wmcInt = stWmcInt;

                  dissolutionVector.push_back(dissStruct);
                  // dissProbStruct(int id_i = 0, double instab_i = 0,
                  //                int posVect_i = 0, double prob_i = 0,
                  //                double prob_01_i = 0) {}
                  dissolutionVectorSize++;
                  break;
                }
              }
            }
          } else {

            posnb = stenb->getInDissolutionVectorPos();
            if (posnb != -1) {
              sumWmcInt += dwmcval;
              // dissolutionVector[posnb].wmc += dwmcval;
              dissolutionVector[posnb].wmcInt = stenb->getWmc();
            }

            // for (int k = 0; k < dissPhaseIDVectSize; k++) {
            //   if (dissPhaseIDVect[k] == nbpid && numLeft[k] > 0){
            //     posnb = stenb->getInDissolutionVectorPos();
            //     sumWmc += dwmcval;
            //     stWmc = dwmcval + dissolutionVector[posnb].wmc;
            //   }
            // }
          }

          // addGrowthSite(ste, nbpid);
          if (ste->getInGrowInterfacePos(nbpid) == -1) {
            // if (chemSys_->getAffinity(nbpid, nbpid) > 10) {//itself affinity!
            addGrowthSite(ste, nbpid);
            // }
          }
          if (i < NUM_NEAREST_NEIGHBORS) {

            for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
                 phaseTmpl++) {
              // if (chemSys_->isGrowthTemplate(phaseTmpl, nbpid)) {
              if (growthTemplate_[phaseTmpl][nbpid]) {
                // addGrowthSite(ste, phaseTmpl);
                if (ste->getInGrowInterfacePos(phaseTmpl) == -1) {
                  addGrowthSite(ste, phaseTmpl);
                }
              }
            }
          }
        } else if (nbpid == ELECTROLYTEID) {
          nbid = stenb->getId();
          growth_local = site_[nbid].getGrowthPhases();
          grLocSize = growth_local.size();
          for (ii = 0; ii < grLocSize; ii++) {
            phaseid_exist = false;
            for (jj = 0; jj < NN_NNN; jj++) {
              nbpid = site_[nbid].nb(jj)->getMicroPhaseId();
              if (growth_local[ii] == nbpid) {
                phaseid_exist = true;
                break;
              }
              if (jj < NUM_NEAREST_NEIGHBORS) {
                // if (chemSys_->isGrowthTemplate(growth_local[ii], nbpid)) {
                if (growthTemplate_[growth_local[ii]][nbpid]) {
                  phaseid_exist = true;
                  break;
                }
              }
            }
            if (phaseid_exist == false) {
              removeGrowthSite_diss(stenb, growth_local[ii]);
            }
          }
        }
      }

      // update dissolutionVector & involved sites
      site_[dissolutionVector[isitePos].id].setInDissolutionVectorPos(-1);
      if (isitePos != dissolutionVectorSize - 1) {
        dissolutionVector[isitePos] =
            dissolutionVector[dissolutionVectorSize - 1];
        site_[dissolutionVector[isitePos].id].setInDissolutionVectorPos(
            isitePos);
      }
      dissolutionVector.pop_back();
      dissolutionVectorSize--;

      numLeftTot--;
      numChangeTot++;

      numLeft[posVect]--;
      numChange[posVect]++;

    } catch (out_of_range &oor) {
      EOBException ex("Lattice", "dissolvePhase", "site_", site_.size(), i);
      ex.printException();

      cout << endl << "Lattice::dissolvePhase error" << endl;
      cout << endl
           << "totalTRC trc_g bcl numLeftTot numChangeTot  :  " << totalTRC
           << "   " << trc_d << "   " << bcl << "   " << numLeftTot << "   "
           << numChangeTot << endl;
      cout << endl
           << "steId pid dissolutionVectorSize :  " << ste->getId() << "   "
           << pid << "   " << dissolutionVectorSize << endl;

      for (i = 0; i < dissPhaseIDVectSize; i++) {
        phaseID = dissPhaseIDVect[i];
        isite = interface_[phaseID].getDissolutionSites();
        dim_isite[i] = isite.size();
        cout << "        phaseid count_ dim_isite numleft numchange  :  "
             << phaseID << "   " << "   " << count_[phaseID] << "   "
             << dim_isite[i] << "   " << numLeft[i] << "   " << numChange[i]
             << endl;
      }
      cout << "stop program" << endl;
      exit(1);
    }

    // dissolution probabilities based on wmc
    for (int i = 0; i < dissPhaseIDVectSize; i++) {
      if (numLeft[i] == 0) {
        if (writeFirst[i]) {
          continue;
        } else {
          writeFirst[i] = true;

          phaseID = dissPhaseIDVect[i];
          dim_isite[i] = dissolutionInterfaceSize_[phaseID];
          /*
          cout << "    Lattice::dissolvePhase DISS_END for i = " << i
               << "   totalTRC/trc_d/bcl " << totalTRC << "/" << trc_d << "/"
               << bcl << endl;
          cout << "      DISS_END dissPhaseIDVectSize = " << dissPhaseIDVectSize
               << "   dissolutionVectorSize = " << dissolutionVectorSize
               << "   numLeftTot = " << numLeftTot
               << "   numChangeTot = " << numChangeTot << endl;
          cout << "        DISS_END phaseid count_ dim_isite numleft numchange "
                  ":  "
               << setw(3) << dissPhaseIDVect[i] << "   " << setw(8) << right <<
          count_[dissPhaseIDVect[i]]
               << "   " << setw(8) << dim_isite[i] << "   "
               << setw(8) << numLeft[i] << "   " << setw(8) << numChange[i] <<
          endl; cout.flush();
          */

          dissolutionVectorSize = dissolutionVector.size();
          for (int ij = 0; ij < dissolutionVectorSize; ij++) {
            site_[dissolutionVector[ij].id].setInDissolutionVectorPos(-1);
          }
          dissolutionVector.clear();
          sumWmcInt = 0.;
          posProbVect = 0;
          for (int i_i = 0; i_i < dissPhaseIDVectSize; i_i++) {
            if (numLeft[i_i] > 0) {
              phaseID = dissPhaseIDVect[i_i];
              isite = interface_[phaseID].getDissolutionSites();
              dim_isite[i_i] = isite.size();
              for (jj = 0; jj < dim_isite[i_i]; jj++) {

                stId = isite[jj].getId();
                stWmcInt = site_[stId].getWmc();
                sumWmcInt += stWmcInt;
                dissStruct.id = stId;
                dissStruct.posVect = i_i;
                dissStruct.wmcInt = stWmcInt;

                dissolutionVector.push_back(dissStruct);
                // dissProbStruct(int id_i = 0, double instab_i = 0, int
                // posVect_i = 0, double prob_i = 0, double prob_01_i = 0) {}

                site_[stId].setInDissolutionVectorPos(posProbVect);
                posProbVect++;
              }
            }
          }
          dissolutionVectorSize = dissolutionVector.size();
        }
      }
    }

    if (callGEM > -1) {
      dissolutionVectorSize = dissolutionVector.size();
      for (int i = 0; i < dissolutionVectorSize; i++) {
        site_[dissolutionVector[i].id].setInDissolutionVectorPos(-1);
      }
      break;
    }

    // dissolutionVectorSize = dissolutionVector.size();
    // cout << "     *** out totalTRC trc_d bcl :  " << totalTRC
    //      << "   " << trc_d << "   " << bcl << endl;
  }

  numadded_D = numChangeTot;
  vector<int> numleft;
  numleft.clear();
  callGEM = 0;
  for (int i = 0; i < dissPhaseIDVectSize; i++) {
    if (numLeft[i] > 0)
      callGEM++;
    numleft.push_back(numLeft[i]);
  }
  numleft.push_back(callGEM);
  return (numleft);
}

void Lattice::addDissolutionSite(Site *ste, int pid) {
  // try {
  ste->setInDissInterfacePos(dissolutionInterfaceSize_[pid]);
  interface_[pid].addDissolutionSite(ste);
  dissolutionInterfaceSize_[pid]++;
  ste->clearGrowth();
  // } catch (out_of_range &oor) {
  //   cout << endl << "EOB Lattice::addDissolutionSite pid = " << pid << " =>
  //   exit" << endl; exit(1);
  // }
}

void Lattice::addGrowthSite(Site *ste, int pid) {
  // try {
  interface_[pid].addGrowthSite(ste);
  ste->setInGrowInterfacePos(pid, growthInterfaceSize_[pid]);
  growthInterfaceSize_[pid]++;
  ste->addGrowthPhaseId(pid);
  // } catch (out_of_range &oor) {
  //   cout << endl << "EOB Lattice::addGrowthSite pid = " << pid << " => exit"
  //   << endl; exit(1);
  // }
}

void Lattice::removeDissolutionSite(Site *ste0, int pid) {
  // try {
  int pos1 = dissolutionInterfaceSize_[pid] - 1;
  int sid = interface_[pid].getDissolutionSitesId(pos1);
  Site *ste1 = &site_[sid];
  int pos0 = ste0->getInDissInterfacePos();
  interface_[pid].removeDissolutionSite(pos0, pos1);
  // if (pos0 != pos1)
  ste1->setInDissInterfacePos(pos0);
  ste0->setInDissInterfacePos(-1);
  dissolutionInterfaceSize_[pid]--;
  // } catch (out_of_range &oor) {
  //   cout << endl << "EOB Lattice::removeDissolutionSite pid = " << pid << "
  //   => exit" << endl; exit(1);
  // }
}

void Lattice::removeGrowthSite_nucleation(Site *ste0) {
  // try {
  int pos0, pos1, grPhId;
  Site *ste1;
  int sid;
  vector<int> plist = ste0->getGrowthPhases();
  int plsize = plist.size();
  for (int j = 0; j < plsize; j++) {
    grPhId = plist[j];
    pos1 = growthInterfaceSize_[grPhId] - 1;
    sid = interface_[grPhId].getGrowthSitesId(pos1);
    ste1 = &site_[sid];
    pos0 = ste0->getInGrowInterfacePos(grPhId);
    interface_[grPhId].removeGrowthSite(pos0, pos1);
    // if (pos0 != pos1)
    ste1->setInGrowInterfacePos(grPhId, pos0);
    ste0->setInGrowInterfacePos(grPhId, -1);

    growthInterfaceSize_[grPhId]--;
  }
  ste0->clearGrowth();
  // } catch (out_of_range &oor) {
  //   cout << endl << "EOB Lattice::removeGrowthSite_nucleation pid = " << pid
  //   << " => exit" << endl; exit(1);
  // }
}

void Lattice::removeGrowthSite_grow(Site *ste0, int grPhId) {
  // try {
  int pos0, pos1, sid;
  // Site *ste1;

  pos0 = ste0->getInGrowInterfacePos(grPhId);
  pos1 = growthInterfaceSize_[grPhId] - 1;
  sid = interface_[grPhId].getGrowthSitesId(pos1);

  interface_[grPhId].removeGrowthSite(pos0, pos1);
  // if (pos0 != pos1)
  site_[sid].setInGrowInterfacePos(grPhId, pos0);
  ste0->setInGrowInterfacePos(grPhId, -1);

  growthInterfaceSize_[grPhId]--;
  // } catch (out_of_range &oor) {
  //   cout << endl << "EOB Lattice::removeGrowthSite_grow grPhId = " << grPhId
  //   << " => exit" << endl; exit(1);
  // }
}

void Lattice::removeGrowthSite_diss(Site *ste0, int pid) {
  // try {
  int pos1 = growthInterfaceSize_[pid] - 1;
  int sid = interface_[pid].getGrowthSitesId(pos1);
  Site *ste1 = &site_[sid];
  int pos0 = ste0->getInGrowInterfacePos(pid);
  interface_[pid].removeGrowthSite(pos0, pos1);
  // if (pos0 != pos1)
  ste1->setInGrowInterfacePos(pid, pos0);
  ste0->setInGrowInterfacePos(pid, -1);

  ste0->removeGrowthSite(pid);

  growthInterfaceSize_[pid]--;
  // } catch (out_of_range &oor) {
  //   cout << endl << "EOB Lattice::removeGrowthSite_diss pid = " << pid << "
  //   => exit" << endl; exit(1);
  // }
}

double Lattice::changeSaturationState(double aqFracToEmpty) {

  double aqFracEmptied = 0.0;
  if (aqFracToEmpty > 0.0) {
    aqFracEmptied = emptyPorosity(aqFracToEmpty);
  } else if (aqFracToEmpty < 0.0) {
    aqFracEmptied = -fillPorosity(-aqFracToEmpty);
  }
  return (aqFracEmptied);
}

double Lattice::emptyPorosity(double aqFracToEmpty) {
  // Start removing any saturated pore voxels if they exist
  // How many saturated sites should be removed?
  double dNumsites = static_cast<double>(numSites_);
  int numToEmpty = static_cast<int>((aqFracToEmpty * numSites_) + 0.5);
  int actualFracEmptied = emptyVoxelPorosity(numToEmpty);
  double dFracEmptied = static_cast<double>(actualFracEmptied);
  double aqFracEmptied = dFracEmptied / dNumsites;

  // What is the remaining volume fraction to remove?
  aqFracToEmpty -= aqFracEmptied;
  if (aqFracToEmpty > 0.0) {
    calculatePoreSizeDistribution();
    double subVoxEmptied = emptySubVoxelPorosity(aqFracToEmpty);
    aqFracEmptied += subVoxEmptied;
    aqFracToEmpty -= subVoxEmptied;
  }
  return (aqFracEmptied);
}

double Lattice::emptySubVoxelPorosity(double aqFracToEmpty) {

  double aqFracEmptied = 0.0;
  double volfrac_filled = 0.0;
  int masterPoreSizeDistSize = masterPoreSizeDist_.size();
  for (int i = masterPoreSizeDistSize - 1; (i >= 0) && (aqFracToEmpty > 0.0);
       --i) {

    // volfrac_filled is the volume fraction of size-i saturated porosity on
    // a total microstructure volume basis
    volfrac_filled =
        masterPoreSizeDist_[i].volfrac * masterPoreSizeDist_[i].volfracsat;

    // Deduct this amount or aqFracToEmpty, whichever is less
    double newvfracfilled = volfrac_filled - aqFracToEmpty;
    if (newvfracfilled < 0.0)
      newvfracfilled = 0.0;
    masterPoreSizeDist_[i].volfracsat =
        newvfracfilled / masterPoreSizeDist_[i].volfrac;

    aqFracEmptied = (volfrac_filled - newvfracfilled);
    aqFracToEmpty -= aqFracEmptied;
  }
  return (aqFracEmptied);
}

double Lattice::fillPorosity(double aqFracToFill) {
  // Start with the smallest unsaturated sub-voxel pores
  calculatePoreSizeDistribution();
  double dNumsites = static_cast<double>(numSites_);
  double subVoxFilled = fillSubVoxelPorosity(aqFracToFill);
  aqFracToFill -= subVoxFilled;
  double aqFracFilled = subVoxFilled;

  // How many void voxels should be changed to electrolyte?
  if (aqFracToFill > 0.0) {
    int aqNumToFill = static_cast<int>((aqFracToFill * numSites_) + 0.5);
    int aqNumFilled = fillVoxelPorosity(aqNumToFill);
    double dAqNumFilled = static_cast<double>(aqNumFilled);
    aqFracFilled += (dAqNumFilled / dNumsites);
    aqFracToFill -= aqFracFilled;
  }

  return (aqFracFilled);
}

double Lattice::fillSubVoxelPorosity(double aqFracToFill) {

  double aqFracFilled = 0.0;
  double volfrac_empty = 0.0;
  int masterPoreSizeDistSize = masterPoreSizeDist_.size();
  for (int i = masterPoreSizeDistSize - 1; (i >= 0) && (aqFracToFill > 0.0);
       --i) {

    // volfrac_empty is the volume fraction of size-i empty porosity on
    // a total microstructure volume basis
    volfrac_empty = masterPoreSizeDist_[i].volfrac *
                    (1.0 - masterPoreSizeDist_[i].volfracsat);

    // Deduct this amount or aqFracToFill from the empty porosity, whichever is
    // less
    double newvfracempty = volfrac_empty - aqFracToFill;
    if (newvfracempty < 0.0)
      newvfracempty = 0.0;
    masterPoreSizeDist_[i].volfracsat =
        1.0 - (newvfracempty / masterPoreSizeDist_[i].volfrac);

    aqFracFilled = (volfrac_empty - newvfracempty);
    aqFracToFill -= aqFracFilled;
  }
  return (aqFracFilled);
}

int Lattice::emptyVoxelPorosity(int numToEmpty) {
  int maxsearchsize = 3;
  Site *stenb;

  ///
  /// Finding all potential VOID sites.
  ///
  /// @todo Consider removing some of the standard output, or setting a flag for
  /// it.
  ///

  // cout << "    Lattice::emptyPorosity - check for cyc = " << cyc
  //      << " :      numToEmpty = " << numToEmpty;
  // cout.flush();

  vector<int> distVect =
      findDomainSizeDistribution(ELECTROLYTEID, numToEmpty, maxsearchsize, 0);
  int distVectSize = distVect.size();

  // cout << "      distVect.size() = " << distVectSize << endl;
  // cout.flush();

  ///
  /// We want to empty the sites with the largest pore count
  ///

  // if (distVectSize < numToEmpty) {
  //   cout << endl
  //        << "    Lattice::emptyPorosity - Ran out of water in the system for "
  //           "cyc = "
  //        << cyc << "   distVect.size() < numToEmpty : " << distVectSize << " < "
  //        << numsites << endl;
  //   cout << "    Lattice::emptyPorosity -> normal end of the program" << endl;
  //   // exit(1);
  //   bool is_Error = false;
  //   throw MicrostructureException("Lattice", "emptyPorosity",
  //                                 "ran out of water in the system - normal end",
  //                                 is_Error);
  // }

  int numemptied = 0;
  int siteID;
  vector<int> grVect;
  int dimGrVect;
  Site *ste0, *ste1;
  int sid, pid, pos0, pos1;
  for (int ii = 0; ii < distVectSize; ii++) {
    siteID = distVect[ii];
    ste0 = &site_[siteID];
    setMicroPhaseId(siteID, VOIDID);
    grVect = site_[siteID].getGrowthPhases();
    dimGrVect = grVect.size();
    for (int i = 0; i < dimGrVect; i++) {
      pid = grVect[i];
      pos0 = ste0->getInGrowInterfacePos(pid);
      pos1 = growthInterfaceSize_[pid] - 1;
      sid = interface_[pid].getGrowthSitesId(pos1);
      ste1 = &site_[sid];
      interface_[pid].removeGrowthSite(pos0, pos1);
      ste1->setInGrowInterfacePos(pid, pos0);
      ste0->setInGrowInterfacePos(pid, -1);
      growthInterfaceSize_[pid]--;
    }
    site_[siteID].setWmc0(voidIntPorosity_);
    site_[siteID].dWmc(-electrolyteIntPorosity_);
    for (int i = 0; i < NN_NNN; i++) {
      stenb = site_[siteID].nb(i);
      stenb->dWmc(-electrolyteIntPorosity_);
      if ((stenb->getWmc() == 0) && (stenb->getMicroPhaseId() > ELECTROLYTEID))
        removeDissolutionSite(stenb, stenb->getMicroPhaseId());
    }
    site_[siteID].clearGrowth();
    numemptied++;
  }

  return (numemptied);
}

int Lattice::fillVoxelPorosity(int numToFill) {

  int maxsearchsize = 10;

  ///
  /// Finding all potential VOID sites.
  ///
  /// @todo Consider removing some of the standard output, or setting a flag for
  /// it.
  ///

  if (verbose_) {
    cout << "Lattice::fillPorosity Finding and sorting all "
         << "potential water sites ... ";
    cout.flush();
  }

  // cout << "    Lattice::fillPorosity - check for cyc = " << cyc
  //      << " :      numToFill = " << numToFill;
  // cout.flush();

  vector<int> distVect =
      findDomainSizeDistribution(VOIDID, numToFill, maxsearchsize, 1);
  int distVectSize = distVect.size();

  // cout << "      distVect.size() = " << distVectSize << endl;
  // cout.flush();

  ///
  /// We want to fill the sites with the smallest pore count
  ///

  if (distVectSize < numToFill) {
    cout << endl
         << "    Lattice::fillPorosity - not enough voids in the system to "
            "be filled <-> distVect.size() < numToFill : " << distVectSize
         << " < " << numToFill << endl;
    // exit(1);
    cout << "    Lattice::fillPorosity -> normal end of the program" << endl;
    bool is_Error = false;
    throw MicrostructureException(
        "Lattice", "fillPorosity",
        "not enough voidsr in the system - normal end", is_Error);
  }
  int numfilled = 0;
  int siteID;
  int nbpid;
  Site *stenb;
  for (int ii = 0; ii < distVectSize; ii++) {
    siteID = distVect[ii];
    setMicroPhaseId(siteID, ELECTROLYTEID);
    site_[siteID].setWmc0(electrolyteIntPorosity_);
    site_[siteID].dWmc(electrolyteIntPorosity_);
    site_[siteID].clearGrowth();
    for (int i = 0; i < NN_NNN; i++) {
      stenb = site_[siteID].nb(i);
      // wmcIni = stenb->getWmc();
      stenb->dWmc(electrolyteIntPorosity_);
      nbpid = stenb->getMicroPhaseId();
      if (nbpid > ELECTROLYTEID) {
        if (stenb->getInDissInterfacePos() == -1)
          addDissolutionSite(stenb, nbpid);
        if (site_[siteID].getInGrowInterfacePos(nbpid) == -1)
          addGrowthSite(&site_[siteID], nbpid);
        if (i < NUM_NEAREST_NEIGHBORS) {
          for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
               phaseTmpl++) {
            // if (chemSys_->isGrowthTemplate(phaseTmpl, nbpid)) {
            if (growthTemplate_[phaseTmpl][nbpid]) {
              if (site_[siteID].getInGrowInterfacePos(phaseTmpl) == -1) {
                addGrowthSite(&site_[siteID], phaseTmpl);
              }
            }
          }
        }
      }
    }
    numfilled++;
  }

  return (numfilled);
}

int Lattice::countBox(int boxsize, unsigned int siteid) {
  // string msg;
  int boxhalf = boxsize / 2;
  int nfound = 0;
  int ix, iy, iz, hx, hy, hz;

  int qxlo = site_[siteid].getX() - boxhalf;
  int qxhi = site_[siteid].getX() + boxhalf;
  int qylo = site_[siteid].getY() - boxhalf;
  int qyhi = site_[siteid].getY() + boxhalf;
  int qzlo = site_[siteid].getZ() - boxhalf;
  int qzhi = site_[siteid].getZ() + boxhalf;

  ///
  /// Make sure that the box bounds are all sensible given the
  /// boundary conditions
  ///

  qxlo += checkBC(qxlo, xdim_);
  qxhi += checkBC(qxhi, xdim_);
  qylo += checkBC(qylo, ydim_);
  qyhi += checkBC(qyhi, ydim_);
  qzlo += checkBC(qzlo, zdim_);
  qzhi += checkBC(qzhi, zdim_);

  ///
  /// Count the number of water or void sites in the box
  ///

  for (iz = qzlo; iz <= qzhi; iz++) {
    hz = iz + checkBC(iz, zdim_);
    for (iy = qylo; iy <= qyhi; iy++) {
      hy = iy + checkBC(iy, ydim_);
      for (ix = qxlo; ix <= qxhi; ix++) {
        hx = ix + checkBC(ix, xdim_);
        if (site_[getIndex(hx, hy, hz)].getMicroPhaseId() == ELECTROLYTEID ||
            site_[getIndex(hx, hy, hz)].getMicroPhaseId() == VOIDID) {
          nfound++;
        }
      }
    }
  }
  return (nfound);
}

void Lattice::setResolution(const double res) {
  ///
  /// Make sure that resolution is a valid value
  ///

  string msg;
  if (res <= 1.0e-9) {
    cout << endl;
    msg = "Lattice resolution <= 1 nm";
    throw DataException("Lattice", "setResolution", msg);
  }

  if (verbose_) {
    cout << "Lattice::setResolution Changing lattice resolution from ";
    cout << resolution_ << " to " << res << endl;
    cout.flush();
  }

  // This will now be used in kinetic models to scale surface to volume ratio
  resolution_ = res; // in meters

  areaPerFace_ = res * res;             // in m2 units
  volumePerVoxel_ = areaPerFace_ * res; // in m3 units
  return;
}

void Lattice::setNeighborhoodSA(void) {
  int xp, yp, zp;
  // double dist;

  int xc, yc, zc;
  int ind;

  for (int n = 0; n < numSites_; n++) {
    xc = site_[n].getX();
    yc = site_[n].getY();
    zc = site_[n].getZ();

    for (int k = -1; k <= 1; k++) {
      zp = zc + k;
      for (int j = -1; j <= 1; j++) {
        yp = yc + j;
        for (int i = -1; i <= 1; i++) {
          xp = xc + i;

          ind = getIndex(xp, yp, zp);

          site_[n].setNbSA(ind);
        }
      }
    }
  }
}

vector<int> Lattice::getNeighborhood(const int sitenum, const int size) {
  int xp, yp, zp;
  // double dist;

  vector<int> nh;

  int xc = site_[sitenum].getX();
  int yc = site_[sitenum].getY();
  int zc = site_[sitenum].getZ();

  if (size == 0) {
    nh.push_back(getIndex(xc, yc, zc));
    return nh;
  }

  for (int k = -size; k <= size; k++) {
    zp = zc + k;
    for (int j = -size; j <= size; j++) {
      yp = yc + j;
      for (int i = -size; i <= size; i++) {
        xp = xc + i;
        /*
        dist = static_cast<double>(((xc - xp) * (xc - xp)));
        dist += static_cast<double>(((yc - yp) * (yc - yp)));
        dist += static_cast<double>(((zc - zp) * (zc - zp)));

        if ((sqrt(dist) - 0.5) <= static_cast<double>(size)) {
            nh.push_back(getIndex(xp,yp,zp));
        }
        */
        nh.push_back(getIndex(xp, yp, zp));
      }
    }
  }

  return nh;
}

void Lattice::getNeighborhood(const int sitenum, vector<int> & ids,
                              vector<int> & phs, int & numW,
                              int & totPor) {

  ids = site_[sitenum].getNbSA();
  int size = ids.size(); // 27!
  int phId;
  // cout << endl << ">>>>> Lattice::getNeighborhood - ids.size() = " << size << endl;

  for (int i = 0; i < size; i++) {
    phId = site_[ids[i]].getMicroPhaseId();
    if (phId == ELECTROLYTEID)
      numW++;
    phs.push_back(phId);
    totPor += site_[ids[i]].getWmc0();
  }
}

int Lattice::getIndex(int ix, int iy, int iz) const {
  if (ix < 0) {
    if (BC != 1) {
      ix += xdim_;
    } else {
      ix = 0;
    }
  } else if (ix >= xdim_) {
    if (BC != 1) {
      ix -= xdim_;
    } else {
      ix = xdim_ - 1;
    }
  }

  if (iy < 0) {
    if (BC != 2) {
      iy += ydim_;
    } else {
      iy = 0;
    }
  } else if (iy >= ydim_) {
    if (BC != 2) {
      iy -= ydim_;
    } else {
      iy = ydim_ - 1;
    }
  }

  if (iz < 0) {
    if (BC != 3) {
      iz += zdim_;
    } else {
      iz = 0;
    }
  } else if (iz >= zdim_) {
    if (BC != 3) {
      iz -= zdim_;
    } else {
      iz = zdim_ - 1;
    }
  }

  return (ix + (xdim_ * iy) + ((xdim_ * ydim_) * iz));
}

int Lattice::changeMicrostructure(double time, const int simtype,
                                  vector<int> &vectPhNumDiff,
                                  vector<int> &vectPhIdDiff,
                                  vector<string> &vectPhNameDiff, int recalls,
                                  bool &stopPrg, int cyc) {

  // int i;
  // int numadded, numadded_actual;
  // unsigned int tpid;
  int cursites, newsites;
  int wcursites, wnewsites;
  // double td, tvol, tmass;
  int needRecallGEM = 0;
  vector<double> vol_next, vfrac_next;
  vector<int> netsites(numMicroPhases_, 0);
  vector<string> phasenames;

  vector<int> vectPhNumRecall(numMicroPhases_, 0);

  static int totalTRC, normalTRC, totalRecalls;
  totalTRC++;
  if (recalls == 0) {
    normalTRC++;
  } else {
    totalRecalls++;

    for (int i = 0; i < static_cast<int>(vectPhIdDiff.size()); i++){
      if (chemSys_->getIsKinetic()[vectPhIdDiff[i]])
        vectPhNumRecall[vectPhIdDiff[i]] = vectPhNumDiff[i];
    }

    vectPhNumDiff.clear();
    vectPhIdDiff.clear();
    vectPhNameDiff.clear();
  }

  // checkSite(8); cout << endl << " exit Lattice::changeMicrostructure " <<
  // endl;

  microPhasePorosityInt_.clear();
  microPhasePorosityInt_ = chemSys_->getMicroPhasePorosityInt();

  ///
  /// @todo This function is very large; consider breaking it into small pieces
  /// for ease of maintenance and readability.
  ///
  /// Assign time of this state to the lattice time_ variable
  ///

  time_ = time;

  ///
  /// The next block, if uncommented, reads prior expansion
  /// strain data from a file and loads it into the appropriate
  /// class members.
  ///

  ///
  /// Zero out the amount of water to add to the microstructure.
  ///
  /// @todo Find out why; this class variable is not used again.
  ///

  waterChange_ = 0.0;

  ///
  /// Next, determine how many sites of each phase need to be assigned
  /// The ChemicalSystem object calculates the new volume fractions, so we
  /// just load them up here from the ChemicalSystem object.
  ///

  vol_next = chemSys_->getMicroPhaseVolume();
  phasenames = chemSys_->getMicroPhaseName();
  int volNextSize = vol_next.size();

  // double waterVol_0 = vol_next[1];
  // cout << endl << "waterVol_0 = " << waterVol_0 << endl;

  if (verbose_) {
    cout << endl
         << "Lattice::changeMicrostructure Before adjustMicrostructureVolumes "
            "cyc = "
         << cyc << endl;
    for (int iii = 0; iii < volNextSize; ++iii) {
      cout << "Lattice::changeMicrostructure    Volume of " << phasenames[iii]
           << " = " << vol_next[iii] << " m3"
           << "  &  count_ = " << count_[iii] << endl;
    }
    cout.flush();
  }

  adjustMicrostructureVolumes(vol_next, volNextSize, cyc);

  if (verbose_) {
    cout << endl
         << "Lattice::changeMicrostructure After adjustMicrostructureVolumes "
            "cyc = "
         << cyc << endl;
    for (int iii = 0; iii < volNextSize; ++iii) {
      cout << "Lattice::changeMicrostructure    Volume of " << phasenames[iii]
           << " = " << vol_next[iii] << " m3"
           << "  &  count_ = " << count_[iii] << endl;
    }
  }

  vfrac_next = vol_next;

  adjustMicrostructureVolFracs(phasenames, vol_next, vfrac_next, volNextSize);

  /*
  double waterDensity = waterMolarMass_ / waterMolarVol_ / 1.0e6; // g/cm3
  double waterTotMass_0 = (waterVol_0/initialMicrostructureVolume_) *
  waterDensity * 100 / initSolidMass_; double waterTotMass = vfrac_next[1] *
  waterDensity * 100 / initSolidMass_; double waterTotMoles = waterTotMass /
  waterMolarMass_; // mol double waterTotMoles_0 = waterTotMass_0 /
  waterMolarMass_; // mol double waterDCMolesChemSys =
  chemSys_->getDCMoles(waterDCId_);
  // if (abs(waterDCMolesChemSys - waterTotMoles) >= 1.e-6) {
    cout << endl << "T1 waterDCMolesChemSys - waterTotMoles = "
         << waterDCMolesChemSys - waterTotMoles_0 << "   cyc = " << cyc << endl;
    cout << "   waterDCMolesChemSys = " << waterDCMolesChemSys << endl;
    cout << "   waterDCMolesLattice = " << waterTotMoles << endl;
    cout << "   waterTotMoles_0 = " << waterTotMoles_0 << endl;
    cout << "   initialMicrostructureVolume_ = "<< initialMicrostructureVolume_
  << endl;
  //   cout << endl << "exit" << endl;
  //   exit(0);
  // }
  */

  ///
  /// Calculate number of sites of each phase in next state
  ///

  // Remember we added two extra slots onto the end of vfrac_next
  // to hold the volume fraction of voxel-scale space and subvoxel pore space
  // But those two are not needed here anymore

  if (verbose_) {
    cout << endl << "Lattice::changeMicrostructure Calculating volume "
         << "of each phase to be added... cyc = " << cyc << endl;
    for (int i = 0; i < volNextSize; i++) {
      cout << "Lattice::changeMicrostructure **** i = " << i << "  Volume fraction["
           << phasenames[i] << "] in next state should be = " << vfrac_next[i];
      cout << ", or " << static_cast<int>(numSites_ * vfrac_next[i]) << " sites"
           << endl;
    }
    // cout << "Lattice::changeMicrostructure ****Volume fraction[voxel "
    //      << "pores] in next state"
    //      << " should be = " << voxelPoreVolumeFraction_ << endl;
    // cout << "Lattice::changeMicrostructure ****Volume fraction[subvoxel "
    //         "pores] "
    //         "in next state"
    //      << " should be = " << subvoxelPoreVolumeFraction_ << endl;
    cout.flush();
  }

  ///
  /// The next block is executed only if there will eventually be some
  /// sulfate attack during this simulation.
  ///

  if (simtype == SULFATE_ATTACK && time_ > sulfateAttackTime_) {
    // if (simtype == SULFATE_ATTACK) {
    /*
    ///
    ///  Normalize to get volume fractions and compute number
    ///  of sites of each phase needed
    ///
    ///  @todo Find out why we need to do all of this just because
    ///  there will eventually be sulfate attack.  Why not wait until
    ///  sulfate attack simulation actually starts?
    ///

    /// Used to be a bunch of sulfate atack stuff here
    */
    try {
      // Hard wiring for sulfate attack here
      // @todo Generalize to other phases

      ///
      /// Next block gets executed only if we are now simulating
      /// sulfate attack.
      ///

      // if (time_ >= sattackTime_) {

      // netsites.clear();
      // netsites.resize(numMicroPhases_, 0);
      // int growingSize = growingVectSA_.size();

      // for (int ii = 0; ii < sizeGrowingVectSA_; ++ii) {
      for (int i = FIRST_SOLID; i < volNextSize; i++) {
        cursites = count_[i];
        newsites = ceil(numSites_ * vfrac_next[i]);

        netsites[i] = newsites - cursites;
        // // pid[i] = i;
        // // if (i == growing[ii] && isFirst) {
        // if (i == growingVectSA_[ii] && cyc == 1) {
        //   netsites.at(i) = 0;
        //   count_.at(i) = newsites;
        // }
        // if (netsites[i] != 0) {
        //   cout << "  Lattice::changeMicrostructure SA ***for "
        //        << setw(15) << left << phasenames[i]
        //        << " => i: " << setw(3) << right << i
        //        << "   count_[i]: " << setw(8) << right << count_[i]
        //        << "   vfrac_next[i]: " << right << vfrac_next[i]
        //        << "   netsites in this state: " << setw(9) << right <<
        //        netsites[i]
        //        << "   cursites: " << setw(9) << right << cursites
        //        << "   newsites: " << setw(9) << right << newsites
        //        << "   -> phaseID: " << setw(3) << right
        //        << chemSys_->getMicroPhaseId(phasenames[i]) << endl;
        //   cout.flush();
        // }
      }
      // }
      if (cyc == 1) {
        for (int ii = 0; ii < sizeGrowingVectSA_; ++ii) {
          netsites.at(growingVectSA_[ii]) = 0;
          newsites = ceil(numSites_ * vfrac_next[growingVectSA_[ii]]);
          count_.at(ii) = newsites;
        }
      }
      cout << endl << "  Lattice::changeMicrostructure SA -ini-" << endl;
      for (int i = FIRST_SOLID; i < volNextSize; i++) {
        if (netsites[i] != 0) {
          cursites = count_[i];
          newsites = cursites + netsites[i];
          cout << "  Lattice::changeMicrostructure SA ***for " << setw(15)
               << left << phasenames[i] << " => i: " << setw(3) << right << i
               << "   count_[i]: " << setw(8) << right << count_[i]
               << "   vfrac_next[i]: " << right << vfrac_next[i]
               << "   sites to grow/dissolve: " << setw(9) << right
               << netsites[i] << "   cursites: " << setw(9) << right
               << cursites << "   newsites: " << setw(9) << right << newsites
               // << "   -> phaseID: " << setw(3) << right
               // << chemSys_->getMicroPhaseId(phasenames[i])
               << endl;
          // cout.flush();
        }
      }
      cout.flush();
      if (verbose_) {
        cout << "Lattice::changeMicrostructure Crystal-pressure "
             << "transform at time_ = " << time_ << endl;
        cout.flush();
      }

      ///
      /// The relevant stress-free molar volume ratios for sulfate
      /// attack phase transformations.
      ///

      vector<int> dissPhaseIDVect;
      vector<int> numSiteDissVect;
      vector<string> dissPhNameVect;
      vector<double> volumeRatio;

      int numTotSitesToDissolve = 0;

      int growPhId, shrinkid;
      int shrinkingSize;

      // sizeGrowingVectSA_ = 1 !!!
      cout << endl
           << "  Lattice::changeMicrostructure SA: sizeGrowingVectSA_ = "
           << sizeGrowingVectSA_ << endl;
      for (int i = 0; i < sizeGrowingVectSA_; ++i) {
        shrinkingSize = static_cast<int>(shrinking_[i].size());
        growPhId = growingVectSA_[i]; // AFt mPhId
        cout << "     i/growPhId/growPhName/shrinkingSize : " << i << " / "
             << growPhId << " / " << phasenames[growPhId] << " / "
             << shrinking_[i].size() << endl;
        for (int j = 0; j < shrinkingSize; j++) {
          shrinkid = shrinking_[i][j]; // C4ASH12 mPhId
          cout << "        j/shrinkId/shrinkName : " << j << " / " << shrinkid
               << " / " << phasenames[shrinkid] << endl;
        }
        if (netsites[growPhId] > 0) {
          for (int j = 0; j < shrinkingSize; j++) {
            shrinkid = shrinking_[i][j]; // C4ASH12 mPhId
            if (netsites[shrinkid] < 0) {
              dissPhaseIDVect.push_back(shrinkid);
              numSiteDissVect.push_back(-netsites[shrinkid]);
              dissPhNameVect.push_back(phasenames[shrinkid]);
              volumeRatio.push_back(volratios_[i][j]);
              numTotSitesToDissolve += -netsites[shrinkid];
            }
          }
        }
      }

      if (dissPhaseIDVect.size() > 0) { // for a single growid i.e. AFt!!!
        // call transform-dissolve function
        vector<int> correct_netsites;
        correct_netsites.resize(shrinkingSize + 1, 0);
        // correct_netsites[shrinkingSize + 1] -> AFt and all
        // correct_netsites[i] > 0
        int numadded_D = 0;

        correct_netsites = transformPhase(
            growPhId, netsites[growPhId], dissPhaseIDVect, numSiteDissVect,
            dissPhNameVect, volumeRatio, numadded_D, totalTRC);

        // correction of site numbers to be changed
        for (int i = 0; i < shrinkingSize; i++) {
          if (correct_netsites[i] < 0) {
            cout << endl
                 << "  Lattice::changeMicrostructure SA - transformPhase() error : "
                    "correct_netsites[" << i << "] < 0 (" << correct_netsites[i]
                 << " < 0) " << endl;
            cout << "    shrinkingSize = " << shrinkingSize << endl;
            cout << endl << "    exit" << endl;
            exit(0);
          } else {
            netsites[dissPhaseIDVect[i]] = - correct_netsites[i];
          }
        }
        netsites[growPhId] -= correct_netsites[shrinkingSize]; // for AFt

        if (correct_netsites[shrinkingSize + 1] > -1) {
          int phDiff = dissPhaseIDVect[needRecallGEM];
          string nameDiff = dissPhNameVect[needRecallGEM];
          int numDiff = count_[phDiff];

          cout << endl
               << "  Lattice::changeMicrostructure - transformPhase() anormal "
                  "end"
               << endl;
          cout << "    phDiff,nameDiff,numDiff,count_[phDiff] : " << phDiff
               << " , " << nameDiff << " , " << numDiff << " , "
               << count_[phDiff] << endl;
          cout << "    numTotSitesToDissolve = " << numTotSitesToDissolve
               << "  while numTotSitesDissolved = " << numadded_D << endl;
          cout
              << "    => recall GEM after (re)setDCLowerLimit according to the "
                 "system configuration (lattice)"
              << endl;

          return 0;
        } else {
          cout << endl
               << "  Lattice::changeMicrostructure SA => transformPhase()  "
                  "normal end"
               << endl;
          cout << "  Lattice::changeMicrostructure SA => shrinkingSize = "
               << shrinkingSize
               << "   &   correct_netsites[shrinkingSize + 1] = "
               << correct_netsites[shrinkingSize + 1] << endl;
        }
        cout << "  Lattice::changeMicrostructure SA -fin-" << endl << endl;
        for (int i = FIRST_SOLID; i < numMicroPhases_; i++) { // from i = 2 !!!
          if (netsites[i] != 0) {
            cursites = count_[i];
            newsites = cursites + netsites[i];
            cout << "  Lattice::changeMicrostructure ***for " << setw(15)
                 << left << phasenames[i] << " => i: " << setw(3) << right << i
                 << "   count_[i]: " << setw(8) << right << count_[i]
                 << "   vfrac_next[i]: " << right << vfrac_next[i]
                 << "   sites to grow/dissolve: " << setw(9) << right
                 << netsites[i] << "   cursites: " << setw(9) << right
                 << cursites << "   newsites: " << setw(9) << right << newsites
                 // << "   -> phaseID: " << setw(3) << right
                 // << chemSys_->getMicroPhaseId(phasenames[i])
                 << endl;
            // cout.flush();
          }
        }
        cout.flush();

      } else {
        cout << endl
             << "  Lattice::changeMicrostructure Sulfate Attack => "
                "dissPhaseIDVect.size() = "
             << dissPhaseIDVect.size() << "   => normal dissolution" << endl
             << endl;
      }

      //***********************************************************************

      /*
        @param alphaseid is the microstructure phase id of the Al-bearing phase
      to dissolve
        @param netsitesAlphaseid is the number of Al-bearing sites to dissolve
        @param ettrid is the microstructure phase id of ettringite
        @param netsitesEttrid is the number of ettringite sites to grow
        @return vector (na,ne) where na is the number of Al-bearing sites
      actually changed, and ne is the number of ettringite sites actually grown

      // vector<int> transform(int alphaseid, int netsitesAlphaseid, int ettrid,
      //                       int netsitesEttrid, double volumeratio);

              vector<int> numchanged;
              numchanged.clear();
              // int growid = growing[ii];
              // int shrinkid;
              double volrat;
              int shrinkingsize = shrinking[ii].size();

              for (int iii = 0; iii < shrinkingsize; ++iii) {
                numchanged.resize(2, 0);
                shrinkid = shrinking[ii][iii];
                volrat = volratios[ii][iii];
                if ((netsites.at(shrinkid) < 0) && (netsites.at(growid) > 0)) {
                  numchanged = transform(shrinkid, netsites.at(shrinkid),
      growid, netsites.at(growid), volrat);

                  netsites.at(shrinkid) += numchanged[0];
                  netsites.at(growid) -= numchanged[1];
                  if (verbose_) {
                    cout << "Lattice::changeMicrostructure netsites.at(" <<
      shrinkid
                         << ") is: " << netsites.at(shrinkid) << endl;
                    cout << "Lattice::changeMicrostructure netsites.at(" <<
      growid
                         << ") is: " << netsites.at(growid) << endl;
                    cout.flush();
                  }
                }
              }
      */
      // **************************************************************************
      // }
      // }
    }

    catch (out_of_range &oor) {
      throw EOBException(
          "Lattice", "changeMicrostructure",
          "phasenames or count_ or pid or netsites or vfrac_next",
          phasenames.size(), -100);
    }

  } else {
    ///
    ///  simtype != SULFATE_ATTACK
    ///    or
    ///  simtype == SULFATE_ATTACK && time <= sulfateAttackTime_
    ///
    ///  Normalize to get volume fractions and compute number
    ///  of sites of each phase needed.
    ///

    cout << endl
         << "  Lattice::changeMicrostructure - cyc = " << cyc
         << "  =>  initial values for count_[0]/count_[1] = " << count_[0]
         << " / " << count_[1] << endl;
    for (int i = FIRST_SOLID; i < volNextSize; i++) { // from i = 2 !!!
      if (vfrac_next[i] < 0) {
        cout << endl
             << "Lattice::changeMicrostructure  error - vfrac_next[i] < 0 for "
                "cyc = "
             << cyc << " :" << endl;
        for (int i = 1; i < numMicroPhases_; i++) {
          cout << "   " << i << " : phName/count_/vfrac_next: " << phasenames[i]
               << " / " << count_[i] << " / " << vfrac_next[i] << " / " << endl;
        }
        cout << endl << "end program" << endl;
        bool is_Error = false;
        throw MicrostructureException("Lattice", "changeMicrostructure",
                                      "vfrac_next[i] < 0", is_Error);
        // exit(0);
      }
      cursites = count_[i];

      if (vectPhNumRecall[i] > 0) {
        newsites = vectPhNumRecall[i];
      } else {
        newsites = ceil(numSites_ * vfrac_next[i]);
      // newsites = ceil((numSites_ - count_[0])  * vfrac_next[i]);
      }

      netsites[i] = newsites - cursites;
      if (netsites[i] != 0) {
        cout << "  Lattice::changeMicrostructure ***for " << setw(15) << left
             << phasenames[i] << " => i: " << setw(3) << right << i
             << "   count_[i]: " << setw(8) << right << count_[i]
             << "   vfrac_next[i]: " << right << vfrac_next[i]
             << "   sites to grow/dissolve: " << setw(9) << right
             << netsites[i] << "   cursites: " << setw(9) << right
             << cursites << "   newsites: " << setw(9) << right << newsites
             // << "   -> phaseID: " << setw(3) << right
             // << chemSys_->getMicroPhaseId(phasenames[i])
             << endl;
        // cout.flush();
      }
    }
    cout.flush();
  }

  ///
  /// Sort netsites in ascending order, except we will handle
  /// void space differently.

  int netsitesSize_;
  netsitesSize_ = netsites.size();
  vector<int> dissPhaseIDVect, growPhaseIDVect;
  vector<int> numSiteDissVect, numSiteGrowVect;
  vector<string> dissPhNameVect, growPhNameVect;

  int numTotSitesToDissolve = 0, numTotSitesToGrow = 0;

  for (int i = FIRST_SOLID; i < netsitesSize_; i++) {
    if (netsites[i] < 0) {
      dissPhaseIDVect.push_back(i);
      numSiteDissVect.push_back(-netsites[i]);
      dissPhNameVect.push_back(phasenames[i]);
      numTotSitesToDissolve += -netsites[i];
    } else if (netsites[i] > 0) {
      growPhaseIDVect.push_back(i);
      numSiteGrowVect.push_back(netsites[i]);
      growPhNameVect.push_back(phasenames[i]);
      numTotSitesToGrow += netsites[i];
    }
  }

  ///
  /// The next loop starts at FIRST_SOLID because we exclude void and water
  /// phases
  ///
  /// @todo Consider making the starting index more general
  ///

  // cout << "  Lattice::changeMicrostructure netsites.size() = " <<
  // netsitesSize_
  //      << " -> normalTRC/recalls (totalTRC/normalTRC/totalRecalls): " <<
  //      normalTRC
  //      << " / " << recalls << "   (" << totalTRC << " / " << normalTRC << " /
  //      "
  //      << totalRecalls << ")" << endl;
  cout << "  Lattice::changeMicrostructure - cyc = " << cyc
       << "  =>  number of phases to dissolve : " << setw(3) << right
       << dissPhaseIDVect.size() << endl;
  cout << "  Lattice::changeMicrostructure - cyc = " << cyc
       << "  =>  number of phases to grow     : " << setw(3) << right
       << growPhaseIDVect.size() << endl;

  int numadded_D = 0;
  int numadded_G = 0;
  vector<int> nucleatedPhases;
  // nucleatedPhases.clear();
  // nucleatedPhases.resize(growPhaseIDVect.size(), 0);

  bool nucEffectOK = false;
  try {
    if (dissPhaseIDVect.size() > 0) {
      needRecallGEM = 0;
      vector<int> numLeftDiss;
      numLeftDiss.clear();
      numLeftDiss = dissolvePhase(dissPhaseIDVect, numSiteDissVect,
                                  dissPhNameVect, numadded_D, totalTRC);
      int numLeftDissSize = numLeftDiss.size();
      needRecallGEM = numLeftDiss[numLeftDissSize - 1];

      if (needRecallGEM > 0) {
        cout << endl
             << "  Lattice::changeMicrostructure - anormal end for "
             << needRecallGEM << " phases" << endl;
        for (int i = 0; i < numLeftDissSize - 1; i++) {
          if (numLeftDiss[i] > 0) {
            vectPhIdDiff.push_back(dissPhaseIDVect[i]);
            vectPhNameDiff.push_back(dissPhNameVect[i]);
            vectPhNumDiff.push_back(count_[dissPhaseIDVect[i]]);
            cout << "    phDiff,nameDiff,numDiff,count_[phDiff] : "
                 << dissPhaseIDVect[i] << " , " << dissPhNameVect[i] << " , "
                 << numLeftDiss[i] << " , " << count_[dissPhaseIDVect[i]]
                 << endl;
            cout.flush();
          }
        }

        cout << "    numTotSitesToDissolve = " << numTotSitesToDissolve
             << "  while numTotSitesDissolved = " << numadded_D << endl;
        cout << "    => recall GEM after (re)setDCLowerLimit according to the "
                "system configuration (lattice)"
             << endl;

        return 0;
      }
    }

    if (growPhaseIDVect.size() > 0) {
      try {

        nucleatedPhases = growPhase(growPhaseIDVect, numSiteGrowVect,
                                    growPhNameVect, numadded_G, nucEffectOK, totalTRC);
      } catch (MicrostructureException mex) {
        cout << endl
             << "  Lattice::changeMicrostructure - MicroEx from growPhase - "
                "cyc = "
             << cyc << endl;
        throw mex;
      } catch (out_of_range &oor) {
        EOBException ex("Lattice", "changeMicrostructure", "after growth", 1,
                        0);
        ex.printException();
        cout << endl
             << "Lattice::changeMicrostructure -after growth excp- cyc = "
             << cyc << endl;
        cout.flush();
        exit(1);
      }
    }

    if ((numadded_D != numTotSitesToDissolve) ||
        ((numadded_G < numTotSitesToGrow) && !nucEffectOK)) {
      cout << endl << "      Lattice::changeMicrostructure error: " << endl;
      cout << "(numadded_D != numTotSitesToDissolve) || (numadded_G != "
              "numTotSitesToGrow)"
           << endl;
      cout << endl
           << "numadded_D = " << numadded_D
           << "       numTotSitesToDissolve  = " << numTotSitesToDissolve
           << endl;
      cout << endl
           << "numadded_G = " << numadded_G
           << "       numTotSitesToGrow  = " << numTotSitesToGrow
           << "       nucEffectOK        = " << nucEffectOK << endl;
      cout << endl
           << "   cyc / totalTRC / normalTRC : " << cyc << " / " << totalTRC
           << " / " << normalTRC << endl;
      cout << endl
           << "   recalls / totalRecalls : " << recalls << " / " << totalRecalls
           << endl;
      cout << endl << "   nucleatedPhases :" << endl;
      for (int i = 0; i < static_cast<int>(growPhaseIDVect.size()); i++) {
        cout << " i: " << i << "   phID: " << growPhaseIDVect[i]
             << "   phName: " << growPhNameVect[i]
             << "   toGrow: " << numSiteGrowVect[i]
             << "   grownByNucl: " << nucleatedPhases[i] << endl;
      }
      cout << endl << "stop program" << endl;
      exit(1);
    }
  } catch (out_of_range &oor) {
    throw EOBException("Lattice", "changeMicrostructure",
                       "dissPhaseIDVect.size() or dissPhaseIDVect.size()", 0,
                       -1000);
  }

  if (!nucEffectOK) {
    // nucEffectOK is false => program run further

    cursites = count_[VOIDID];
    newsites = static_cast<int>((numSites_ * vfrac_next[VOIDID]) + 0.5);
    double dcursites = static_cast<double>(cursites);
    double dnumsites = static_cast<double>(numSites_);
    double dnewsites = static_cast<double>(newsites);
    double dcurvoidfrac = dcursites / dnumsites;
    double dnewvoidfrac = dnewsites / dnumsites;
    double aqFracToEmpty = dnewvoidfrac - dcurvoidfrac;
    wcursites = count_[ELECTROLYTEID];
    wnewsites = wcursites - (newsites - cursites);

    int numEmptyFill = newsites - cursites;
    // int numEmpty = 0, numFill = 0;
    // numEmpty = emptyPorosity(newsites - cursites, cyc);
    cout << endl
         << "  Lattice::changeMicrostructure - cyc = " << cyc
         << " :  numEmptyFill = " << numEmptyFill
         << "   aqFracToEmpty = " << aqFracToEmpty
         << "   count_[0]/count_[1] = " << count_[0] << " / " << count_[1]
            << endl;
    // cout.flush();

    double aqFracEmptied = changeSaturationState(aqFracToEmpty);

    cout << "   =>   aqFracEmptied = " << aqFracEmptied
         << "   count_[0]/count_[1] = " << count_[0] << " / " << count_[1]
            << endl;

    if (verbose_) {
      cout << "Lattice::changeMicrostructure ***netsites[" << phasenames[VOIDID]
              << "] in this state = " << (newsites - cursites)
              << "; cursites = " << cursites << " and newsites = " << newsites
              << endl;
      cout << "Lattice::changeMicrostructure ***netsites["
           << phasenames[ELECTROLYTEID]
              << "] in this state = " << (wnewsites - wcursites)
              << "; cursites = " << wcursites << " and newsites = " << wnewsites
              << endl
              << endl;

      // When creating void from water, we should
      // update the target volume fraction of water even though
      // it is not used in any further calculations at this point

      cout << "Lattice::changeMicrostructure Target CAPIILARY WATER "
           << "volume fraction IS " << vfrac_next[ELECTROLYTEID] << endl;

      // vfrac_next[ELECTROLYTEID] -= (static_cast<double>(newsites -
      // cursites)/static_cast<double>(numSites_)); cout << "But WILL BE " <<
      // vfrac_next[ELECTROLYTEID] << " after creating void space" << endl;

      cout << "Lattice::changeMicrostructure Number VOXEL-SCALE VOXELS "
           << "actually emptied was:  " << numEmptyFill << endl;

      ///
      /// Report on target and actual mass fractions
      ///

      cout << "Lattice::changeMicrostructure "
           << "*******************************" << endl;
      cout.flush();
    }

    int totcount = 0;
    for (int i = 0; i < volNextSize; i++) {
      volumeFraction_[i] =
          static_cast<double>(count_[i]) / static_cast<double>(numSites_);
      totcount += count_[i];
      if (volumeFraction_[i] < 0) {
        cout << endl
             << "Lattice::changeMicrostructure  error - volumeFraction_[i] < 0 "
                "for cyc = "
             << cyc << " :" << endl;
        for (int i = 1; i < numMicroPhases_; i++) {
          cout << "   " << i << " : phName/count_/vfrac_next/volumeFraction_ : "
               << phasenames[i] << " / " << count_[i] << " / " << vfrac_next[i]
                  << " / " << volumeFraction_[i] << endl;
        }
        cout << endl << "end program" << endl;
        bool is_Error = false;
        throw MicrostructureException("Lattice", "changeMicrostructure",
                                      "volumeFraction_[i] < 0", is_Error);
        // exit(0);
      }
      if (verbose_) {
        cout << "  Lattice::changeMicrostructure Phase " << i
             << " Target volume fraction was " << vfrac_next[i]
                << " and actual is " << volumeFraction_[i] << ", and " << totcount
                << " of " << site_.size() << " sites claimed so far" << endl;
        cout.flush();
      }
    }
    if (totcount != numSites_) {
      cout << endl
           << "Lattice::changeMicrostructure error => totcount != numSites_ : "
           << totcount << " != " << numSites_ << endl;
      cout << "time,simtype : " << time << " , " << simtype << " , " << endl;
      cout << "recalls  : " << recalls << endl;
      cout << "cyc  : " << cyc << endl;
      cout << "stop program" << endl;
      exit(1);
    }

    if (abs(aqFracEmptied - aqFracToEmpty) > 1.0e-8) {
      cout << "    Lattice::changeMicrostructure - Ran out of water ";
      cout << "    Lattice::changeMicrostructure in the system for cyc = " << cyc
           << endl;
      cout << "        void frac change needed = " << aqFracToEmpty << endl;
      cout << "        void frac change accomplished = " << aqFracEmptied << endl;

      double water_volume = chemSys_->getMicroPhaseVolume(ELECTROLYTEID);
      double water_volfrac = water_volume / initialMicrostructureVolume_;
      cout << "        electrolyte volume fraction left = " << water_volfrac
           << endl;
      cout << "    Lattice::changeMicrostructure -> normal end of the program"
           << endl;
      cout << endl;
      bool is_Error = false;
      throw MicrostructureException("Lattice", "changeMicrostructure",
                                    "no more water in system", is_Error);
    }

    ///  This is a local variable and the value is never used.
    ///
    ///  @todo Why not eliminate this line completely?
    ///

    // double surfa = getSurfaceArea(chemSys_->getMicroPhaseId(CSHMicroName));

    /*
    int totPorTest = 0; // test!
    for (int i = 0; i < numSites_; i++) {
      totPorTest = site_[i].getWmc0();
      for (int j = 0; j < NN_NNN; j++) {
        // stenb = site_[i].nb(j);
        // totPorTest += stenb->getWmc0();
        totPorTest += site_[i].nb(j)->getWmc0();
      }
      if (totPorTest != site_[i].getWmc()) {
        cout << endl << ">>>>> Lattice::getNeighborhood - error for cyc = " << cyc
             << "  &  siteId = " << i << " : totPorTest != site_[" << i << "].getWmc() <=> "
             << totPorTest << " != " << site_[i].getWmc() << endl;
        cout << "       phaseId = " << site_[i].getMicroPhaseId() << endl;
        cout << "       wmc0 = " << site_[i].getWmc0() << endl;
        cout << endl << ">>>>> exit <<<<<" << endl;
        exit(0);
      }
    }
    */

    if (volumeFraction_[ELECTROLYTEID] <= 0.0) {

      cout << endl
           << "  Lattice::changeMicrostructure - volumeFraction_[ELECTROLYTEID] <= 0.0 : "
           << endl;
      cout << "    volumeFraction_[ELECTROLYTEID] = " << volumeFraction_[ELECTROLYTEID]
              << endl;
      cout << "    no VOXEL-SCALE electrolyte in the system" << endl;

    }

    cout << endl
         << "  Lattice::changeMicrostructure -normal end- cyc = " << cyc << endl;

    return 1;
  } else {
    // nucEffectOK is true => program must stop after writing a last set of output files
    // for SA must write .damage files too!
    stopPrg = true;
    return 1;
  }
}

void Lattice::adjustMicrostructureVolumes(vector<double> &vol, int volSize, int cyc) {
  int i = 0;

#ifdef DEBUG
  cout << "Lattice::adjustMicrostructureVolumes" << endl;
  cout.flush();
#endif
  subvoxelWaterVolume_ = 0;
  voxelWaterVolume_ = 0;
  waterVolume_ = vol[ELECTROLYTEID];

  voidVolume_ = 0;

  double solidvolume = 0;

  // Find the total system volume according to GEMS, in m3
  // units.  The individual microstructure phase volumes
  // have already been adjusted for subvoxel porosity in the
  // ChemicalSystem::calculateState function

  // The lattice has a fixed volume given by a fixed number
  // of voxels multiplied by the volume per voxel

  // The voxel-scale water volume is now calculated, during
  // which the subvoxel pore volume is calculated as well

  solidVolumeWithPores_ = 0.0;
  // for (i = 0; i < volSize; ++i) {
  //   if (i != ELECTROLYTEID && i != VOIDID)
  //     solidVolumeWithPores_ += vol[i];
  // }

  // if (solidVolumeWithPores_ <= 0.0)
  //   throw DataException("Lattice", "adjustMicrostructureVolumes",
  //                       "totvolume is NOT positive");

  // The current microstructure volume as predicted by GEMS
  // The initial microstructure volume is calculated the first
  // time that Lattice::changeMicrostructure is called
  // We currently insist remain the volume of the system, so if
  // the microstructure volume deviates from the system volume,
  // we add or subtract voxel-scale porosity to keep them equal

  // Calculate the subvoxel pore volume
  subvoxelPoreVolume_ = 0.0;
  // for (i = 0; i < volSize; ++i) {
  //   if (i != ELECTROLYTEID && i != VOIDID) {
  //     subvoxelPoreVolume_ += (vol[i] * chemSys_->getMicroPhasePorosity(i));
  //   }
  // }

  for (i = FIRST_SOLID; i < volSize; ++i) {
      solidVolumeWithPores_ += vol[i];
      // subvoxelPoreVolume_ += (vol[i] * chemSys_->getMicroPhasePorosity(i));
      subvoxelPoreVolume_ += (vol[i] * microPhasePorosityInt_[i]) / convFactDbl2IntPor_;
  }

  if (solidVolumeWithPores_ <= 0.0)
    throw DataException("Lattice", "adjustMicrostructureVolumes",
                        "solidVolumeWithPores_ is NOT positive");

  solidvolume = solidVolumeWithPores_ - subvoxelPoreVolume_;
  // microstructureVolume_ = chemSys_->getMicroVolume();

  // nonSolidVolume_: total porosity
  // nonSolidVolume_: pores at subvoxel level + voxels occupied by voids + voxels occupied by electrolyte
  if (initialMicrostructureVolume_ > solidvolume) {
    if (initialMicrostructureVolume_ > solidVolumeWithPores_) {
      nonSolidVolume_ = initialMicrostructureVolume_ - solidvolume;
    } else {
      nonSolidVolume_ = subvoxelPoreVolume_;
    }

  } else {
    // initialMicrostructureVolume_ <= solidvolume
    nonSolidVolume_ = subvoxelPoreVolume_;
  }

  // nonSolidVolume_ IS voxelVoidVolume_ + voxelWaterVolume_ +
  // subvoxelPoreVolume_ voxelPoreVolume_ = nonSolidVolume_ -
  // subvoxelPoreVolume_;

  // voxelPoreVolume_: voxels occupied by voids + voxels occupied by electrolyte
  if (initialMicrostructureVolume_ > solidVolumeWithPores_) {
    voxelPoreVolume_ = initialMicrostructureVolume_ - solidVolumeWithPores_;
  } else {
    voxelPoreVolume_ = 0;
  }

  // calcVoxelWaterVolume
  waterVolume_ = vol[ELECTROLYTEID];

  // Up to this point, vol[ELECTROLYTEID] is the total volume
  // of pore solution in the system regardless of whether it is
  // voxel or subvoxel
  // First need to trim off any extra water that lies outside the system

  // voidVolume_ : total empty porosity
  // voids at voxel level and empty pores at subvoxel level
  voidVolume_ = nonSolidVolume_ - waterVolume_;

  if (voidVolume_ < 0.0)
    waterVolume_ = nonSolidVolume_;

  voxelWaterVolume_ = waterVolume_ - subvoxelPoreVolume_;

  if (voxelWaterVolume_ < 0.0) {
    subvoxelWaterVolume_ = waterVolume_;
    voxelWaterVolume_ = 0;
  } else {
    subvoxelWaterVolume_ = subvoxelPoreVolume_;
  }

  voxelVoidVolume_ = voxelPoreVolume_ - voxelWaterVolume_;

  if (chemSys_->isSaturated()) { // System is saturated - test
    bool stPrg = false;
    if (abs(voxelVoidVolume_) > 1.e-15){
      cout << endl << endl
           << "   Lattice::adjustMicrostructureVolumes - error_1 - cyc = "
           << cyc << "   &   waterVolume_ = " << waterVolume_ << endl;
      cout << endl
           << "   voxelVoidVolume_ = " << voxelVoidVolume_
           << "   voxelPoreVolume_ = " << voxelPoreVolume_
           << "   voxelWaterVolume_ = " << voxelWaterVolume_
           << "   subvoxelWaterVolume_ = " << subvoxelWaterVolume_
           << "   subvoxelPoreVolume_ = " << subvoxelPoreVolume_
           << endl;
      stPrg = true;
    }
    voxelVoidVolume_ = 0.0; // check!

    // cout << endl << "   Lattice::adjustMicrostructureVolumes - cyc = " << cyc
    //      << "   subvoxelPoreVolume_ = " << subvoxelPoreVolume_
    //      << "   waterVolume_ = " <<  waterVolume_
    //      << "   abs = " << abs(waterVolume_ - subvoxelPoreVolume_) << endl;

    if ((waterVolume_ < subvoxelPoreVolume_) &&
        (abs(waterVolume_ - subvoxelPoreVolume_) > 1.e-15) &&
        (subvoxelPoreVolume_ > 0)) {
      cout << endl << endl
           << "   Lattice::adjustMicrostructureVolumes - error_2 - cyc = "
           << cyc << "   &   waterVolume_ = " << waterVolume_ << endl;
      cout << endl
           << "   voxelVoidVolume_ = " << voxelVoidVolume_
           << "   voxelPoreVolume_ = " << voxelPoreVolume_
           << "   voxelWaterVolume_ = " << voxelWaterVolume_
           << "   subvoxelWaterVolume_ = " << subvoxelWaterVolume_
           << "   subvoxelPoreVolume_ = " << subvoxelPoreVolume_
           << endl;
      stPrg = true;
    }
    if (stPrg) {
      cout << endl << "stop" << endl;
      exit(0);
    }
  }

  if (voxelVoidVolume_ < 0.0) {
    cout << endl << "   Lattice::adjustMicrostructureVolumes - cyc = " << cyc
         << " - test: voxelVoidVolume_ = " << voxelVoidVolume_ << endl;
    voxelVoidVolume_ = 0.0;
  }

  // cout << endl
  //      << "  Lattice::adjustMicrostructureVolumes - cyc = " << cyc
  //      << " test: voxelVoidVolume_/waterVolume_/subvoxelPoreVolume_/"
  //                 "voxelWaterVolume_/subvoxelWaterVolume_/voxelPoreVolume_ = "
  //      << voxelVoidVolume_ << " / " << waterVolume_ << " / " << subvoxelPoreVolume_ << " / "
  //      << voxelWaterVolume_ << " / " << subvoxelWaterVolume_ << " / " << voxelPoreVolume_
  //      << endl;

  // if (chemSys_->isSaturated()) { // System is saturated
  //   voxelVoidVolume_ = 0.0;
  //   if (waterVolume_ <  subvoxelPoreVolume_) {
  //     waterVolume_ = subvoxelPoreVolume_;
  //     voxelWaterVolume_ = 0;
  //   } else {
  //     voxelWaterVolume_ = waterVolume_ - subvoxelPoreVolume_;
  //   }
  //   // voxelWaterVolume_ = waterVolume_ - subvoxelPoreVolume_;
  //   subvoxelWaterVolume_ = subvoxelPoreVolume_;
  // }

  /// From here to the end of this time cycle, vol[ELECTROLYTEID] is the
  /// volume of voxel-scale pore water only, not the total volume of water.
  /// We can always recover the total volume of water by asking ChemicalSystem
  /// for it.

  voxelPoreVolume_ = voxelVoidVolume_ + voxelWaterVolume_;
  vol[ELECTROLYTEID] = voxelWaterVolume_;
  vol[VOIDID] = voxelVoidVolume_;

  /// If this has been coded correctly, all we have done is determine how
  /// much of the microstructure water should be assigned to voxel-scale pores
  /// And since only voxel-scale pore water is visible in the microstructure,
  /// we assign the voxel-scale pore water to ELECTROLYTEID and neglect any
  /// subvoxel electrolyte contribution

  if (verbose_) {
    cout << "Lattice::adjustMicrostructureVolumes" << endl;
    cout << "Lattice::adjustMicrostructureVolumesRESULTS:" << endl;
    cout << "Lattice::adjustMicrostructureVolumesAll water volume = "
         << waterVolume_ << " m3" << endl;
    cout << "Lattice::adjustMicrostructureVolumes All void volume = "
         << voidVolume_ << endl;
    cout << "Lattice::adjustMicrostructureVolumes Voxel water volume = "
         << voxelWaterVolume_ << " m3" << endl;
    cout << "Lattice::adjustMicrostructureVolumes Voxel void volume = "
         << voxelVoidVolume_ << " m3" << endl;
    cout << "Lattice::adjustMicrostructureVolumes Subvoxel water volume = "
         << subvoxelWaterVolume_ << " m3" << endl;
    cout << "Lattice::adjustMicrostructureVolumes Subvoxel pore volume = "
         << subvoxelPoreVolume_ << " m3" << endl;
    cout << "Lattice::adjustMicrostructureVolumes" << endl;
    cout.flush();
  }

  // cout << endl
  //      << "   Lattice::adjustMicrostructureVolumes - cyc = " << cyc
  //      << endl;
  // cout << "     waterVolume_/capPoreVol/capVoidVol/capWaterVol"
  //         "/subvoxelWaterVolume_ :   "
  //      << waterVolume_ << " / " << voxelPoreVolume_ << " / "
  //      << voxelVoidVolume_ << " / " << voxelWaterVolume_<< " / "
  //      << subvoxelWaterVolume_ << endl;

  ///
  /// End of manual adjustment
  ///

  return;
}

void Lattice::adjustMicrostructureVolFracs(vector<string> &names,
                                           const vector<double> vol,
                                           vector<double> &vfrac, int volSize) {
  int i = 0;

#ifdef DEBUG
  cout << "Lattice::adjustMicrostructureVolFracs" << endl;
  cout.flush();
#endif

  // Remember there are now two extra slots at the end of vfrac, just
  // like at the end of vol. These two extra slots hold the voxel-scale
  // pore space and the subvoxel pore space

  // vfrac.clear();
  // vfrac.resize(vol.size(), 0.0);

  // This model tries to reconcile two different volumes:
  //     (1) The total system volume according to GEMS (solid + liquid).
  //         This volume can change during the simulation.  If the
  //         GEMS volume is less than the initial volume of (solid + liquid),
  //         we make up the difference with void space.
  //     (2) The fixed RVE volume which is equal to the number of
  //         voxels multiplied by the volume per voxel
  //
  // We can either let the RVE change volume, which it generally does
  // prior to setting, or we can force it to remain constant, which
  // it approximately does after setting.
  //
  // As of 2023 April 9, we choose to keep the RVE volume fixed, so
  // the only way to reconcile its volume to the GEMS volume is to
  // ensure that the VOLUME FRACTIONS are the same in both.
  //
  // Microstructure volume consists of the sum of (1) the *apparent
  // volumes* of all the solid phases (solid plus subvoxel pores)
  // and (2) the voxel-scale pore volume

#ifdef DEBUG
  cout << "Lattice::adjustMicrostructureVolFracs - Calculated "
       << "total microstructure volume is " << chemSys_->getInitMicroVolume()
       << endl;
  cout.flush();
#endif

  // Calculate volume fractions based on total GEMS adjusted volume
  for (i = 0; i < volSize; ++i) {
    vfrac[i] = vol[i] / initialMicrostructureVolume_;
  }

  if (verbose_) {
    for (i = 0; i < volSize; ++i) {
      cout << "Lattice::adjustMicrostructureVolFracs volume fraction of "
           << names[i] << " should be " << vfrac[i] << ", ("
           << vol[i] << "/" << initialMicrostructureVolume_ // totmicvolume
           << ") and volume fraction NOW is "
           << static_cast<double>(count_[i]) / static_cast<double>(numSites_)
           << endl;
      cout.flush();
    }
  }

  // Calculate volume fraction of subvoxel porosity and
  // voxel-scale porosity whether saturated or not
  voxelPoreVolumeFraction_ =
      voxelPoreVolume_ / initialMicrostructureVolume_; // totmicvolume;
  subvoxelPoreVolumeFraction_ =
      subvoxelPoreVolume_ / initialMicrostructureVolume_; // totmicvolume;

  return;
}

// void Lattice::calcSubvoxelPoreVolume(vector<double> &vol) {

  // Find the total system volume according to GEMS, in m3
  // units.  The individual microstructure phase volumes
  // have already // been adjusted for subvoxel porosity in the
  // ChemicalSystem::calculateState function

  // The lattice has a fixed volume given by a fixed number
  // of voxels multiplied by the volume per voxel

  // This will hold the subvoxel pore volume (m3)

//   subvoxelPoreVolume_ = 0.0;
//   int size = vol.size();
  // double phi; // Holds the subvoxel porosity of a microstructurephase
//   for (int i = 0; i < size; ++i) {
//     if (i != ELECTROLYTEID && i != VOIDID) {
//       subvoxelPoreVolume_ += (vol[i] * chemSys_->getMicroPhasePorosity(i));
//     }
//   }

  // The total amount of non-solid space in the microstructure
// }

// void Lattice::calcSolidVolumeWithPores(vector<double> &vol) {

  // Find the total system volume according to GEMS, in m3
  // units.  The individual microstructure phase volumes
  // have already // been adjusted for subvoxel porosity in the
  // ChemicalSystem::calculateState function

  // The lattice has a fixed volume given by a fixed number
  // of voxels multiplied by the volume per voxel

  // This will hold the subvoxel pore volume (m3)

//   solidVolumeWithPores_ = 0.0;
//   int size = vol.size();
//   for (int i = 0; i < size; ++i) {
//     if (i != ELECTROLYTEID && i != VOIDID) {
//       solidVolumeWithPores_ += vol[i];
//     }
//   }
// }

// void Lattice::calcVoxelWaterVolume(vector<double> &vol) {
//   calcSubvoxelPoreVolume(vol);
//   voxelWaterVolume_ = vol[ELECTROLYTEID] - subvoxelPoreVolume_;
//   if (voxelWaterVolume_ < 0.0)
//     voxelWaterVolume_ = 0.0;
//   if (chemSys_->isSaturated()) {
//     voxelWaterVolume_ = vol[ELECTROLYTEID] - subvoxelPoreVolume_;
//   }
// }

// void Lattice::calcVoxelVoidVolume() {
//   voxelPoreVolume_ = nonSolidVolume_ - subvoxelPoreVolume_;
//   voxelVoidVolume_ = voxelPoreVolume_ - voxelWaterVolume_;
// }

void Lattice::calculatePoreSizeDistribution(void) {
  // Define some commonly used local variables

  int i; // vector indices
  // STEP 1:
  // The microstructure has a PSD
  // that is made up of the contributions of the PSDs
  // of each microstructure phase.
  //
  // Calculate the contribution that each phase makes
  // to the microstructure's overall PSD

  vector<vector<struct PoreSizeData>> phasePoreSizeDist =
      getPhasePoreSizeDistributions();

  // STEP 2:
  // Determine the maximum pore diameter (nm) from among all pore
  // size distributions defined
  // double maxUpperDiameter = getMaxPhasePoreDiameter(phasePoreSizeDist);

  // STEP 3:
  // Create a histogram of the combined pore size distributions
  // of all phases, taking account of the phase volume fractions also,
  // using pre-defined diameter bins
  vector<double> histogramDiameters = setPSDiameters();
  calcMasterPoreSizeDist(histogramDiameters, phasePoreSizeDist);

  // cout << "Lattice::calculatePoreSizeDistribution  resolution_ = "
  //      << resolution_ << " m" << endl;
  // cout << "Lattice::calculatePoreSizeDistribution  maxUpperDiameter = "
  //      << maxUpperDiameter << " nm" << endl;
  // cout << "Lattice::calculatePoreSizeDistribution Histogram:" << endl;
  // for (i = 0; i < masterPoreSizeDist_.size(); ++i) {
  //   cout << "diam = " << masterPoreSizeDist_[i].diam
  //        << ", volfrac = " << masterPoreSizeDist_[i].volfrac
  //        << ", volfracsat = " << masterPoreSizeDist_[i].volfracsat << endl;
  //   cout.flush();
  // }
  // cout << endl;
  // cout.flush();

  // At this point we have a complete pore volume distribution
  // for the microstructure.  We next need to determine
  // the saturation state.

  // When we arrive at this function, we know the following
  // information:
  //
  // Volume fraction of saturated voxel-scale pores on a total
  // microstructure volume basis = volumeFraction_[ELECTROLYTEID]
  //
  // Volume fraction of unsaturated voxel-scale pores
  // on a total microstructure volume basis
  // = volumeFraction_[VOIDID]
  //
  // Volume fraction of voxel-scale pores on a total
  // microstructure volume basis = voxelporevolumeFraction_
  //
  // Volume fraction of subvoxel pores on a total
  // microstructure volume basis = subvoxelporevolumeFraction_
  //
  // We do NOT yet know the fraction of subvoxel pores
  // that are saturated

  // water_volume is the absolute volume (m3) of aqueous
  // solution according to GEMS, whether in voxel-scale
  // pores, subvoxel pores, or squeezed outside the
  // system due to lack of porosity to contain it.

  double water_volume = chemSys_->getMicroPhaseVolume(ELECTROLYTEID);

  // microvol is the absolute volume (m3) of all the defined
  // microstructure phases, including subvoxel pores in
  // solid phases

  // double microvol = chemSys_->getMicroVolume();

  // microvol is the absolute initial volume (m3) of all the defined
  // microstructure phases, including subvoxel pores in
  // solid phases

  // double initmicrovol = chemSys_->getInitMicroVolume();

  // This is the volume fraction of liquid water whether
  // in voxel-scale or subvoxel porosity, on a total
  // microstructure volume basis

  // We may have had to push some of the system's
  // voxel-scale water outside the boundary of the
  // microstructure due to lack of space, so maybe
  // water_volume is larger than the actual volume
  // fraction of water available within the microstructure
  // We can check for that now, though.

  // double excesswater =
  //    water_volume - voxelPoreVolume_ - subvoxelPoreVolume_;

  // If excess water > 0 then we have a completely
  // saturated system and it is easy to write out
  // the water partitioning among pore sizes

  double water_volfrac = water_volume / initialMicrostructureVolume_;

  // This is the total porosity including voxel-scale
  // pore volume fraction and subvoxel pore volume
  // fraction, all on a total microstructure volume
  // basis

  // double pore_volfrac =
  //     voxelporevolumeFraction_ + subvoxelporevolumeFraction_;

  // Fill sub-voxel porosity with water for now
  double volfrac_avail;
  double volfrac_filled;
  int masterPoreSizeDistSize = masterPoreSizeDist_.size();

  bool isSealed = false;
  if (!(chemSys_->isSaturated()))
    isSealed = true;
  voxelPoreVolumeFractionSaturated_ = subvoxelPoreVolumeFractionSaturated_ =
      0.0;
  for (i = 0; i < masterPoreSizeDistSize; ++i) {
    volfrac_avail = masterPoreSizeDist_[i].volfrac;
    if (volfrac_avail > 0.0) {
      volfrac_filled = water_volfrac / volfrac_avail;
      if (volfrac_filled > 1.0)
        volfrac_filled = 1.0;
    } else {
      volfrac_filled = 0.0;
    }
    if (isSealed) {
      masterPoreSizeDist_[i].volfracsat = volfrac_filled;
      water_volfrac -= volfrac_avail;
      if (water_volfrac < 0.0)
        water_volfrac = 0.0;
    } else {
      masterPoreSizeDist_[i].volfracsat = volfrac_filled;
    }
    if (masterPoreSizeDist_[i].diam >= 1.0e9 * resolution_) {
      voxelPoreVolumeFractionSaturated_ +=
          (masterPoreSizeDist_[i].volfrac * masterPoreSizeDist_[i].volfracsat);
    } else {
      subvoxelPoreVolumeFractionSaturated_ +=
          (masterPoreSizeDist_[i].volfrac * masterPoreSizeDist_[i].volfracsat);
    }
  }

  return;
}

vector<vector<struct PoreSizeData>>
    Lattice::getPhasePoreSizeDistributions(void) {

  // STEP 1:
  // Get the input PSD of each microstructure phase
  // including void and electrolyte voxels
  vector<vector<struct PoreSizeData>> poreSizeDist =
      chemSys_->getPoreSizeDistribution();

  // Chemical system does not know the voxel resolution, so we
  // must set the diameter of electrolyte and void voxels manually
  poreSizeDist[VOIDID][0].diam = poreSizeDist[ELECTROLYTEID][0].diam =
      1.0e9 * resolution_; // in nanometer units

  // STEP 2: Get the internal porosity of each phase (dimensionless)
  // which is the internal volume fraction of pores within a phase
  vector<double> pore_volfrac = getPoreVolumeFractions();

  // STEP 3:
  // Partition the total porosity of each phase
  // among its different pore diameters
  //
  /// @note The PSDs must be indexed
  /// the same as the THAMES phase id numbers, or else change this
  /// code below.

  int numPhasesWithPSD = poreSizeDist.size();
  int numPSDRowsForThisPhase;
  for (int i = 0; i < numPhasesWithPSD; ++i) {
    numPSDRowsForThisPhase = poreSizeDist[i].size();
    for (int j = 0; j < numPSDRowsForThisPhase; ++j) {
      poreSizeDist[i][j].volfrac *= pore_volfrac[i];
    }
  }

  return (poreSizeDist);
}

double Lattice::getMaxPhasePoreDiameter(
    const vector<vector<struct PoreSizeData>> phasePoreSizeDist) {

  int numPSDRowsForThisPhase = 0;
  double upperDiameterForThisPhase = 0.0;
  double maxUpperDiameter = 0.0;
  for (int i = 0; i < static_cast<int>(phasePoreSizeDist.size()); ++i) {
    numPSDRowsForThisPhase = phasePoreSizeDist[i].size();
    if (numPSDRowsForThisPhase != 0) {
      upperDiameterForThisPhase =
          phasePoreSizeDist[i][numPSDRowsForThisPhase - 1].diam;
      if (upperDiameterForThisPhase > maxUpperDiameter)
        maxUpperDiameter = upperDiameterForThisPhase;
    }
  }

  return (maxUpperDiameter);
}

vector<double> Lattice::getPoreVolumeFractions(void) {
  // Get the volume fraction of the microstructure occupied
  // by internal pores of each phase
  // int numMicroPhases = volumeFraction_.size();
  vector<double> pore_volfrac(numMicroPhases_, 0.0);

  // The variable phi will hold the "subvoxel" porosity of a phase
  double phi = 0.0;
  double d_convFactDbl2IntPor = static_cast<double>(convFactDbl2IntPor_);

  for (int i = 0; i < numMicroPhases_; ++i) {
    // phi = chemSys_->getMicroPhasePorosity(i);
    phi = microPhasePorosityInt_[i] / d_convFactDbl2IntPor;
    pore_volfrac[i] = volumeFraction_[i] * phi;
  }
  return (pore_volfrac);
}

vector<double> Lattice::setPSDiameters(void) {
  vector<double> binnedDiams;
  binnedDiams.clear();
  double upperDiameter = 1.0e9 * resolution_;

  for (double xit = 0.3; pow(10.0, xit) < upperDiameter; xit += 0.05) {
    binnedDiams.push_back(pow(10.0, xit));
  }
  binnedDiams.push_back(upperDiameter);
  return (binnedDiams);
}

void Lattice::calcMasterPoreSizeDist(const vector<double> histogramDiameters,
         const vector<vector<struct PoreSizeData>> phasePoreSizeDist) {

  int numHistogramDiameters = histogramDiameters.size();
  masterPoreSizeDist_.clear();

  struct PoreSizeData dpsv;
  dpsv.diam = dpsv.volfrac = 0.0;
  dpsv.volfracsat = 1.0;

  int numPhasesWithPSD = phasePoreSizeDist.size();

  double diameter = 0.0;
  int i, j, k;
  int numPSDRowsForThisPhase = 0;
  vector<int> stoppedAt(numPhasesWithPSD, 0);
  subvoxelPoreVolumeFraction_ = voxelPoreVolumeFraction_ = 0.0;

  for (k = 0; k < numHistogramDiameters; ++k) {
    diameter = histogramDiameters[k];
    dpsv.diam = diameter;
    for (j = 0; j < numPhasesWithPSD; ++j) {
      numPSDRowsForThisPhase = phasePoreSizeDist[j].size();
      if (numPSDRowsForThisPhase > 0) {
        for (i = stoppedAt[j]; (phasePoreSizeDist[j][i].diam <= diameter) &&
                               (i < numPSDRowsForThisPhase);
             ++i) {
          dpsv.volfrac += phasePoreSizeDist[j][i].volfrac;
        }
        stoppedAt[j] = i;
      }
    }
    if (diameter < 1.0e9 * resolution_) {
      subvoxelPoreVolumeFraction_ += dpsv.volfrac;
    } else {
      voxelPoreVolumeFraction_ += dpsv.volfrac;
    }
    masterPoreSizeDist_.push_back(dpsv);
    dpsv.volfrac = 0.0;
  }

  return;
}

void Lattice::writePoreSizeDistribution(const double curtime,
                                        const string timeString) {

  // At this point we have a complete pore volume distribution
  // for the microstructure.  We next need to determine
  // the saturation state.

  // When we arrive at this function, we know the following
  // information:
  //
  // Volume fraction of saturated voxel-scale pores on a total
  // microstructure volume basis = volumeFraction_[ELECTROLYTEID]
  //
  // Volume fraction of unsaturated voxel-scale pores
  // on a total microstructure volume basis
  // = volumeFraction_[VOIDID]
  //
  // Volume fraction of voxel-scale pores on a total
  // microstructure volume basis = voxelPoreVolumeFraction_
  //
  // Volume fraction of subvoxel pores on a total
  // microstructure volume basis = subvoxelPoreVolumeFraction_
  //
  // We do NOT yet know the fraction of subvoxel pores
  // that are saturated

  // double water_volfrac = waterVolume_ / initialMicrostructureVolume_;

  // This is the total porosity including voxel-scale
  // pore volume fraction and subvoxel pore volume
  // fraction, all on a total microstructure volume
  // basis

  double pore_volfrac = voxelPoreVolumeFraction_ + subvoxelPoreVolumeFraction_;

  string ofileName(jobRoot_);

  ostringstream ostrT;
  ostrT << setprecision(3) << temperature_;
  string tempstr(ostrT.str());

  ofileName = ofileName + "_PoreSizeDistribution." + timeString + "."
                        + tempstr + "K.csv";

  ofstream out(ofileName.c_str());

  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writePoreSizeDistribution", ofileName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  // Write the header

  if (verbose_) {
    cout << "Time = " << curtime << " h" << endl;
    cout << "Voxel-scale pore volume fraction (>= " << 1.0e9 * resolution_
         << " nm) = " << voxelPoreVolumeFraction_ << endl;
    cout << "Saturated voxel-scale pore volume fraction = "
         << voxelPoreVolumeFractionSaturated_ << endl;
    cout << "Empty voxel-scale void volume fraction = "
         << voxelPoreVolumeFraction_ - voxelPoreVolumeFractionSaturated_
         << endl;
    cout << "Subvoxel volume fraction (< " << 1.0e9 * resolution_
         << " nm) = " << subvoxelPoreVolumeFraction_ << endl;
    cout << "Saturated subvoxel volume fraction (< " << 1.0e9 * resolution_
         << " nm) = " << subvoxelPoreVolumeFractionSaturated_ << endl;
    cout << "Empty subvoxel-scale void volume fraction = "
         << subvoxelPoreVolumeFraction_ - subvoxelPoreVolumeFractionSaturated_
         << endl;
    cout << "Total pore volume fraction = " << pore_volfrac << endl;
    cout << "Total void volume fraction = " << volumeFraction_[VOIDID] << endl;
    cout << "Pore size saturation data:" << endl;
    cout << "Diameter (nm),Volume Fraction,Fraction Saturated" << endl;
    cout << "Masterporevolume size = " << masterPoreSizeDist_.size() << endl;
    cout.flush();
  }

  out << "Time = " << curtime << " h,," << endl;
  out << "Voxel-scale pore volume fraction (>= " << 1.0e9 * resolution_
      << " nm) = " << voxelPoreVolumeFraction_ << ",," << endl;
  out << "Saturated voxel-scale pore volume fraction = "
      << voxelPoreVolumeFractionSaturated_ << ",," << endl;
  out << "Empty voxel-scale void volume fraction = "
      << voxelPoreVolumeFraction_ - voxelPoreVolumeFractionSaturated_ << ",,"
      << endl;
  out << "Subvoxel volume fraction (< " << 1.0e9 * resolution_
      << " nm) = " << subvoxelPoreVolumeFraction_ << ",," << endl;
  out << "Saturated subvoxel volume fraction (< " << 1.0e9 * resolution_
      << " nm) = " << subvoxelPoreVolumeFractionSaturated_ << ",," << endl;
  out << "Empty subvoxel-scale void volume fraction = "
      << subvoxelPoreVolumeFraction_ - subvoxelPoreVolumeFractionSaturated_
      << ",," << endl;
  out << "Total pore volume fraction = " << pore_volfrac << ",," << endl;
  out << "Total void volume fraction = " << volumeFraction_[VOIDID] << ",,"
      << endl;
  out << "Pore size saturation data:" << ",," << endl;
  out << "Masterporevolume size = " << masterPoreSizeDist_.size() << ",,"
      << endl;
  out << "Diameter (nm),Volume Fraction,Fraction Saturated" << endl;

  int masterPoreSizeDistSize = masterPoreSizeDist_.size();
  for (int i = 0; i < masterPoreSizeDistSize - 1; i++) {
    if (masterPoreSizeDist_[i].volfrac > 0.0) {
      out << masterPoreSizeDist_[i].diam << ","
          << masterPoreSizeDist_[i].volfrac << ","
          << masterPoreSizeDist_[i].volfracsat << endl;
    }
  }

  out << ">=" << masterPoreSizeDist_[masterPoreSizeDistSize - 1].diam << ","
      << masterPoreSizeDist_[masterPoreSizeDistSize - 1].volfrac << ","
      << masterPoreSizeDist_[masterPoreSizeDistSize - 1].volfracsat << endl;
  out.flush();
  out.close();

  return;
}

void Lattice::writeMicroColors() {
  string ofileName(jobRoot_);
  ofileName.append("_Colors.csv");

  ofstream out(ofileName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writeLattice", ofileName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  // int microPhaseId;
  string microPhaseName;
  // int red, green, blue;
  vector<int> colors;
  out << numMicroPhases_ << endl;
  for (int i = 0; i < numMicroPhases_; i++) {
    microPhaseName = chemSys_->getMicroPhaseName(i);
    // colors = chemSys_->getColor(microPhaseId);
    colors = chemSys_->getRGB(i);
    out << i << " " << microPhaseName << " " << colors[0] << " "
        << colors[1] << " " << colors[2] << endl;
  }
  out.flush();
  out.close();
  return;
}

void Lattice::writeLattice(const string timeString) {
  string ofileName(jobRoot_);

  ostringstream ostrT;
  ostrT << setprecision(3) << temperature_;
  string tempstr(ostrT.str());

  ofileName = ofileName + "." + timeString + "." + tempstr + "K.img";

  ofstream out(ofileName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writeLattice", ofileName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  // Write image header information first

  out << VERSIONSTRING << " " << thamesVersion_ << endl;
  out << XSIZESTRING << " " << xdim_ << endl;
  out << YSIZESTRING << " " << ydim_ << endl;
  out << ZSIZESTRING << " " << zdim_ << endl;
  out << IMGRESSTRING << " " << resolution_ << endl;

  // int index;
  // for (int k = 0; k < zdim_; k++) {
  //   for (int j = 0; j < ydim_; j++) {
  //     for (int i = 0; i < xdim_; i++) {
  //       index = getIndex(i, j, k);
  //       out << site_[index].getMicroPhaseId() << endl;
  //     }
  //   }
  // }
  for (int i = 0; i < numSites_; i++) {
    out << site_[i].getMicroPhaseId() << endl;
  }
  out.close();
}

void Lattice::writeNewLattice(int newZdim) {
  string ofileName("newInputImg");
  ostringstream ostr_newZdim;
  ostr_newZdim << setprecision(3) << newZdim;
  string newZstr(ostr_newZdim.str());
  ofileName = ofileName + "_newZdim." + newZstr + ".img";

  ofstream out(ofileName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writeLattice", ofileName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  // Write image header information first

  out << VERSIONSTRING << " " << thamesVersion_ << endl;
  out << XSIZESTRING << " " << xdim_ << endl;
  out << YSIZESTRING << " " << ydim_ << endl;
  out << ZSIZESTRING << " " << newZdim << endl;
  out << IMGRESSTRING << " " << 1 << endl;

  int numSites = xdim_ * ydim_ * newZdim;
  for (int i = 0; i < numSites; i++) {
    out << site_[i].getMicroPhaseId() << endl;
  }
  out.close();
}

void Lattice::writeLatticeXYZ(const double curtime,
                              const string timeString) {
  string ofileName(jobRoot_);

  ostringstream ostr1;
  ostr1 << setfill('0') << setw(6)
        << static_cast<int>((curtime * 60.0) + 0.5); // minutes
  string timestr(ostr1.str());

  ostringstream ostrT;
  ostrT << setprecision(3) << temperature_;
  string tempstr(ostrT.str());

  ofileName = ofileName + "_allSites." + timeString + "." + tempstr + "K.xyz";

  if (verbose_) {
    cout << "    In Lattice::writeLatticeXYZ, curtime = " << curtime
         << ", ofileName = " << ofileName << endl;
    cout.flush();
  }

  ofstream out(ofileName.c_str());

  // Write the file headers
  out << numSites_ << endl; // Number of voxels to visualize
  out << "Lattice=\"" << static_cast<float>(xdim_) << " 0.0 0.0 0.0 "
      << static_cast<float>(ydim_) << " 0.0 0.0 0.0 " << static_cast<float>(zdim_)
      << "\" ";
  out << "Time=" << timestr << " ";
  // out << "Properties=pos:R:3:color:R:3:transparency:R:1 " << endl;
  out << "Properties=phaseID:I:1:element:S:1:pos:R:3:vector_color:R:3:"
         "radius:R:1:transparency:R:1"
      << endl;

  // Loop over all voxels and write out the solid ones

  float x, y, z;
  int mPhId;
  vector<int> colors;
  string symb, mPhName;

  for (int i = 0; i < numSites_; i++) {
    x = site_[i].getX();
    y = site_[i].getY();
    z = site_[i].getZ();
    mPhId = site_[i].getMicroPhaseId();
    symb = cfgElem_[mPhId].symb; // mPhIdgetElemSymb(mPhId);
    mPhName = chemSys_->getMicroPhaseName(mPhId);
    // colors = chemSys_->getColor(mPhId);
    colors = chemSys_->getRGB(mPhId);
    // out << mPhId << "\t" << symb << "\t" << x << "\t" << y << "\t" << z << "\t"
    out << setw(2) << right << mPhId << "  " << setw(15) << left << mPhName
        << "\t" << x << "\t" << y << "\t" << z << "\t"
        << colors[0] << "\t" << colors[1] << "\t" << colors[2] << "\t"
        << particRadius_ << "\t" << "0.0" << endl;
    // out << x << "\t" << y << "\t" << z
    //     << "\t" << colors[0] << "\t" << colors[1] << "\t" << colors[2]
    //     << "\t" << "0.0" << endl;
  }
  out.close();
}

void Lattice::appendXYZ(double curtime) {

  string ofileName(jobRoot_);
  ofileName = ofileName + ".xyz";
  ofstream out;

  if (curtime < 1.0e-8) {
    // Create a new file
    out.open(ofileName.c_str(), ios::out);
  } else {
    // File should already exist so append to it
    out.open(ofileName.c_str(), ios::app);
  }

  // Write the file headers
  out << numSites_ << endl; // Number of voxels to visualize
  out << "Lattice=\"" << static_cast<float>(xdim_) << " 0.0 0.0 0.0 "
      << static_cast<float>(ydim_) << " 0.0 0.0 0.0 " << static_cast<float>(zdim_)
      << "\" ";
  out << "Properties=pos:R:3:color:R:3:transparency:R:1 ";
  out << "Time=" << curtime << endl;

  // Loop over all voxels and write out the solid ones

  float x, y, z;
  int mPhId;
  vector<float> colors;
  // string symb;
  float transparency = 0.7;

  for (int i = 0; i < numSites_; i++) {
    x = site_[i].getX();
    y = site_[i].getY();
    z = site_[i].getZ();
    mPhId = site_[i].getMicroPhaseId();
    // symb = cfgElem_[mPhId].symb; //mPhIdgetElemSymb(mPhId);
    colors = chemSys_->getRGBf(mPhId);

    transparency = (mPhId > 1) ? 0.0 : 0.7;
    out << x << "\t" << y << "\t" << z << "\t" << colors[0] << "\t" << colors[1]
        << "\t" << colors[2] << "\t" << transparency << endl;
  }
  out.close();
}

void Lattice::writeLatticeCFG(const string timeString) {

  string ofileNameCFG(jobRoot_);
  string ofileNameUSR(jobRoot_);

  ostringstream ostrT;
  ostrT << setprecision(3) << temperature_;
  string tempstr(ostrT.str());

  ofileNameCFG = ofileNameCFG + "." + timeString + "." + tempstr + "K.cfg";
  ofileNameUSR = ofileNameUSR + "." + timeString + "." + tempstr + "K.usr";

  ofstream outCFG(ofileNameCFG.c_str());
  ofstream outUSR(ofileNameUSR.c_str());

  int ord;
  double x, y, z;
  double xdim1 = xdim_ - 1;
  double ydim1 = ydim_ - 1;
  double zdim1 = zdim_ - 1;
  vector<int> colors;
  int i, mPhId;
  bool mPhNotExist;

  int numPart = numSites_;

  // write cfg header
  outCFG << "Number of particles = " << numPart << endl;

  outCFG << "A = 3.0 Angstrom (basic length-scale)" << endl;

  outCFG << "H0(1,1) = " << xdim1 << " A" << endl;
  outCFG << "H0(1,2) = 0 A" << endl;
  outCFG << "H0(1,3) = 0 A" << endl;

  outCFG << "H0(2,1) = 0 A" << endl;
  outCFG << "H0(2,2) = " << ydim1 << " A" << endl;
  outCFG << "H0(2,3) = 0 A" << endl;

  outCFG << "H0(3,1) = 0 A" << endl;
  outCFG << "H0(3,2) = 0 A" << endl;
  outCFG << "H0(3,3) = " << zdim1 << " A" << endl;
  outCFG << ".NO_VELOCITY." << endl;
  outCFG << "entry_count = 3" << endl;

  // write cfg coordinates & usr file (colors)
  ord = 0;
  for (mPhId = 0; mPhId < numMicroPhases_; mPhId++) {
    // for (mPhId = FIRST_SOLID; mPhId < numMicroPhases_; mPhId++) {
    mPhNotExist = true;
    // colors = chemSys_->getColor(mPhId);
    colors = chemSys_->getRGB(mPhId);
    for (i = 0; i < numSites_; i++) {
      if (site_[i].getMicroPhaseId() == mPhId) {
        if (mPhNotExist) {
          outCFG << "  " << cfgElem_[mPhId].mass << endl;
          outCFG << cfgElem_[mPhId].symb << endl;
          mPhNotExist = false;
        }
        x = (site_[i].getX()) / xdim1;
        y = (site_[i].getY()) / ydim1;
        z = (site_[i].getZ()) / zdim1;
        outCFG << x << "\t" << y << "\t" << z << endl;
        outUSR << ord << "\t" << colors[0] << "\t" << colors[1] << "\t"
               << colors[2] << "\t0.7051" << endl;
        ord++;
      }
    }
  }
  outCFG.close();
  outUSR.close();

  if (ord != numPart) {
    cout << endl
         << "*************** error in writeLatticeCFG! ******************"
         << endl;
    cout << endl << "                    ord != numSites_" << endl;
    cout << "                                   STOP" << endl;
    exit(1);
  }
}

void Lattice::writeDamageLattice(const string timeString) {
  string ofileName(jobRoot_);
  ostringstream ostrT;
  ostrT << setprecision(3) << temperature_;
  string tempstr(ostrT.str());

  ofileName = ofileName + "." + timeString + "." + tempstr + "K.damage.img";

  cout << endl
       << "  Lattice::writeDamageLattice - ofileName = " << ofileName << endl;

  ofstream out(ofileName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writeLattice", ofileName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  // Write image header information first

  out << VERSIONSTRING << " " << thamesVersion_ << endl;
  out << XSIZESTRING << " " << xdim_ << endl;
  out << YSIZESTRING << " " << ydim_ << endl;
  out << ZSIZESTRING << " " << zdim_ << endl;
  out << IMGRESSTRING << " " << resolution_ << endl;

  for (int i = 0; i < numSites_; i++) {
    if (site_[i].IsDamage()) {
      out << "1" << endl;
    } else {
      out << "0" << endl;
    }
  }

  out.close();
}

void Lattice::writeLatticePNG(const string timeString) {
  int i, j;
  string oppmName(jobRoot_);
  string opngName(jobRoot_);

  vector<double> dumvec;
  vector<int> idumvec;
  vector<vector<int>> image;
  vector<vector<double>> dshade;
  dumvec.resize(ydim_, 0.0);
  idumvec.resize(ydim_, 0);
  dshade.resize(xdim_, dumvec);
  image.resize(xdim_, idumvec);
  bool done;
  int resCallSystem;

  ///
  /// Construct the name of the output file
  ///

  ostringstream ostrT;
  ostrT << setprecision(3) << temperature_;
  string tempstr(ostrT.str());

  oppmName = oppmName + "." + timeString + "." + tempstr + "K.ppm";
  opngName = opngName + "." + timeString + "." + tempstr + "K.png";

  ///
  /// Open the output file
  ///

  ofstream out(oppmName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writeLatticePNG", oppmName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  ///
  /// Write PPM header for full color image
  ///

  out << "P3" << endl;
  out << xdim_ << " " << ydim_ << endl;
  out << COLORSATVAL << endl;

  int slice = zdim_ / 2;
  int nd, izz;
  int sitenum;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      if (depthEffect_) {
        done = false;
        nd = 0;
        izz = slice;
        do {
          sitenum = getIndex(i, j, izz);
          if (nd == 10 || site_[sitenum].getMicroPhaseId() > 1) {
            done = true;
          } else {
            nd++;
            izz++;
            if (izz >= zdim_)
              izz -= zdim_;
          }
        } while (!done);
        sitenum = getIndex(i, j, izz);
        image[i][j] = site_[sitenum].getMicroPhaseId();
        dshade[i][j] = 0.1 * (10.0 - (static_cast<double>(nd)));
      } else {
        sitenum = getIndex(i, j, slice);
        image[i][j] = site_[sitenum].getMicroPhaseId();
        dshade[i][j] = 1.0;
      }
    }
  }

  int red, green, blue;
  vector<int> colors;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      colors = chemSys_->getRGB(image[i][j]);
      red = dshade[i][j] * colors[0] + 0.5;
      green = dshade[i][j] * colors[1] + 0.5;
      blue = dshade[i][j] * colors[2] + 0.5;
      out << red << " " << green << " " << blue << endl;
    }
  }

  out.close();

  ///
  /// Execute system call to convert PPM to PNG.
  ///
  /// @warning This relies on installation of ImageMagick
  ///

  string buff;
  buff = ConvertCommand + " " + oppmName + " " + opngName;
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    // handle the error;
    cout << endl
         << endl
         << "    Lattice.cc - error in writeLatticePNG() : resCallSystem = -1"
         << endl;
    cout << endl << "    STOP program" << endl;
    // throw HandleException ("writeLatticePNG", "Lattice.cc",
    //                 "system(buff.c_str())", "err : resCallSystem = -1");
    exit(1);
  }
  return;
}

void Lattice::writeDamageLatticePNG(const string timeString) {
  int i, j;
  string oppmName(jobRoot_);
  string opngName(jobRoot_);

  vector<double> dumvec;
  vector<unsigned int> idumvec;
  vector<vector<unsigned int>> image;
  vector<vector<double>> dshade;
  dumvec.resize(ydim_, 0.0);
  idumvec.resize(ydim_, 0);
  dshade.resize(xdim_, dumvec);
  image.resize(xdim_, idumvec);
  bool done;
  int resCallSystem;

  ///
  /// Construct the name of the output file
  ///

  ostringstream ostrT;
  ostrT << setprecision(3) << temperature_;
  string tempstr(ostrT.str());
  string buff;

  oppmName = oppmName + "." + timeString + "." + tempstr + "K.damage.ppm";
  opngName = opngName + "." + timeString + "." + tempstr + "K.damage.png";

  ///
  /// Open the output file
  ///

  ofstream out(oppmName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writeLatticePNG", oppmName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  ///
  /// Write PPM header for full color image
  ///

  out << "P3" << endl;
  out << ydim_ << " " << xdim_ << endl;
  out << COLORSATVAL << endl;

  int slice = zdim_ / 2;
  int nd, izz;
  int sitenum;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      done = false;
      nd = 0;
      izz = slice;
      do {
        sitenum = getIndex(i, j, izz);
        if (nd == 10 || site_[sitenum].getMicroPhaseId() > 1) {
          done = true;
        } else {
          nd++;
          izz++;
          if (izz >= zdim_)
            izz -= zdim_;
        }
      } while (!done);
      sitenum = getIndex(i, j, izz);
      if (site_[sitenum].IsDamage()) {
        image[i][j] = 1;
      } else {
        image[i][j] = 5;
      }
      dshade[i][j] = 1.0;
    }
  }

  int red, green, blue;
  vector<int> colors;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      colors = chemSys_->getRGB(image[i][j]);
      red = dshade[i][j] * colors[0] + 0.5;
      green = dshade[i][j] * colors[1] + 0.5;
      blue = dshade[i][j] * colors[2] + 0.5;
      out << red << " " << green << " " << blue << endl;
    }
  }

  out.close();

  ///
  /// Execute system call to convert PPM to PNG.
  ///
  /// @warning This relies on installation of ImageMagick
  ///

  buff = ConvertCommand + " " + oppmName + " " + opngName;
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    // handle the error;
    cout << endl
         << endl
         << "    Lattice.cc - error in writeDamageLatticePNG() : resCallSystem "
            "= -1"
         << endl;
    cout << endl << "    STOP program" << endl;
    // throw HandleException ("writeDamageLatticePNG", "Lattice.cc",
    //                 "system(buff.c_str())", "err : resCallSystem = -1");
    exit(1);
  }
  return;
}

void Lattice::makeMovie() {
  int i, j, k;
  string ofileName(jobRoot_);
  string ofbasename(jobRoot_);
  string ofgifileName(jobRoot_);
  string ofgifbasename(jobRoot_);

  vector<double> dumvec;
  vector<int> idumvec;
  vector<vector<int>> image;
  vector<vector<double>> dshade;
  dumvec.resize(ydim_, 0.0);
  idumvec.resize(ydim_, 0);
  dshade.resize(zdim_, dumvec);
  image.resize(zdim_, idumvec);
  bool done;
  int resCallSystem;

  int slice;
  int nd, izz;
  int sitenum;

  ///
  /// Construct the name of the output file
  ///

  string buff;
  ostringstream ostr1, ostr2, ostr3;
  ostr1 << static_cast<int>(time_ * 10.0); // tenths of an hour
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());

  ///
  /// Loop over number of slices in the x direction, making one image per
  /// slice and appending it to the end of the master file.
  ///

  for (k = 10; k < zdim_; k++) {

    ///
    /// Open the output file.
    ///

    ostr3.clear();
    ostr3 << static_cast<int>(k); // x slice number
    string kstr(ostr3.str());
    ofileName =
        ofbasename + "." + timestr + "." + tempstr + "K." + kstr + ".ppm";
    ofgifileName =
        ofgifbasename + "." + timestr + "." + tempstr + "K." + kstr + ".gif";

    ofstream out(ofileName.c_str());
    if (!out.is_open()) {
      throw FileException("Lattice", "makeMovie", ofileName, "Could not open");
    }

    ///
    /// Write PPM header for full color image.
    ///

    out << "P3" << endl;
    out << xdim_ << " " << ydim_ << endl;
    out << COLORSATVAL << endl;

    slice = k;
    for (j = 0; j < ydim_; j++) {
      for (i = 0; i < xdim_; i++) {
        done = false;
        nd = 0;
        // sitenum = getIndex(i, j, slice);
        izz = slice;
        do {
          sitenum = getIndex(i, j, izz);
          if (nd == 10 || site_[sitenum].getMicroPhaseId() > 1) {
            done = true;
          } else {
            nd++;
            izz++;
            if (izz >= zdim_)
              izz -= zdim_;
          }
        } while (!done);
        sitenum = getIndex(i, j, izz);
        image[i][j] = site_[sitenum].getMicroPhaseId();
        dshade[i][j] = 0.1 * (10.0 - (static_cast<double>(nd)));
      }
    }

    int red, green, blue;
    vector<int> colors;
    for (j = 0; j < ydim_; j++) {
      for (i = 0; i < xdim_; i++) {
        // colors = chemSys_->getColor(image[i][j]);
        colors = chemSys_->getRGB(image[i][j]);
        red = dshade[i][j] * colors[0] + 0.5;
        green = dshade[i][j] * colors[1] + 0.5;
        blue = dshade[i][j] * colors[2] + 0.5;
        out << red << " " << green << " " << blue << endl;
      }
    }
    out.close();

    ///
    /// Execute system call to convert PPM to GIF.
    ///
    /// @warning This relies on installation of ImageMagick
    ///

    buff = ConvertCommand + " " + ofileName + " " + ofgifileName;
    resCallSystem = system(buff.c_str());
    if (resCallSystem == -1) {
      // handle the error;
      cout << endl
           << endl
           << "    Lattice.cc - error(1) in makeMovie() : resCallSystem = -1"
           << endl;
      cout << endl << "    STOP program" << endl;
      // throw HandleException ("makeMovie", "Lattice.cc",
      //                 "system(buff.c_str())", " err_1 : resCallSystem = -1");
      exit(1);
    }
  }

  ///
  /// Execute system call to convert GIF frames to animated GIF.
  ///
  /// @warning This relies on installation of gifsicle
  ///

  buff = "gifsicle --delay=10 " + ofgifbasename + "*.gif > " + ofgifbasename +
         ".movie.gif";
  resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    // handle the error;
    cout << endl
         << endl
         << "    Lattice.cc - error(2) in makeMovie() : resCallSystem = -1"
         << endl;
    cout << endl << "    STOP program" << endl;
    // throw HandleException ("makeMovie", "Lattice.cc",
    //                 "system(buff.c_str())", "err_2 : resCallSystem = -1");
    exit(1);
  }
}

double Lattice::fillAllPorosity(const int cyc) {

  double waterVoidMolesM = 0;
  double waterVoidMolesV = 0;
  double waterVoidMass = 0;

  double waterDensity = waterMolarMass_ / waterMolarVol_ / 1.0e6; // g/cm3

  cout << "       Lattice::fillAllPorosity ini - cyc = " << cyc
       << "  :  count_[VOIDID] = " << count_[VOIDID]
       << "  &  count_[ELECTROLYTEID] = " << count_[ELECTROLYTEID] << endl;

  if (count_[VOIDID] > 0) {
    // double massAddWater = 0;
    int countVoid = 0;
    int countElectrolyte = 0;
    double countWmc0 = 0;
    double volFracVoid = 0;
    double volFracPorosity = 0;
    double volFracWmc0 = 0;
    double volFracEle = 0;
    vector<int> voidVect;

    int pId;
    for (int i = 0; i < numSites_; i++) {
      pId = site_[i].getMicroPhaseId();
      if (pId == VOIDID) {
        voidVect.push_back(i);
        countVoid++;
      } else if (pId == ELECTROLYTEID) {
        countElectrolyte++;
      } else {
        countWmc0 += site_[i].getWmc0();
      }
    }

    if (count_[VOIDID] == countVoid) {
      int siteID;
      int nbpid;
      Site *stenb;
      for (int ii = 0; ii < countVoid; ii++) {
        siteID = voidVect[ii];
        setMicroPhaseId(siteID, ELECTROLYTEID);
        site_[siteID].setWmc0(electrolyteIntPorosity_);
        site_[siteID].dWmc(electrolyteIntPorosity_);
        site_[siteID].clearGrowth();
        for (int i = 0; i < NN_NNN; i++) {
          stenb = site_[siteID].nb(i);
          // wmcIni = stenb->getWmc();
          stenb->dWmc(electrolyteIntPorosity_);
          nbpid = stenb->getMicroPhaseId();
          if (nbpid > ELECTROLYTEID) {
            if (stenb->getInDissInterfacePos() == -1)
              addDissolutionSite(stenb, nbpid);
            if (site_[siteID].getInGrowInterfacePos(nbpid) == -1)
              addGrowthSite(&site_[siteID], nbpid);
            if (i < NUM_NEAREST_NEIGHBORS) {
              for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
                   phaseTmpl++) {
                // if (chemSys_->isGrowthTemplate(phaseTmpl, nbpid)) {
                if (growthTemplate_[phaseTmpl][nbpid]) {
                  if (site_[siteID].getInGrowInterfacePos(phaseTmpl) == -1) {
                    addGrowthSite(&site_[siteID], phaseTmpl);
                  }
                }
              }
            }
          }
        }
      }
    } else {
      cout << endl
           << "Lattice::fillAllPorosity(...) error cyc = " << cyc
           << "  :  count_[VOIDID] = " << count_[VOIDID]
           << "  while  countVoid = " << countVoid << endl;
      cout << endl << "stop program" << endl;
      exit(0);
    }

    volFracWmc0 = countWmc0 / numSites_;
    volFracEle = static_cast<double>(countElectrolyte) / numSites_;
    volFracVoid = static_cast<double>(countVoid) / numSites_;
    volFracPorosity = (countWmc0 + countVoid + countElectrolyte) / numSites_;

    cout << "       Lattice::fillAllPorosity fin - cyc = " << cyc
         << " : volFracWmc0 = " << volFracWmc0
         << "   volFracELECTR = " << volFracEle
         << "   volFracVoid = " << volFracVoid
         << "   volFracPorosity = " << volFracPorosity << endl;

    // add water for subvoxel volume that is empty???

    waterVoidMass = volFracVoid * waterDensity * 100 / initSolidMass_;

    waterVoidMolesM = waterVoidMass / waterMolarMass_; // mol
    waterVoidMolesV =
        volFracVoid * initialMicrostructureVolume_ / waterMolarVol_; // mol
  }

  cout << "       Lattice::fillAllPorosity fin - cyc = " << cyc
       << " : waterVoidMass = " << waterVoidMass
       << "   waterVoidMolesM = " << waterVoidMolesM
       << "   waterVoidMolesV = " << waterVoidMolesV << endl;

  // return waterAddMoles;
  return waterVoidMolesM;
}

vector<int> Lattice::transformPhase(int growPhId, int netsites_growPhId,
                                    vector<int> dissPhaseIDVect,
                                    vector<int> numSiteDissVect,
                                    vector<string> dissPhNameVect,
                                    vector<double> volumeRatio, int &numadded_D,
                                    int totalTRC) {

  // bool damagedTrue = false;

  //*** controll
  int bcl = 0;
  int static trc_t;
  trc_t++;

  int i, jj;

  Site *ste;
  int pid;
  int isitePos, phaseID;
  int stId;
  int posVect;
  double rng, probRNG;

  int dissPhaseIDVectSize = dissPhaseIDVect.size();
  vector<int> numChange(dissPhaseIDVectSize, 0);
  vector<int> dim_isite(dissPhaseIDVectSize, 0);
  vector<int> numLeft = numSiteDissVect; // numtotake
  vector<Isite> isite;

  cout << endl << "    Lattice::transformPhase totalTRC = " << totalTRC << " =>" << endl;
  cout << "      phaseId/phaseName/sitesToBeDissolved :" << endl;
  int numChangeTot = 0;
  int numLeftTot = 0;
  for (i = 0; i < dissPhaseIDVectSize; i++) {
    cout << "        "
         << setw(3) << dissPhaseIDVect[i] << " / "
         << setw(15) << left << dissPhNameVect[i] << " / "
         << setw(8) << right << numSiteDissVect[i] << endl;
    numLeftTot += numSiteDissVect[i];
  }
  cout << "      total sitesToBeDissolved  : " << setw(8) << right
       << numLeftTot << endl;
  cout << "      total sitesToBeGrown (AFt): " << setw(8) << right
       << netsites_growPhId << endl;

  if (numLeftTot == 0) {
    cout << "Lattice::transformPhase error numLeftTot = 0" << endl;
    cout << "   totalTRC/trc_t/bcl :  "
         << "   " << totalTRC << "/" << trc_t << "/" << bcl << endl;
    cout << "stop program" << endl;
    exit(0);
  }

  // transform probabilities : equal probability
  vector<structDissVect> dissolutionVector;
  structDissVect dissStruct;
  int posProbVect = 0;
  // double sumWmc = 0; // numLeftTot
  long int sumWmcInt = 0;
  for (i = 0; i < dissPhaseIDVectSize; i++) {
    phaseID = dissPhaseIDVect[i];
    isite = interface_[phaseID].getDissolutionSites();
    dim_isite[i] = isite.size();
    for (jj = 0; jj < dim_isite[i]; jj++) {

      stId = isite[jj].getId();
      // stWmc =  site_[stId].getWmc();
      // sumWmc += stWmc;
      dissStruct.id = stId;
      dissStruct.posVect = i;
      // dissStruct.wmc = stWmc;
      dissStruct.wmcInt = 1; // equal probability!

      dissolutionVector.push_back(dissStruct);
      // dissProbStruct(int id_i = 0, double instab_i = 0, int posVect_i = 0,
      // double prob_i = 0, double prob_01_i = 0) {}

      site_[stId].setInDissolutionVectorPos(posProbVect);
      posProbVect++;
    }
  }
  int dissolutionVectorSize = dissolutionVector.size();
  sumWmcInt = dissolutionVectorSize;

  cout << "    Lattice::transformPhase DISS_INI totalTRC/trc_t/bcl/sumWmcInt : "
       << totalTRC << " / " << trc_t << " / " << bcl << " / " << sumWmcInt << endl;
  cout << "      DISS_INI dissPhaseIDVectSize = " << dissPhaseIDVectSize
       << "   dissolutionVectorSize = " << dissolutionVectorSize
       << "   numLeftTot = " << numLeftTot
       << "   numChangeTot = " << numChangeTot << endl;
  for (i = 0; i < dissPhaseIDVectSize; i++) {
    phaseID = dissPhaseIDVect[i];
    // isite = interface_[phaseID].getDissolutionSites();
    // dim_isite = isite.size();
    cout << "        DISS_INI for i = " << setw(3) << i
         << "  => phaseID phaseName count_ dim_isite numleft numchange  :  "
         << setw(3) << phaseID << "   " << setw(15) << left << dissPhNameVect[i]
         << "   " << setw(8) << right << count_[phaseID] << "   " << setw(8)
         << dim_isite[i] << "   " << setw(8) << numLeft[i] << "   " << setw(8)
         << numChange[i] << endl;
  }
  cout << "        WAIT..." << endl;
  cout.flush();

  // int isitePosError = 0;

  int callGEM = -1;

  vector<double> expval;
  expval.clear();
  expval.resize(3, 0.0);

  vector<int> coordin;
  coordin.clear();
  coordin.resize(3, 0);

  double alreadygrown = 0.0;
  int max = -10000;

  vector<Site *> porousneighbor, waterneighbor;
  int waterneighborSize;
  Site *stenbtp;

  // string fileName(jobRoot_ + "_alsubvol.dat");
  vector<int> alnbSiteId;
  vector<int> alnbPhId;
  // int stenb_mPhId;
  // int size;
  // Site *alstenb;

  int numWater, numPorous;
  // double totalPorosity;
  int totalPorosity;

  double volumeratio;
  // int numGrowth = 0;

  int smallerThanMax = 0;
  int greaterThanMax = 0;

  int countExpPos = 0;
  int countExpNeg = 0;

  double sumWmcDbl;

  while (numLeftTot > 0 &&              // monosulphate
         (dissolutionVectorSize > 0) && // monosulphate dissolution interface
         (alreadygrown < netsites_growPhId)) { // AFt

    try {
      bcl++;

      // equal probability!
      rng = callRNG();
      sumWmcDbl = static_cast<double>(sumWmcInt);

      // sumWmc = 0.;
      // for (int i = 0; i < dissolutionVectorSize; i++) {
      //   sumWmc +=  dissolutionVector[i].instab;
      // }
      probRNG = dissolutionVector[0].wmcInt / sumWmcDbl;
      // probRNG_0 = 1./ sumWmc;
      // probRNG = probRNG_0;
      if (rng <= probRNG) {
        isitePos = 0;
      } else {
        for (isitePos = 1; isitePos < dissolutionVectorSize; isitePos++) {
          probRNG += (dissolutionVector[isitePos].wmcInt / sumWmcDbl);
          if (rng <= probRNG)
            break;
        }
      }

      if (isitePos >= dissolutionVectorSize) {
        // isitePosError++;
        cout << endl
             << "transformPhase:     *** isitePosError: "
                "bcl/isitePos/rng/probRNG = "
             << bcl << " / " << isitePos << " / " << rng << " / " << probRNG
             << endl;
        cout << endl << "    exit" << endl;
        exit(1);
      } else {
        long int sumWmcT = 0;
        for (int i = 0; i < dissolutionVectorSize; i++) {
          sumWmcT += dissolutionVector[isitePos].wmcInt;
        }
        if (abs(sumWmcInt - sumWmcT) > 0) {
          cout << endl
               << "transformPhase:     *** sumError: "
                  "bcl/sumWmc/sumWmcT/(sumWmcInt - sumWmcT) = "
               << bcl << " / " << sumWmcInt << " / " << sumWmcT << " / "
               << (sumWmcInt - sumWmcT) << endl;
          cout << endl << "    exit" << endl;
          exit(1);
        }
      }

      // pos to convert

      ste = &site_[dissolutionVector[isitePos].id];
      pid = ste->getMicroPhaseId(); // intrebare pid diff phaseid ???
      if (pid != 14) { // for the moment - only Monosulfate
        cout << endl
             << ">>> for isitePos = " << isitePos
             << "  =>  pid = " << " (!= 14!!!) for bcl = " << bcl << endl;
        cout << " >stop<" << endl;
        exit(0);
      }

      posVect =
          dissolutionVector[isitePos].posVect; // pos pid in dissPhaseIDVect
      volumeratio = volumeRatio[posVect];
      max = static_cast<int>(volumeratio);

      sumWmcInt -= 1; // ste->getWmc();

      porousneighbor.clear();
      waterneighbor.clear();
      for (int j = 0; j < NUM_NEAREST_NEIGHBORS; j++) {
        stenbtp = ste->nb(j);
        if (stenbtp->getMicroPhaseId() == ELECTROLYTEID) {
          waterneighbor.push_back(stenbtp);
        } else if (chemSys_->isPorous(stenbtp->getMicroPhaseId())) {
          porousneighbor.push_back(stenbtp);
        }
      }
      waterneighborSize = waterneighbor.size();

      if (ste->getInDissInterfacePos() == -1) {
        cout << endl
             << "    Lattice::transformPhase error: "
                "ste->getInDissInterfacePos() = -1"
             << endl;
        cout << "    Lattice::transformPhase error: steId/pid/posVect/isitePos  "
             << ste->getId() << "/" << pid << "/" << posVect << "/" << isitePos << endl;
        cout << "    Lattice::transformPhase error: totalTRC/trc_t/bcl  "
             << totalTRC << "/" << trc_t << "/" << bcl << endl;
        cout << "    Lattice::transformPhase error: exit" << endl;
        exit(0);
      }

      if ((waterneighborSize + 1) <= max) { // count the site itself

        smallerThanMax++;

        /// Expansion should occur
        ///
        /// 1. Take the subvolume centered on aluminate site that will dissolve
        ///

        alnbSiteId.clear();
        alnbPhId.clear();
        numWater = numPorous = 0;
        totalPorosity = 0;

        // string fileName(jobRoot_ + "_alsubvol.dat"); // check!

        // cout << "Lattice::transformPhase fileName = " << fileName << endl;

        // alnbSiteId =
        //   writeSubVolume(fileName, ste, 1); // all 26 neighbors + ste itself

        // cout << "Lattice::transformPhase alnbSiteId.size() = " <<
        // alnbSiteId.size() << endl;

        // alnbSiteId = getNeighborhood(ste->getId(), 1);

        getNeighborhood(ste->getId(), alnbSiteId, alnbPhId, numWater, totalPorosity);

        // ste->setDamage();
        // size = alnbSiteId.size();
        // for (int nb = 0; nb < size; nb++) {
          // alstenb = &site_[alnbSiteId[nb]];
          // stenb_mPhId = alstenb->getMicroPhaseId();
          // alnbPhId.push_back(stenb_mPhId);
          // if (stenb_mPhId == ELECTROLYTEID) {
          // if (alnbPhId[nb] == ELECTROLYTEID) {
          //   numWater++;
          // }
          // else if (chemSys_->isPorous(stenb_mPhId)) {
          //   // numPorous++;
          //   alstenb->setDamage();
          // } else if (chemSys_->isWeak(stenb_mPhId)) {
          //   alstenb->setDamage();
          // }
          // totalPorosity += alstenb->getWmc0(); // including ste
          // totalPorosity += site_[alnbSiteId[nb]].getWmc0();
        // }

        ///
        /// 2. Calculate the effective bulk modulus of this subvolume
        ///

        // cout << "Lattice::transformPhase bf-FEsolver_" << endl;

        // double subbulk = FEsolver_->getBulkModulus(fileName);

        // cout << "Lattice::transformPhase af-FEsolver_" << endl;

        // double subbulk = FEsolver_->getBulkModulus(alnbPhId);

        vector<int> *p_alnbPhId = &alnbPhId;
        double subbulk = FEsolver_->getBulkModulus(p_alnbPhId);

        subbulk = subbulk * 1.0e3; // convert GPa to MPa

        double subsolidbulk =
            subbulk * ((1.0 + numWater / 27.0) / (1.0 - numWater / 27.0));

        ///
        /// 3. Calculate crystallization strain in this sub volume;
        ///    porevolfrac is the volume fraction of pore space occupied by
        ///    crystal
        ///
        ///    @todo generalize porous phase porosities in the block below
        ///    instead of 0.25
        ///

        // double porevolfrac = 0.0;
        // if (numWater != 0 || numPorous != 0) {
        //   // porevolfrac = static_cast<double>(volumeratio) / (numWater + (numPorous *
        //   0.25)); porevolfrac = volumeratio / totalPorosity;
        // } else {
        //   porevolfrac = 1.0;
        // }
        double porevolfrac = 1.;
        if (totalPorosity > 0)
          // porevolfrac = volumeratio / (totalPorosity / static_cast<double>(convFactDbl2IntPor_));
          porevolfrac = (volumeratio / totalPorosity) * convFactDbl2IntPor_;

        //  This is hard-wired right now
        //  @todo generalize crystallization pressure to more phases

        // double exp = solut_->calculateCrystalStrain(SI_[growingid],
        // porevolfrac,
        //                                             subbulk, subsolidbulk);
        double exp = chemSys_->calculateCrystalStrain(growPhId, porevolfrac,
                                                      subbulk, subsolidbulk);

        ///
        /// 4. Apply expansion strain on each voxel in this sub volume
        ///

        if (exp > 0) {
          countExpPos++;
          applyExpansion(alnbSiteId, exp); //*
        } else {
          countExpNeg++;
        }

        //************

        transformSolSol(ste, pid, growPhId, totalTRC); // sol to sol (nonporous)

        // update dissolutionVector & involved sites
        // site_[dissolutionVector[isitePos].id].setInDissolutionVectorPos(-1);
        if (isitePos != dissolutionVectorSize - 1) {
          dissolutionVector[isitePos] =
              dissolutionVector[dissolutionVectorSize - 1];
          site_[dissolutionVector[isitePos].id].setInDissolutionVectorPos(
              isitePos);
        }
        dissolutionVector.pop_back();
        dissolutionVectorSize--;

        vector<int> extractSteIdVect;
        for (int i = 0; i < waterneighborSize; i++) {
          extractSteIdVect.clear();
          extractSteIdVect = transformLiqSol(waterneighbor[i], growPhId, totalTRC); // liq to sol
          if (extractSteIdVect.size() > 0) {
            int pid;
            string pname;
            int size = extractSteIdVect.size();
            // cout << endl << "        >>>> forLoop -> bcl/extractSteIdVect.size   : " // for tests
            //      << setw(8) << right << bcl << "  "
            //      << setw(3) << right << extractSteIdVect.size() << endl;
            for (int i = 0; i < size; i++) {
              pid = site_[extractSteIdVect[i]].getMicroPhaseId();
              pname = chemSys_->getMicroPhaseName(pid);
              // cout << "          i/extractSteIdVect[i]/microPhaseId/microPhaseName : "
              //      << setw(3) << right << i << "  "
              //      << setw(10) << right << extractSteIdVect[i] << "  "
              //      << setw(3) << right << pid << "  " << setw(15) << left << pname
              //      << endl;
              // if (pid == 14) { // dissPhaseIDVect[0] now
              if (pid == dissPhaseIDVect[0]) { // dissPhaseIDVect[0] = 14 now!!!
                // cout << "            * - extract siteId " << extractSteIdVect[i]
                //      << " from dissolutionVector !!!" << endl;
                int posDissVect = site_[extractSteIdVect[i]].getInDissolutionVectorPos();
                site_[extractSteIdVect[i]].setInDissolutionVectorPos(-1);
                if (posDissVect != dissolutionVectorSize - 1) {
                  dissolutionVector[posDissVect] =
                      dissolutionVector[dissolutionVectorSize - 1];
                  site_[dissolutionVector[posDissVect].id].setInDissolutionVectorPos(
                        posDissVect);
                }
                dissolutionVector.pop_back();
                dissolutionVectorSize--;
                sumWmcInt -= 1;
              } // else { // for tests
                // cout << "            * - nothing to do for siteId " << extractSteIdVect[i]
                //      << endl;
              // }
            }
          }

          alreadygrown++;
        }

      } else { // if ((waterneighborSize + 1) <= max) { // count the site itself

        greaterThanMax++;

        transformSolSol(ste, pid, growPhId, totalTRC); // sol to sol

        // update dissolutionVector & involved sites
        // site_[dissolutionVector[isitePos].id].setInDissolutionVectorPos(-1);
        if (isitePos != dissolutionVectorSize - 1) {
          dissolutionVector[isitePos] =
              dissolutionVector[dissolutionVectorSize - 1];
          site_[dissolutionVector[isitePos].id].setInDissolutionVectorPos(
              isitePos);
        }
        dissolutionVector.pop_back();
        dissolutionVectorSize--;

        // same_0 ... same_1

        double thresh = 0.0, g = 0.0;
        thresh = volumeratio - max;
        g = callRNG();
        int upperindex = (g < thresh) ? max : max - 1;
        int contor = 0;
        int ij;
        vector<int> extractSteIdVect;

        while (contor < upperindex) {
          g = callRNG();
          for (ij = 0; ij < waterneighborSize; ij++) {
            if (g < (ij + 1.0) / waterneighborSize) {
              break;
            }
          }

          // change phaseId for site ij
          extractSteIdVect.clear();
          extractSteIdVect = transformLiqSol(waterneighbor[ij], growPhId, totalTRC); // liq to sol
          if (extractSteIdVect.size() > 0) {
            int pid;
            string pname;
            int size = extractSteIdVect.size();
            // cout << endl << "        >>>> whileLoop -> bcl/extractSteIdVect.size : " // for tests
            //      << setw(8) << right << bcl << "  "
            //      << setw(3) << right << extractSteIdVect.size() << endl;
            for (int i = 0; i < size; i++) {
              pid = site_[extractSteIdVect[i]].getMicroPhaseId();
              pname = chemSys_->getMicroPhaseName(pid);
              // cout << "          i/extractSteIdVect[i]/microPhaseId/microPhaseName : "
              //      << setw(3) << right << i << "  "
              //      << setw(10) << right << extractSteIdVect[i] << "  "
              //      << setw(3) << right << pid << "  " << setw(15) << left << pname
              //      << endl;
              // if (pid == 14) { // dissPhaseIDVect[0] now
              if (pid == dissPhaseIDVect[0]) { // dissPhaseIDVect[0] = 14 now!!!
                // cout << "            * - extract siteId " << extractSteIdVect[i]
                //      << " from dissolutionVector !!!" << endl;
                int posDissVect = site_[extractSteIdVect[i]].getInDissolutionVectorPos();
                site_[extractSteIdVect[i]].setInDissolutionVectorPos(-1);
                if (posDissVect != dissolutionVectorSize - 1) {
                  dissolutionVector[posDissVect] =
                      dissolutionVector[dissolutionVectorSize - 1];
                  site_[dissolutionVector[posDissVect].id].setInDissolutionVectorPos(
                        posDissVect);
                }
                dissolutionVector.pop_back();
                dissolutionVectorSize--;
                sumWmcInt -= 1;
              } // else { // for tests
                // cout << "            * - nothing to do for siteId " << extractSteIdVect[i]
                //      << endl;
              // }
            }
          }

          // extract site ij from waterneighbor
          waterneighbor[ij] = waterneighbor[waterneighborSize - 1];
          waterneighborSize--;

          contor++;
          alreadygrown++;
        }
      }

      numLeftTot--;
      numChangeTot++;

      numLeft[posVect]--;
      numChange[posVect]++;

      alreadygrown++;

    } catch (out_of_range &oor) {
      EOBException ex("Lattice", "transformPhase", "site_", site_.size(), i);
      ex.printException();

      cout << endl << "Lattice::transformPhase error" << endl;
      cout << endl
           << "totalTRC trc_t bcl numLeftTot numChangeTot  :  " << totalTRC
           << "   " << trc_t << "   " << bcl << "   " << numLeftTot << "   "
           << numChangeTot << endl;
      cout << endl
           << "steId pid dissolutionVectorSize :  " << ste->getId() << "   "
           << pid << "   " << dissolutionVectorSize << endl;
      cout.flush();
      for (i = 0; i < dissPhaseIDVectSize; i++) {
        phaseID = dissPhaseIDVect[i];
        isite = interface_[phaseID].getDissolutionSites();
        dim_isite[i] = isite.size();
        cout << "        phaseid count_ dim_isite numleft numchange  :  "
             << phaseID << "   " << "   " << count_[phaseID] << "   "
             << dim_isite[i] << "   " << numLeft[i] << "   " << numChange[i]
             << endl;
      }
      cout << "stop program" << endl;
      exit(1);
    }
    // cout << "     *** out totalTRC trc_d bcl :  " << totalTRC
    //      << "   " << trc_d << "   " << bcl << endl;
  }

  if (dissolutionVectorSize != static_cast<int>(dissolutionVector.size())) {
    cout << endl
         << "Lattice::transformPhase error => dissolutionVectorSize != "
            "dissolutionVector.size() : "
         << dissolutionVectorSize << " != " << dissolutionVector.size() << endl;
    cout << endl
         << "totalTRC trc_t bcl numLeftTot numChangeTot  :  " << totalTRC
         << "   " << trc_t << "   " << bcl << "   " << numLeftTot << "   "
         << numChangeTot << endl;
    cout << "stop program" << endl;
    exit(1);
  } else {
    cout << endl
         << "    Lattice::transformPhase => dissolutionVectorSize = "
            "dissolutionVector.size() = "
         << dissolutionVectorSize << endl;

    for (int i = 0; i < dissolutionVectorSize; i++) {
      site_[dissolutionVector[i].id].setInDissolutionVectorPos(-1);
    }
  }

  cout << endl
       << "    Lattice::transformPhase => totalTRC trc_t bcl numLeftTot numChangeTot  :  "
       << totalTRC << "   " << trc_t << "   " << bcl << "   " << numLeftTot << "   "
       << numChangeTot << endl;
  cout << "    Lattice::transformPhase : max = " << max
       << " & ((waterneighborSize + 1) <= max)? => smallerThanMax & greaterThanMax = "
       << smallerThanMax << "  &  " << greaterThanMax
       << endl;
  cout << "    Lattice::transformPhase : for smallerThanMax"
          " exp = chemSys_->calculateCrystalStrain  => countExpPos = "
       << countExpPos << " (applyExpansion!)  &  countExpNeg = "
       << countExpNeg << endl;

  numadded_D = numChangeTot;

  vector<int> numLeftFin = numLeft;
  numLeftFin.push_back(alreadygrown);
  numLeftFin.push_back(callGEM);

  return (numLeftFin);
} // transformPhase() - end

void Lattice::transformSolSol(Site *ste, int oldPhId, int newPhId,
                              int totalTRC) {

  //*** for control
  int static trc_cT;
  trc_cT++;

  int nbid, nbpid, nb_pid;
  vector<int> growth_local;
  int grLocSize;
  // double wmcIni, wmcEnd, dwmcval;
  int wmcIni, wmcEnd, dwmcval;
  Site *stenb;
  bool phaseid_exist;
  vector<int> inGrowInterfacePos;
  int pos;
  int affInt;

  int steId = ste->getId();

  if (verbose_) {
    cout << endl
         << "    Lattice::transformSolSol() INI totalTRC/trc_cT  "
         << totalTRC << "/" << trc_cT << "  : steId = " << setw(3) << steId
         << "  =>  oldPhId -> newPhId :  " << setw(3) << oldPhId << "  ->  "
         << setw(3) << left << newPhId << endl;
    cout.flush();
  }

  // lattice update for ste site
  wmcIni = ste->getWmc0(); // normally wmcIni = 1
  // sumWmc -= 1; //ste->getWmc();

  removeDissolutionSite(ste, oldPhId);
  setMicroPhaseId(ste, newPhId);

  /// Weighted mean curvature (wmc) is changed by the difference
  /// between the growing phase's porosity and the template's porosity.
  ///

  //wmcEnd = chemSys_->getMicroPhasePorosityInt(newPhId); // normally wmcEnd = 1
  wmcEnd = microPhasePorosityInt_[newPhId];
  ste->setWmc0(wmcEnd);

  dwmcval = wmcEnd - wmcIni; // normally dwmcval = 0
  ste->dWmc(dwmcval);

  double steWmc = ste->getWmc();
  if (steWmc > 0.0) {
    addDissolutionSite(ste, newPhId);
  }

  for (int i = 0; i < NN_NNN; i++) {
    stenb = ste->nb(i);
    stenb->dWmc(dwmcval);
  }

  for (int i = 0; i < NN_NNN; i++) {
    stenb = ste->nb(i);
    nbpid = stenb->getMicroPhaseId();

    if (nbpid == ELECTROLYTEID) {
      nbid = stenb->getId();
      growth_local = site_[nbid].getGrowthPhases();
      grLocSize = growth_local.size();
      for (int ii = 0; ii < grLocSize; ii++) {
        phaseid_exist = false;
        for (int jj = 0; jj < NN_NNN; jj++) {
          nb_pid = site_[nbid].nb(jj)->getMicroPhaseId();
          if (growth_local[ii] == nb_pid) {
            phaseid_exist = true;
            break;
          }
          if (jj < NUM_NEAREST_NEIGHBORS) {
            // if (chemSys_->isGrowthTemplate(growth_local[ii], nb_pid)) {
            if (growthTemplate_[growth_local[ii]][nb_pid]) {
              phaseid_exist = true;
              break;
            }
          }
        }
        if (phaseid_exist == false) {
          removeGrowthSite_diss(stenb, growth_local[ii]);
        }
      }
    }
  }

  // update dissolutionVector & involved site
  site_[steId].setInDissolutionVectorPos(-1);

  for (int j = 0; j < NN_NNN; j++) {
    stenb = ste->nb(j);
    nbpid = stenb->getMicroPhaseId();
    if (nbpid == ELECTROLYTEID) {
      inGrowInterfacePos = stenb->getInGrowInterfacePosVector();
      for (int k = FIRST_SOLID; k < numMicroPhases_; k++) {
        pos = inGrowInterfacePos[k]; // instenb->getInGrowInterfacePos(k);
        if (pos != -1) {
          // aff = chemSys_->getAffinity(k, newPhId) -
          //       chemSys_->getAffinity(k, oldPhId);
          // interface_[k].updateAffinity(pos, aff);
          affInt = affinityInt_[k][newPhId] - affinityInt_[k][oldPhId];
          interface_[k].updateAffinityInt(pos, affInt);
        } else {
          if (k == newPhId)
            addGrowthSite(stenb, newPhId);
        }
      }
    }
  }
  for (int j = 0; j < NUM_NEAREST_NEIGHBORS; j++) {
    stenb = ste->nb(j);
    nbpid = stenb->getMicroPhaseId();
    if (nbpid == ELECTROLYTEID) {
      for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
           phaseTmpl++) {
        // if (chemSys_->isGrowthTemplate(phaseTmpl, newPhId))
        if (growthTemplate_[phaseTmpl][newPhId])
          if (stenb->getInGrowInterfacePos(phaseTmpl) == -1)
            addGrowthSite(stenb, phaseTmpl);
      }
    }
  }
}

vector<int> Lattice::transformLiqSol(Site *ste, int growPhID, int totalTRC) {

  //*** for controll
  int static trc_gT;
  trc_gT++;

  vector<int> inGrowInterfacePos;
  vector<int> plist;
  int plistSize;

  vector<int> extVect;

  // double wmcIni, wmcEnd, dwmcval;
  // double steWmc, stenbWmc;
  int wmcIni, wmcEnd, dwmcval;
  int steWmc, stenbWmc;
  // double rng;

  int mPhId;
  int pos;
  int aftyInt;

  Site *stenb;

  // if (verbose_) {
  //   int steId = ste->getId();
  //   string growPhName = chemSys_->getMicroPhaseName(growPhID);
  //   cout << endl
  //        << "    Lattice::transformLiqSol() INI totalTRC/trc_gT  "
  //        << totalTRC << "/" << trc_gT << "  :  steId = "
  //        << steId << "  => growPhID = " << setw(3) << growPhID
  //        << "   growPhName = " << setw(15) << left << growPhName
  //        << "   count_ = " << setw(8) << right << count_[growPhID]
  //        << "   growthInterfaceSize_ = " << setw(8)
  //        << growthInterfaceSize_[growPhID] << endl;
  //   cout.flush();
  // }

  mPhId = ste->getMicroPhaseId(); // always ELECTROLYTEID !!

  if (mPhId != ELECTROLYTEID) {
    int steId = ste->getId();
    cout << endl
         << "Lattice::transformLiqSol() error : mPhId != ELECTROLYTEID  =>"
            "  steId/mPhId/growPhID/totalTRC/trc_gT :  "
         << steId << " / " << mPhId << " / " << growPhID << "   " << totalTRC
         << "/" << trc_gT << endl;
    cout << "STOP" << endl;
    exit(0);
  }

  wmcIni = ste->getWmc0(); // always 1 !!!

  plist = ste->getGrowthPhases();
  plistSize = plist.size();
  for (int j = 0; j < plistSize; j++) {
    removeGrowthSite_grow(ste, plist[j]);
  }
  ste->clearGrowth();
  setMicroPhaseId(ste, growPhID);

  ///
  /// Weighted mean curvature (wmc) is changed by the difference
  /// between the growing phase's porosity and the template's porosity.

  // double por = chemSys_->getMicroPhasePorosity(growPhID);
  wmcEnd = microPhasePorosityInt_[growPhID];
  ste->setWmc0(wmcEnd);

  dwmcval = wmcEnd - wmcIni;
  ste->dWmc(dwmcval);

  ///
  /// Now that the site has been added, it is eligible for dissolution
  /// later on, so we add it to the list of dissolution sites.
  ///

  steWmc = ste->getWmc();

  if (steWmc > 0) {
    addDissolutionSite(ste, growPhID);
  }

  ///
  /// Update the wmc of the neighboring sites.  This can be done
  /// because the wmc is originally calculated within a box around each
  /// site, so any time the id of a site within that box changes, it
  /// will change the wmc of the site at the box's center.
  ///
  /// NN_NNN = NUM_NEAREST_NEIGHBORS +  NUM_SECONDNEAREST_NEIGHBORS;
  ///

  for (int j = 0; j < NN_NNN; j++) {
    stenb = ste->nb(j);
    stenb->dWmc(dwmcval);
    stenbWmc = stenb->getWmc();
    if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
      inGrowInterfacePos = stenb->getInGrowInterfacePosVector();
      for (int k = FIRST_SOLID; k < numMicroPhases_; k++) {
        pos = inGrowInterfacePos[k]; // instenb->getInGrowInterfacePos(k);
        if (pos != -1) {
          // afty = chemSys_->getAffinity(k, growPhID);
          // interface_[k].updateAffinity(pos, afty);
          aftyInt = affinityInt_[k][growPhID];
          interface_[k].updateAffinityInt(pos, aftyInt);
        } else {
          if (k == growPhID) {
            addGrowthSite(stenb, growPhID);
          }
        }
      }
    } else if ((stenbWmc == 0.0) &&
               (stenb->getMicroPhaseId() > ELECTROLYTEID)) {
      extVect.push_back(stenb->getId());
      removeDissolutionSite(stenb, stenb->getMicroPhaseId());
    }
  }

  for (int j = 0; j < NUM_NEAREST_NEIGHBORS; j++) {
    stenb = ste->nb(j);
    if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
      for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
           phaseTmpl++) {
        // if (chemSys_->isGrowthTemplate(phaseTmpl, growPhID))
        if (growthTemplate_[phaseTmpl][growPhID])
          if (stenb->getInGrowInterfacePos(phaseTmpl) == -1)
            addGrowthSite(stenb, phaseTmpl);
      }
    }
  }

  // if (verbose_) {
  //   int steId = ste->getId();
  //   cout << endl
  //        << "    Lattice::transformLiqSol() END totalTRC/trc_gT  "
  //        << totalTRC << "/" << trc_gT << "  :  steId = "
  //        << steId << "  => growPhID = " << setw(3) << growPhID
  //        << "   growPhName = " << setw(15) << left << growPhName
  //        << "   count_ = " << setw(8) << right << count_[growPhID]
  //        << "   growthInterfaceSize_ = " << setw(8)
  //        << growthInterfaceSize_[growPhID] << endl;
  //   cout.flush();
  // }

  return extVect;
}

void Lattice::createGrowingVectSA() {

  chemSys_->getVectorsSA(growingVectSA_, shrinking_, volratios_);

  int gSize = growingVectSA_.size();
  int sSize = shrinking_.size();
  int vSize = volratios_.size();

  sizeGrowingVectSA_ = gSize;

  if (gSize == 0 || sSize == 0 || vSize == 0) {
    string msg = "gSize = 0 || sSize = 0 || vSize = 0";
    cout << endl
         << "Lattice::createGrowingVectSA() - error: " << msg << endl;
    cout << "     growingVectSA_.size() : gSize = " << gSize << endl;
    cout << "     shrinking_.size()     : sSize = " << sSize << endl;
    cout << "     volratios_.size()     : vSize = " << vSize << endl;
    cout << endl << "end program" << endl;
    throw DataException("Lattice", "createGrowingVectSA", msg);
    // exit(0);
  } else {

    cout << endl << "   Lattice::createGrowingVectSA() - create growingVectSA_, "
                    "shrinking_ and volratios_ vectors:" << endl;
    cout << "     growingVectSA_.size() = " << gSize << endl;
    cout << "     shrinking_.size()     = " << sSize << endl;
    cout << "     volratios_.size()     = " << vSize << endl;
    for (int i = 0; i < gSize; i++) {
      cout << endl << "     growingVectSA_[" << i << "] = " << growingVectSA_[i]
           << "  :  " << chemSys_->getMicroPhaseName(growingVectSA_[i])
           << endl;
      sSize = shrinking_[i].size();
      vSize = volratios_[i].size();
      if (sSize != vSize) {
        ostringstream ostr;
        ostr << i;
        string msg = "sSize != vSize for i = " + ostr.str();

        cout << endl
             << "Lattice::createGrowingVectSA() - error: " << msg << endl;
        cout << "     growingVectSA_.size() : gSize = " << gSize << endl;
        cout << "     shrinking_.size()     : sSize = " << sSize << endl;
        cout << "     volratios_.size()     : vSize = " << vSize << endl;
        cout << endl << "end program" << endl;
        throw DataException("Lattice", "createGrowingVectSA", msg);
        // exit(0);
      }
      for (int j = 0; j < sSize; j++) {
        cout << "       shrinking_[" << i << "," << j << "] = " << shrinking_[i][j]
             << "  :  " << chemSys_->getMicroPhaseName(shrinking_[i][j]) << endl;
        cout << "          volratios_[" << i << "," << j << "] = "
             << setprecision(3) << volratios_[i][j] << endl;
      }
    }
  }
  cout << setprecision(15);

  // cout << "     AFTMicroName       : " << AFTMicroName
  //      << " (id = " << chemSys_->getMicroPhaseId(AFTMicroName) << ") : " << endl;
  // cout << "       MonocarbMicroName  : " << MonocarbMicroName
  //      << " (id = " << chemSys_->getMicroPhaseId(MonocarbMicroName) << ")"
  //      << endl;
  //      shrinking_[i].push_back(chemSys_->getMicroPhaseId(MonocarbMicroName));
  //      volratios_[i].push_back(2.699);
  // cout << "       HemicarbMicroName  : " << HemicarbMicroName
  //      << " (id = " << chemSys_->getMicroPhaseId(HemicarbMicroName) << ")"
  //      << endl;
  //      shrinking_[i].push_back(chemSys_->getMicroPhaseId(HemicarbMicroName));
  //      volratios_[i].push_back(2.485);
  // cout << "       HydrotalcMicroName : " << HydrotalcMicroName
  //      << " (id = " << chemSys_->getMicroPhaseId(HydrotalcMicroName) << ")"
  //      << endl;
  //      shrinking_[i].push_back(chemSys_->getMicroPhaseId(HydrotalcMicroName));
  //      volratios_[i].push_back(3.211);
  // cout << "       MonosulfMicroName  : " << MonosulfMicroName
  //      << " (id = " << chemSys_->getMicroPhaseId(MonosulfMicroName) << ")"
  //      << endl;
  //      shrinking_[i].push_back(chemSys_->getMicroPhaseId(MonosulfMicroName));
  //      volratios_[i].push_back(2.288);
  //
  // cout << "       CSHMicroName       : " << CSHMicroName
  //      << " (id = " << chemSys_->getMicroPhaseId(CSHMicroName) << ")" << endl;
}

vector<int> Lattice::writeSubVolume(string fileName, Site *centerste,
                                    int size) {
  ofstream out(fileName.c_str());

  // out << "Version: 7.0" << endl;
  // out << "X_Size: 3" << endl;
  // out << "Y_Size: 3" << endl;
  // out << "Z_Size: 3" << endl;
  // out << "Image_Resolution: 1" << endl;

  out << VERSIONSTRING << " " << thamesVersion_ << endl;
  out << XSIZESTRING << " " << 3 << endl;
  out << YSIZESTRING << " " << 3 << endl;
  out << ZSIZESTRING << " " << 3 << endl;
  out << IMGRESSTRING << " " << resolution_ << endl;

  vector<int> alnb = getNeighborhood(centerste->getId(), size);
  int alnbsize = alnb.size();

  int phaseid;
  for (int j = 0; j < alnbsize; j++) {
    phaseid = site_[alnb[j]].getMicroPhaseId();
    out << phaseid << endl;
  }
  out.close();

  return alnb;
}

void Lattice::applyExpansion(vector<int> alnb, double exp) {
  Site *ste;
  vector<double> expval(3, exp);
  // vector<int> coordin(3, 0);
  // if (exp > 0.0) {
  int size = alnb.size();
  for (int i = 0; i < size; i++) {
    ste = &site_[alnb[i]];
    if (exp > ste->getExpansionStrain()) {
      ste->setExpansionStrain(exp);
      map<int, vector<double>>::iterator p = expansion_.find(alnb[i]);
      if (p != expansion_.end()) {
        (p->second)[0] = exp;
        (p->second)[1] = exp;
        (p->second)[2] = exp;
      } else {
        expansion_.insert(make_pair(alnb[i], expval));
      }
    }

    // map<int, vector<int>>::iterator pp = expansion_coordin_.find(alnb[i]);
    // if (pp == expansion_coordin_.end()) {
    //   coordin[0] = ste->getX();
    //   coordin[1] = ste->getY();
    //   coordin[2] = ste->getZ();
    //   expansion_coordin_.insert(make_pair(alnb[i], coordin));
    // }
  }
  // }

  return;
}

int Lattice::findDomainSize(int siteid, int maxsize) {

  int boxhalf = maxsize / 2;
  int nfound = 0;

  int phaseid = (site_[siteid]).getMicroPhaseId();
  int qx = (site_[siteid]).getX();
  int qy = (site_[siteid]).getY();
  int qz = (site_[siteid]).getZ();

  int qxlo = qx - boxhalf;
  int qxhi = qx + boxhalf;
  int qylo = qy - boxhalf;
  int qyhi = qy + boxhalf;
  int qzlo = qz - boxhalf;
  int qzhi = qz + boxhalf;

  /***
   *    Count the number of requisite sites in the
   *    3-D cube box using whatever boundaries are specified
   ***/

  for (int ix = qxlo; ix <= qxhi; ix++) {
    for (int iy = qylo; iy <= qyhi; iy++) {
      for (int iz = qzlo; iz <= qzhi; iz++) {

        /// Count if phase id only

        if (site_[getIndex(ix, iy, iz)].getMicroPhaseId() == phaseid) {
          nfound++;
        }
      }
    }
  }

  return nfound;
}

vector<int> Lattice::findDomainSizeDistribution(int phaseid, const int numsites,
                                                int maxsize,
                                                int sortorder = 0) {
  // if sortorder is 0 => sorting in descending order

  int domainsize = 0;
  int ix, iy, iz;
  int qx, qy, qz;
  int qxlo, qxhi, qylo, qyhi, qzlo, qzhi;
  int sizeCateg, pos, stId;
  double rng;

  int boxhalf = maxsize / 2;
  int dim = maxsize;
  if (maxsize % 2 == 0)
    dim++;
  int totSize = static_cast<int>(pow(dim, 3)) + 1;

  vector<int> init;
  vector<vector<int>> siteDomainSizeDistribution(totSize, init);
  vector<int> selectedSites;

  for (int i = 0; i < numSites_; i++) {
    if (site_[i].getMicroPhaseId() == phaseid) {

      qx = (site_[i]).getX();
      qy = (site_[i]).getY();
      qz = (site_[i]).getZ();

      qxlo = qx - boxhalf;
      qxhi = qx + boxhalf;
      qylo = qy - boxhalf;
      qyhi = qy + boxhalf;
      qzlo = qz - boxhalf;
      qzhi = qz + boxhalf;

      /***
       *    Count the number of requisite sites in the
       *    3-D cube box using whatever boundaries are specified
       ***/

      domainsize = 0;
      for (ix = qxlo; ix <= qxhi; ix++) {
        for (iy = qylo; iy <= qyhi; iy++) {
          for (iz = qzlo; iz <= qzhi; iz++) {

            /// Count if phase id only

            if (site_[getIndex(ix, iy, iz)].getMicroPhaseId() == phaseid) {
              domainsize++;
            }
          }
        }
      }

      siteDomainSizeDistribution[domainsize].push_back(i);
    }
  }

  int j;
  int contor = 0;
  int diff = numsites;
  if (sortorder == 0) { // emptyPorosity
    j = totSize;
    while (diff > 0) {
      sizeCateg = siteDomainSizeDistribution[j - 1].size();
      if (diff < sizeCateg) {
        for (int k = 0; k < diff; k++) {
          rng = callRNG();
          pos = static_cast<int>(rng * (sizeCateg - 1));
          stId = siteDomainSizeDistribution[j - 1][pos];
          selectedSites.push_back(stId);
          siteDomainSizeDistribution[j - 1][pos] =
              siteDomainSizeDistribution[j - 1][sizeCateg - 1];

          sizeCateg--;
          contor++;
        }
      } else {
        if (sizeCateg > 0) {
          selectedSites.insert(selectedSites.end(),
                               siteDomainSizeDistribution[j - 1].begin(),
                               siteDomainSizeDistribution[j - 1].end());
          contor += sizeCateg;
        }
        j--;
      }
      diff = numsites - contor;
      if (j == 0)
        break;
    }
  } else { // fillPorosity
    j = 1;
    while (diff > 0) {
      sizeCateg = siteDomainSizeDistribution[j].size();
      if (diff < sizeCateg) {
        for (int k = 0; k < diff; k++) {
          rng = callRNG();
          pos = static_cast<int>(rng * (sizeCateg - 1));
          stId = siteDomainSizeDistribution[j][pos];
          selectedSites.push_back(stId);
          siteDomainSizeDistribution[j][pos] =
              siteDomainSizeDistribution[j][sizeCateg - 1];

          sizeCateg--;
          contor++;
        }
      } else {
        if (sizeCateg > 0) {
          selectedSites.insert(selectedSites.end(),
                               siteDomainSizeDistribution[j].begin(),
                               siteDomainSizeDistribution[j].end());
          contor += sizeCateg;
        }
        j++;
      }
      diff = numsites - contor;
      if (j > totSize)
        break;
    }
  }

  return (selectedSites);
}

void Lattice::findIsolatedClusters(void) {
  // voxels without contact with electrolyte i.e. voxels having low
  // probability to dissolve in this step these voxels will not be
  // send to GEMS computing vfrac we use canDissolve vector instead
  // count_ vector
  vector<int> canDissolve;
  int canDissLast;
  bool siteDiss;
  int i, j, jj;

  for (i = 0; i < numMicroPhases_; i++) {
    canDissolve.push_back(0);
    if (i >= FIRST_SOLID) {
      canDissLast = -1;
      while (canDissLast != canDissolve[i]) {
        canDissLast = canDissolve[i];
        for (j = 0; j < numSites_; j++) {
          if (site_[j].getMicroPhaseId() == i) {
            if (site_[j].getVisit() == 0) {
              if (site_[j].getWmc() > 0) {
                canDissolve[i]++;
                site_[j].setVisit(1);
                for (jj = 0; jj < NN_NNN; jj++) {
                  if (site_[j].nb(jj)->getMicroPhaseId() == i &&
                      site_[j].nb(jj)->getVisit() == 0) {
                    site_[j].nb(jj)->setVisit(1);
                    canDissolve[i]++;
                  }
                }
              } else {
                siteDiss = false;
                for (jj = 0; jj < NN_NNN; jj++) {
                  if (site_[j].nb(jj)->getMicroPhaseId() == i &&
                      site_[j].nb(jj)->getVisit() == 1) {
                    canDissolve[i]++;
                    site_[j].setVisit(1);
                    siteDiss = true;
                    break;
                  }
                }
                if (siteDiss) {
                  for (jj = 0; jj < NN_NNN; jj++) {
                    if (site_[j].nb(jj)->getMicroPhaseId() == i &&
                        site_[j].nb(jj)->getVisit() == 0) {
                      site_[j].nb(jj)->setVisit(1);
                      canDissolve[i]++;
                    }
                  }
                }
              }
            }
          }
        } // for
      } // while
    } // if
  } // for

  for (j = 0; j < numSites_; j++) {
    site_[j].setVisit(0);
  }

  cout << endl
       << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
       << endl;
  for (i = 0; i < numMicroPhases_; i++) {
    cout << "     i = " << i << "      canDissolve[" << i
         << "] = " << canDissolve[i] << endl;
  }
  cout << endl
       << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
       << endl;

  return;
}

void Lattice::populateElementData(void) {
  chemElem elem;

  elem.z = 1;
  elem.symb = "H";
  elem.mass = 1.0080;
  cfgElem_.push_back(elem);

  elem.z = 2;
  elem.symb = "He";
  elem.mass = 4.0026;
  cfgElem_.push_back(elem);

  elem.z = 3;
  elem.symb = "Li";
  elem.mass = 6.94;
  cfgElem_.push_back(elem);

  elem.z = 4;
  elem.symb = "Be";
  elem.mass = 9.0122;
  cfgElem_.push_back(elem);

  elem.z = 5;
  elem.symb = "B";
  elem.mass = 10.81;
  cfgElem_.push_back(elem);

  elem.z = 6;
  elem.symb = "C";
  elem.mass = 12.011;
  cfgElem_.push_back(elem);

  elem.z = 7;
  elem.symb = "N";
  elem.mass = 14.007;
  cfgElem_.push_back(elem);

  elem.z = 8;
  elem.symb = "O";
  elem.mass = 15.999;
  cfgElem_.push_back(elem);

  elem.z = 9;
  elem.symb = "F";
  elem.mass = 18.998;
  cfgElem_.push_back(elem);

  elem.z = 10;
  elem.symb = "Ne";
  elem.mass = 20.180;
  cfgElem_.push_back(elem);

  elem.z = 11;
  elem.symb = "Na";
  elem.mass = 22.990;
  cfgElem_.push_back(elem);

  elem.z = 12;
  elem.symb = "Mg";
  elem.mass = 24.305;
  cfgElem_.push_back(elem);

  elem.z = 13;
  elem.symb = "Al";
  elem.mass = 26.982;
  cfgElem_.push_back(elem);

  elem.z = 14;
  elem.symb = "Si";
  elem.mass = 28.085;
  cfgElem_.push_back(elem);

  elem.z = 15;
  elem.symb = "P";
  elem.mass = 30.974;
  cfgElem_.push_back(elem);

  elem.z = 16;
  elem.symb = "S";
  elem.mass = 32.06;
  cfgElem_.push_back(elem);

  elem.z = 17;
  elem.symb = "Cl";
  elem.mass = 35.45;
  cfgElem_.push_back(elem);

  elem.z = 18;
  elem.symb = "Ar";
  elem.mass = 39.95;
  cfgElem_.push_back(elem);

  elem.z = 19;
  elem.symb = "K";
  elem.mass = 39.098;
  cfgElem_.push_back(elem);

  elem.z = 20;
  elem.symb = "Ca";
  elem.mass = 40.078;
  cfgElem_.push_back(elem);

  elem.z = 21;
  elem.symb = "Sc";
  elem.mass = 44.956;
  cfgElem_.push_back(elem);

  elem.z = 22;
  elem.symb = "Ti";
  elem.mass = 47.867;
  cfgElem_.push_back(elem);

  elem.z = 23;
  elem.symb = "V";
  elem.mass = 50.942;
  cfgElem_.push_back(elem);

  elem.z = 24;
  elem.symb = "Cr";
  elem.mass = 51.996;
  cfgElem_.push_back(elem);

  elem.z = 25;
  elem.symb = "Mn";
  elem.mass = 54.938;
  cfgElem_.push_back(elem);

  elem.z = 26;
  elem.symb = "Fe";
  elem.mass = 55.845;
  cfgElem_.push_back(elem);

  elem.z = 27;
  elem.symb = "Co";
  elem.mass = 58.933;
  cfgElem_.push_back(elem);

  elem.z = 28;
  elem.symb = "Ni";
  elem.mass = 58.693;
  cfgElem_.push_back(elem);

  elem.z = 29;
  elem.symb = "Cu";
  elem.mass = 63.546;
  cfgElem_.push_back(elem);

  elem.z = 30;
  elem.symb = "Zn";
  elem.mass = 65.38;
  cfgElem_.push_back(elem);

  elem.z = 31;
  elem.symb = "Ga";
  elem.mass = 69.723;
  cfgElem_.push_back(elem);

  elem.z = 32;
  elem.symb = "Ge";
  elem.mass = 72.630;
  cfgElem_.push_back(elem);

  elem.z = 33;
  elem.symb = "As";
  elem.mass = 74.922;
  cfgElem_.push_back(elem);

  elem.z = 34;
  elem.symb = "Se";
  elem.mass = 78.971;
  cfgElem_.push_back(elem);

  elem.z = 35;
  elem.symb = "Br";
  elem.mass = 79.904;
  cfgElem_.push_back(elem);

  elem.z = 36;
  elem.symb = "Kr";
  elem.mass = 83.798;
  cfgElem_.push_back(elem);

  elem.z = 37;
  elem.symb = "Rb";
  elem.mass = 85.468;
  cfgElem_.push_back(elem);

  elem.z = 38;
  elem.symb = "Sr";
  elem.mass = 87.62;
  cfgElem_.push_back(elem);

  elem.z = 39;
  elem.symb = "Y";
  elem.mass = 88.906;
  cfgElem_.push_back(elem);

  elem.z = 40;
  elem.symb = "Zr";
  elem.mass = 91.224;
  cfgElem_.push_back(elem);

  elem.z = 41;
  elem.symb = "Nb";
  elem.mass = 92.906;
  cfgElem_.push_back(elem);

  elem.z = 42;
  elem.symb = "Mo";
  elem.mass = 95.95;
  cfgElem_.push_back(elem);

  elem.z = 43;
  elem.symb = "Tc";
  elem.mass = 97.0;
  cfgElem_.push_back(elem);

  elem.z = 44;
  elem.symb = "Ru";
  elem.mass = 101.07;
  cfgElem_.push_back(elem);

  elem.z = 45;
  elem.symb = "Rh";
  elem.mass = 102.91;
  cfgElem_.push_back(elem);

  elem.z = 46;
  elem.symb = "Pd";
  elem.mass = 106.42;
  cfgElem_.push_back(elem);

  elem.z = 47;
  elem.symb = "Ag";
  elem.mass = 107.87;
  cfgElem_.push_back(elem);

  elem.z = 48;
  elem.symb = "Cd";
  elem.mass = 112.41;
  cfgElem_.push_back(elem);

  elem.z = 49;
  elem.symb = "In";
  elem.mass = 114.82;
  cfgElem_.push_back(elem);

  elem.z = 50;
  elem.symb = "Sn";
  elem.mass = 118.71;
  cfgElem_.push_back(elem);

  elem.z = 51;
  elem.symb = "Sb";
  elem.mass = 121.76;
  cfgElem_.push_back(elem);

  elem.z = 52;
  elem.symb = "Te";
  elem.mass = 127.60;
  cfgElem_.push_back(elem);

  elem.z = 53;
  elem.symb = "I";
  elem.mass = 126.90;
  cfgElem_.push_back(elem);

  elem.z = 54;
  elem.symb = "Xe";
  elem.mass = 131.29;
  cfgElem_.push_back(elem);

  elem.z = 55;
  elem.symb = "Cs";
  elem.mass = 132.91;
  cfgElem_.push_back(elem);

  elem.z = 56;
  elem.symb = "Ba";
  elem.mass = 137.33;
  cfgElem_.push_back(elem);

  elem.z = 57;
  elem.symb = "La";
  elem.mass = 138.91;
  cfgElem_.push_back(elem);

  elem.z = 58;
  elem.symb = "Ce";
  elem.mass = 140.12;
  cfgElem_.push_back(elem);

  elem.z = 59;
  elem.symb = "Pr";
  elem.mass = 140.91;
  cfgElem_.push_back(elem);

  elem.z = 60;
  elem.symb = "Nd";
  elem.mass = 144.24;
  cfgElem_.push_back(elem);

  elem.z = 61;
  elem.symb = "Pm";
  elem.mass = 145.0;
  cfgElem_.push_back(elem);

  elem.z = 62;
  elem.symb = "Sm";
  elem.mass = 150.36;
  cfgElem_.push_back(elem);

  elem.z = 63;
  elem.symb = "Eu";
  elem.mass = 151.96;
  cfgElem_.push_back(elem);

  elem.z = 64;
  elem.symb = "Gd";
  elem.mass = 157.25;
  cfgElem_.push_back(elem);

  elem.z = 65;
  elem.symb = "Tb";
  elem.mass = 158.93;
  cfgElem_.push_back(elem);

  elem.z = 66;
  elem.symb = "Dy";
  elem.mass = 162.50;
  cfgElem_.push_back(elem);

  elem.z = 67;
  elem.symb = "Ho";
  elem.mass = 164.93;
  cfgElem_.push_back(elem);

  elem.z = 68;
  elem.symb = "Er";
  elem.mass = 167.26;
  cfgElem_.push_back(elem);

  elem.z = 69;
  elem.symb = "Tm";
  elem.mass = 168.93;
  cfgElem_.push_back(elem);

  elem.z = 70;
  elem.symb = "Yb";
  elem.mass = 173.05;
  cfgElem_.push_back(elem);

  elem.z = 71;
  elem.symb = "Lu";
  elem.mass = 174.97;
  cfgElem_.push_back(elem);

  elem.z = 72;
  elem.symb = "Hf";
  elem.mass = 178.49;
  cfgElem_.push_back(elem);

  elem.z = 73;
  elem.symb = "Ta";
  elem.mass = 180.95;
  cfgElem_.push_back(elem);

  elem.z = 74;
  elem.symb = "W";
  elem.mass = 183.84;
  cfgElem_.push_back(elem);

  elem.z = 75;
  elem.symb = "Re";
  elem.mass = 186.21;
  cfgElem_.push_back(elem);

  elem.z = 76;
  elem.symb = "Os";
  elem.mass = 190.23;
  cfgElem_.push_back(elem);

  elem.z = 77;
  elem.symb = "Ir";
  elem.mass = 192.22;
  cfgElem_.push_back(elem);

  elem.z = 78;
  elem.symb = "Pt";
  elem.mass = 195.08;
  cfgElem_.push_back(elem);

  elem.z = 79;
  elem.symb = "Au";
  elem.mass = 196.97;
  cfgElem_.push_back(elem);

  elem.z = 80;
  elem.symb = "Hg";
  elem.mass = 200.59;
  cfgElem_.push_back(elem);

  elem.z = 81;
  elem.symb = "Tl";
  elem.mass = 204.38;
  cfgElem_.push_back(elem);

  elem.z = 82;
  elem.symb = "Pb";
  elem.mass = 207.2;
  cfgElem_.push_back(elem);

  elem.z = 83;
  elem.symb = "Bi";
  elem.mass = 208.98;
  cfgElem_.push_back(elem);

  elem.z = 84;
  elem.symb = "Po";
  elem.mass = 209.0;
  cfgElem_.push_back(elem);

  elem.z = 85;
  elem.symb = "At";
  elem.mass = 210.0;
  cfgElem_.push_back(elem);

  elem.z = 86;
  elem.symb = "Rn";
  elem.mass = 222.0;
  cfgElem_.push_back(elem);

  elem.z = 87;
  elem.symb = "Fr";
  elem.mass = 223.0;
  cfgElem_.push_back(elem);

  elem.z = 88;
  elem.symb = "Ra";
  elem.mass = 226.0;
  cfgElem_.push_back(elem);

  elem.z = 89;
  elem.symb = "Ac";
  elem.mass = 227.0;
  cfgElem_.push_back(elem);

  elem.z = 90;
  elem.symb = "Th";
  elem.mass = 232.04;
  cfgElem_.push_back(elem);

  elem.z = 91;
  elem.symb = "Pa";
  elem.mass = 231.04;
  cfgElem_.push_back(elem);

  elem.z = 92;
  elem.symb = "U";
  elem.mass = 238.03;
  cfgElem_.push_back(elem);

  elem.z = 93;
  elem.symb = "Np";
  elem.mass = 237.0;
  cfgElem_.push_back(elem);

  elem.z = 94;
  elem.symb = "Pu";
  elem.mass = 244.0;
  cfgElem_.push_back(elem);

  elem.z = 95;
  elem.symb = "Am";
  elem.mass = 243.0;
  cfgElem_.push_back(elem);

  elem.z = 96;
  elem.symb = "Cm";
  elem.mass = 247.0;
  cfgElem_.push_back(elem);

  elem.z = 97;
  elem.symb = "Bk";
  elem.mass = 247.0;
  cfgElem_.push_back(elem);

  elem.z = 98;
  elem.symb = "Cf";
  elem.mass = 251.0;
  cfgElem_.push_back(elem);

  elem.z = 99;
  elem.symb = "Es";
  elem.mass = 252.0;
  cfgElem_.push_back(elem);

  elem.z = 100;
  elem.symb = "Fm";
  elem.mass = 257.0;
  cfgElem_.push_back(elem);

  elem.z = 101;
  elem.symb = "Md";
  elem.mass = 258.0;
  cfgElem_.push_back(elem);

  elem.z = 102;
  elem.symb = "No";
  elem.mass = 259.0;
  cfgElem_.push_back(elem);

  elem.z = 103;
  elem.symb = "Lr";
  elem.mass = 266.0;
  cfgElem_.push_back(elem);

  elem.z = 104;
  elem.symb = "Rf";
  elem.mass = 267.0;
  cfgElem_.push_back(elem);

  elem.z = 105;
  elem.symb = "Db";
  elem.mass = 268.0;
  cfgElem_.push_back(elem);

  elem.z = 106;
  elem.symb = "Sg";
  elem.mass = 269.0;
  cfgElem_.push_back(elem);

  elem.z = 107;
  elem.symb = "Bh";
  elem.mass = 270.0;
  cfgElem_.push_back(elem);

  elem.z = 108;
  elem.symb = "Hs";
  elem.mass = 269.0;
  cfgElem_.push_back(elem);

  elem.z = 109;
  elem.symb = "Mt";
  elem.mass = 278.0;
  cfgElem_.push_back(elem);

  elem.z = 110;
  elem.symb = "Ds";
  elem.mass = 281.0;
  cfgElem_.push_back(elem);

  elem.z = 111;
  elem.symb = "Rg";
  elem.mass = 282.0;
  cfgElem_.push_back(elem);

  elem.z = 112;
  elem.symb = "Cn";
  elem.mass = 285.0;
  cfgElem_.push_back(elem);

  elem.z = 113;
  elem.symb = "Nh";
  elem.mass = 286.0;
  cfgElem_.push_back(elem);

  elem.z = 114;
  elem.symb = "Fl";
  elem.mass = 289.0;
  cfgElem_.push_back(elem);

  elem.z = 115;
  elem.symb = "Mc";
  elem.mass = 290.0;
  cfgElem_.push_back(elem);

  elem.z = 116;
  elem.symb = "Lv";
  elem.mass = 293.0;
  cfgElem_.push_back(elem);

  elem.z = 117;
  elem.symb = "Ts";
  elem.mass = 294.0;
  cfgElem_.push_back(elem);

  elem.z = 118;
  elem.symb = "Og";
  elem.mass = 294.0;
  cfgElem_.push_back(elem);

  return;
}

void Lattice::checkSite(int stId) {
  // int phId = site_[stId].getMicroPhaseId();
  cout << endl << " Lattice::checkSite( " << stId << " ):" << endl;
  cout << "    phaseId : " << site_[stId].getMicroPhaseId() << endl;
  cout << "    inDissInterfacePos_ : " << site_[stId].getInDissInterfacePos()
       << endl;
  if (site_[stId].getInDissInterfacePos() != -1) {
    cout << "     in dissInterface on pos inDissInterfacePos_ : "
         << interface_[site_[stId].getMicroPhaseId()].getDissolutionSitesId(
                site_[stId].getInDissInterfacePos())
         << endl;
  }

  std::vector<int> growth = site_[stId].getGrowthPhases();
  int size = growth.size();
  int k;
  cout << endl << " growth_.size() : " << size << endl;
  for (k = 0; k < size; k++) {
    cout << "       k = " << k << "   growth_[k] = " << growth[k] << endl;
  }
  cout << endl << " inGrowInterfacePos_ : " << endl;
  for (k = 0; k < numMicroPhases_; k++) {
    cout << "       k = " << k << "   site_[stId].getInGrowInterfacePos(k) = "
         << site_[stId].getInGrowInterfacePos(k) << endl;
  }

  cout << endl
       << "     in growInterfaces on pos inGrowInterfacePos_ :" << endl;
  for (k = 0; k < numMicroPhases_; k++) {
    cout << "       k_ = " << k << endl;
    cout.flush();
    if (site_[stId].getInGrowInterfacePos(k) > -1) {
      size = growthInterfaceSize_[k];
      cout << "       k = " << k
           << "   pos = " << site_[stId].getInGrowInterfacePos(k)
           << "   size = " << size << endl;
      cout.flush();
      cout << "            siteId in grInt = "
           << interface_[k].getGrowthSitesId(
                  site_[stId].getInGrowInterfacePos(k))
           << endl;
      cout.flush();
    }
  }
}

void Lattice::addSeedCSHQ(bool seedMassCEM, bool seedMassC3S,
                          bool seedMassC2S, double massFraction) {

  double rng;
  int cshqId = chemSys_->getMicroPhaseId("CSHQ");
  double density = 0, porosity = 0, avMolarMass = 0.0, avMolarVol = 0.0;

  std::vector<double> DCporosities = chemSys_->getMicroPhaseDCPorosities(cshqId);

  double numMoles = 0;


  double conc;
  double molarVol, molarMass;
  cout << endl << "  Lattice::addSeedCSHQ - microPhaseId = "
       << cshqId << " => DC components :" << endl;

  std::vector<int> compDC = chemSys_->getMicroPhaseDCMembers(cshqId);
  int sizeCompDC = compDC.size();
  for (int i = 0; i < sizeCompDC; i++) {
    numMoles = chemSys_->getDCMoles(compDC[i]);
    molarVol = chemSys_->getDCMolarVolume(compDC[i]);
    molarMass = chemSys_->getDCMolarMass(compDC[i]);

    conc = chemSys_->getDCConcentration(compDC[i]);
    avMolarMass += (conc * molarMass);    // g/mol
    avMolarVol += (conc * molarVol);      // m3/mol
    cout << "    " << i << "\tDCId: " << compDC[i]
         << "    DCName = " << setw(12) << left << chemSys_->getDCName(compDC[i])
         << "    Molar mass = " << molarMass
         << " g/mol     Molar volume = " << molarVol
         << " m3/mol     porosity = " <<  DCporosities[i] << endl;
  }
  density = avMolarMass / avMolarVol / 1.0e6; // g/cm3
  porosity = chemSys_->getMicroPhasePorosity(cshqId);
  cout << endl << "  Lattice::addSeedCSHQ - CSHQ : "
       << " avMolarMass = " << avMolarMass
       << " g/mol    avMolarVol = " << avMolarVol
       << " m3/mol    density = " << density
       << " g/cm3    porosity = " << porosity
       << endl;

  microPhasePorosityInt_[cshqId] = chemSys_->getMicroPhasePorosityInt(cshqId);

  std::vector<std::string> microPhaseName = chemSys_->getMicroPhaseName();
  std::vector<double> microPhaseMass = chemSys_->getMicroPhaseMass();

  double totalInitMass = 0.0;
  double seedMassCSHQ = 0.0;


  if (seedMassCEM == true) {
    for (int i = FIRST_SOLID; i < numMicroPhases_; i++) {
      if (chemSys_->isCementComponent(i))
        totalInitMass += microPhaseMass[i];
    }
  } else if (seedMassC3S == true && seedMassC2S == true) {
    for (int i = FIRST_SOLID; i < numMicroPhases_; i++) {
      if (microPhaseName[i] == "Alite" || microPhaseName[i] == "Belite")
        totalInitMass += microPhaseMass[i];
    }
  } else if (seedMassC3S == true) {
    for (int i = FIRST_SOLID; i < numMicroPhases_; i++) {
      if (microPhaseName[i] == "Alite") {
        totalInitMass = microPhaseMass[i];
        break;
      }
    }
  } else if (seedMassC2S == true) {
    for (int i = FIRST_SOLID; i < numMicroPhases_; i++) {
      if (microPhaseName[i] == "Belite") {
        totalInitMass = microPhaseMass[i];
        break;
      }
    }
  } else {
    // error!
    cout << endl << "  Lattice::addSeedCSHQ - ERROR : " << endl;
    cout << "    you decided to create some CSHQ seeds but all next variables are false :"
         << endl;
    cout << "      seedMassCEM, seedMassC3S, seedMassC2S" << endl;
    cout << endl << "    => at least one of them must be true!" << endl;
    cout << endl << "  exit" << endl;
    exit(0);
  }

  if (massFraction <= 0.0) {
    // error!
    cout << endl << "  Lattice::addSeedCSHQ - ERROR : " << endl;
    cout << "    massFraction <= 0.0 !!  :  massFraction = " << massFraction
         << "   while massFraction must be positive!!!" << endl;
    cout << endl << "  exit" << endl;
    exit(0);
  }

  seedMassCSHQ = totalInitMass * massFraction;
  cout << endl << "  Lattice::addSeedCSHQ : totalInitMass = "
       << totalInitMass << "  &  massFraction = " << massFraction
       << "  =>  seedMassCSHQ = " << seedMassCSHQ << endl;

  cout << endl << "  Lattice::addSeedCSHQ - VOIDS/ELECTROLYTE : count_[1] = "
       << count_[1] << "  &  count_[2] = " <<  count_[2] << endl;

  double numMolesCSHQ = 0;
  double seedSolidVolCSHQ = 0;
  if (avMolarMass > 0) {
    double sum = 0;
    numMolesCSHQ = seedMassCSHQ / avMolarMass;
    cout << endl << "  Lattice::addSeedCSHQ : numMolesCSHQ = "
         << numMolesCSHQ << " :" << endl;
    for (int i = 0; i < sizeCompDC; i++) {
      conc = chemSys_->getDCConcentration(compDC[i]);
      numMoles = numMolesCSHQ * conc;
      cout << "    i = " << i << "    DCId = " << compDC[i]
           << "    DCName = " << setw(12) << left << chemSys_->getDCName(compDC[i])
           << "    numMoles = " << numMoles << endl;
      chemSys_->setDCMoles(compDC[i], numMoles);
      sum += numMoles;
    }
    cout << "  Lattice::addSeedCSHQ : sum(numMoles)  = " << sum << endl;
  } else {
    // error!
    cout << endl
         << "  Lattice::addSeedCSHQ - ERROR :" << endl;
    cout << "    not positive average molar mass of CSHQ!! : avMolarMass = "
         << avMolarMass << endl;
    cout << endl << "  exit" << endl;
    exit(0);
  }

  // in Lattice:
  seedSolidVolCSHQ = numMolesCSHQ * avMolarVol; // [mol] * [m3/mol]
  chemSys_->setMicroPhaseVolume(cshqId, seedSolidVolCSHQ);
  cout << endl << "  Lattice::addSeedCSHQ : seedSolidVolCSHQ = "
       << seedSolidVolCSHQ << " m3"<< endl;

  double varWaterMoles = seedSolidVolCSHQ / waterMolarVol_;
  cout << "  Lattice::addSeedCSHQ : varWaterMoles    = " << varWaterMoles << endl;
  double waterMoles = chemSys_->getDCMoles(waterDCId_);
  cout << "  Lattice::addSeedCSHQ : waterMoles       = " << waterMoles << " (old)" << endl;
  waterMoles -= varWaterMoles;
  cout << "  Lattice::addSeedCSHQ : waterMoles       = " << waterMoles << " (new)" << endl;
  chemSys_->setDCMoles(waterDCId_, waterMoles);
  chemSys_->setMicroPhaseVolume(ELECTROLYTEID, waterMoles * waterMolarVol_);

  double seedApVolCSHQ = seedSolidVolCSHQ / (1 - porosity);
  cout << endl << "  Lattice::addSeedCSHQ : seedApVolCSHQ    = " << seedApVolCSHQ
       << " m3" << endl;

  int numVoxelsApVolCSHQ = seedApVolCSHQ * numSites_ / initialMicrostructureVolume_;
  cout << "  Lattice::addSeedCSHQ  =>  numVoxelsApVolCSHQ = " << numVoxelsApVolCSHQ
       << " !!!" << endl;

  // cout << endl << "  Lattice::addSeedCSHQ - DCMoles:" << endl;
  // int numDCs = chemSys_->getNumDCs();
  // for (int i = 0; i < numDCs; i++) {
  //   cout << "    i = " << i << "   DCName = " << chemSys_->getDCName(i)
  //        << "   DCMoles = " << chemSys_->getDCMoles(i) << endl;
  // }

  // ******************************************************************
  //replace numVoxelsApVolCSHQ sites of electrolyte with CSHQ - random
  vector<int> waterNucSites;
  int numRemaining = numVoxelsApVolCSHQ;

  cout << endl
       << "  Lattice::addSeedCSHQ phaseID = " << cshqId
       << "   namePhase = CSHQ"  << endl;
  cout << "    sites to add        : numVoxelsApVolCSHQ = "
       << numVoxelsApVolCSHQ << endl;
  cout << "    sites in the system : count_[" << cshqId
       << "] = " << count_[cshqId] << endl;
  cout << "    void sites          : count_[0] = " << count_[0] << endl;
  cout << "    electrolyte sites   : count_[1] = " << count_[1] << endl;

  // make a list of all saturated sites capable of hosting nucleation
  waterNucSites.clear();
  for (int k = 0; k < numSites_; ++k) {
    if (site_[k].getMicroPhaseId() == ELECTROLYTEID) {
      waterNucSites.push_back(k);
    }
  }
  int sizeWS = waterNucSites.size();
  cout << endl
     << "  Lattice::addSeedCSHQ  - eligible water sites : sizeWS = "
     << sizeWS << endl;
  if (numVoxelsApVolCSHQ > sizeWS) {
    cout << endl
         << "      Lattice::addSeedCSHQ : requested number of seeds is"
            " larger than the number of all possible nucleation sites (sizeWSs)!"
         << endl;
    cout << "     numVoxelsApVolCSHQ > all possible nuc sites : "
         << numVoxelsApVolCSHQ << " > " << sizeWS << endl;
  }

  int siteId, pos ;
  while ((sizeWS > 0) && (numRemaining > 0)) {
    rng = callRNG();
    pos = static_cast<int>(rng * sizeWS);
    if (pos == sizeWS)
      pos--;

    siteId = waterNucSites[pos];
    site_[siteId].setMicroPhaseId(cshqId);
    count_[ELECTROLYTEID]--;
    count_[cshqId]++;

    numRemaining--;

    waterNucSites[pos] = waterNucSites[sizeWS - 1];
    waterNucSites.pop_back();

    sizeWS--;
  }

  if (numRemaining > 0) {
    cout << endl
         << "      Lattice::addSeedCSHQ - END ERROR : requested numSeedsCSHQ = "
         << numVoxelsApVolCSHQ << " while number of added seeds = "
         << numVoxelsApVolCSHQ - numRemaining << endl;
    cout << endl
         << "      exit" << endl;
    exit(0);

  } else {

    double totmicvol = 0.0;
    int totCount_ = 0;
    int mPhId;
    cout << endl
         << "  Lattice::addSeedCSHQ - volumes & count_ : " << endl;
    for (int i = 0; i < numMicroPhases_; i++) {
      mPhId = chemSys_->getMicroPhaseId(i);
      if (mPhId != VOIDID) {
        cout << "    mPhId = " << i << "\tvolume = " << chemSys_->getMicroPhaseVolume(mPhId)
             << "\t" << chemSys_->getMicroPhaseName(i)
             << " / count_ = " << count_[i] << endl;
        totmicvol += chemSys_->getMicroPhaseVolume(mPhId);
        totCount_ += count_[i];
      }
    }
    cout << "  Lattice::addSeedCSHQ - totCount_/totmicvol : "
         << totCount_ << " / " << totmicvol << endl;

    cout << endl
         << "  Lattice::addSeedCSHQ - END OK : a number of numSeedsCSHQ = "
         << numVoxelsApVolCSHQ << " CSHQ sites have been nucleated (seeds) !" << endl;
  }

  // cout << endl << "exit Lattice::addSeedCSHQ" << endl;
  // exit(0);
}

void Lattice::calcOneFaceAreaPerHundredGramSolid(void) {

  // Explanation of the normalization factor
  //
  // The surface area initially calculated is just in units of voxel faces
  // We want it to ultimately have units of m2 per 100 g of solid
  // So we must divide surface area by a factor with units of 100 g per m2 per
  // voxel
  //
  // areaPerFace_ has units of m2 per voxel face
  // surfaceArea_ has units of voxel faces per microstructure
  // volumePerVoxel_ has units of m3 per voxel
  // numSites_ has units of voxels
  // initSolidMass_ has units of g per cm3 of whole microstructure
  // 1.0e6 has units of cm3 per m3
  // So new surfaceArea values have units of
  // m2*(m3/cm3)(cm3/g)*(1/voxel)*(voxel/m3)
  // So the normalization below would be m2 / g of initial solid
  // But all our normalizations are on a basis of 100 g
  // So we further multiply by 100 to get m2 per 100 g of initial solid

  // double oneFaceAreaPerHundredGramSolid;

  if (initSolidMass_ > 0.0) {
    double oneFaceAreaPerMicrostructure =
        areaPerFace_ /
        (static_cast<double>(numSites_)); // m2/face/microstructure
    double oneFaceAreaPerMeterCubed = oneFaceAreaPerMicrostructure /
                                      volumePerVoxel_; // m2/face/
                                                       // (m3 of microstructure)

    double oneFaceAreaPerCmCubed =
        1.0e-6 * oneFaceAreaPerMeterCubed; // m2/face/
                                           // (cm3 of microstructure)
    double oneFaceAreaPerGramSolid =
        oneFaceAreaPerCmCubed / initSolidMass_; // m2/face/
                                                // (g of solid)
    oneFaceAreaPerHundredGramSolid_ =
        100.0 * oneFaceAreaPerGramSolid; // m2/face/
                                         // (100 g of solid)
  } else {
    string msg = "Divide by zero error:  initSolidMass_ = 0";
    throw FloatException("Lattice", "calcOneFaceAreaPerHundredGramSolid", msg);
  }

  cout << endl
       << "  Lattice::calcOneFaceAreaPerHundredGramSolid : "
          "oneFaceAreaPerHundredGramSolid_ = "
       << oneFaceAreaPerHundredGramSolid_ << endl;
}


void Lattice::calcSurfaceAreaInit(int phId) {
  // in m2 / (100g of solid)

  Site *ste, *stenb;

  double d_convFactDbl2IntPor = static_cast<double> (convFactDbl2IntPor_);
  for (int i = 0; i < numSites_; i++) {
    ste = &site_[i];
    if (ste->getMicroPhaseId() == phId) {
      for (int j = 0; j < NUM_NEAREST_NEIGHBORS; j++) {
        stenb = ste->nb(j);
        surfaceArea_[phId] += (microPhasePorosityInt_[stenb->getMicroPhaseId()] /
            d_convFactDbl2IntPor);
      }
    }
  }
  surfaceArea_[phId] *= oneFaceAreaPerHundredGramSolid_;

  return;
}

void Lattice::calcSurfaceArea(int phId) {
  // in m2 / (100g of solid)

  Site *stenb;
  int stId;

  vector<Isite> isite = interface_[phId].getDissolutionSites();
  int size = isite.size();
  double d_convFactDbl2IntPor = static_cast<double> (convFactDbl2IntPor_);

  for (int i = 0; i < size; i++) {
    stId = isite[i].getId();
    // for (int j = 0; j < NN_NNN; j++) {
    for (int j = 0; j < NUM_NEAREST_NEIGHBORS; j++) {
      stenb = site_[stId].nb(j); //ste->nb(j);
      surfaceArea_[phId] += (microPhasePorosityInt_[stenb->getMicroPhaseId()] /
          d_convFactDbl2IntPor);
    }
  }

  // surfaceArea_[phId] = (numTotEle + numTotCSHQ * porCSHQ) * oneFaceAreaPerHundredGramSolid_;
  surfaceArea_[phId] *= oneFaceAreaPerHundredGramSolid_;

  return;
}
