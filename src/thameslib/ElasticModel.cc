/**
@file ElasticModel.cc
@brief Method definitions for the ElasticModel base class

*/
#include "ElasticModel.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

ElasticModel::ElasticModel(int nx, int ny, int nz, int dim, ChemicalSystem *cs,
                           int npoints, const bool verbose, const bool warning)
    : chemSys_(cs) {
  ///
  /// Assign the dimensions of the finite element (FE) mesh
  ///

#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  nxy_ = nx_ * ny_;
  ns_ = nx_ * ny_ * nz_;

  if (verbose_) {
    std::clog << "ElasticModel::ElasticModel Constructor nx_ = " << nx_
              << " ny_ = " << ny_ << " nz_ = " << nz_ << " nxy_ = " << nxy_
              << " ns_ = " << ns_ << endl;
    std::clog.flush();
  }

  ///
  /// Initialize the prescribed stresses and strains
  ///

  strxx_ = stryy_ = strzz_ = 0.0;
  strxz_ = stryz_ = strxy_ = 0.0;
  sxx_ = syy_ = szz_ = 0.0;
  sxz_ = syz_ = sxy_ = 0.0;

  ///
  /// Initialize the stresses and strains (arrays) for each element
  ///

  elestress_.clear();
  elestrain_.clear();
  elestress_.resize(ns_);
  elestrain_.resize(ns_);

  for (int i = 0; i < ns_; i++) {
    elestress_[i].resize(6, 0.0);
    elestrain_[i].resize(6, 0.0);
  }

  ///
  /// Initialize the array of strain energies for each element
  ///

  strainengy_.clear();
  strainengy_.resize(ns_, 0.0);

  ///
  /// `nphase_` is the number of phases being considered in the problem.
  /// The values of `pix(m)` will run from 1 to `nphase_`.
  ///

  // nphase_ = nphase;
  nphase_ = chemSys_->getNumMicroPhases();
  numDCs_ = chemSys_->getNumDCs();
  waterDCId_ = chemSys_->getDCId(WaterDCName); // "H2O@"

  double molarVolume;
  for (int i = 0; i < numDCs_; i++) {
    molarVolume = chemSys_->getDCMolarVolume(i);
    convFactDCs_.push_back(molarVolume * 1.0e9);
  }

  vector<int> mPhDCcomp;
  mPhDCcomp.clear();
  mPhDCcomp_.push_back(mPhDCcomp); // void
  mPhDCcomp_.push_back(mPhDCcomp); // water
  for (int mPhId = FIRST_SOLID; mPhId < nphase_; mPhId++) {
    mPhDCcomp.clear();
    mPhDCcomp = chemSys_->getMicroPhaseDCMembers(mPhId);
    mPhDCcomp_.push_back(mPhDCcomp);
  }

  ///
  /// Establish the stopping criterion for convergence (proportional to
  /// the number of elements, `ns_`
  ///

  gtest_ = (1.0e-8) * ns_;
  gg_ = 0.0;

  ///
  /// The parameter `phasemod_[i][j]` is the `bulk[i][0]` and `shear[i][1]`
  /// moduli of the i'th phase. These can be input in terms of Young's moduli
  /// E[i][0] and Poisson's ratio nu[i][1]. The program, then changes them to
  /// bulk and shear moduli. For anisotropic elastic material, one can directly
  /// input the elastic moduli tensor cmod in subroutine femat, and skip this
  /// part.
  ///
  /// @warning If you wish to input in terms of bulk (0) and shear (1), then
  /// make sure to comment out the following loop.

  phasemod_.clear();
  phasemod_.resize(nphase_);
  for (int ijk = 0; ijk < nphase_; ijk++) {
    phasemod_[ijk].resize(2, 0.0);
  }

  ///
  /// Initialize elastic modulus variable `cmod_`
  ///

  cmod_.clear();
  cmod_.resize(nphase_);
  for (int ijk = 0; ijk < nphase_; ijk++) {
    cmod_[ijk].resize(6);
    for (int j = 0; j < 6; j++) {
      cmod_[ijk][j].resize(6, 0.0);
    }
  }

  ElasModul();

  ///
  /// Initialize the neighbor table for each element
  ///

  ib_.clear();
  ib_.resize(dim);
  for (int m = 0; m < dim; m++) {
    ib_[m].resize(27, 0);
  }

  BuildNeighbor();

  ///
  /// `npoints_` is the number of microstructres to use.
  ///

  npoints_ = npoints;

  ///
  /// Initialize `pix_[m]` and volume fraction of each phase `prob_[]`.
  ///

  pix_.clear();
  pix_.resize(ns_, 0);

  prob_.clear();
  prob_.resize(nphase_, 0.0);

  ///
  /// Initialize displacement at each node `u_`.
  ///

  u_.clear();
  u_.resize(dim);
  for (int m = 0; m < dim; m++) {
    u_[m].resize(3, 0.0);
  }

  ///
  /// Initialize :
  /// finite element stiffness matrices, `dk_`, the energy constant, `C_`,
  /// and the linear term coefficient vector, `b_`, required for
  /// computing the energy.

  // cmod_.clear();
  // cmod_.resize(nphase_);
  // for (int ijk = 0; ijk < nphase_; ijk++) {
  //   cmod_[ijk].resize(6);
  //   for (int j = 0; j < 6; j++) {
  //     cmod_[ijk][j].resize(6, 0.0);
  //   }
  // }

  dk_.clear();
  dk_.resize(nphase_);
  for (int ijk = 0; ijk < nphase_; ijk++) {
    dk_[ijk].resize(8);
    for (int i = 0; i < 8; i++) {
      dk_[ijk][i].resize(3);
      for (int k = 0; k < 3; k++) {
        dk_[ijk][i][k].resize(8);
        for (int j = 0; j < 8; j++) {
          dk_[ijk][i][k][j].resize(3, 0.0);
        }
      }
    }
  }

  initStiffness(); // dk_

  b_.clear();
  b_.resize(dim);
  for (int m = 0; m < dim; m++) {
    b_[m].resize(3, 0.0);
  }

  C_ = 0.0;

  ///
  /// Initialize the gradient at each element, `gb_[ns_][3]`,
  /// the auxiliary conjugate gradient vector variable at each element,
  /// `h_[ns_][3]`, and the other conjugate gradient vector variable at each
  /// element, `Ah_[ns_][3]`.
  ///

  gb_.clear();
  h_.clear();
  Ah_.clear();

  gb_.resize(dim);
  h_.resize(dim);
  Ah_.resize(dim);
  for (int m = 0; m < dim; m++) {
    gb_[m].resize(3, 0.0);
    h_[m].resize(3, 0.0);
    Ah_[m].resize(3, 0.0);
  }

  return;
}

void ElasticModel::BuildNeighbor() {

  ///
  /// First construct the 27 neighbor table in terms of delta i, delta j, and
  /// delta k information. (see Table 3 in manual)
  ///

  int in[27], jn[27], kn[27];
  in[0] = 0;
  in[1] = 1;
  in[2] = 1;
  in[3] = 1;
  in[4] = 0;
  in[5] = (-1);
  in[6] = (-1);
  in[7] = (-1);

  jn[0] = 1;
  jn[1] = 1;
  jn[2] = 0;
  jn[3] = (-1);
  jn[4] = (-1);
  jn[5] = (-1);
  jn[6] = 0;
  jn[7] = 1;

  for (int n = 0; n < 8; n++) {
    kn[n] = 0;
    kn[n + 8] = (-1);
    kn[n + 16] = 1;
    in[n + 8] = in[n];
    in[n + 16] = in[n];
    jn[n + 8] = jn[n];
    jn[n + 16] = jn[n];
  }

  in[24] = 0;
  in[25] = 0;
  in[26] = 0;
  jn[24] = 0;
  jn[25] = 0;
  jn[26] = 0;
  kn[24] = (-1);
  kn[25] = 1;
  kn[26] = 0;

  ///
  /// Now construct neighbor table according to 1D labels
  /// Matrix ib_[m][n] gives the 1-d label of the n'th neighbor (n=0,26) of the
  /// node labelled m.
  ///
  /// This set of for loops also manually checks for wrapping due to
  /// periodic boundary conditions in all three directions
  ///

  int m, m1;
  int i1, j1, k1;
  for (int k = 0; k < nz_; k++) {
    for (int j = 0; j < ny_; j++) {
      for (int i = 0; i < nx_; i++) {
        m = nxy_ * k + nx_ * j + i;
        for (int n = 0; n < 27; n++) {
          i1 = i + in[n];
          j1 = j + jn[n];
          k1 = k + kn[n];
          if (i1 < 0)
            i1 += nx_;
          else if (i1 >= nx_)
            i1 -= nx_;
          if (j1 < 0)
            j1 += ny_;
          else if (j1 >= ny_)
            j1 -= ny_;
          if (k1 < 0)
            k1 += nz_;
          else if (k1 >= nz_)
            k1 -= nz_;
          m1 = nxy_ * k1 + nx_ * j1 + i1;
          ib_[m][n] = m1;
        }
      }
    }
  }

  return;
}

// void ElasticModel::ElasModul(string phasemod_fileName, int nphase) {
void ElasticModel::ElasModul(void) {

  ///
  /// Open the file input stream with the specified name.
  /// This file holds the values of the Young's modulus [GPa] and
  /// Poisson's ratio of each phase
  ///

  // ifstream in(phasemod_fileName.c_str());
  // if (!in) {
  //   std::clog << "ElasticModel::ElasModul can't open the file: " <<
  //   phasemod_fileName
  //        << endl;
  //  std::clog.flush();
  //  exit(1);
  // } else {
  string mPhName;
  for (int i = 0; i < nphase_; i++) {
    // in >> phaseid;
    // in >> buff;
    // in >> elsmodul;
    // phasemod_[phaseid][0] = elsmodul;
    // in >> elsmodul;
    // phasemod_[phaseid][1] = elsmodul;
    mPhName = chemSys_->getMicroPhaseName(i);
    // bulk modulus K
    phasemod_[i][0] = chemSys_->getElasticModuliComp(mPhName).K;
    // shear modulus G
    phasemod_[i][1] = chemSys_->getElasticModuliComp(mPhName).G;
  }
  // }

  ///
  /// The program uses bulk modulus (0) and shear modulus (1), so transform
  /// Young's modulus (E) and Poisoon's ratio (n) to them.
  /// (k,G) -> (E,n) : E = 9KG/(3K + G)  &  n = (3K - 2G)/(2(3K + G))
  /// (E,n) -> (k,G) : K = E/(3(1 - 2n)) &  G = E/(2(1 + n))
  ///

  // for (int i = 0; i < nphase_; i++) {
  //   double save = phasemod_[i][0];
  //   phasemod_[i][0] = (phasemod_[i][0] / 3.0) / (1 - 2 * phasemod_[i][1]);
  //   phasemod_[i][1] = (save / 2.0) / (1.0 + phasemod_[i][1]);
  // }

  ///
  /// Set up the elastic modulus tensor for each phase, which uses
  /// engineering notation for the components as described in the
  /// documentation for the `cmod_` member.
  ///
  /// The ck and cmu matrices are used to multiply by the (scalar)
  /// bulk and shear moduli values, `phasemod_[0]` and `phasemod_[1]`,
  /// respectively.
  ///

  double ck[6][6], cmu[6][6];
  ck[0][0] = 1.0;
  ck[0][1] = 1.0;
  ck[0][2] = 1.0;
  ck[0][3] = 0.0;
  ck[0][4] = 0.0;
  ck[0][5] = 0.0;
  ck[1][0] = 1.0;
  ck[1][1] = 1.0;
  ck[1][2] = 1.0;
  ck[1][3] = 0.0;
  ck[1][4] = 0.0;
  ck[1][5] = 0.0;
  ck[2][0] = 1.0;
  ck[2][1] = 1.0;
  ck[2][2] = 1.0;
  ck[2][3] = 0.0;
  ck[2][4] = 0.0;
  ck[2][5] = 0.0;
  ck[3][0] = 0.0;
  ck[3][1] = 0.0;
  ck[3][2] = 0.0;
  ck[3][3] = 0.0;
  ck[3][4] = 0.0;
  ck[3][5] = 0.0;
  ck[4][0] = 0.0;
  ck[4][1] = 0.0;
  ck[4][2] = 0.0;
  ck[4][3] = 0.0;
  ck[4][4] = 0.0;
  ck[4][5] = 0.0;
  ck[5][0] = 0.0;
  ck[5][1] = 0.0;
  ck[5][2] = 0.0;
  ck[5][3] = 0.0;
  ck[5][4] = 0.0;
  ck[5][5] = 0.0;

  cmu[0][0] = 4.0 / 3.0;
  cmu[0][1] = (-2.0) / 3.0;
  cmu[0][2] = (-2.0) / 3.0;
  cmu[0][3] = 0.0;
  cmu[0][4] = 0.0;
  cmu[0][5] = 0.0;
  cmu[1][0] = (-2.0) / 3.0;
  cmu[1][1] = 4.0 / 3.0;
  cmu[1][2] = (-2.0) / 3.0;
  cmu[1][3] = 0.0;
  cmu[1][4] = 0.0;
  cmu[1][5] = 0.0;
  cmu[2][0] = (-2.0) / 3.0;
  cmu[2][1] = (-2.0) / 3.0;
  cmu[2][2] = 4.0 / 3.0;
  cmu[2][3] = 0.0;
  cmu[2][4] = 0.0;
  cmu[2][5] = 0.0;
  cmu[3][0] = 0.0;
  cmu[3][1] = 0.0;
  cmu[3][2] = 0.0;
  cmu[3][3] = 1.0;
  cmu[3][4] = 0.0;
  cmu[3][5] = 0.0;
  cmu[4][0] = 0.0;
  cmu[4][1] = 0.0;
  cmu[4][2] = 0.0;
  cmu[4][3] = 0.0;
  cmu[4][4] = 1.0;
  cmu[4][5] = 0.0;
  cmu[5][0] = 0.0;
  cmu[5][1] = 0.0;
  cmu[5][2] = 0.0;
  cmu[5][3] = 0.0;
  cmu[5][4] = 0.0;
  cmu[5][5] = 1.0;

  ///
  /// Construct `cmod_` tensor for each phase by matrix multiplication
  ///
  /// \f{equation}
  ///     c_{ij} = K c^k_{ij} + \mu c^{\mu}_{ij}
  /// \f}
  ///

  for (int ijk = 0; ijk < nphase_; ijk++) {
    for (int j = 0; j < 6; j++) {
      for (int i = 0; i < 6; i++) {
        cmod_[ijk][i][j] =
            phasemod_[ijk][0] * ck[i][j] + phasemod_[ijk][1] * cmu[i][j];
      }
    }
  }
  return;
}

void ElasticModel::initStiffness(void) {
  double dndx[8], dndy[8], dndz[8];
  double g[3][3][3];
  double es[6][8][3];
  double x, y, z;

  ///
  /// (User) NOTE: complete elastic modulus matrix is used, so an anisotropic
  /// matrix could be directly input at any point, since program is written to
  /// use a general elastic moduli tensor, but is only explicitly implemented
  /// for isotropic materials.
  ///

  ///
  /// Initialize stiffness matrices dk_
  ///

  ///
  /// Set up elastic moduli matrices for each kind of element
  ///

  ///
  /// Set up Simpson's integration rule weight vector
  ///

  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        int nm = 0;
        if (i == 1)
          nm += 1;
        if (j == 1)
          nm += 1;
        if (k == 1)
          nm += 1;
        g[i][j][k] = pow(4.0, static_cast<double>(nm));
      }
    }
  }

  ///
  /// Loop over the nphase kinds of pixels and Simpson's rule quadrature points
  /// in order to compute the stiffness matrices. Stiffness matrices of
  /// trilinear finite elements are quadratic in x, y, and z, so that Simpson's
  /// rule quadrature gives exact results.
  ///

  for (int ijk = 0; ijk < nphase_; ijk++) {
    for (int k = 0; k < 3; k++) {
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
          x = i / 2.0;
          y = j / 2.0;
          z = k / 2.0;

          ///
          /// dndx means the negative derivative, with respect to x,
          /// of the shape matrix
          /// N (see Manual, Sec. 2.2), dndy, and dndz are similar.
          ///

          dndx[0] = (-(1.0 - y)) * (1.0 - z);
          dndx[1] = (1.0 - y) * (1.0 - z);
          dndx[2] = y * (1.0 - z);
          dndx[3] = (-y) * (1.0 - z);
          dndx[4] = (-(1.0 - y)) * z;
          dndx[5] = (1.0 - y) * z;
          dndx[6] = y * z;
          dndx[7] = (-y) * z;
          dndy[0] = (-(1.0 - x)) * (1.0 - z);
          dndy[1] = (-x) * (1.0 - z);
          dndy[2] = x * (1.0 - z);
          dndy[3] = (1.0 - x) * (1.0 - z);
          dndy[4] = (-(1.0 - x)) * z;
          dndy[5] = (-x) * z;
          dndy[6] = x * z;
          dndy[7] = (1.0 - x) * z;
          dndz[0] = (-(1.0 - x)) * (1.0 - y);
          dndz[1] = (-x) * (1.0 - y);
          dndz[2] = (-x) * y;
          dndz[3] = (-(1.0 - x)) * y;
          dndz[4] = (1.0 - x) * (1.0 - y);
          dndz[5] = x * (1.0 - y);
          dndz[6] = x * y;
          dndz[7] = (1.0 - x) * y;

          ///
          /// Now build strain matrix
          ///

          for (int n1 = 0; n1 < 6; n1++) {
            for (int n2 = 0; n2 < 8; n2++) {
              for (int n3 = 0; n3 < 3; n3++) {
                es[n1][n2][n3] = 0.0;
              }
            }
          }

          for (int n = 0; n < 8; n++) {
            es[0][n][0] = dndx[n];
            es[1][n][1] = dndy[n];
            es[2][n][2] = dndz[n];
            es[3][n][0] = dndz[n];
            es[3][n][2] = dndx[n];
            es[4][n][1] = dndz[n];
            es[4][n][2] = dndy[n];
            es[5][n][0] = dndy[n];
            es[5][n][1] = dndx[n];
          }

          ///
          /// Matrix multiply to determine value at (x,y,z), multiply by proper
          /// weight, and sum into dk, the stiffness matrix.
          ///

          for (int mm = 0; mm < 3; mm++) {
            for (int nn = 0; nn < 3; nn++) {
              for (int ii = 0; ii < 8; ii++) {
                for (int jj = 0; jj < 8; jj++) {
                  double sum = 0.0;
                  for (int kk = 0; kk < 6; kk++) {
                    for (int ll = 0; ll < 6; ll++) {
                      sum +=
                          es[kk][ii][mm] * cmod_[ijk][kk][ll] * es[ll][jj][nn];
                    }
                  }
                  dk_[ijk][ii][mm][jj][nn] += g[i][j][k] * sum / 216.0;
                }
              }
            }
          }
        }
      }
    }
  }
}

// void ElasticModel::ppixel(string fileName, int nphase) {
void ElasticModel::ppixel(string fileName) {
  ///
  /// If you want to set up a test image inside the program, instead of
  /// reading it in from a file, this should be done inside this method
  ///

  string buff, version;
  double resolution;
  int m;

  ///
  /// Open and read the input file stream with the microstructure data
  ///

  ifstream in(fileName.c_str());

  if (!in) {

    std::clog << "can't open the file: " << fileName << endl;
    cerr << "can't open the file: " << fileName << endl;
    exit(1);

  } else {

    /*
    in >> buff;
    if (buff == VERSIONSTRING) {
      in >> version;
      in >> buff;
      if (buff == XSIZESTRING) {
        in >> nx_;
        in >> buff;
        in >> ny_;
        in >> buff;
        in >> nz_;
      }
      in >> buff;
      if (buff == IMGRESSTRING) {
        in >> resolution;
      }
    } else {

      ///
      /// Image file does not have the expected header.
      /// Assume default size to 100 and resolution to 1.0
      /// micrometers.
      ///

      version = "2.0";
      resolution = 1.0;
      nx_ = ny_ = nz_ = 100;
    }
    */
    in >> buff; // VERSIONSTRING
    in >> version;
    in >> buff; // XSIZESTRING
    in >> nx_;
    in >> buff; // YSIZESTRING
    in >> ny_;
    in >> buff; // ZSIZESTRING
    in >> nz_;
    in >> buff; // IMGRESSTRING
    in >> resolution;

    ///
    /// Each line of the microstructure file contains the phase id
    /// to assign, and the microstructure file must be written in
    /// the same order as the finite elements are populated
    ///

    // ns_ = nx_ * ny_ * nz_;
    // nxy_ = nx_ * ny_;
    for (int k = 0; k < nz_; k++) {
      for (int j = 0; j < ny_; j++) {
        for (int i = 0; i < nx_; i++) {
          m = nxy_ * k + nx_ * j + i;
          in >> pix_[m];
        }
      }
    }

    in.close();

    // cout << endl << "ini ppixel vector (6):" << endl;
    // for (int i = 0; i < ns_; i++) {
    //   cout << pix_[i] << endl;
    // }
    // cout << "end ppixel vector (6):" << endl;
    // exit(0);

    ///
    /// Check for wrong phase labels--less than 1 or greater than nphase_.
    /// Note that nothing is done about the error; it is just reported
    /// to standard out.
    ///
    /// @todo See if it makes sense to do value checking as exception handling
    ///

    // for (int m = 0; m < ns_; m++) {
    //   if (pix_[m] < 0) {
    //     std::clog << "Phase label in pix < 0 --- error at " << m << endl;
    //   } else if (pix_[m] >= nphase_) {
    //     std::clog << "Phase label in pix >= nphase_ --- error at " << m <<
    //     endl;
    //   }
    // }
  }

  return;
}

// void ElasticModel::ppixel(string fileName, int nphase) {
void ElasticModel::ppixel(vector<int> vectPhId) {

  ///
  /// Each line of the microstructure file contains the phase id
  /// to assign, and the microstructure file must be written in
  /// the same order as the finite elements are populated
  ///

  int m;
  int kji = -1;
  // ns_ = nx_ * ny_ * nz_;
  // nxy_ = nx_ * ny_;
  for (int k = 0; k < nz_; k++) {
    for (int j = 0; j < ny_; j++) {
      for (int i = 0; i < nx_; i++) {
        m = nxy_ * k + nx_ * j + i;
        kji++;
        pix_[m] = vectPhId[kji];
      }
    }
  }

  return;
}

void ElasticModel::ppixel(vector<int> *p_vectPhId) {

  ///
  /// Each line of the microstructure file contains the phase id
  /// to assign, and the microstructure file must be written in
  /// the same order as the finite elements are populated
  ///

  for (int i = 0; i < ns_; i++) {
    pix_[i] = p_vectPhId->at(i);
  }

  return;
}

void ElasticModel::getAvgStrainengy() {
  ///
  /// Clear and initialize the array that will hold the
  /// volume averaged strain energy of each phase in the microstructure
  ///

  avgStrainengy_.clear();
  avgStrainengy_.resize(nphase_, 0.0);

  ///
  /// Sum the strain energy in each phase by one sweep through the mesh
  ///

  for (int m = 0; m < ns_; m++) {
    avgStrainengy_[pix_[m]] += strainengy_[m];
  }

  ///
  /// Normalize by the volume of the phase
  ///

  for (int i = 0; i < nphase_; i++) {
    avgStrainengy_[i] = avgStrainengy_[i] / (prob_[i] * ns_);
  }

  ///
  /// The variable `strainenergy` (not `strainenergy_`) is defined
  /// in the StrainEnergy.h header file.  It is used in the
  /// GEM-IPM node class to facilitate communication of phase
  /// strain energies back and forth between the two, because
  /// strain energy adds to the Gibbs energy of a phase and therefore
  /// should influence equilibrium calculations.
  ///
  /// The block below relates the COMPUTED average strain energy density,
  /// `avgStrainengy_`, with units of GJ/m3, to the strain
  /// energy PER MOLE in the corresponding GEM dependent components that
  /// make up that phase, `strainenergy`, in units of J/mol.
  /// To make the conversion, we must multiply avgStrainengy_ by the molar
  /// volume (units of m3/mol) and then multiply by e9 to convert GJ to J.
  ///
  /// To take the example of C3S, it is a microstructure phase (id = 2),
  /// with one associated GEM DC (also called C3S, with GEM id 115 and
  /// a molar volume of 7.318e-5 m3/mol.  Multiplying that molar volume by
  /// e9 to convert from GJ to J gives an overall converstion factor of
  /// 7.318e4 J m3 / GJ mol.
  ///
  /// @note This is only true for one particular incarnation of GEMS
  /// and one particular set of ICs, DCs, etc.
  ///
  /// @remarks It would be better if these index values
  /// and molar volumes were not hard-coded like this.
  ///
  /// @todo Generalize the indices and molar volumes
  ///

  // strainenergy.clear();
  // strainenergy.resize(156, 0.0);
  strainenergy.clear();
  strainenergy.resize(numDCs_, 0.0);

  // double molarVolume;
  // vector<double> convFactDCs;
  // for (int i = 0; i < numDCs_; i++) {
  //   molarVolume = chemSys_->getDCMolarVolume(i);
  //   convFactDCs_.push_back(molarVolume * 1.e9);
  // }

  strainenergy[waterDCId_] =
      avgStrainengy_[ELECTROLYTEID] * convFactDCs_[waterDCId_];

  int size;
  int DCId;
  for (int mPhId = FIRST_SOLID; mPhId < nphase_; mPhId++) {
    // mPhDCcomp = chemSys_->getMicroPhaseDCMembers(mPhId);
    size = mPhDCcomp_[mPhId].size();
    for (int j = 0; j < size; j++) {
      DCId = mPhDCcomp_[mPhId][j];
      strainenergy[DCId] = avgStrainengy_[mPhId] * convFactDCs_[DCId];
    }
  }

  /*
  ///
  /// DC component making up the H2O microstructure phase
  ///

  strainenergy[70] = avgStrainengy_[1] * (1.8068E4); // Water

    ...

  strainenergy[120] = avgStrainengy_[18] * (1.93985E5); // CAH10

  ///
  /// DC component making up the LIME microstructure phase
  ///

  strainenergy[131] = avgStrainengy_[19] * (1.6764E4); // Free lime
  */

  return;
}

void ElasticModel::writeStress(string &root, double time, int index) {
  if (index >= 0 && index < 6) {
    double min, max;
    min = max = 0.0;

    ///
    /// Create and initialize the local rgb vector
    ///

    vector<int> color;
    color.clear();
    color.resize(3, 0);

    ///
    /// Specify the file name and open the output stream
    ///

    ostringstream ostr;
    ostr << static_cast<int>(time * 60.0);
    string timestr(ostr.str());
    string ofileName(root);
    string ofpngname(root);
    if (index == 0) {
      ofileName = ofileName + "." + "stress-xx." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "stress-xx." + timestr + ".png";
    } else if (index == 1) {
      ofileName = ofileName + "." + "stress-yy." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "stress-yy." + timestr + ".png";
    } else if (index == 2) {
      ofileName = ofileName + "." + "stress-zz." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "stress-zz." + timestr + ".png";
    } else if (index == 3) {
      ofileName = ofileName + "." + "stress-xz." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "stress-xz." + timestr + ".png";
    } else if (index == 4) {
      ofileName = ofileName + "." + "stress-yz." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "stress-yz." + timestr + ".png";
    } else if (index == 5) {
      ofileName = ofileName + "." + "stress-xy." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "stress-xy." + timestr + ".png";
    }
    ofstream out(ofileName.c_str());
    if (!out.is_open()) {
      std::clog << ofileName << "could not open." << endl;
      exit(1);
    }

    ///
    /// Write the PPM header for a full color image
    /// Currently only uses cyan of different intensities
    ///

    out << "P3" << endl;
    out << nx_ << " " << nz_ << endl;
    out << "255" << endl;

    int slice = nx_ / 2;
    int m;
    for (int j = 0; j < ny_; j++) {
      for (int k = 0; k < nz_; k++) {
        m = nxy_ * k + nx_ * j + slice;
        if (min > elestress_[m][index])
          min = elestress_[m][index];
        if (max < elestress_[m][index])
          max = elestress_[m][index];
      }
    }

    if (verbose_) {
      std::clog << "ElasticModel::writeStress minimum stress-" << index
                << " is: " << min << endl;
      std::clog << "ElasticModel::writeStress maximum stress-" << index
                << " is: " << max << endl;
      std::clog.flush();
    }

    for (int k = 0; k < nz_; k++) {
      for (int j = 0; j < nz_; j++) {
        m = nxy_ * k + nx_ * j + slice;
        color[1] = static_cast<int>(
            ((elestress_[m][index] - min) / (max - min)) * 255);
        color[2] = static_cast<int>(
            ((elestress_[m][index] - min) / (max - min)) * 255);
        out << color[0] << " " << color[1] << " " << color[2] << endl;
      }
    }

    out.close();

    ///
    /// PPM file is finished.  Now convert to PNG using ImageMagick convert
    /// command via a system call (not recommended).
    ///

    // string buff = "convert " + ofileName + " " + ofpngname;
    string buff = ConvertCommand + " " + ofileName + " " + ofpngname;
    int resCallSystem = system(buff.c_str());
    if (resCallSystem == -1) {
      // handle the error;
      std::clog
          << endl
          << endl
          << "    ElasticModel.cc - error in writeStress() : resCallSystem = -1"
          << endl;
      std::clog << endl << "    STOP program" << endl;
      // throw HandleException ("writeStress", "ElasticModel.cc",
      //                "system(buff.c_str())", "resCallSystem = -1");
      exit(1);
    }

    return;

  } else {

    std::clog << "index out of range. should be between 0 and 6." << endl;
    exit(1);
  }
}

void ElasticModel::writeStrain(string &root, double time, int index) {
  if (index >= 0 && index < 6) {

    double min, max;
    min = max = 0.0;

    ///
    /// Create and initialize the local rgb vector
    ///

    vector<int> color;
    color.clear();
    color.resize(3, 0);

    ///
    /// Specify the file name and open the output stream
    ///

    ostringstream ostr;
    ostr << static_cast<int>(time * 60.0);
    string timestr(ostr.str());
    string ofileName(root);
    string ofpngname(root);
    if (index == 0) {
      ofileName = ofileName + "." + "strain-xx." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "strain-xx." + timestr + ".png";
    } else if (index == 1) {
      ofileName = ofileName + "." + "strain-yy." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "strain-yy." + timestr + ".png";
    } else if (index == 2) {
      ofileName = ofileName + "." + "strain-zz." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "strain-zz." + timestr + ".png";
    } else if (index == 3) {
      ofileName = ofileName + "." + "strain-xz." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "strain-xz." + timestr + ".png";
    } else if (index == 4) {
      ofileName = ofileName + "." + "strain-yz." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "strain-yz." + timestr + ".png";
    } else if (index == 5) {
      ofileName = ofileName + "." + "strain-xy." + timestr + ".ppm";
      ofpngname = ofpngname + "." + "strain-xy." + timestr + ".png";
    }
    ofstream out(ofileName.c_str());
    if (!out.is_open()) {
      std::clog << ofileName << "could not open." << endl;
      exit(1);
    }

    ///
    /// Write the PPM header for a full color image
    /// Currently only uses cyan of different intensities
    ///

    out << "P3" << endl;
    out << nx_ << " " << nz_ << endl;
    out << "255" << endl;

    int slice = nx_ / 2;
    int m;
    for (int j = 0; j < ny_; j++) {
      for (int k = 0; k < nz_; k++) {
        m = nxy_ * k + nx_ * j + slice;
        if (min > elestrain_[m][index])
          min = elestrain_[m][index];
        if (max < elestrain_[m][index])
          max = elestrain_[m][index];
      }
    }
    for (int k = 0; k < nz_; k++) {
      for (int j = 0; j < nz_; j++) {
        m = nxy_ * k + nx_ * j + slice;
        color[1] = static_cast<int>(
            ((elestrain_[m][index] - min) / (max - min)) * 255);
        color[2] = static_cast<int>(
            ((elestrain_[m][index] - min) / (max - min)) * 255);
        out << color[0] << " " << color[1] << " " << color[2] << endl;
      }
    }

    out.close();

    ///
    /// PPM file is finished.  Now convert to PNG using ImageMagick convert
    /// command via a system call (not recommended).
    ///

    // string buff = "convert " + ofileName + " " + ofpngname;
    string buff = ConvertCommand + " " + ofileName + " " + ofpngname;
    int resCallSystem = system(buff.c_str());
    if (resCallSystem == -1) {
      // handle the error;
      std::clog
          << endl
          << endl
          << "    ElasticModel.cc - error in writeStrain() : resCallSystem = -1"
          << endl;
      std::clog << endl << "    STOP program" << endl;
      // throw HandleException ("writeStrain", "ElasticModel.cc",
      //                "system(buff.c_str())", "resCallSystem = -1");
      exit(1);
    }
    return;

  } else {

    std::clog << "index out of range. should be between 0 to 5." << endl;
    exit(1);
  }
}

void ElasticModel::writeDisp(string &root, string timeString) {

#ifdef DEBUG
  std::clog << "ElasticModel::writeDisp" << endl;
  std::clog.flush();
#endif

  ///
  /// Specify the file name and open the output stream
  ///

  // ostringstream ostr;
  // ostr << static_cast<int>(time * 60.0);
  // string timestr(ostr.str());
  string ofileName(root);
  // ofileName = ofileName + "." + "disp." + timestr + ".dat";
  ofileName = ofileName + "." + timeString + ".disp.dat";
  ofstream out(ofileName.c_str());

  if (!out.is_open()) {
    std::clog << ofileName << "could not open." << endl;
    exit(1);
  }

  // int m;
  // for (int k = 0; k < nz_; k++) {
  //   for (int j = 0; j < ny_; j++) {
  //     for (int i = 0; i < nx_; i++) {
  //       m = nxy_ * k + nx_ * j + i;
  //       out << u_[m][0] << "    " << u_[m][1] << "    " << u_[m][2] << endl;
  //     }
  //   }
  // }

  for (int m = 0; m < ns_; m++) {
    out << u_[m][0] << "    " << u_[m][1] << "    " << u_[m][2] << endl;
  }

  out.close();

  return;
}

void ElasticModel::writeStrainEngy(string &root, double time) {
#ifdef DEBUG
  std::clog << "ElasticModel::writeStrainEngy" << endl;
  std::clog.flush();
#endif

  double min, max;
  min = max = 0.0;

  ///
  /// Create and initialize the local rgb vector
  ///

  vector<int> color;
  color.clear();
  color.resize(3, 0);

  ///
  /// Specify the file name and open the output stream
  ///

  ostringstream ostr;
  ostr << static_cast<int>(time * 60.0);
  string timestr(ostr.str());
  string ofileName(root);
  string ofpngname(root);
  ofileName = ofileName + "." + "strainengy." + timestr + ".ppm";
  ofpngname = ofpngname + "." + "strainengy." + timestr + ".png";
  ofstream out(ofileName.c_str());
  if (!out.is_open()) {
    std::clog << ofileName << "could not open." << endl;
    exit(1);
  }

  ///
  /// Write the PPM header for a full color image
  /// Currently only uses cyan of different intensities
  ///

  out << "P3" << endl;
  out << nx_ << " " << nz_ << endl;
  out << "255" << endl;

  int m;
  int slice = nx_ / 2;
  for (int j = 0; j < ny_; j++) {
    for (int k = 0; k < nz_; k++) {
      m = nxy_ * k + nx_ * j + slice;
      if (min > strainengy_[m])
        min = strainengy_[m];
      if (max < strainengy_[m])
        max = strainengy_[m];
    }
  }

  for (int k = 0; k < nz_; k++) {
    for (int j = 0; j < nz_; j++) {
      m = nxy_ * k + nx_ * j + slice;
      color[1] = static_cast<int>(((strainengy_[m] - min) / (max - min)) * 255);
      color[2] = static_cast<int>(((strainengy_[m] - min) / (max - min)) * 255);
      out << color[0] << " " << color[1] << " " << color[2] << endl;
    }
  }

  out.close();

  ///
  /// PPM file is finished.  Now convert to PNG using ImageMagick convert
  /// command via a system call (not recommended).
  ///

  // string buff = "convert " + ofileName + " " + ofpngname;
  string buff = ConvertCommand + " " + ofileName + " " + ofpngname;
  int resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    // handle the error;
    std::clog
        << endl
        << endl
        << "    ElasticModel.cc - error in writeStrainEngy() : resCallSystem = "
           "-1"
        << endl;
    std::clog << endl << "    STOP program" << endl;
    // throw HandleException ("writeStrainEngy", "ElasticModel.cc",
    //                "system(buff.c_str())", "resCallSystem = -1");
    exit(1);
  }
  return;
}

void ElasticModel::assig(void) {
  for (int i = 0; i < nphase_; i++) {
    prob_[i] = 0.0;
  }

  for (int m = 0; m < ns_; m++) {
    prob_[pix_[m]] += 1;
  }

  double ns_dbl = ns_;
  for (int i = 0; i < nphase_; i++) {
    prob_[i] = prob_[i] / ns_dbl;
  }

  return;
}
