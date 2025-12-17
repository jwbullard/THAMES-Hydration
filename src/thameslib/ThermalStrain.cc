/**
@file ThermalStrain.cc
@brief Defines methods for the ThermalStrain class.

*/
#include "ThermalStrain.h"

using std::cout;
using std::endl;
using std::map;
using std::string;
using std::vector;

ThermalStrain::ThermalStrain(int nx, int ny, int nz, int dim,
                             ChemicalSystem *cs, const bool hasAggregateSlab,
                             const bool verbose, const bool warning)
    : ElasticModel(nx, ny, nz, dim, cs, hasAggregateSlab, verbose, warning) {

#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  if (verbose_) {
    std::clog << "ThermalStrain::ThermalStrain" << endl;
    std::clog.flush();
  }

  kmax_ = 80; // 3;
  std::clog
      << endl
      << "ThermalStrain::ThermalStrain - "
         "the number of relaxation steps for elastic computation :  kmax_ = "
      << kmax_ << endl;

  isFirst_ = true;

  Y_ = 0.0;

  localgg_ = 0.0;

  boxsize_ = 27;
  boxnum_ = static_cast<int>(
      pow(static_cast<double>(static_cast<int>(boxsize_ / 2) * 2 + 1), 3.0));
  localgtest_ = 1.0e-30 * boxnum_;

  if (verbose_) {
    std::clog << "ThermslStrain::ThermalStrain localgtest_ is: " << localgtest_
              << endl;
    std::clog.flush();
  }

  eigen_.clear();
  eigen_.resize(ns_);
  for (int n = 0; n < ns_; n++) {
    eigen_[n].resize(6, 0.0);
  }

  b0_.clear();
  b1_.clear();
  b2_.clear();
  b3_.clear();
  b4_.clear();
  b5_.clear();
  b0_.resize(dim);
  b1_.resize(dim);
  b2_.resize(dim);
  b3_.resize(dim);
  b4_.resize(dim);
  b5_.resize(dim);
  for (int n = 0; n < dim; n++) {
    b0_[n].resize(3, 0.0);
    b1_[n].resize(3, 0.0);
    b2_[n].resize(3, 0.0);
    b3_[n].resize(3, 0.0);
    b4_[n].resize(3, 0.0);
    b5_[n].resize(3, 0.0);
  }

  exp_.clear();

  T_.clear();
  T_.resize(dim);
  for (int n = 0; n < dim; n++) {
    T_[n].resize(3, 0.0);
  }

  ///
  /// Indexes for CSH, ettringite, and hydrotalcite are hard-wired into
  /// the code.
  ///
  /// @todo Generalize this so that is not fragile to changes in the indexing
  /// of phases.
  ///

  // write a new method!

  tstrength_.clear();
  tstrength_.resize(nphase_, 5.0); // [MPa]

  // tstrength_[12] = 1.0; // CSH
  // tstrength_[15] = 1.0; // ETTR
  // tstrength_[17] = 1.0; // HYDROTALC
  int mPhId;
  mPhId = chemSys_->getMicroPhaseId_SA("CSHQ");
  if (mPhId != -1) {
    tstrength_[mPhId] = 1.0; // CSHQ [MPa]
  } else {
    std::clog << "chemSys_->getMicroPhaseId_SA(\"CSHQ\") = -1 !!!" << endl;
  }

  mPhId = chemSys_->getMicroPhaseId_SA("AFt");
  if (mPhId != -1) {
    tstrength_[mPhId] = 1.0; // AFt [MPa]
  } else {
    std::clog << "chemSys_->getMicroPhaseId_SA(\"AFt\") = -1 !!!" << endl;
  }

  mPhId = chemSys_->getMicroPhaseId_SA("Hydrotalcite");
  if (mPhId != -1) {
    tstrength_[mPhId] = 1.0; // Hydrotalcite [MPa]
  } else {
    std::clog << "chemSys_->getMicroPhaseId_SA(\"Hydrotalcite) = -1 !!!"
              << endl;
  }

  mPhId = chemSys_->getMicroPhaseId_SA("Alite");
  if (mPhId != -1) {
    tstrength_[mPhId] = 10.0; // C3S - Alite [MPa]
  } else {
    std::clog << "chemSys_->getMicroPhaseId_SA(\"Alite\") = -1 !!!" << endl;
  }

  mPhId = chemSys_->getMicroPhaseId_SA("Belite");
  if (mPhId != -1) {
    tstrength_[mPhId] = 10.0; // C2S - Belite [MPa]
  } else {
    std::clog << "chemSys_->getMicroPhaseId_SA(\"Belite\") = -1 !!!" << endl;
  }

  mPhId = chemSys_->getMicroPhaseId_SA("Aluminate");
  if (mPhId != -1) {
    tstrength_[mPhId] = 10.0; // C3A - Aluminate [MPa]
  } else {
    std::clog << "chemSys_->getMicroPhaseId_SA(\"Aluminate\") = -1 !!!" << endl;
  }

  mPhId = chemSys_->getMicroPhaseId_SA("Ferrite");
  if (mPhId != -1) {
    tstrength_[mPhId] = 10.0; // C4AF - Ferrite [MPa]
  } else {
    std::clog << "chemSys_->getMicroPhaseId_SA(\"Ferrite\") = -1 !!!" << endl;
  }

  zcon_.clear();
  zcon_.resize(2);
  for (int i = 0; i < 2; i++) {
    zcon_[i].resize(3);
    for (int k = 0; k < 3; k++) {
      zcon_[i][k].resize(2);
      for (int j = 0; j < 2; j++) {
        zcon_[i][k][j].resize(3, 0.0);
      }
    }
  }

  ss_.clear();
  ss_.resize(ns_);
  for (int i = 0; i < ns_; i++) {
    ss_[i].resize(8);
    for (int j = 0; j < 8; j++) {
      ss_[i][j].resize(3, 0.0);
    }
  }
}

// void ThermalStrain::femat(int nx, int ny, int nz, int ns, int nphase,
//                           int iskip) {
void ThermalStrain::femat(int iskip) {
  double dndx[8], dndy[8], dndz[8];
  // double g[3][3][3];
  double es[6][8][3], delta[8][3];
  int is[8];
  // double x, y, z;
  int nx1 = nx_ - 1;
  int ny1 = ny_ - 1;
  int nz1 = nz_ - 1;
  double nx_dbl = nx_;
  double ny_dbl = ny_;
  double nz_dbl = nz_;
  int m;
  double sum = 0;
  double exx, eyy, ezz, exz, eyz, exy;

  /*
  ///
  /// Generate dk, zcon, T and Y on first pass. After that they are constant,
  /// since they are independent of the macrostrains. Only b gets upgraded
  /// as the macrostrains change.
  ///

  if (iskip == 0) {

    ///
    /// Initialize stiffness matrices.
    ///

    for (int m = 0; m < nphase_; m++) {
      for (int l = 0; l < 3; l++) {
        for (int k = 0; k < 3; k++) {
          for (int j = 0; j < 8; j++) {
            for (int i = 0; i < 8; i++) {
              dk_[m][i][k][j][l] = 0.0;
            }
          }
        }
      }
    }

    ///
    /// Initialize zcon matrix (gives C term for arbitrary macrostrains)
    ///

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        for (int mi = 0; mi < 3; mi++) {
          for (int mj = 0; mj < 3; mj++) {
            zcon_[i][mi][j][mj] = 0.0;
          }
        }
      }
    }

    ///
    /// @note An anisotropic elastic moduli tensor could be input at this point,
    /// bypassing this part, which assumes isotropic elasticity, so that there
    /// are only two independent numbers making up the elastic moduli tensor,
    /// the bulk modulus K and the shear modulus G.
    ///

    // ElasModul(phasemod_fileName_, nphase_);
    ElasModul();

    ///
    /// Set up Simpson's integration rule weight vector
    ///

    for (int k = 0; k < 3; k++) {
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
          int nm = 0;
          if (i == 1)
            nm++;
          if (j == 1)
            nm++;
          if (k == 1)
            nm++;
          g[i][j][k] = pow(4.0, static_cast<double>(nm));
        }
      }
    }

    ///
    /// Loop over the nphase kinds of pixels and Simpson's rule quadrature
    /// points in order to compute the stiffness matrices. Stiffness matrices of
    /// trilinear finite elements are quadratic in x, y, and z, so that
    /// Simpson's rule quadrature gives exact results.
    ///

    for (int ijk = 0; ijk < nphase_; ijk++) {
      for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            x = i / 2.0;
            y = j / 2.0;
            z = k / 2.0;

            ///
            /// dndx means the negative derivative with respect to x, of the
            /// shape matrix N (see manual, Sec. 2.2), dndy and dndz are
            /// similar.
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
            /// Matrix multiply to determine value at (x,y,z),
            /// multiply by proper weight,
            /// and sum into dk, the stiffness matrix.
            ///

            for (int mm = 0; mm < 3; mm++) {
              for (int nn = 0; nn < 3; nn++) {
                for (int ii = 0; ii < 8; ii++) {
                  for (int jj = 0; jj < 8; jj++) {
                    double sum = 0.0;
                    for (int kk = 0; kk < 6; kk++) {
                      for (int ll = 0; ll < 6; ll++) {
                        sum += es[kk][ii][mm] * cmod_[ijk][kk][ll] *
                               es[ll][jj][nn];
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
    */

  if (iskip == 0) {
    ///
    /// Initialize zcon matrix (gives C term for arbitrary macrostrains)
    ///

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        for (int mi = 0; mi < 3; mi++) {
          for (int mj = 0; mj < 3; mj++) {
            zcon_[i][mi][j][mj] = 0.0;
          }
        }
      }
    }

    ///
    /// Now compute the ss matrices, which give the thermal strain terms for
    /// the i'th phase, single pixel.
    ///

    dndx[0] = -0.25;
    dndx[1] = 0.25;
    dndx[2] = 0.25;
    dndx[3] = -0.25;
    dndx[4] = -0.25;
    dndx[5] = 0.25;
    dndx[6] = 0.25;
    dndx[7] = -0.25;
    dndy[0] = -0.25;
    dndy[1] = -0.25;
    dndy[2] = 0.25;
    dndy[3] = 0.25;
    dndy[4] = -0.25;
    dndy[5] = -0.25;
    dndy[6] = 0.25;
    dndy[7] = 0.25;
    dndz[0] = -0.25;
    dndz[1] = -0.25;
    dndz[2] = -0.25;
    dndz[3] = -0.25;
    dndz[4] = 0.25;
    dndz[5] = 0.25;
    dndz[6] = 0.25;
    dndz[7] = 0.25;

    ///
    /// Build average strain matrix.
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

    for (int mmm = 0; mmm < ns_; mmm++) {
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          sum = 0.0;
          for (int nm = 0; nm < 6; nm++) {
            for (int n = 0; n < 6; n++) {
              sum += cmod_[pix_[mmm]][n][nm] * es[n][mm][nn] * eigen_[mmm][nm];
            }
          }
          ss_[mmm][mm][nn] = sum;
        }
      }
    }

    ///
    /// Call method to generate `zcon_`
    /// `zcon_` is a (2,3) X (2,3) matrix.
    ///

    // constfunc(ns_, nx_, ny_, nz_);
    constfunc(); // depends on pix_

    ///
    /// Set up linear term, T_, for thermal energy. It does not depend on the
    /// macrostrains or displacement, so there is no need to update it as the
    /// macrostrains change. T_ is built up out of the ss_ matrices.
    ///

    int nss = ns_ + 2;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < nss; m++) {
        T_[m][m3] = 0.0;
      }
    }

    ///
    /// For all cases, the correspondence between 0-7 finite element node
    /// (corner) labels  and 0-26 neighbor labels is (see Table 4 in manual):
    ///
    ///   - 0 = ib_[m][26]
    ///   - 1 = ib_[m][2]
    ///   - 2 = ib_[m][1]
    ///   - 3 = ib_[m][0]
    ///   - 4 = ib_[m][25]
    ///   - 5 = ib_[m][18]
    ///   - 6 = ib_[m][17]
    ///   - 7 = ib_[m][16]
    ///

    is[0] = 26;
    is[1] = 2;
    is[2] = 1;
    is[3] = 0;
    is[4] = 25;
    is[5] = 18;
    is[6] = 17;
    is[7] = 16;

    ///
    /// Do all points, but no macrostrain terms
    ///
    /// @note Factor of 2 on linear thermal term is cancelled
    /// by factor of 1/2 out in front of total energy term.
    ///

    for (int k = 0; k < nz_; k++) {
      for (int j = 0; j < ny_; j++) {
        for (int i = 0; i < nx_; i++) {
          m = nxy_ * k + nx_ * j + i;
          for (int mm = 0; mm < 8; mm++) {
            for (int nn = 0; nn < 3; nn++) {
              T_[ib_[m][is[mm]]][nn] -= ss_[m][mm][nn];
            }
          }
        }
      }
    }

    ///
    /// Pick up and sum in all terms multiplying macrostrains.
    ///

    for (int ipp = 0; ipp < 2; ipp++) {
      for (int jpp = 0; jpp < 3; jpp++) {
        exx = 0.0;
        eyy = 0.0;
        ezz = 0.0;
        exz = 0.0;
        eyz = 0.0;
        exy = 0.0;
        if ((ipp == 0) && (jpp == 0))
          exx = 1.0;
        if ((ipp == 0) && (jpp == 1))
          eyy = 1.0;
        if ((ipp == 0) && (jpp == 2))
          ezz = 1.0;
        if ((ipp == 1) && (jpp == 0))
          exz = 1.0;
        if ((ipp == 1) && (jpp == 1))
          eyz = 1.0;
        if ((ipp == 1) && (jpp == 2))
          exy = 1.0;

        // x = nx - 1, face
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 1) || (i8 == 2) || (i8 == 5) || (i8 == 6)) {
              delta[i8][0] = exx * nx_dbl;
              delta[i8][1] = exy * nx_dbl;
              delta[i8][2] = exz * nx_dbl;
            }
          }
        }
        // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
        for (int j = 0; j < ny1; j++) {
          for (int k = 0; k < nz1; k++) {
            m = nxy_ * k + nx_ * j + nx1;
            for (int nn = 0; nn < 3; nn++) {
              for (int mm = 0; mm < 8; mm++) {
                T_[ns_ + ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
              }
            }
          }
        }

        // y = (ny - 1), face
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 2) || (i8 == 3) || (i8 == 6) || (i8 == 7)) {
              delta[i8][0] = exy * ny_dbl;
              delta[i8][1] = eyy * ny_dbl;
              delta[i8][2] = eyz * ny_dbl;
            }
          }
        }
        // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
        for (int i = 0; i < nx1; i++) {
          for (int k = 0; k < nz1; k++) {
            m = nxy_ * k + nx_ * ny1 + i;
            for (int nn = 0; nn < 3; nn++) {
              for (int mm = 0; mm < 8; mm++) {
                T_[ns_ + ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
              }
            }
          }
        }

        // z = (nz - 1), face
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 4) || (i8 == 5) || (i8 == 6) || (i8 == 7)) {
              delta[i8][0] = exz * nz_dbl;
              delta[i8][1] = eyz * nz_dbl;
              delta[i8][2] = ezz * nz_dbl;
            }
          }
        }
        // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
        for (int i = 0; i < nx1; i++) {
          for (int j = 0; j < ny1; j++) {
            m = nxy_ * nz1 + nx_ * j + i;
            for (int nn = 0; nn < 3; nn++) {
              for (int mm = 0; mm < 8; mm++) {
                T_[ns_ + ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
              }
            }
          }
        }

        // x = (nx - 1), y = (ny - 1), edge
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 1) || (i8 == 5)) {
              delta[i8][0] = exx * nx_dbl;
              delta[i8][1] = exy * nx_dbl;
              delta[i8][2] = exz * nx_dbl;
            }
            if ((i8 == 3) || (i8 == 7)) {
              delta[i8][0] = exy * ny_dbl;
              delta[i8][1] = eyy * ny_dbl;
              delta[i8][2] = eyz * ny_dbl;
            }
            if ((i8 == 2) || (i8 == 6)) {
              delta[i8][0] = exy * ny_dbl + exx * nx_dbl;
              delta[i8][1] = eyy * ny_dbl + exy * nx_dbl;
              delta[i8][2] = eyz * ny_dbl + exz * nx_dbl;
            }
          }
        }
        // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
        for (int k = 0; k < nz1; k++) {
          m = nxy_ * k + nx_ * ny1 + nx1;
          for (int nn = 0; nn < 3; nn++) {
            for (int mm = 0; mm < 8; mm++) {
              T_[ns_ + ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
            }
          }
        }

        // x = (nx - 1), z = (nz - 1), edge
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 1) || (i8 == 2)) {
              delta[i8][0] = exx * nx_dbl;
              delta[i8][1] = exy * nx_dbl;
              delta[i8][2] = exz * nx_dbl;
            }
            if ((i8 == 4) || (i8 == 7)) {
              delta[i8][0] = exz * nz_dbl;
              delta[i8][1] = eyz * nz_dbl;
              delta[i8][2] = ezz * nz_dbl;
            }
            if ((i8 == 5) || (i8 == 6)) {
              delta[i8][0] = exz * nz_dbl + exx * nx_dbl;
              delta[i8][1] = eyz * nz_dbl + exy * nx_dbl;
              delta[i8][2] = ezz * nz_dbl + exz * nx_dbl;
            }
          }
        }
        // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
        for (int j = 0; j < ny1; j++) {
          m = nxy_ * nz1 + nx_ * j + nx1;
          for (int nn = 0; nn < 3; nn++) {
            for (int mm = 0; mm < 8; mm++) {
              T_[ns_ + ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
            }
          }
        }

        // y = (ny - 1), z = (nz - 1), edge
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if ((i8 == 4) || (i8 == 5)) {
              delta[i8][0] = exz * nz_dbl;
              delta[i8][1] = eyz * nz_dbl;
              delta[i8][2] = ezz * nz_dbl;
            }
            if ((i8 == 2) || (i8 == 3)) {
              delta[i8][0] = exy * ny_dbl;
              delta[i8][1] = eyy * ny_dbl;
              delta[i8][2] = eyz * ny_dbl;
            }
            if ((i8 == 6) || (i8 == 7)) {
              delta[i8][0] = exy * ny_dbl + exz * nz_dbl;
              delta[i8][1] = eyy * ny_dbl + eyz * nz_dbl;
              delta[i8][2] = eyz * ny_dbl + ezz * nz_dbl;
            }
          }
        }
        // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
        for (int i = 0; i < nx1; i++) {
          m = nxy_ * nz1 + nx_ * ny1 + i;
          for (int nn = 0; nn < 3; nn++) {
            for (int mm = 0; mm < 8; mm++) {
              T_[ns_ + ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
            }
          }
        }

        // x = (nx - 1), y = (ny - 1), z = (nz - 1), corner
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i8 = 0; i8 < 8; i8++) {
            delta[i8][i3] = 0.0;
            if (i8 == 1) {
              delta[i8][0] = exx * nx_dbl;
              delta[i8][1] = exy * nx_dbl;
              delta[i8][2] = exz * nx_dbl;
            }
            if (i8 == 3) {
              delta[i8][0] = exy * ny_dbl;
              delta[i8][1] = eyy * ny_dbl;
              delta[i8][2] = eyz * ny_dbl;
            }
            if (i8 == 4) {
              delta[i8][0] = exz * nz_dbl;
              delta[i8][1] = eyz * nz_dbl;
              delta[i8][2] = ezz * nz_dbl;
            }
            if (i8 == 7) {
              delta[i8][0] = exy * ny_dbl + exz * nz_dbl;
              delta[i8][1] = eyy * ny_dbl + eyz * nz_dbl;
              delta[i8][2] = eyz * ny_dbl + ezz * nz_dbl;
            }
            if (i8 == 5) {
              delta[i8][0] = exx * nx_dbl + exz * nz_dbl;
              delta[i8][1] = exy * nx_dbl + eyz * nz_dbl;
              delta[i8][2] = exz * nx_dbl + ezz * nz_dbl;
            }
            if (i8 == 2) {
              delta[i8][0] = exx * nx_dbl + exy * ny_dbl;
              delta[i8][1] = exy * nx_dbl + eyy * ny_dbl;
              delta[i8][2] = exz * nx_dbl + eyz * ny_dbl;
            }
            if (i8 == 6) {
              delta[i8][0] = exx * nx_dbl + exy * ny_dbl + exz * nz_dbl;
              delta[i8][1] = exy * nx_dbl + eyy * ny_dbl + eyz * nz_dbl;
              delta[i8][2] = exz * nx_dbl + eyz * ny_dbl + ezz * nz_dbl;
            }
          }
        }
        // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
        m = nxy_ * nz1 + nx_ * ny1 + nx1;
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            T_[ns_ + ipp][jpp] -= ss_[m][mm][nn] * delta[mm][nn];
          }
        }
      }
    }

    ///
    /// Compute Y, the 0.5 (eigenstrain<sub>i</sub>) C<sub>ij</sub>
    /// (eigenstrain<sub>j</sub>) energy, doesn't ever change with macrostrain
    /// or displacements.
    ///

    Y_ = 0.0;
    for (int m = 0; m < ns_; m++) {
      for (int n = 0; n < 6; n++) {
        for (int nn = 0; nn < 6; nn++) {
          Y_ += 0.5 * eigen_[m][n] * cmod_[pix_[m]][n][nn] * eigen_[m][nn];
        }
      }
    }
  } // if (iskip == 0) {

  ///
  /// Use auxiliary variables (exx, etc.) instead of u_[] variable,
  /// for convenience, and to make the following code easier to read.
  ///

  exx = u_[ns_][0];
  eyy = u_[ns_][1];
  ezz = u_[ns_][2];
  exz = u_[ns_ + 1][0];
  eyz = u_[ns_ + 1][1];
  exy = u_[ns_ + 1][2];

  ///
  /// Set up vector for linear term that comes from periodic boundary
  /// conditions. Notation and conventions same as for `T_` term.
  /// This is done using the stiffness matrices, and the periodic terms in
  /// the macrostrains. It is easier to set up `b_` this way than to
  /// analytically write out all the terms involved.
  ///

  for (int m3 = 0; m3 < 3; m3++) {
    for (int m = 0; m < ns_; m++) {
      b_[m][m3] = 0.0;
    }
  }

  ///
  /// For all cases, the correspondence between 0-7 finite element node (corner)
  /// labels  and 0-26 neighbor labels is (see Table 4 in manual):
  ///
  ///   - 0 = ib_[m][26]
  ///   - 1 = ib_[m][2]
  ///   - 2 = ib_[m][1]
  ///   - 3 = ib_[m][0]
  ///   - 4 = ib_[m][25]
  ///   - 5 = ib_[m][18]
  ///   - 6 = ib_[m][17]
  ///   - 7 = ib_[m][16]
  ///

  is[0] = 26;
  is[1] = 2;
  is[2] = 1;
  is[3] = 0;
  is[4] = 25;
  is[5] = 18;
  is[6] = 17;
  is[7] = 16;

  C_ = 0.0;

  // x = nx - 1, face
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 1) || (i8 == 2) || (i8 == 5) || (i8 == 6)) {
        delta[i8][0] = exx * nx_dbl;
        delta[i8][1] = exy * nx_dbl;
        delta[i8][2] = exz * nx_dbl;
      }
    }
  }
  // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
  for (int j = 0; j < ny1; j++) {
    for (int k = 0; k < nz1; k++) {
      m = nxy_ * k + nx_ * j + nx1;
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          b_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
  }

  // y = ny - 1, face
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 2) || (i8 == 3) || (i8 == 6) || (i8 == 7)) {
        delta[i8][0] = exy * ny_dbl;
        delta[i8][1] = eyy * ny_dbl;
        delta[i8][2] = eyz * ny_dbl;
      }
    }
  }
  // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
  for (int i = 0; i < nx1; i++) {
    for (int k = 0; k < nz1; k++) {
      m = nxy_ * k + nx_ * ny1 + i;
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          b_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
  }

  // z = nz - 1, face
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 4) || (i8 == 5) || (i8 == 6) || (i8 == 7)) {
        delta[i8][0] = exz * nz_dbl;
        delta[i8][1] = eyz * nz_dbl;
        delta[i8][2] = ezz * nz_dbl;
      }
    }
  }
  // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
  for (int i = 0; i < nx1; i++) {
    for (int j = 0; j < ny1; j++) {
      m = nxy_ * nz1 + nx_ * j + i;
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          b_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
  }

  // x = nx - 1, y = ny - 1, edge
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 1) || (i8 == 5)) {
        delta[i8][0] = exx * nx_dbl;
        delta[i8][1] = exy * nx_dbl;
        delta[i8][2] = exz * nx_dbl;
      }
      if ((i8 == 3) || (i8 == 7)) {
        delta[i8][0] = exy * ny_dbl;
        delta[i8][1] = eyy * ny_dbl;
        delta[i8][2] = eyz * ny_dbl;
      }
      if ((i8 == 2) || (i8 == 6)) {
        delta[i8][0] = exy * ny_dbl + exx * nx_dbl;
        delta[i8][1] = eyy * ny_dbl + exy * nx_dbl;
        delta[i8][2] = eyz * ny_dbl + exz * nx_dbl;
      }
    }
  }
  // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
  for (int k = 0; k < nz1; k++) {
    m = nxy_ * k + nx_ * ny1 + nx1;
    for (int nn = 0; nn < 3; nn++) {
      for (int mm = 0; mm < 8; mm++) {
        sum = 0.0;
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m8 = 0; m8 < 8; m8++) {
            sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
          }
        }
        b_[ib_[m][is[mm]]][nn] += sum;
      }
    }
  }

  // x = nx - 1, z = nz - 1, edge
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 1) || (i8 == 2)) {
        delta[i8][0] = exx * nx_dbl;
        delta[i8][1] = exy * nx_dbl;
        delta[i8][2] = exz * nx_dbl;
      }
      if ((i8 == 4) || (i8 == 7)) {
        delta[i8][0] = exz * nz_dbl;
        delta[i8][1] = eyz * nz_dbl;
        delta[i8][2] = ezz * nz_dbl;
      }
      if ((i8 == 5) || (i8 == 6)) {
        delta[i8][0] = exz * nz_dbl + exx * nx_dbl;
        delta[i8][1] = eyz * nz_dbl + exy * nx_dbl;
        delta[i8][2] = ezz * nz_dbl + exz * nx_dbl;
      }
    }
  }
  // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
  for (int j = 0; j < ny1; j++) {
    m = nxy_ * nz1 + nx_ * j + nx1;
    for (int nn = 0; nn < 3; nn++) {
      for (int mm = 0; mm < 8; mm++) {
        sum = 0.0;
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m8 = 0; m8 < 8; m8++) {
            sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
          }
        }
        b_[ib_[m][is[mm]]][nn] += sum;
      }
    }
  }

  // y = ny - 1, z = nz - 1, edge
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 4) || (i8 == 5)) {
        delta[i8][0] = exz * nz_dbl;
        delta[i8][1] = eyz * nz_dbl;
        delta[i8][2] = ezz * nz_dbl;
      }
      if ((i8 == 2) || (i8 == 3)) {
        delta[i8][0] = exy * ny_dbl;
        delta[i8][1] = eyy * ny_dbl;
        delta[i8][2] = eyz * ny_dbl;
      }
      if ((i8 == 6) || (i8 == 7)) {
        delta[i8][0] = exy * ny_dbl + exz * nz_dbl;
        delta[i8][1] = eyy * ny_dbl + eyz * nz_dbl;
        delta[i8][2] = eyz * ny_dbl + ezz * nz_dbl;
      }
    }
  }
  // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
  for (int i = 0; i < nx1; i++) {
    m = nxy_ * nz1 + nx_ * ny1 + i;
    for (int nn = 0; nn < 3; nn++) {
      for (int mm = 0; mm < 8; mm++) {
        sum = 0.0;
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m8 = 0; m8 < 8; m8++) {
            sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
          }
        }
        b_[ib_[m][is[mm]]][nn] += sum;
      }
    }
  }

  // x = nx - 1, y = ny - 1, z = nz - 1, corner
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if (i8 == 1) {
        delta[i8][0] = exx * nx_dbl;
        delta[i8][1] = exy * nx_dbl;
        delta[i8][2] = exz * nx_dbl;
      }
      if (i8 == 3) {
        delta[i8][0] = exy * ny_dbl;
        delta[i8][1] = eyy * ny_dbl;
        delta[i8][2] = eyz * ny_dbl;
      }
      if (i8 == 4) {
        delta[i8][0] = exz * nz_dbl;
        delta[i8][1] = eyz * nz_dbl;
        delta[i8][2] = ezz * nz_dbl;
      }
      if (i8 == 7) {
        delta[i8][0] = exy * ny_dbl + exz * nz_dbl;
        delta[i8][1] = eyy * ny_dbl + eyz * nz_dbl;
        delta[i8][2] = eyz * ny_dbl + ezz * nz_dbl;
      }
      if (i8 == 5) {
        delta[i8][0] = exx * nx_dbl + exz * nz_dbl;
        delta[i8][1] = exy * nx_dbl + eyz * nz_dbl;
        delta[i8][2] = exz * nx_dbl + ezz * nz_dbl;
      }
      if (i8 == 2) {
        delta[i8][0] = exx * nx_dbl + exy * ny_dbl;
        delta[i8][1] = exy * nx_dbl + eyy * ny_dbl;
        delta[i8][2] = exz * nx_dbl + eyz * ny_dbl;
      }
      if (i8 == 6) {
        delta[i8][0] = exx * nx_dbl + exy * ny_dbl + exz * nz_dbl;
        delta[i8][1] = exy * nx_dbl + eyy * ny_dbl + eyz * nz_dbl;
        delta[i8][2] = exz * nx_dbl + eyz * ny_dbl + ezz * nz_dbl;
      }
    }
  }

  // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
  m = nxy_ * nz1 + nx_ * ny1 + nx1;
  for (int nn = 0; nn < 3; nn++) {
    for (int mm = 0; mm < 8; mm++) {
      sum = 0.0;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m8 = 0; m8 < 8; m8++) {
          sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
        }
      }
      b_[ib_[m][is[mm]]][nn] += sum;
    }
  }

  return;
}

// void ThermalStrain::bgrad(int nx, int ny, int nz, double exx,
//                           double eyy, double ezz, double exz, double eyz,
//                           double exy) {
void ThermalStrain::bgrad(double exx, double eyy, double ezz, double exz,
                          double eyz, double exy) {
  int is[8];
  double delta[8][3];
  int nx1 = nx_ - 1;
  int ny1 = ny_ - 1;
  int nz1 = nz_ - 1;
  double nx_dbl = nx_;
  double ny_dbl = ny_;
  double nz_dbl = nz_;
  int m;
  double sum;

  ///
  /// exx, eyy, ezz, exz, eyz, exy are the artificial macrostrains used to
  /// get the gradient terms (appropriate combinations of 1's and 0's).

  ///
  /// Set up vector for linear term.
  ///

  is[0] = 26;
  is[1] = 2;
  is[2] = 1;
  is[3] = 0;
  is[4] = 25;
  is[5] = 18;
  is[6] = 17;
  is[7] = 16;

  // x = nx - 1, face
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 1) || (i8 == 2) || (i8 == 5) || (i8 == 6)) {
        delta[i8][0] = exx * nx_dbl;
        delta[i8][1] = exy * nx_dbl;
        delta[i8][2] = exz * nx_dbl;
      }
    }
  }
  // nx1 = nx_ - 1, ny1 = ny_ - 1, nz1 = nz_ - 1
  for (int j = 0; j < ny1; j++) {
    for (int k = 0; k < nz1; k++) {
      m = nxy_ * k + nx_ * j + nx1;
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          if (exx == 1.0)
            b0_[ib_[m][is[mm]]][nn] += sum;
          if (eyy == 1.0)
            b1_[ib_[m][is[mm]]][nn] += sum;
          if (ezz == 1.0)
            b2_[ib_[m][is[mm]]][nn] += sum;
          if (exz == 1.0)
            b3_[ib_[m][is[mm]]][nn] += sum;
          if (eyz == 1.0)
            b4_[ib_[m][is[mm]]][nn] += sum;
          if (exy == 1.0)
            b5_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
  }

  // y = ny - 1, face
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 2) || (i8 == 3) || (i8 == 6) || (i8 == 7)) {
        delta[i8][0] = exy * ny_dbl;
        delta[i8][1] = eyy * ny_dbl;
        delta[i8][2] = eyz * ny_dbl;
      }
    }
  }
  // nx1 = nx_ - 1, ny1 = ny_ - 1, nz1 = nz_ - 1
  for (int i = 0; i < nx1; i++) {
    for (int k = 0; k < nz1; k++) {
      m = nxy_ * k + nx_ * ny1 + i;
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          if (exx == 1.0)
            b0_[ib_[m][is[mm]]][nn] += sum;
          if (eyy == 1.0)
            b1_[ib_[m][is[mm]]][nn] += sum;
          if (ezz == 1.0)
            b2_[ib_[m][is[mm]]][nn] += sum;
          if (exz == 1.0)
            b3_[ib_[m][is[mm]]][nn] += sum;
          if (eyz == 1.0)
            b4_[ib_[m][is[mm]]][nn] += sum;
          if (exy == 1.0)
            b5_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
  }

  // z = nz - 1, face
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 4) || (i8 == 5) || (i8 == 6) || (i8 == 7)) {
        delta[i8][0] = exz * nz_dbl;
        delta[i8][1] = eyz * nz_dbl;
        delta[i8][2] = ezz * nz_dbl;
      }
    }
  }
  // nx1 = nx_ - 1, ny1 = ny_ - 1, nz1 = nz_ - 1
  for (int i = 0; i < nx1; i++) {
    for (int j = 0; j < ny1; j++) {
      m = nxy_ * nz1 + nx_ * j + i;
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          sum = 0.0;
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
            }
          }
          if (exx == 1.0)
            b0_[ib_[m][is[mm]]][nn] += sum;
          if (eyy == 1.0)
            b1_[ib_[m][is[mm]]][nn] += sum;
          if (ezz == 1.0)
            b2_[ib_[m][is[mm]]][nn] += sum;
          if (exz == 1.0)
            b3_[ib_[m][is[mm]]][nn] += sum;
          if (eyz == 1.0)
            b4_[ib_[m][is[mm]]][nn] += sum;
          if (exy == 1.0)
            b5_[ib_[m][is[mm]]][nn] += sum;
        }
      }
    }
  }

  // x = nx - 1, y = ny - 1, edge
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 1) || (i8 == 5)) {
        delta[i8][0] = exx * nx_dbl;
        delta[i8][1] = exy * nx_dbl;
        delta[i8][2] = exz * nx_dbl;
      }
      if ((i8 == 3) || (i8 == 7)) {
        delta[i8][0] = exy * ny_dbl;
        delta[i8][1] = eyy * ny_dbl;
        delta[i8][2] = eyz * ny_dbl;
      }
      if ((i8 == 2) || (i8 == 6)) {
        delta[i8][0] = exy * ny_dbl + exx * nx_dbl;
        delta[i8][1] = eyy * ny_dbl + exy * nx_dbl;
        delta[i8][2] = eyz * ny_dbl + exz * nx_dbl;
      }
    }
  }
  // nx1 = nx_ - 1, ny1 = ny_ - 1, nz1 = nz_ - 1
  for (int k = 0; k < nz1; k++) {
    m = nxy_ * k + nx_ * ny1 + nx1;
    for (int nn = 0; nn < 3; nn++) {
      for (int mm = 0; mm < 8; mm++) {
        sum = 0.0;
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m8 = 0; m8 < 8; m8++) {
            sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
          }
        }
        if (exx == 1.0)
          b0_[ib_[m][is[mm]]][nn] += sum;
        if (eyy == 1.0)
          b1_[ib_[m][is[mm]]][nn] += sum;
        if (ezz == 1.0)
          b2_[ib_[m][is[mm]]][nn] += sum;
        if (exz == 1.0)
          b3_[ib_[m][is[mm]]][nn] += sum;
        if (eyz == 1.0)
          b4_[ib_[m][is[mm]]][nn] += sum;
        if (exy == 1.0)
          b5_[ib_[m][is[mm]]][nn] += sum;
      }
    }
  }

  // x = nx - 1, z = nz - 1, edge
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 1) || (i8 == 2)) {
        delta[i8][0] = exx * nx_dbl;
        delta[i8][1] = exy * nx_dbl;
        delta[i8][2] = exz * nx_dbl;
      }
      if ((i8 == 4) || (i8 == 7)) {
        delta[i8][0] = exz * nz_dbl;
        delta[i8][1] = eyz * nz_dbl;
        delta[i8][2] = ezz * nz_dbl;
      }
      if ((i8 == 5) || (i8 == 6)) {
        delta[i8][0] = exz * nz_dbl + exx * nx_dbl;
        delta[i8][1] = eyz * nz_dbl + exy * nx_dbl;
        delta[i8][2] = ezz * nz_dbl + exz * nx_dbl;
      }
    }
  }
  // nx1 = nx_ - 1, ny1 = ny_ - 1, nz1 = nz_ - 1
  for (int j = 0; j < ny1; j++) {
    m = nxy_ * nz1 + nx_ * j + nx1;
    for (int nn = 0; nn < 3; nn++) {
      for (int mm = 0; mm < 8; mm++) {
        sum = 0.0;
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m8 = 0; m8 < 8; m8++) {
            sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
          }
        }
        if (exx == 1.0)
          b0_[ib_[m][is[mm]]][nn] += sum;
        if (eyy == 1.0)
          b1_[ib_[m][is[mm]]][nn] += sum;
        if (ezz == 1.0)
          b2_[ib_[m][is[mm]]][nn] += sum;
        if (exz == 1.0)
          b3_[ib_[m][is[mm]]][nn] += sum;
        if (eyz == 1.0)
          b4_[ib_[m][is[mm]]][nn] += sum;
        if (exy == 1.0)
          b5_[ib_[m][is[mm]]][nn] += sum;
      }
    }
  }

  // y = ny - 1, z = nz - 1, edge
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 4) || (i8 == 5)) {
        delta[i8][0] = exz * nz_dbl;
        delta[i8][1] = eyz * nz_dbl;
        delta[i8][2] = ezz * nz_dbl;
      }
      if ((i8 == 2) || (i8 == 3)) {
        delta[i8][0] = exy * ny_dbl;
        delta[i8][1] = eyy * ny_dbl;
        delta[i8][2] = eyz * ny_dbl;
      }
      if ((i8 == 6) || (i8 == 7)) {
        delta[i8][0] = exy * ny_dbl + exz * nz_dbl;
        delta[i8][1] = eyy * ny_dbl + eyz * nz_dbl;
        delta[i8][2] = eyz * ny_dbl + ezz * nz_dbl;
      }
    }
  }
  // nx1 = nx_ - 1, ny1 = ny_ - 1, nz1 = nz_ - 1
  for (int i = 0; i < nx1; i++) {
    m = nxy_ * nz1 + nx_ * ny1 + i;
    for (int nn = 0; nn < 3; nn++) {
      for (int mm = 0; mm < 8; mm++) {
        sum = 0.0;
        for (int m3 = 0; m3 < 3; m3++) {
          for (int m8 = 0; m8 < 8; m8++) {
            sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
          }
        }
        if (exx == 1.0)
          b0_[ib_[m][is[mm]]][nn] += sum;
        if (eyy == 1.0)
          b1_[ib_[m][is[mm]]][nn] += sum;
        if (ezz == 1.0)
          b2_[ib_[m][is[mm]]][nn] += sum;
        if (exz == 1.0)
          b3_[ib_[m][is[mm]]][nn] += sum;
        if (eyz == 1.0)
          b4_[ib_[m][is[mm]]][nn] += sum;
        if (exy == 1.0)
          b5_[ib_[m][is[mm]]][nn] += sum;
      }
    }
  }

  // x = nx - 1, y = ny - 1, z = nz - 1, corner
  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if (i8 == 1) {
        delta[i8][0] = exx * nx_dbl;
        delta[i8][1] = exy * nx_dbl;
        delta[i8][2] = exz * nx_dbl;
      }
      if (i8 == 3) {
        delta[i8][0] = exy * ny_dbl;
        delta[i8][1] = eyy * ny_dbl;
        delta[i8][2] = eyz * ny_dbl;
      }
      if (i8 == 4) {
        delta[i8][0] = exz * nz_dbl;
        delta[i8][1] = eyz * nz_dbl;
        delta[i8][2] = ezz * nz_dbl;
      }
      if (i8 == 7) {
        delta[i8][0] = exy * ny_dbl + exz * nz_dbl;
        delta[i8][1] = eyy * ny_dbl + eyz * nz_dbl;
        delta[i8][2] = eyz * ny_dbl + ezz * nz_dbl;
      }
      if (i8 == 5) {
        delta[i8][0] = exx * nx_dbl + exz * nz_dbl;
        delta[i8][1] = exy * nx_dbl + eyz * nz_dbl;
        delta[i8][2] = exz * nx_dbl + ezz * nz_dbl;
      }
      if (i8 == 2) {
        delta[i8][0] = exx * nx_dbl + exy * ny_dbl;
        delta[i8][1] = exy * nx_dbl + eyy * ny_dbl;
        delta[i8][2] = exz * nx_dbl + eyz * ny_dbl;
      }
      if (i8 == 6) {
        delta[i8][0] = exx * nx_dbl + exy * ny_dbl + exz * nz_dbl;
        delta[i8][1] = exy * nx_dbl + eyy * ny_dbl + eyz * nz_dbl;
        delta[i8][2] = exz * nx_dbl + eyz * ny_dbl + ezz * nz_dbl;
      }
    }
  }
  // nx1 = nx_ - 1, ny1 = ny_ - 1, nz1 = nz_ - 1
  m = nxy_ * nz1 + nx_ * ny1 + nx1;
  for (int nn = 0; nn < 3; nn++) {
    for (int mm = 0; mm < 8; mm++) {
      sum = 0.0;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m8 = 0; m8 < 8; m8++) {
          sum += delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn];
        }
      }
      if (exx == 1.0)
        b0_[ib_[m][is[mm]]][nn] += sum;
      if (eyy == 1.0)
        b1_[ib_[m][is[mm]]][nn] += sum;
      if (ezz == 1.0)
        b2_[ib_[m][is[mm]]][nn] += sum;
      if (exz == 1.0)
        b3_[ib_[m][is[mm]]][nn] += sum;
      if (eyz == 1.0)
        b4_[ib_[m][is[mm]]][nn] += sum;
      if (exy == 1.0)
        b5_[ib_[m][is[mm]]][nn] += sum;
    }
  }

  return;
}

// void ThermalStrain::constfunc(int ns, int nx, int ny, int nz) {
void ThermalStrain::constfunc() {
  double delta[8][3];
  double pp[6][6], s[6][6];
  double econ = 0;

  ///
  /// Routine to set up 6 x 6 matrix for energy term involving macrostrains
  /// only, pulled out of femat.
  ///
  /// Idea is to evaluate the quadratic term in the macrostrains repeatedly for
  /// choices of strain like exx = 1, exy = 1, all others = 0, build up 21
  /// choices, then recombine to get matrix elements by themselves.
  ///

  int nx1 = nx_ - 1;
  int ny1 = ny_ - 1;
  int nz1 = nz_ - 1;
  double nx_dbl = nx_;
  double ny_dbl = ny_;
  double nz_dbl = nz_;
  int m;

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      s[i][j] = 0.0;
      pp[i][j] = 0.0;
    }
  }

  double exx = 0, eyy = 0, ezz = 0, exz = 0, eyz = 0, exy = 0;

  for (int ii = 0; ii < 6; ii++) {
    for (int jj = ii; jj < 6; jj++) {
      econ = 0.0;
      exx = eyy = ezz = 0.0;
      exz = eyz = exy = 0.0;
      if ((ii == 0) && (jj == 0))
        exx = 1.0;
      if ((ii == 1) && (jj == 1))
        eyy = 1.0;
      if ((ii == 2) && (jj == 2))
        ezz = 1.0;
      if ((ii == 3) && (jj == 3))
        exz = 1.0;
      if ((ii == 4) && (jj == 4))
        eyz = 1.0;
      if ((ii == 5) && (jj == 5))
        exy = 1.0;
      if ((ii == 0) && (jj == 1)) {
        exx = 1.0;
        eyy = 1.0;
      }
      if ((ii == 0) && (jj == 2)) {
        exx = 1.0;
        ezz = 1.0;
      }
      if ((ii == 0) && (jj == 3)) {
        exx = 1.0;
        exz = 1.0;
      }
      if ((ii == 0) && (jj == 4)) {
        exx = 1.0;
        eyz = 1.0;
      }
      if ((ii == 0) && (jj == 5)) {
        exx = 1.0;
        exy = 1.0;
      }
      if ((ii == 1) && (jj == 2)) {
        eyy = 1.0;
        ezz = 1.0;
      }
      if ((ii == 1) && (jj == 3)) {
        eyy = 1.0;
        exz = 1.0;
      }
      if ((ii == 1) && (jj == 4)) {
        eyy = 1.0;
        eyz = 1.0;
      }
      if ((ii == 1) && (jj == 5)) {
        eyy = 1.0;
        exy = 1.0;
      }
      if ((ii == 2) && (jj == 3)) {
        ezz = 1.0;
        exz = 1.0;
      }
      if ((ii == 2) && (jj == 4)) {
        ezz = 1.0;
        eyz = 1.0;
      }
      if ((ii == 2) && (jj == 5)) {
        ezz = 1.0;
        exy = 1.0;
      }
      if ((ii == 3) && (jj == 4)) {
        exz = 1.0;
        eyz = 1.0;
      }
      if ((ii == 3) && (jj == 5)) {
        exz = 1.0;
        exy = 1.0;
      }
      if ((ii == 4) && (jj == 5)) {
        eyz = 1.0;
        exy = 1.0;
      }
      // x = nx - 1, face
      for (int i3 = 0; i3 < 3; i3++) {
        for (int i8 = 0; i8 < 8; i8++) {
          delta[i8][i3] = 0.0;
          if ((i8 == 1) || (i8 == 2) || (i8 == 5) || (i8 == 6)) {
            delta[i8][0] = exx * nx_dbl;
            delta[i8][1] = exy * nx_dbl;
            delta[i8][2] = exz * nx_dbl;
          }
        }
      }
      // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
      for (int j = 0; j < ny1; j++) {
        for (int k = 0; k < nz1; k++) {
          m = nxy_ * k + nx_ * j + nx1;
          for (int nn = 0; nn < 3; nn++) {
            for (int mm = 0; mm < 8; mm++) {
              for (int m3 = 0; m3 < 3; m3++) {
                for (int m8 = 0; m8 < 8; m8++) {
                  econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                          delta[mm][nn];
                }
              }
            }
          }
        }
      }

      // y = ny - 1, face
      for (int i3 = 0; i3 < 3; i3++) {
        for (int i8 = 0; i8 < 8; i8++) {
          delta[i8][i3] = 0.0;
          if ((i8 == 2) || (i8 == 3) || (i8 == 6) || (i8 == 7)) {
            delta[i8][0] = exy * ny_dbl;
            delta[i8][1] = eyy * ny_dbl;
            delta[i8][2] = eyz * ny_dbl;
          }
        }
      }
      // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
      for (int i = 0; i < nx1; i++) {
        for (int k = 0; k < nz1; k++) {
          m = nxy_ * k + nx_ * ny1 + i;
          for (int nn = 0; nn < 3; nn++) {
            for (int mm = 0; mm < 8; mm++) {
              for (int m3 = 0; m3 < 3; m3++) {
                for (int m8 = 0; m8 < 8; m8++) {
                  econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                          delta[mm][nn];
                }
              }
            }
          }
        }
      }

      // z = nz - 1, face
      for (int i3 = 0; i3 < 3; i3++) {
        for (int i8 = 0; i8 < 8; i8++) {
          delta[i8][i3] = 0.0;
          if ((i8 == 4) || (i8 == 5) || (i8 == 6) || (i8 == 7)) {
            delta[i8][0] = exz * nz_dbl;
            delta[i8][1] = eyz * nz_dbl;
            delta[i8][2] = ezz * nz_dbl;
          }
        }
      }
      // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
      for (int i = 0; i < nx1; i++) {
        for (int j = 0; j < ny1; j++) {
          m = nxy_ * nz1 + nx_ * j + i;
          for (int nn = 0; nn < 3; nn++) {
            for (int mm = 0; mm < 8; mm++) {
              for (int m3 = 0; m3 < 3; m3++) {
                for (int m8 = 0; m8 < 8; m8++) {
                  econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                          delta[mm][nn];
                }
              }
            }
          }
        }
      }

      // x = nx - 1, y = ny - 1, edge
      for (int i3 = 0; i3 < 3; i3++) {
        for (int i8 = 0; i8 < 8; i8++) {
          delta[i8][i3] = 0.0;
          if ((i8 == 1) || (i8 == 5)) {
            delta[i8][0] = exx * nx_dbl;
            delta[i8][1] = exy * nx_dbl;
            delta[i8][2] = exz * nx_dbl;
          }
          if ((i8 == 3) || (i8 == 7)) {
            delta[i8][0] = exy * ny_dbl;
            delta[i8][1] = eyy * ny_dbl;
            delta[i8][2] = eyz * ny_dbl;
          }
          if ((i8 == 2) || (i8 == 6)) {
            delta[i8][0] = exy * ny_dbl + exx * nx_dbl;
            delta[i8][1] = eyy * ny_dbl + exy * nx_dbl;
            delta[i8][2] = eyz * ny_dbl + exz * nx_dbl;
          }
        }
      }
      // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
      for (int k = 0; k < nz1; k++) {
        m = nxy_ * k + nx_ * ny1 + nx1;
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            for (int m3 = 0; m3 < 3; m3++) {
              for (int m8 = 0; m8 < 8; m8++) {
                econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                        delta[mm][nn];
              }
            }
          }
        }
      }

      // x = nx - 1, z = nz - 1, edge
      for (int i3 = 0; i3 < 3; i3++) {
        for (int i8 = 0; i8 < 8; i8++) {
          delta[i8][i3] = 0.0;
          if ((i8 == 1) || (i8 == 2)) {
            delta[i8][0] = exx * nx_dbl;
            delta[i8][1] = exy * nx_dbl;
            delta[i8][2] = exz * nx_dbl;
          }
          if ((i8 == 4) || (i8 == 7)) {
            delta[i8][0] = exz * nz_dbl;
            delta[i8][1] = eyz * nz_dbl;
            delta[i8][2] = ezz * nz_dbl;
          }
          if ((i8 == 5) || (i8 == 6)) {
            delta[i8][0] = exz * nz_dbl + exx * nx_dbl;
            delta[i8][1] = eyz * nz_dbl + exy * nx_dbl;
            delta[i8][2] = ezz * nz_dbl + exz * nx_dbl;
          }
        }
      }
      // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
      for (int j = 0; j < ny1; j++) {
        m = nxy_ * nz1 + nx_ * j + nx1;
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            for (int m3 = 0; m3 < 3; m3++) {
              for (int m8 = 0; m8 < 8; m8++) {
                econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                        delta[mm][nn];
              }
            }
          }
        }
      }

      // y = ny - 1, z = nz - 1, edge
      for (int i3 = 0; i3 < 3; i3++) {
        for (int i8 = 0; i8 < 8; i8++) {
          delta[i8][i3] = 0.0;
          if ((i8 == 4) || (i8 == 5)) {
            delta[i8][0] = exz * nz_dbl;
            delta[i8][1] = eyz * nz_dbl;
            delta[i8][2] = ezz * nz_dbl;
          }
          if ((i8 == 2) || (i8 == 3)) {
            delta[i8][0] = exy * ny_dbl;
            delta[i8][1] = eyy * ny_dbl;
            delta[i8][2] = eyz * ny_dbl;
          }
          if ((i8 == 6) || (i8 == 7)) {
            delta[i8][0] = exy * ny_dbl + exz * nz_dbl;
            delta[i8][1] = eyy * ny_dbl + eyz * nz_dbl;
            delta[i8][2] = eyz * ny_dbl + ezz * nz_dbl;
          }
        }
      }
      // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
      for (int i = 0; i < nx1; i++) {
        m = nxy_ * nz1 + nx_ * ny1 + i;
        for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 8; mm++) {
            for (int m3 = 0; m3 < 3; m3++) {
              for (int m8 = 0; m8 < 8; m8++) {
                econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                        delta[mm][nn];
              }
            }
          }
        }
      }

      // x = nx - 1, y = ny - 1, z = nz - 1, corner
      for (int i3 = 0; i3 < 3; i3++) {
        for (int i8 = 0; i8 < 8; i8++) {
          delta[i8][i3] = 0.0;
          if (i8 == 1) {
            delta[i8][0] = exx * nx_dbl;
            delta[i8][1] = exy * nx_dbl;
            delta[i8][2] = exz * nx_dbl;
          }
          if (i8 == 3) {
            delta[i8][0] = exy * ny_dbl;
            delta[i8][1] = eyy * ny_dbl;
            delta[i8][2] = eyz * ny_dbl;
          }
          if (i8 == 4) {
            delta[i8][0] = exz * nz_dbl;
            delta[i8][1] = eyz * nz_dbl;
            delta[i8][2] = ezz * nz_dbl;
          }
          if (i8 == 7) {
            delta[i8][0] = exy * ny_dbl + exz * nz_dbl;
            delta[i8][1] = eyy * ny_dbl + eyz * nz_dbl;
            delta[i8][2] = eyz * ny_dbl + ezz * nz_dbl;
          }
          if (i8 == 5) {
            delta[i8][0] = exx * nx_dbl + exz * nz_dbl;
            delta[i8][1] = exy * nx_dbl + eyz * nz_dbl;
            delta[i8][2] = exz * nx_dbl + ezz * nz_dbl;
          }
          if (i8 == 2) {
            delta[i8][0] = exx * nx_dbl + exy * ny_dbl;
            delta[i8][1] = exy * nx_dbl + eyy * ny_dbl;
            delta[i8][2] = exz * nx_dbl + eyz * ny_dbl;
          }
          if (i8 == 6) {
            delta[i8][0] = exx * nx_dbl + exy * ny_dbl + exz * nz_dbl;
            delta[i8][1] = exy * nx_dbl + eyy * ny_dbl + eyz * nz_dbl;
            delta[i8][2] = exz * nx_dbl + eyz * ny_dbl + ezz * nz_dbl;
          }
        }
      }
      // nx1 = nx - 1, ny1 = ny - 1, nz1 = nz - 1
      m = nxy_ * nz1 + nx_ * ny1 + nx1;
      for (int nn = 0; nn < 3; nn++) {
        for (int mm = 0; mm < 8; mm++) {
          for (int m3 = 0; m3 < 3; m3++) {
            for (int m8 = 0; m8 < 8; m8++) {
              econ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                      delta[mm][nn];
            }
          }
        }
      }
      pp[ii][jj] = econ * 2;
    }
  }

  for (int i = 0; i < 6; i++) {
    for (int j = i; j < 6; j++) {
      if (i == j) {
        s[i][j] = pp[i][j];
      } else {
        s[i][j] = pp[i][j] - pp[i][i] - pp[j][j];
      }
    }
  }
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      pp[i][j] = 0.5 * (s[i][j] + s[j][i]);
    }
  }

  ///
  /// Map pp[i][j] into zcon_[2][3][2][3]
  ///

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int mi = 0; mi < 3; mi++) {
        for (int mj = 0; mj < 3; mj++) {
          int ii, jj;
          if (i == 0)
            ii = i + mi;
          if (i == 1)
            ii = i + mi + 2;
          if (j == 0)
            jj = j + mj;
          if (j == 1)
            jj = j + mj + 2;
          zcon_[i][mi][j][mj] = pp[ii][jj];
          if (verbose_) {
            std::clog << "ThermalStrain::constfunc i, mi, j, mj, zcon_[" << i
                      << "][" << mi << "][" << j << "][" << mj << "]" << i
                      << " " << mi << " " << j << " " << mj << " "
                      << zcon_[i][mi][j][mj] << endl;
          }
        }
      }
    }
  }

  return;
}

double ThermalStrain::energy() {
  int nss = ns_ + 2;
  double utot = 0;
  double sum = 0;

  ///
  /// Do global matrix multiply via small stiffness matrices, Ah_ = A * h_
  /// The long statement below correctly brings in all the terms from the
  /// global matrix A using only the small stiffness matrices.
  ///

  for (int m3 = 0; m3 < 3; m3++) {
    for (int m = 0; m < nss; m++) {
      gb_[m][m3] = 0.0;
    }
  }

  for (int j = 0; j < 3; j++) {
    for (int n = 0; n < 3; n++) {
      for (int m = 0; m < ns_; m++) {
        gb_[m][j] += u_[ib_[m][0]][n] * (dk_[pix_[ib_[m][26]]][0][j][3][n] +
                                         dk_[pix_[ib_[m][6]]][1][j][2][n] +
                                         dk_[pix_[ib_[m][24]]][4][j][7][n] +
                                         dk_[pix_[ib_[m][14]]][5][j][6][n]) +
                     u_[ib_[m][1]][n] * (dk_[pix_[ib_[m][26]]][0][j][2][n] +
                                         dk_[pix_[ib_[m][24]]][4][j][6][n]) +
                     u_[ib_[m][2]][n] * (dk_[pix_[ib_[m][26]]][0][j][1][n] +
                                         dk_[pix_[ib_[m][4]]][3][j][2][n] +
                                         dk_[pix_[ib_[m][12]]][7][j][6][n] +
                                         dk_[pix_[ib_[m][24]]][4][j][5][n]) +
                     u_[ib_[m][3]][n] * (dk_[pix_[ib_[m][4]]][3][j][1][n] +
                                         dk_[pix_[ib_[m][12]]][7][j][5][n]) +
                     u_[ib_[m][4]][n] * (dk_[pix_[ib_[m][5]]][2][j][1][n] +
                                         dk_[pix_[ib_[m][4]]][3][j][0][n] +
                                         dk_[pix_[ib_[m][13]]][6][j][5][n] +
                                         dk_[pix_[ib_[m][12]]][7][j][4][n]) +
                     u_[ib_[m][5]][n] * (dk_[pix_[ib_[m][5]]][2][j][0][n] +
                                         dk_[pix_[ib_[m][13]]][6][j][4][n]) +
                     u_[ib_[m][6]][n] * (dk_[pix_[ib_[m][5]]][2][j][3][n] +
                                         dk_[pix_[ib_[m][6]]][1][j][0][n] +
                                         dk_[pix_[ib_[m][13]]][6][j][7][n] +
                                         dk_[pix_[ib_[m][14]]][5][j][4][n]) +
                     u_[ib_[m][7]][n] * (dk_[pix_[ib_[m][6]]][1][j][3][n] +
                                         dk_[pix_[ib_[m][14]]][5][j][7][n]) +
                     u_[ib_[m][8]][n] * (dk_[pix_[ib_[m][24]]][4][j][3][n] +
                                         dk_[pix_[ib_[m][14]]][5][j][2][n]) +
                     u_[ib_[m][9]][n] * (dk_[pix_[ib_[m][24]]][4][j][2][n]) +
                     u_[ib_[m][10]][n] * (dk_[pix_[ib_[m][12]]][7][j][2][n] +
                                          dk_[pix_[ib_[m][24]]][4][j][1][n]) +
                     u_[ib_[m][11]][n] * (dk_[pix_[ib_[m][12]]][7][j][1][n]) +
                     u_[ib_[m][12]][n] * (dk_[pix_[ib_[m][12]]][7][j][0][n] +
                                          dk_[pix_[ib_[m][13]]][6][j][1][n]) +
                     u_[ib_[m][13]][n] * (dk_[pix_[ib_[m][13]]][6][j][0][n]) +
                     u_[ib_[m][14]][n] * (dk_[pix_[ib_[m][13]]][6][j][3][n] +
                                          dk_[pix_[ib_[m][14]]][5][j][0][n]) +
                     u_[ib_[m][15]][n] * (dk_[pix_[ib_[m][14]]][5][j][3][n]) +
                     u_[ib_[m][16]][n] * (dk_[pix_[ib_[m][26]]][0][j][7][n] +
                                          dk_[pix_[ib_[m][6]]][1][j][6][n]) +
                     u_[ib_[m][17]][n] * (dk_[pix_[ib_[m][26]]][0][j][6][n]) +
                     u_[ib_[m][18]][n] * (dk_[pix_[ib_[m][26]]][0][j][5][n] +
                                          dk_[pix_[ib_[m][4]]][3][j][6][n]) +
                     u_[ib_[m][19]][n] * (dk_[pix_[ib_[m][4]]][3][j][5][n]) +
                     u_[ib_[m][20]][n] * (dk_[pix_[ib_[m][4]]][3][j][4][n] +
                                          dk_[pix_[ib_[m][5]]][2][j][5][n]) +
                     u_[ib_[m][21]][n] * (dk_[pix_[ib_[m][5]]][2][j][4][n]) +
                     u_[ib_[m][22]][n] * (dk_[pix_[ib_[m][5]]][2][j][7][n] +
                                          dk_[pix_[ib_[m][6]]][1][j][4][n]) +
                     u_[ib_[m][23]][n] * (dk_[pix_[ib_[m][6]]][1][j][7][n]) +
                     u_[ib_[m][24]][n] * (dk_[pix_[ib_[m][13]]][6][j][2][n] +
                                          dk_[pix_[ib_[m][12]]][7][j][3][n] +
                                          dk_[pix_[ib_[m][14]]][5][j][1][n] +
                                          dk_[pix_[ib_[m][24]]][4][j][0][n]) +
                     u_[ib_[m][25]][n] * (dk_[pix_[ib_[m][5]]][2][j][6][n] +
                                          dk_[pix_[ib_[m][4]]][3][j][7][n] +
                                          dk_[pix_[ib_[m][26]]][0][j][4][n] +
                                          dk_[pix_[ib_[m][6]]][1][j][5][n]) +
                     u_[ib_[m][26]][n] * (dk_[pix_[ib_[m][26]]][0][j][0][n] +
                                          dk_[pix_[ib_[m][6]]][1][j][1][n] +
                                          dk_[pix_[ib_[m][5]]][2][j][2][n] +
                                          dk_[pix_[ib_[m][4]]][3][j][3][n] +
                                          dk_[pix_[ib_[m][24]]][4][j][4][n] +
                                          dk_[pix_[ib_[m][14]]][5][j][5][n] +
                                          dk_[pix_[ib_[m][13]]][6][j][6][n] +
                                          dk_[pix_[ib_[m][12]]][7][j][7][n]);
      }
    }
  }

  utot = 0.0;
  for (int m3 = 0; m3 < 3; m3++) {
    for (int m = 0; m < ns_; m++) {
      utot += 0.5 * u_[m][m3] * gb_[m][m3] + b_[m][m3] * u_[m][m3];

      ///
      /// This is the gradient of energy with respect to normal displacements.
      ///

      gb_[m][m3] += b_[m][m3];
    }
  }

  ///
  /// Compute "constant" macrostrain energy term.
  ///

  C_ = 0.0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int mi = 0; mi < 3; mi++) {
        for (int mj = 0; mj < 3; mj++) {
          C_ += 0.5 * u_[ns_ + i][mi] * zcon_[i][mi][j][mj] * u_[ns_ + j][mj];
        }
      }
    }
  }

  utot += C_;

  ///
  /// Add in constant term from thermal energy, Y.
  ///

  utot += Y_;

  ///
  /// Add in linear term in thermal energy.
  ///

  for (int m3 = 0; m3 < 3; m3++) {
    for (int m = 0; m < nss; m++) {
      utot += T_[m][m3] * u_[m][m3];
    }
  }

  ///
  /// Compute the gradient with respect to macrostrains;
  /// put in piece from first derivative of zcon quadratic term.
  ///

  for (int i = 0; i < 2; i++) {
    for (int mi = 0; mi < 3; mi++) {
      sum = 0.0;
      for (int j = 0; j < 2; j++) {
        for (int mj = 0; mj < 3; mj++) {
          sum += zcon_[i][mi][j][mj] * u_[ns_ + j][mj];
        }
      }
      gb_[ns_ + i][mi] = sum;
    }
  }

  ///
  /// Add in piece of gradient, for displacements as well as macrostrains,
  /// that come from linear term in thermal energy.
  ///

  for (int m3 = 0; m3 < 3; m3++) {
    for (int m = 0; m < nss; m++) {
      gb_[m][m3] += T_[m][m3];
    }
  }

  ///
  /// Generate part that comes from b*u term; do this by calling
  /// `b_` generation with appropriate macrostrain set to 1 to
  /// get that partial derivative. Just use bgrad (taken from femat),
  /// skip `dk_` and `zcon_` part.
  ///

  for (int ii = 0; ii < 6; ii++) {
    sum = 0.0;
    if (ii == 0) {
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns_; m++) {
          sum += u_[m][m3] * b0_[m][m3];
        }
      }
      gb_[ns_][0] += sum;
    }
    if (ii == 1) {
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns_; m++) {
          sum += u_[m][m3] * b1_[m][m3];
        }
      }
      gb_[ns_][1] += sum;
    }
    if (ii == 2) {
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns_; m++) {
          sum += u_[m][m3] * b2_[m][m3];
        }
      }
      gb_[ns_][2] += sum;
    }
    if (ii == 3) {
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns_; m++) {
          sum += u_[m][m3] * b3_[m][m3];
        }
      }
      gb_[ns_ + 1][0] += sum;
    }
    if (ii == 4) {
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns_; m++) {
          sum += u_[m][m3] * b4_[m][m3];
        }
      }
      gb_[ns_ + 1][1] += sum;
    }
    if (ii == 5) {
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < ns_; m++) {
          sum += u_[m][m3] * b5_[m][m3];
        }
      }
      gb_[ns_ + 1][2] += sum;
    }
  }

  return utot;
}

// int ThermalStrain::dembx(int ns, double gg, int ldemb, int kkk) {
int ThermalStrain::dembx(int ldemb, int kkk) {
  double lambda, gamma;
  double hAh, gglast;
  int Lstep;
  int nss = ns_ + 2;

  ///
  /// Initialize the conjugate direction vector on first call to dembx only.
  /// For calls to dembx after the first, we want to continue using the value
  /// of h_ determined in the previous call.
  ///

  if (kkk == 0) {
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < nss; m++) {
        h_[m][m3] = gb_[m][m3];
      }
    }
  }

  ///
  /// Lstep counts the number of conjugate gradient steps taken
  /// in each call to dembx.
  ///

  Lstep = 0;

  for (int ijk = 0; ijk < ldemb; ijk++) {
    Lstep++;

    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < nss; m++) {
        Ah_[m][m3] = 0.0;
      }
    }

    ///
    /// Do global matrix multiply via small stiffness matrices, Ah_ = A * h_
    /// The long statement below correctly brings in all the terms from the
    /// global matrix A using only the small stiffness matrices.
    ///

    for (int j = 0; j < 3; j++) {
      for (int n = 0; n < 3; n++) {
        for (int m = 0; m < ns_; m++) {
          Ah_[m][j] +=
              h_[ib_[m][0]][n] * (dk_[pix_[ib_[m][26]]][0][j][3][n] +
                                  dk_[pix_[ib_[m][6]]][1][j][2][n] +
                                  dk_[pix_[ib_[m][24]]][4][j][7][n] +
                                  dk_[pix_[ib_[m][14]]][5][j][6][n]) +
              h_[ib_[m][1]][n] * (dk_[pix_[ib_[m][26]]][0][j][2][n] +
                                  dk_[pix_[ib_[m][24]]][4][j][6][n]) +
              h_[ib_[m][2]][n] * (dk_[pix_[ib_[m][26]]][0][j][1][n] +
                                  dk_[pix_[ib_[m][4]]][3][j][2][n] +
                                  dk_[pix_[ib_[m][12]]][7][j][6][n] +
                                  dk_[pix_[ib_[m][24]]][4][j][5][n]) +
              h_[ib_[m][3]][n] * (dk_[pix_[ib_[m][4]]][3][j][1][n] +
                                  dk_[pix_[ib_[m][12]]][7][j][5][n]) +
              h_[ib_[m][4]][n] * (dk_[pix_[ib_[m][5]]][2][j][1][n] +
                                  dk_[pix_[ib_[m][4]]][3][j][0][n] +
                                  dk_[pix_[ib_[m][13]]][6][j][5][n] +
                                  dk_[pix_[ib_[m][12]]][7][j][4][n]) +
              h_[ib_[m][5]][n] * (dk_[pix_[ib_[m][5]]][2][j][0][n] +
                                  dk_[pix_[ib_[m][13]]][6][j][4][n]) +
              h_[ib_[m][6]][n] * (dk_[pix_[ib_[m][5]]][2][j][3][n] +
                                  dk_[pix_[ib_[m][6]]][1][j][0][n] +
                                  dk_[pix_[ib_[m][13]]][6][j][7][n] +
                                  dk_[pix_[ib_[m][14]]][5][j][4][n]) +
              h_[ib_[m][7]][n] * (dk_[pix_[ib_[m][6]]][1][j][3][n] +
                                  dk_[pix_[ib_[m][14]]][5][j][7][n]) +
              h_[ib_[m][8]][n] * (dk_[pix_[ib_[m][24]]][4][j][3][n] +
                                  dk_[pix_[ib_[m][14]]][5][j][2][n]) +
              h_[ib_[m][9]][n] * (dk_[pix_[ib_[m][24]]][4][j][2][n]) +
              h_[ib_[m][10]][n] * (dk_[pix_[ib_[m][12]]][7][j][2][n] +
                                   dk_[pix_[ib_[m][24]]][4][j][1][n]) +
              h_[ib_[m][11]][n] * (dk_[pix_[ib_[m][12]]][7][j][1][n]) +
              +h_[ib_[m][12]][n] * (dk_[pix_[ib_[m][12]]][7][j][0][n] +
                                    dk_[pix_[ib_[m][13]]][6][j][1][n]) +
              h_[ib_[m][13]][n] * (dk_[pix_[ib_[m][13]]][6][j][0][n]) +
              h_[ib_[m][14]][n] * (dk_[pix_[ib_[m][13]]][6][j][3][n] +
                                   dk_[pix_[ib_[m][14]]][5][j][0][n]) +
              h_[ib_[m][15]][n] * (dk_[pix_[ib_[m][14]]][5][j][3][n]) +
              h_[ib_[m][16]][n] * (dk_[pix_[ib_[m][26]]][0][j][7][n] +
                                   dk_[pix_[ib_[m][6]]][1][j][6][n]) +
              h_[ib_[m][17]][n] * (dk_[pix_[ib_[m][26]]][0][j][6][n]) +
              h_[ib_[m][18]][n] * (dk_[pix_[ib_[m][26]]][0][j][5][n] +
                                   dk_[pix_[ib_[m][4]]][3][j][6][n]) +
              h_[ib_[m][19]][n] * (dk_[pix_[ib_[m][4]]][3][j][5][n]) +
              h_[ib_[m][20]][n] * (dk_[pix_[ib_[m][4]]][3][j][4][n] +
                                   dk_[pix_[ib_[m][5]]][2][j][5][n]) +
              h_[ib_[m][21]][n] * (dk_[pix_[ib_[m][5]]][2][j][4][n]) +
              h_[ib_[m][22]][n] * (dk_[pix_[ib_[m][5]]][2][j][7][n] +
                                   dk_[pix_[ib_[m][6]]][1][j][4][n]) +
              h_[ib_[m][23]][n] * (dk_[pix_[ib_[m][6]]][1][j][7][n]) +
              h_[ib_[m][24]][n] * (dk_[pix_[ib_[m][13]]][6][j][2][n] +
                                   dk_[pix_[ib_[m][12]]][7][j][3][n] +
                                   dk_[pix_[ib_[m][14]]][5][j][1][n] +
                                   dk_[pix_[ib_[m][24]]][4][j][0][n]) +
              h_[ib_[m][25]][n] * (dk_[pix_[ib_[m][5]]][2][j][6][n] +
                                   dk_[pix_[ib_[m][4]]][3][j][7][n] +
                                   dk_[pix_[ib_[m][26]]][0][j][4][n] +
                                   dk_[pix_[ib_[m][6]]][1][j][5][n]) +
              h_[ib_[m][26]][n] * (dk_[pix_[ib_[m][26]]][0][j][0][n] +
                                   dk_[pix_[ib_[m][6]]][1][j][1][n] +
                                   dk_[pix_[ib_[m][5]]][2][j][2][n] +
                                   dk_[pix_[ib_[m][4]]][3][j][3][n] +
                                   dk_[pix_[ib_[m][24]]][4][j][4][n] +
                                   dk_[pix_[ib_[m][14]]][5][j][5][n] +
                                   dk_[pix_[ib_[m][13]]][6][j][6][n] +
                                   dk_[pix_[ib_[m][12]]][7][j][7][n]);
        }
      }
    }

    ///
    /// The above accurately gives the second derivative matrix with respect to
    /// nodal displacements, but fails to give the 2nd derivative terms that
    /// include the macrostrains [du d(strain) and d(strain)d(strain)].
    /// Use repeated calls to bgrad to generate mixed 2nd derivatives terms,
    /// plus use zcon_ in order to correct the matrix multiply and correctly
    /// bring in macrostrain terms (see manual, Sec. 2.4).
    ///

    for (int ii = 0; ii < 6; ii++) {

      ///
      /// Fill in terms from matrix multiply right hand sides, 1 to ns_.
      ///

      for (int m = 0; m < ns_; m++) {
        for (int m1 = 0; m1 < 3; m1++) {
          if (ii == 0)
            Ah_[m][m1] += b0_[m][m1] * h_[ns_][0];
          if (ii == 1)
            Ah_[m][m1] += b1_[m][m1] * h_[ns_][1];
          if (ii == 2)
            Ah_[m][m1] += b2_[m][m1] * h_[ns_][2];
          if (ii == 3)
            Ah_[m][m1] += b3_[m][m1] * h_[ns_ + 1][0];
          if (ii == 4)
            Ah_[m][m1] += b4_[m][m1] * h_[ns_ + 1][1];
          if (ii == 5)
            Ah_[m][m1] += b5_[m][m1] * h_[ns_ + 1][2];
        }
      }

      // now do across bottom, 1 to ns_.
      for (int m = 0; m < ns_; m++) {
        if (ii == 0)
          Ah_[ns_][0] += b0_[m][0] * h_[m][0] + b0_[m][1] * h_[m][1] +
                         b0_[m][2] * h_[m][2];
        if (ii == 1)
          Ah_[ns_][1] += b1_[m][0] * h_[m][0] + b1_[m][1] * h_[m][1] +
                         b1_[m][2] * h_[m][2];
        if (ii == 2)
          Ah_[ns_][2] += b2_[m][0] * h_[m][0] + b2_[m][1] * h_[m][1] +
                         b2_[m][2] * h_[m][2];
        if (ii == 3)
          Ah_[ns_ + 1][0] += b3_[m][0] * h_[m][0] + b3_[m][1] * h_[m][1] +
                             b3_[m][2] * h_[m][2];
        if (ii == 4)
          Ah_[ns_ + 1][1] += b4_[m][0] * h_[m][0] + b4_[m][1] * h_[m][1] +
                             b4_[m][2] * h_[m][2];
        if (ii == 5)
          Ah_[ns_ + 1][2] += b5_[m][0] * h_[m][0] + b5_[m][1] * h_[m][1] +
                             b5_[m][2] * h_[m][2];
      }
      // now do right hand corner terms, ns to (ns+1).
      for (int m = 0; m < 2; m++) {
        for (int m1 = 0; m1 < 3; m1++) {
          if (ii == 0)
            Ah_[ns_][0] += zcon_[0][0][m][m1] * h_[ns_ + m][m1];
          if (ii == 1)
            Ah_[ns_][1] += zcon_[0][1][m][m1] * h_[ns_ + m][m1];
          if (ii == 2)
            Ah_[ns_][2] += zcon_[0][2][m][m1] * h_[ns_ + m][m1];
          if (ii == 3)
            Ah_[ns_ + 1][0] += zcon_[1][0][m][m1] * h_[ns_ + m][m1];
          if (ii == 4)
            Ah_[ns_ + 1][1] += zcon_[1][1][m][m1] * h_[ns_ + m][m1];
          if (ii == 5)
            Ah_[ns_ + 1][2] += zcon_[1][2][m][m1] * h_[ns_ + m][m1];
        }
      }
    }

    hAh = 0.0;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < nss; m++) {
        hAh += h_[m][m3] * Ah_[m][m3];
      }
    }

    lambda = gg_ / hAh;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < nss; m++) {
        u_[m][m3] -= lambda * h_[m][m3];
        gb_[m][m3] -= lambda * Ah_[m][m3];
      }
    }

    /*
    string outfilename = "displacement.dat";
    ofstream out1(outfilename.c_str(),ios::app);
    out1 << "after one step of relaxation..." << endl;
    for (int k = 0; k < nz_; k++) {
      int m = nx_ * ny_ * k + nx_ * 50 + 50;
      out1 << u_[m][0] << endl;
    }
    out1 << "macrostrain:" << endl;
    out1 << u_[ns_][0] << endl;
    out1 << u_[ns_][1] << endl;
    out1 << u_[ns_][2] << endl;
    out1 << u_[ns_+1][0] << endl;
    out1 << u_[ns_+1][1] << endl;
    out1 << u_[ns_+1][2] << endl;
    */

    gglast = gg_;
    gg_ = 0.0;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < nss; m++) {
        gg_ += gb_[m][m3] * gb_[m][m3];
      }
    }
    if (gg_ < gtest_) {
      return Lstep;
    } else {
      gamma = gg_ / gglast;
      for (int m3 = 0; m3 < 3; m3++) {
        for (int m = 0; m < nss; m++) {
          h_[m][m3] = gb_[m][m3] + gamma * h_[m][m3];
        }
      }
    }
  }

  return Lstep;
}

void ThermalStrain::relax(int cyc, int kmax) {
  int ldemb = 150, ltot = 0;
  // int iskip;
  double utot;
  int Lstep;
  int nss = ns_ + 2;

  // iskip = 1;
  // femat(nx_, ny_, nz_, ns_, nphase_, iskip);
  femat(1);

  // utot = energy(nx_, ny_, nz_, ns_);
  utot = energy();

  gg_ = 0.0;
  for (int m3 = 0; m3 < 3; m3++) {
    for (int m = 0; m < nss; m++) {
      gg_ += gb_[m][m3] * gb_[m][m3];
    }
  }

  if (verbose_) {
    std::clog << "ThermalStrain::relax gg_ = " << gg_ << endl;
    std::clog.flush();
  }

  /*
  string outfilename = "displacement.dat";
  ofstream out(outfilename.c_str(),ios::app);
  out << "displacement got from previous step:" << endl;
  for (int k = 0; k < nz_; k++) {
    int m = nx_ * ny_ * k + nx_ * 50 + 50;
    out << u_[m][0] << endl;
  }
  out << "macrostrain:" << endl;
  out << u_[ns_][0] << endl;
  out << u_[ns_][1] << endl;
  out << u_[ns_][2] << endl;
  out << u_[ns_+1][0] << endl;
  out << u_[ns_+1][1] << endl;
  out << u_[ns_+1][2] << endl;
  */

  std::clog << "    ThermalStrain::relax - ini - cyc = " << cyc
            << "   kmax_ = " << kmax << "   gtest_ = " << gtest_
            << "   =>   WAIT..." << endl;
  int kkk;
  for (kkk = 0; kkk < kmax; kkk++) {
    ///
    /// Call dembx to implement conjugate gradient routine.
    ///

    // Lstep = dembx(ns_, gg_, ldemb, kkk);
    Lstep = dembx(ldemb, kkk);
    ltot += Lstep;

    ///
    /// Call energy to compute energy after dembx call. If gg_ < gtest_,
    /// this will be the final energy. If gg_ is still larger than gtest_, then
    /// this will give an intermediate energy with which to check how the
    /// relaxation process is coming along. The call to energy does not change
    /// the gradient or the value of gg_. Need to first call femat to update the
    /// vector b_, as the value of the components of b_ depend on the
    /// macrostrains.
    ///

    // iskip = 1;
    // femat(nx_, ny_, nz_, ns_, nphase_, iskip);
    femat(1);

    utot = energy();

    if (verbose_) {
      std::clog << "ThermalStrain::relax energy = " << utot << " gg_ = " << gg_
                << " ltot = " << ltot << endl;
      std::clog.flush();
    }

    ///
    /// If relaxation process is finished, jump out of loop.
    ///

    if (gg_ >= gtest_) {
      ///
      /// Output stresses, strains, and macrostrains as an additional aid in
      /// judging how well the relaxation process is proceeding.
      ///

      stress();

      if (verbose_) {
        std::clog << "ThermalStrain::relax stresses: xx, "
                  << "yy, zz, xz, yz, xy " << strxx_ << " " << stryy_ << " "
                  << strzz_ << " " << strxz_ << " " << stryz_ << " " << strxy_
                  << endl;
        std::clog << "ThermalStrain::relax strains: xx, "
                  << "yy, zz, xz, yz, xy " << sxx_ << " " << syy_ << " " << szz_
                  << " " << sxz_ << " " << syz_ << " " << sxy_ << endl;
        std::clog << "ThermalStrain::relax macrostrains "
                  << "in same order " << u_[ns_][0] << " " << u_[ns_][1] << " "
                  << u_[ns_][2] << " " << u_[ns_ + 1][0] << " "
                  << u_[ns_ + 1][1] << " " << u_[ns_ + 1][2] << " " << endl;
        std::clog << "ThermalStrain::relax avg = "
                  << (u_[ns_][0] + u_[ns_][1] + u_[ns_][2]) / 3.0 << endl;
        std::clog << "ThermalStrain::relax energy = " << utot
                  << " gg_ = " << gg_ << " ltot = " << ltot << endl;
        std::clog.flush();
      }

    } else {
      break;
    }
  }

  stress();

  if (verbose_) {
    std::clog << "ThermalStrain::relax stresses: xx, "
              << "yy, zz, xz, yz, xy " << strxx_ << " " << stryy_ << " "
              << strzz_ << " " << strxz_ << " " << stryz_ << " " << strxy_
              << endl;
    std::clog << "ThermalStrain::relax strains: xx, "
              << "yy, zz, xz, yz, xy " << sxx_ << " " << syy_ << " " << szz_
              << " " << sxz_ << " " << syz_ << " " << sxy_ << endl;

    std::clog << "ThermalStrain::relax macrostrains in same order "
              << u_[ns_][0] << " " << u_[ns_][1] << " " << u_[ns_][2] << " "
              << u_[ns_ + 1][0] << " " << u_[ns_ + 1][1] << " "
              << u_[ns_ + 1][2] << " " << endl;
    std::clog << "ThermalStrain::relax avg = "
              << (u_[ns_][0] + u_[ns_][1] + u_[ns_][2]) / 3.0 << endl;
    std::clog.flush();
  }

  std::ostringstream ostr2;
  ostr2 << ltot;
  // string step(ostr2.str());

  // writeStress(time,step,0);
  // writeStrainEngy(time,step);

  std::ofstream outstr("macrostrains.dat");
  outstr << u_[ns_][0];
  outstr.close();

  std::clog << "    ThermalStrain::relax - fin - cyc = " << cyc
            << "   kkk = " << kkk << "   =>   gg_ = " << gg_
            << "   =>   gg <= gtest_ ! (gtest_ = " << gtest_
            << ") &  macrostrains.dat file has been updated" << endl;

  return;
}

// void ThermalStrain::stress(int nx, int ny, int nz, int ns) {
void ThermalStrain::stress() {
  double dndx[8], dndy[8], dndz[8];
  double es[6][8][3];
  double uu[8][3];
  double str11, str22, str33, str13, str23, str12;
  double s11, s22, s33, s13, s23, s12;
  double exx, eyy, ezz, exz, eyz, exy;

  int nx1 = nx_ - 1;
  int ny1 = ny_ - 1;
  int nz1 = nz_ - 1;
  double nx_dbl = nx_;
  double ny_dbl = ny_;
  double nz_dbl = nz_;
  int m;

  exx = u_[ns_][0];
  eyy = u_[ns_][1];
  ezz = u_[ns_][2];
  exz = u_[ns_ + 1][0];
  eyz = u_[ns_ + 1][1];
  exy = u_[ns_ + 1][2];

  ///
  /// Set up single element strain matrix.
  ///

  dndx[0] = (-0.25);
  dndx[1] = 0.25;
  dndx[2] = 0.25;
  dndx[3] = (-0.25);
  dndx[4] = (-0.25);
  dndx[5] = 0.25;
  dndx[6] = 0.25;
  dndx[7] = (-0.25);
  dndy[0] = (-0.25);
  dndy[1] = (-0.25);
  dndy[2] = 0.25;
  dndy[3] = 0.25;
  dndy[4] = (-0.25);
  dndy[5] = (-0.25);
  dndy[6] = 0.25;
  dndy[7] = 0.25;
  dndz[0] = (-0.25);
  dndz[1] = (-0.25);
  dndz[2] = (-0.25);
  dndz[3] = (-0.25);
  dndz[4] = 0.25;
  dndz[5] = 0.25;
  dndz[6] = 0.25;
  dndz[7] = 0.25;

  ///
  /// Build average strain matrix, follows code in femat,
  /// but for average strain over an element, not the strain at a point.
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
  /// Compute average stresses and strains in each pixel.
  ///

  sxx_ = 0.0;
  syy_ = 0.0;
  szz_ = 0.0;
  sxz_ = 0.0;
  syz_ = 0.0;
  sxy_ = 0.0;
  strxx_ = 0.0;
  stryy_ = 0.0;
  strzz_ = 0.0;
  strxz_ = 0.0;
  stryz_ = 0.0;
  strxy_ = 0.0;

  double exx_nx_dbl = exx * nx_dbl;
  double exy_nx_dbl = exy * nx_dbl;
  double exz_nx_dbl = exz * nx_dbl;

  double exy_ny_dbl = exy * ny_dbl;
  double eyy_ny_dbl = eyy * ny_dbl;
  double eyz_ny_dbl = eyz * ny_dbl;

  double exz_nz_dbl = exz * nz_dbl;
  double eyz_nz_dbl = eyz * nz_dbl;
  double ezz_nz_dbl = ezz * nz_dbl;

  int i, j, k, nxy_k, nx_j;
  for (k = 0; k < nz_; k++) {
    nxy_k = nxy_ * k;
    for (j = 0; j < ny_; j++) {
      nx_j = nxy_k + nx_ * j;
      for (i = 0; i < nx_; i++) {
        // m = nxy_ * k + nx_ * j + i;
        m = nx_j + i;

        ///
        /// Load in elements of 8-vector using periodic boundary conditions.
        ///

        for (int mm = 0; mm < 3; mm++) {
          uu[0][mm] = u_[m][mm];
          uu[1][mm] = u_[ib_[m][2]][mm];
          uu[2][mm] = u_[ib_[m][1]][mm];
          uu[3][mm] = u_[ib_[m][0]][mm];
          uu[4][mm] = u_[ib_[m][25]][mm];
          uu[5][mm] = u_[ib_[m][18]][mm];
          uu[6][mm] = u_[ib_[m][17]][mm];
          uu[7][mm] = u_[ib_[m][16]][mm];
        }

        ///
        /// Correct for periodic boundary conditions, some displacements are
        /// wrong for an element on a periodic boundary. Since they come from an
        /// opposite face, need to put in applied strain to correct them.
        ///

        // nx1 = nx - 1
        if (i == nx1) {
          uu[1][0] += exx_nx_dbl; // exx_nx_dbl = exx * nx_dbl;
          uu[1][1] += exy_nx_dbl;
          uu[1][2] += exz_nx_dbl;
          uu[2][0] += exx_nx_dbl;
          uu[2][1] += exy_nx_dbl;
          uu[2][2] += exz_nx_dbl;
          uu[5][0] += exx_nx_dbl;
          uu[5][1] += exy_nx_dbl;
          uu[5][2] += exz_nx_dbl;
          uu[6][0] += exx_nx_dbl;
          uu[6][1] += exy_nx_dbl;
          uu[6][2] += exz_nx_dbl;
        }
        // ny1 = ny - 1
        if (j == ny1) {
          uu[2][0] += exy_ny_dbl; // exy_ny_dbl = exy * ny_dbl;
          uu[2][1] += eyy_ny_dbl;
          uu[2][2] += eyz_ny_dbl;
          uu[3][0] += exy_ny_dbl;
          uu[3][1] += eyy_ny_dbl;
          uu[3][2] += eyz_ny_dbl;
          uu[6][0] += exy_ny_dbl;
          uu[6][1] += eyy_ny_dbl;
          uu[6][2] += eyz_ny_dbl;
          uu[7][0] += exy_ny_dbl;
          uu[7][1] += eyy_ny_dbl;
          uu[7][2] += eyz_ny_dbl;
        }
        // nz1 = nz - 1
        if (k == nz1) {
          uu[4][0] += exz_nz_dbl; // exz_nz_dbl = exz * nz_dbl;
          uu[4][1] += eyz_nz_dbl;
          uu[4][2] += ezz_nz_dbl;
          uu[5][0] += exz_nz_dbl;
          uu[5][1] += eyz_nz_dbl;
          uu[5][2] += ezz_nz_dbl;
          uu[6][0] += exz_nz_dbl;
          uu[6][1] += eyz_nz_dbl;
          uu[6][2] += ezz_nz_dbl;
          uu[7][0] += exz_nz_dbl;
          uu[7][1] += eyz_nz_dbl;
          uu[7][2] += ezz_nz_dbl;
        }

        ///
        /// Stresses and strains in an element.
        ///

        str11 = 0.0;
        str22 = 0.0;
        str33 = 0.0;
        str13 = 0.0;
        str23 = 0.0;
        str12 = 0.0;
        s11 = 0.0;
        s22 = 0.0;
        s33 = 0.0;
        s13 = 0.0;
        s23 = 0.0;
        s12 = 0.0;

        ///
        /// Compute average stress and strain tensor in each element
        /// First put thermal strain-induced stresses into stress tensor
        ///

        for (int n = 0; n < 6; n++) {
          str11 -= cmod_[pix_[m]][0][n] * eigen_[m][n];
          str22 -= cmod_[pix_[m]][1][n] * eigen_[m][n];
          str33 -= cmod_[pix_[m]][2][n] * eigen_[m][n];
          str13 -= cmod_[pix_[m]][3][n] * eigen_[m][n];
          str23 -= cmod_[pix_[m]][4][n] * eigen_[m][n];
          str12 -= cmod_[pix_[m]][5][n] * eigen_[m][n];
        }
        for (int n3 = 0; n3 < 3; n3++) {
          for (int n8 = 0; n8 < 8; n8++) {
            s11 += es[0][n8][n3] * uu[n8][n3];
            s22 += es[1][n8][n3] * uu[n8][n3];
            s33 += es[2][n8][n3] * uu[n8][n3];
            s13 += es[3][n8][n3] * uu[n8][n3];
            s23 += es[4][n8][n3] * uu[n8][n3];
            s12 += es[5][n8][n3] * uu[n8][n3];
            for (int n = 0; n < 6; n++) {

              ///
              /// Compute stresses in each element that include both non-thermal
              /// and thermal strains.
              ///

              str11 += cmod_[pix_[m]][0][n] * es[n][n8][n3] * uu[n8][n3];
              str22 += cmod_[pix_[m]][1][n] * es[n][n8][n3] * uu[n8][n3];
              str33 += cmod_[pix_[m]][2][n] * es[n][n8][n3] * uu[n8][n3];
              str13 += cmod_[pix_[m]][3][n] * es[n][n8][n3] * uu[n8][n3];
              str23 += cmod_[pix_[m]][4][n] * es[n][n8][n3] * uu[n8][n3];
              str12 += cmod_[pix_[m]][5][n] * es[n][n8][n3] * uu[n8][n3];
            }
          }
        }
        elestress_[m][0] = str11;
        elestress_[m][1] = str22;
        elestress_[m][2] = str33;
        elestress_[m][3] = str13;
        elestress_[m][4] = str23;
        elestress_[m][5] = str12;
        elestrain_[m][0] = s11;
        elestrain_[m][1] = s22;
        elestrain_[m][2] = s33;
        elestrain_[m][3] = s13;
        elestrain_[m][4] = s23;
        elestrain_[m][5] = s12;

        ///
        /// Compute the strain energy for each element.
        ///

        for (int ii = 0; ii < 6; ii++) {
          for (int jj = 0; jj < 6; jj++) {
            strainengy_[m] = 0.5 * elestrain_[m][ii] * cmod_[pix_[m]][ii][jj] *
                             elestrain_[m][jj];
          }
        }

        ///
        /// Sum local stresses and strains into global stresses and strains.
        ///

        strxx_ += str11;
        stryy_ += str22;
        strzz_ += str33;
        strxz_ += str13;
        stryz_ += str23;
        strxy_ += str12;
        sxx_ += s11;
        syy_ += s22;
        szz_ += s33;
        sxz_ += s13;
        syz_ += s23;
        sxy_ += s12;
      }
    }
  }

  ///
  /// Volume average global stresses and strain.
  ///

  double ns_dbl = ns_;
  strxx_ = strxx_ / ns_dbl;
  stryy_ = stryy_ / ns_dbl;
  strzz_ = strzz_ / ns_dbl;
  strxz_ = strxz_ / ns_dbl;
  stryz_ = stryz_ / ns_dbl;
  strxy_ = strxy_ / ns_dbl;

  sxx_ = sxx_ / ns_dbl;
  syy_ = syy_ / ns_dbl;
  szz_ = szz_ / ns_dbl;
  sxz_ = sxz_ / ns_dbl;
  syz_ = syz_ / ns_dbl;
  sxy_ = sxy_ / ns_dbl;

  ///
  /// Calculate the volume-averaged strain energy for each phase/DC
  /// (called ONCE after all element strain energies are computed)
  ///
  getAvgStrainengy();

  return;
}

void ThermalStrain::calc(int cyc, double time, std::vector<int> *p_vectPhId,
                         double exx, double eyy, double ezz, double exz,
                         double eyz, double exy) {
  int kmax = kmax_; // 60;
  // int iskip;
  // double utot;
  int m;

  ///
  /// Read in a microstructure in function setMicrostructure, and set up pix_[m]
  /// with the appropriate phase assignments.
  ///

  std::clog << endl
            << "  ThermalStrain::Calc - ini - cyc = " << cyc
            << "   time = " << time << endl;

  // setMicrostructure(fileName);
  setMicrostructure(p_vectPhId);

  ///
  /// Count and output the volume fractions of the different phases.
  ///

  assig();

#ifdef DEBUG
  for (int i = 0; i < nphase_; i++) {
    std::clog << "ThermalStrain::Calc Phase " << i
              << " bulk = " << phasemod_[i][0] << " shear = " << phasemod_[i][1]
              << endl;
  }

  for (int i = 0; i < nphase_; i++) {
    std::clog << "ThermalStrain::Calc Volume fraction of phase " << i << "  is "
              << prob_[i] << endl;
  }
  std::clog.flush();

  ///
  /// Output thermal strains for each phase.
  ///

  std::clog << "Thermal Strains" << endl;
  for (int i = 0; i < ns_; i++) {
    bool flag = false;
    for (int j = 0; j < 6; j++) {
      if (eigen_[i][j] != 0.0)
        flag = true;
    }
    if (flag)
      std::clog << "eigen of site[" << i << "]: " << eigen_[i][0] << "   "
                << eigen_[i][1] << "   " << eigen_[i][2] << "   "
                << eigen_[i][3] << "   " << eigen_[i][4] << "   "
                << eigen_[i][5] << "   " << "phase is: " << pix_[i] << endl;
  }
#endif

  if (isFirst_) {
#ifdef DEBUG
    std::clog << "ThermalStrain::Calc This is the first time "
              << "Calc to be called, "
              << "so u_ should be initialized." << endl;
#endif

    ///
    /// @note Set applied strains. Actual shear strain applied in is `exy_`,
    /// `exz_`, and `eyz_` as given in the statements below. The engineering
    /// shear strain, by which the shear modulus is usually defined, is twice
    /// these values.
    ///

    u_[ns_][0] = exx;
    u_[ns_][1] = eyy;
    u_[ns_][2] = ezz;
    u_[ns_ + 1][0] = exz;
    u_[ns_ + 1][1] = eyz;
    u_[ns_ + 1][2] = exy;

#ifdef DEBUG
    std::clog << "ThermalStrain::Calc Applied engineering strains" << endl;
    std::clog << "ThermalStrain::Calc exx,   eyy,   ezz,   exz,   eyz,   exy"
              << endl;
    std::clog << "ThermalStrain::Calc " << exx << "   " << eyy << "   " << ezz
              << "   " << exz << "   " << eyz << "   " << exy << endl;
    std::clog.flush();
#endif

    ///
    /// Apply homogeneous macroscopic strain as the initial condition
    /// to displacement variables.
    ///

    int i, j, k, nxy_k, nx_j;
    for (k = 0; k < nz_; k++) {
      nxy_k = nxy_ * k;
      for (j = 0; j < ny_; j++) {
        nx_j = nxy_k + nx_ * j;
        for (i = 0; i < nx_; i++) {
          // m = nxy_ * k + nx_ * j + i;
          m = nx_j + i;
          // double x = static_cast<double>(i);
          // double y = static_cast<double>(j);
          // double z = static_cast<double>(k);
          // u_[m][0] = x * u_[ns_][0] + y * u_[ns_ + 1][2] + z * u_[ns_ +
          // 1][0]; u_[m][1] = x * u_[ns_ + 1][2] + y * u_[ns_][1] + z * u_[ns_
          // + 1][1]; u_[m][2] = x * u_[ns_ + 1][0] + y * u_[ns_ + 1][1] + z *
          // u_[ns_][2];
          u_[m][0] = i * u_[ns_][0] + j * u_[ns_ + 1][2] + k * u_[ns_ + 1][0];
          u_[m][1] = i * u_[ns_ + 1][2] + j * u_[ns_][1] + k * u_[ns_ + 1][1];
          u_[m][2] = i * u_[ns_ + 1][0] + j * u_[ns_ + 1][1] + k * u_[ns_][2];
        }
      }
    }
    if (time > 0) {
#ifdef DEBUG
      std::clog << "ThermalStrain::Calc Trying to open "
                << "displacement.dat file to initialize u_" << endl;
      std::clog.flush();
#endif

      std::ifstream in("displacement.dat");
      if (!in) {
#ifdef DEBUG
        std::clog << "ThermalStrain::Calc Cannot find "
                  << "displacement.dat file." << endl;
        std::clog.flush();
#endif
      } else {
        double buff;
        // for (int k = 0; k < nz_; k++) {
        //   for (int j = 0; j < ny_; j++) {
        //     for (int i = 0; i < nx_; i++) {
        //       int m = nx_ * ny_ * k + nx_ * j + i;
        //       in >> buff;
        //       u_[m][0] = buff;
        //       in >> buff;
        //       u_[m][1] = buff;
        //       in >> buff;
        //       u_[m][2] = buff;
        //     }
        //   }
        // }
        for (int m = 0; m < ns_; m++) {
          in >> buff;
          u_[m][0] = buff;
          in >> buff;
          u_[m][1] = buff;
          in >> buff;
          u_[m][2] = buff;
        }
        in >> buff;
        u_[ns_][0] = buff;
        in >> buff;
        u_[ns_][1] = buff;
        in >> buff;
        u_[ns_][2] = buff;
        in >> buff;
        u_[ns_ + 1][0] = buff;
        in >> buff;
        u_[ns_ + 1][1] = buff;
        in >> buff;
        u_[ns_ + 1][2] = buff;
      }
    }

    isFirst_ = false;
  }

  ///
  /// Set up the finite element stiffness matrices, the constant, `C_`,
  /// the vector, `b_`, required for the energy. `b_` and `C_` depend on the
  /// macrostrains.
  /// When they are updated, the values of `b_` and `C_` are updated too via
  /// calling subroutine femat. Only compute the thermal strain terms the first
  /// time femat is calles, (iskip = 0) as they are unaffected by later changes
  /// (iskip = 1) in displacements and macrostrains. Compute initial value of
  /// gradient `gb_` and `gg_ = gb_ * gb_`.
  ///

  int iskip = 0;
  femat(iskip);
  /// Note: We need to call the energy function to calculate
  /// the global gradient `gb_`, but we don't need the return
  /// value from the function so we can just call it without assigning
  /// its return value to anything!

  energy();
  gg_ = 0.0;
  int nss = ns_ + 2;
  for (int m3 = 0; m3 < 3; m3++) {
    for (int m = 0; m < nss; m++) {
      gg_ += gb_[m][m3] * gb_[m][m3];
    }
  }

#ifdef DEBUG
  std::clog << "ThermalStrain::Calc energy = " << utot << " gg_ = " << gg_
            << endl;
  std::clog.flush();
#endif

  for (int m1 = 0; m1 < 3; m1++) {
    for (int m = 0; m < nss; m++) {
      b0_[m][m1] = 0.0;
      b1_[m][m1] = 0.0;
      b2_[m][m1] = 0.0;
      b3_[m][m1] = 0.0;
      b4_[m][m1] = 0.0;
      b5_[m][m1] = 0.0;
    }
  }

  ///
  /// Since `b0`, `b1`, `b2`, `b3`, `b4`, `b5` always keep constant,
  /// `bgrad` is only needed to be called for 6 times to get the values of them.
  ///

  double e11, e22, e33;
  double e13, e23, e12;
  for (int ii = 0; ii < 6; ii++) {
    e11 = e22 = e33 = 0.0;
    e13 = e23 = e12 = 0.0;
    if (ii == 0) {
      e11 = 1.0;
      bgrad(e11, e22, e33, e13, e23, e12);
    }
    if (ii == 1) {
      e22 = 1.0;
      bgrad(e11, e22, e33, e13, e23, e12);
    }
    if (ii == 2) {
      e33 = 1.0;
      bgrad(e11, e22, e33, e13, e23, e12);
    }
    if (ii == 3) {
      e13 = 1.0;
      bgrad(e11, e22, e33, e13, e23, e12);
    }
    if (ii == 4) {
      e23 = 1.0;
      bgrad(e11, e22, e33, e13, e23, e12);
    }
    if (ii == 5) {
      e12 = 1.0;
      bgrad(e11, e22, e33, e13, e23, e12);
    }
  }

  ///
  /// This is the main relaxation loop.
  /// @note `kmax` is the maximum number of times dembx will be called, with
  /// `ldemb` conjugate gradient steps performed during each call. The total
  /// number of conjugate gradient steps allowed for a given elastic computation
  /// is `kmax * ldemb`.
  ///

  int xlo, xhi, ylo, yhi, zlo, zhi;
  int halfbox = boxsize_ / 2;
  int expindex;
  vector<int> expcor;
  int x, y, z;

  for (int count = 0; count < 0; count++) {
#ifdef DEBUG
    std::clog << "ThermalStrain::Calc boxsize_ is: " << boxsize_ << endl;
    std::clog.flush();
#endif

    for (map<int, vector<int>>::iterator it = exp_.begin(); it != exp_.end();
         it++) {
      localgg_ = 0.0;
      expindex = it->first;
      expcor = it->second;
      xlo = expcor[0] - halfbox;
      xhi = expcor[0] + halfbox;
      ylo = expcor[1] - halfbox;
      yhi = expcor[1] + halfbox;
      zlo = expcor[2] - halfbox;
      zhi = expcor[2] + halfbox;

      femat(1);

      /// Note: We need to call the energy function to calculate
      /// the global gradient `gb_`, but we don't need the return
      /// value from the function so we can just call it without assigning
      /// its return value to anything!

      energy();

      for (int m3 = 0; m3 < 3; m3++) {
        for (int k = zlo; k <= zhi; k++) {
          for (int j = ylo; j <= yhi; j++) {
            for (int i = xlo; i <= xhi; i++) {
              x = i;
              y = j;
              z = k;
              if (x < 0)
                x += nx_;
              if (x >= nx_)
                x -= nx_;
              if (y < 0)
                y += ny_;
              if (y >= ny_)
                y -= ny_;
              if (z < 0)
                z += nz_;
              if (z >= nz_)
                z -= nz_;
              m = nxy_ * z + nx_ * y + x;
              localgg_ += gb_[m][m3] * gb_[m][m3];
            }
          }
        }
      }
#ifdef DEBUG
      std::clog << "ThermalStrain::Calc localgg_ is: " << localgg_ << endl;
      std::clog.flush();
#endif

      // localRelax(boxsize_, expcor[0], expcor[1], expcor[2], expindex);
      localRelax(xlo, xhi, ylo, yhi, zlo, zhi, expcor[0], expcor[1], expcor[2],
                 expindex);
    }
  }

  relax(cyc, kmax);

  string outfilename = "displacement.dat";
  std::ofstream out(outfilename.c_str());
  // for (int k = 0; k < nz_; k++) {
  //   for (int j = 0; j < ny_; j++) {
  //     for (int i = 0; i < nx_; i++) {
  //       m = nxy_ * k + nx_ * j + i;
  //       out << u_[m][0] << " " << u_[m][1] << " " << u_[m][2] << endl;
  //     }
  //   }
  // }
  for (int m = 0; m < ns_; m++) {
    out << u_[m][0] << " " << u_[m][1] << " " << u_[m][2] << endl;
  }
  out << u_[ns_][0] << " " << u_[ns_][1] << " " << u_[ns_][2] << endl;
  out << u_[ns_ + 1][0] << " " << u_[ns_ + 1][1] << " " << u_[ns_ + 1][2]
      << endl;
  out.close();

  std::clog << "  ThermalStrain::Calc - fin - cyc = " << cyc
            << "   time = " << time
            << "   =>   displacement.dat file has been updated" << endl;
  std::clog.flush();

  return;
}

// int ThermalStrain::localDembx(int boxsize, int x, int y, int z, int
// localldemb,
//                               int kkk) {
int ThermalStrain::localDembx(int xlo, int xhi, int ylo, int yhi, int zlo,
                              int zhi, int localldemb, int kkk) {
  int nss = ns_ + 2;
  double lambda, gamma;
  double hAh, gglast;
  int Lstep;
  // int halfbox = static_cast<int>(boxsize / 2);
  bool flag =
      false; // determines whether bgrad function needs to be called or not

  // int xlo, xhi, ylo, yhi, zlo, zhi;
  // xlo = x - halfbox;
  // xhi = x + halfbox;
  // ylo = y - halfbox;
  // yhi = y + halfbox;
  // zlo = z - halfbox;
  // zhi = z + halfbox;
  int x, y, z;
  int m;

  if (kkk == 0) {
    for (int k = zlo; k <= zhi; k++) {
      for (int j = ylo; j <= yhi; j++) {
        for (int i = xlo; i <= xhi; i++) {
          for (int m3 = 0; m3 < 3; m3++) {
            x = i;
            y = j;
            z = k;
            if (x < 0) {
              x += nx_;
              flag = true;
            }
            if (x >= nx_) {
              x -= nx_;
              flag = true;
            }
            if (y < 0) {
              y += ny_;
              flag = true;
            }
            if (y >= ny_) {
              y -= ny_;
              flag = true;
            }
            if (z < 0) {
              z += nz_;
              flag = true;
            }
            if (z >= nz_) {
              z -= nz_;
              flag = true;
            }
            m = nxy_ * z + nx_ * y + x;
            h_[m][m3] = gb_[m][m3];
          }
        }
      }
    }
  }

  ///
  /// Lstep counts the number of conjugate gradient steps taken
  /// in each call to dembx.
  ///

  Lstep = 0;

  for (int ijk = 0; ijk < localldemb; ijk++) {
    Lstep++;

    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < nss; m++) {
        Ah_[m][m3] = 0.0;
      }
    }

    ///
    /// Do global matrix multiply via small stiffness matrices, `Ah_ = A * h_`
    /// The long statement below correctly brings in all the terms from the
    /// global matrix A using only the small stiffness matrices.
    ///

    for (int j = 0; j < 3; j++) {
      for (int n = 0; n < 3; n++) {
        for (int kk = zlo; kk <= zhi; kk++) {
          for (int jj = ylo; jj <= yhi; jj++) {
            for (int ii = xlo; ii <= xhi; ii++) {
              x = ii;
              y = jj;
              z = kk;
              if (x < 0)
                x += nx_;
              if (x >= nx_)
                x -= nx_;
              if (y < 0)
                y += ny_;
              if (y >= ny_)
                y -= ny_;
              if (z < 0)
                z += nz_;
              if (z >= nz_)
                z -= nz_;
              m = nxy_ * z + nx_ * y + x;
              Ah_[m][j] +=
                  h_[ib_[m][0]][n] * (dk_[pix_[ib_[m][26]]][0][j][3][n] +
                                      dk_[pix_[ib_[m][6]]][1][j][2][n] +
                                      dk_[pix_[ib_[m][24]]][4][j][7][n] +
                                      dk_[pix_[ib_[m][14]]][5][j][6][n]) +
                  h_[ib_[m][1]][n] * (dk_[pix_[ib_[m][26]]][0][j][2][n] +
                                      dk_[pix_[ib_[m][24]]][4][j][6][n]) +
                  h_[ib_[m][2]][n] * (dk_[pix_[ib_[m][26]]][0][j][1][n] +
                                      dk_[pix_[ib_[m][4]]][3][j][2][n] +
                                      dk_[pix_[ib_[m][12]]][7][j][6][n] +
                                      dk_[pix_[ib_[m][24]]][4][j][5][n]) +
                  h_[ib_[m][3]][n] * (dk_[pix_[ib_[m][4]]][3][j][1][n] +
                                      dk_[pix_[ib_[m][12]]][7][j][5][n]) +
                  h_[ib_[m][4]][n] * (dk_[pix_[ib_[m][5]]][2][j][1][n] +
                                      dk_[pix_[ib_[m][4]]][3][j][0][n] +
                                      dk_[pix_[ib_[m][13]]][6][j][5][n] +
                                      dk_[pix_[ib_[m][12]]][7][j][4][n]) +
                  h_[ib_[m][5]][n] * (dk_[pix_[ib_[m][5]]][2][j][0][n] +
                                      dk_[pix_[ib_[m][13]]][6][j][4][n]) +
                  h_[ib_[m][6]][n] * (dk_[pix_[ib_[m][5]]][2][j][3][n] +
                                      dk_[pix_[ib_[m][6]]][1][j][0][n] +
                                      dk_[pix_[ib_[m][13]]][6][j][7][n] +
                                      dk_[pix_[ib_[m][14]]][5][j][4][n]) +
                  h_[ib_[m][7]][n] * (dk_[pix_[ib_[m][6]]][1][j][3][n] +
                                      dk_[pix_[ib_[m][14]]][5][j][7][n]) +
                  h_[ib_[m][8]][n] * (dk_[pix_[ib_[m][24]]][4][j][3][n] +
                                      dk_[pix_[ib_[m][14]]][5][j][2][n]) +
                  h_[ib_[m][9]][n] * (dk_[pix_[ib_[m][24]]][4][j][2][n]) +
                  h_[ib_[m][10]][n] * (dk_[pix_[ib_[m][12]]][7][j][2][n] +
                                       dk_[pix_[ib_[m][24]]][4][j][1][n]) +
                  h_[ib_[m][11]][n] * (dk_[pix_[ib_[m][12]]][7][j][1][n]) +
                  +h_[ib_[m][12]][n] * (dk_[pix_[ib_[m][12]]][7][j][0][n] +
                                        dk_[pix_[ib_[m][13]]][6][j][1][n]) +
                  h_[ib_[m][13]][n] * (dk_[pix_[ib_[m][13]]][6][j][0][n]) +
                  h_[ib_[m][14]][n] * (dk_[pix_[ib_[m][13]]][6][j][3][n] +
                                       dk_[pix_[ib_[m][14]]][5][j][0][n]) +
                  h_[ib_[m][15]][n] * (dk_[pix_[ib_[m][14]]][5][j][3][n]) +
                  h_[ib_[m][16]][n] * (dk_[pix_[ib_[m][26]]][0][j][7][n] +
                                       dk_[pix_[ib_[m][6]]][1][j][6][n]) +
                  h_[ib_[m][17]][n] * (dk_[pix_[ib_[m][26]]][0][j][6][n]) +
                  h_[ib_[m][18]][n] * (dk_[pix_[ib_[m][26]]][0][j][5][n] +
                                       dk_[pix_[ib_[m][4]]][3][j][6][n]) +
                  h_[ib_[m][19]][n] * (dk_[pix_[ib_[m][4]]][3][j][5][n]) +
                  h_[ib_[m][20]][n] * (dk_[pix_[ib_[m][4]]][3][j][4][n] +
                                       dk_[pix_[ib_[m][5]]][2][j][5][n]) +
                  h_[ib_[m][21]][n] * (dk_[pix_[ib_[m][5]]][2][j][4][n]) +
                  h_[ib_[m][22]][n] * (dk_[pix_[ib_[m][5]]][2][j][7][n] +
                                       dk_[pix_[ib_[m][6]]][1][j][4][n]) +
                  h_[ib_[m][23]][n] * (dk_[pix_[ib_[m][6]]][1][j][7][n]) +
                  h_[ib_[m][24]][n] * (dk_[pix_[ib_[m][13]]][6][j][2][n] +
                                       dk_[pix_[ib_[m][12]]][7][j][3][n] +
                                       dk_[pix_[ib_[m][14]]][5][j][1][n] +
                                       dk_[pix_[ib_[m][24]]][4][j][0][n]) +
                  h_[ib_[m][25]][n] * (dk_[pix_[ib_[m][5]]][2][j][6][n] +
                                       dk_[pix_[ib_[m][4]]][3][j][7][n] +
                                       dk_[pix_[ib_[m][26]]][0][j][4][n] +
                                       dk_[pix_[ib_[m][6]]][1][j][5][n]) +
                  h_[ib_[m][26]][n] * (dk_[pix_[ib_[m][26]]][0][j][0][n] +
                                       dk_[pix_[ib_[m][6]]][1][j][1][n] +
                                       dk_[pix_[ib_[m][5]]][2][j][2][n] +
                                       dk_[pix_[ib_[m][4]]][3][j][3][n] +
                                       dk_[pix_[ib_[m][24]]][4][j][4][n] +
                                       dk_[pix_[ib_[m][14]]][5][j][5][n] +
                                       dk_[pix_[ib_[m][13]]][6][j][6][n] +
                                       dk_[pix_[ib_[m][12]]][7][j][7][n]);
            }
          }
        }
      }
    }

    ///
    /// The above accurately gives the second derivative matrix with respect to
    /// nodal displacements, but fails to give the 2nd derivative terms that
    /// include the macrostrains [du d(strain) and d(strain)d(strain)].
    /// Use repeated calls to bgrad to generate mixed 2nd derivatives terms,
    /// plus use zcon_ in order to correct the matrix multiply and correctly
    /// bring in macrostrain terms (see manual, Sec. 2.4).
    ///

    for (int ii = 0; ii < 6; ii++) {
      if (flag) {

        ///
        /// Fill in terms from matrix multiply right hand sides, 1 to ns_.
        ///

        for (int m1 = 0; m1 < 3; m1++) {
          for (int k = zlo; k <= zhi; k++) {
            for (int j = ylo; j <= yhi; j++) {
              for (int i = xlo; i <= xhi; i++) {
                x = i;
                y = j;
                z = k;
                if (x < 0)
                  x += nx_;
                if (x >= nx_)
                  x -= nx_;
                if (y < 0)
                  y += ny_;
                if (y >= ny_)
                  y -= ny_;
                if (z < 0)
                  z += nz_;
                if (z >= nz_)
                  z -= nz_;
                m = nxy_ * z + nx_ * y + x;
                if (ii == 0)
                  Ah_[m][m1] += b0_[m][m1] * h_[ns_][0];
                if (ii == 1)
                  Ah_[m][m1] += b1_[m][m1] * h_[ns_][1];
                if (ii == 2)
                  Ah_[m][m1] += b2_[m][m1] * h_[ns_][2];
                if (ii == 3)
                  Ah_[m][m1] += b3_[m][m1] * h_[ns_ + 1][0];
                if (ii == 4)
                  Ah_[m][m1] += b4_[m][m1] * h_[ns_ + 1][1];
                if (ii == 5)
                  Ah_[m][m1] += b5_[m][m1] * h_[ns_ + 1][2];
              }
            }
          }
        }
        // now do across bottom, 1 to ns_.
        for (int k = zlo; k <= zhi; k++) {
          for (int j = ylo; j <= yhi; j++) {
            for (int i = xlo; i <= xhi; i++) {
              x = i;
              y = j;
              z = k;
              if (x < 0)
                x += nx_;
              if (x >= nx_)
                x -= nx_;
              if (y < 0)
                y += ny_;
              if (y >= ny_)
                y -= ny_;
              if (z < 0)
                z += nz_;
              if (z >= nz_)
                z -= nz_;
              m = nxy_ * z + nx_ * y + x;
              if (ii == 0)
                Ah_[ns_][0] += b0_[m][0] * h_[m][0] + b0_[m][1] * h_[m][1] +
                               b0_[m][2] * h_[m][2];
              if (ii == 1)
                Ah_[ns_][1] += b1_[m][0] * h_[m][0] + b1_[m][1] * h_[m][1] +
                               b1_[m][2] * h_[m][2];
              if (ii == 2)
                Ah_[ns_][2] += b2_[m][0] * h_[m][0] + b2_[m][1] * h_[m][1] +
                               b2_[m][2] * h_[m][2];
              if (ii == 3)
                Ah_[ns_ + 1][0] += b3_[m][0] * h_[m][0] + b3_[m][1] * h_[m][1] +
                                   b3_[m][2] * h_[m][2];
              if (ii == 4)
                Ah_[ns_ + 1][1] += b4_[m][0] * h_[m][0] + b4_[m][1] * h_[m][1] +
                                   b4_[m][2] * h_[m][2];
              if (ii == 5)
                Ah_[ns_ + 1][2] += b5_[m][0] * h_[m][0] + b5_[m][1] * h_[m][1] +
                                   b5_[m][2] * h_[m][2];
            }
          }
        }
      }
      // now do right hand corner terms, ns to (ns+1).
      for (int m = 0; m < 2; m++) {
        for (int m1 = 0; m1 < 3; m1++) {
          if (ii == 0)
            Ah_[ns_][0] += zcon_[0][0][m][m1] * h_[ns_ + m][m1];
          if (ii == 1)
            Ah_[ns_][1] += zcon_[0][1][m][m1] * h_[ns_ + m][m1];
          if (ii == 2)
            Ah_[ns_][2] += zcon_[0][2][m][m1] * h_[ns_ + m][m1];
          if (ii == 3)
            Ah_[ns_ + 1][0] += zcon_[1][0][m][m1] * h_[ns_ + m][m1];
          if (ii == 4)
            Ah_[ns_ + 1][1] += zcon_[1][1][m][m1] * h_[ns_ + m][m1];
          if (ii == 5)
            Ah_[ns_ + 1][2] += zcon_[1][2][m][m1] * h_[ns_ + m][m1];
        }
      }
    }

    hAh = 0.0;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int k = zlo; k <= zhi; k++) {
        for (int j = ylo; j <= yhi; j++) {
          for (int i = xlo; i <= xhi; i++) {
            x = i;
            y = j;
            z = k;
            if (x < 0)
              x += nx_;
            if (x >= nx_)
              x -= nx_;
            if (y < 0)
              y += ny_;
            if (y >= ny_)
              y -= ny_;
            if (z < 0)
              z += nz_;
            if (z >= nz_)
              z -= nz_;
            m = nxy_ * z + nx_ * y + x;
            hAh += h_[m][m3] * Ah_[m][m3];
          }
        }
      }
    }

    lambda = localgg_ / hAh;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int k = zlo; k <= zhi; k++) {
        for (int j = ylo; j <= yhi; j++) {
          for (int i = xlo; i <= xhi; i++) {
            x = i;
            y = j;
            z = k;
            if (x < 0)
              x += nx_;
            if (x >= nx_)
              x -= nx_;
            if (y < 0)
              y += ny_;
            if (y >= ny_)
              y -= ny_;
            if (z < 0)
              z += nz_;
            if (z >= nz_)
              z -= nz_;
            m = nxy_ * z + nx_ * y + x;
            u_[m][m3] -= lambda * h_[m][m3];
            gb_[m][m3] -= lambda * Ah_[m][m3];
          }
        }
      }
    }

    gglast = localgg_;
    localgg_ = 0.0;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int k = zlo; k <= zhi; k++) {
        for (int j = ylo; j <= yhi; j++) {
          for (int i = xlo; i <= xhi; i++) {
            x = i;
            y = j;
            z = k;
            if (x < 0)
              x += nx_;
            if (x >= nx_)
              x -= nx_;
            if (y < 0)
              y += ny_;
            if (y >= ny_)
              y -= ny_;
            if (z < 0)
              z += nz_;
            if (z >= nz_)
              z -= nz_;
            m = nxy_ * z + nx_ * y + x;
            localgg_ += gb_[m][m3] * gb_[m][m3];
          }
        }
      }
    }

#ifdef DEBUG
    std::clog << "ThermalStrain::localDembx localgg_ = " << localgg_ << endl;
    std::clog.flush();
#endif

    if (localgg_ < localgtest_) {
      return Lstep;
    } else {
      gamma = localgg_ / gglast;
#ifdef DEBUG
      std::clog << "ThermalStrain::localDembx gamma = " << gamma << endl;
      std::clog.flush();
#endif
      for (int m3 = 0; m3 < 3; m3++) {
        for (int k = zlo; k <= zhi; k++) {
          for (int j = ylo; j <= yhi; j++) {
            for (int i = xlo; i <= xhi; i++) {
              x = i;
              y = j;
              z = k;
              if (x < 0)
                x += nx_;
              if (x >= nx_)
                x -= nx_;
              if (y < 0)
                y += ny_;
              if (y >= ny_)
                y -= ny_;
              if (z < 0)
                z += nz_;
              if (z >= nz_)
                z -= nz_;
              m = nxy_ * z + nx_ * y + x;
              h_[m][m3] = gb_[m][m3] + gamma * h_[m][m3];
            }
          }
        }
      }
    }
  }

  return Lstep;
}

// void ThermalStrain::localRelax(int boxsize, int x, int y, int z, int index) {
void ThermalStrain::localRelax(int xlo, int xhi, int ylo, int yhi, int zlo,
                               int zhi, int x, int y, int z, int index) {

#ifdef DEBUG
  std::clog << "ThermalStrain::localRelax " << endl;
  std::clog << "ThermalStrain::localRelax x, y, z and index are: " << x << " "
            << y << " " << z << " " << index << endl;
  std::clog.flush();
#endif

  int localldemb = 60; //, ltot = 0;
  // int iskip;
  // double utot;
  // int Lstep;

  // int halfbox = (int)(boxsize / 2);

  // Lstep = localDembx(boxsize, x, y, z, localldemb, 0);

  // ltot += Lstep;

  // ltot = localDembx(xlo, xhi, ylo, yhi, zlo, zhi, localldemb, 0);
  localDembx(xlo, xhi, ylo, yhi, zlo, zhi, localldemb, 0);

#ifdef DEBUG
  std::clog << "ThermalStrain::localRelax localgg_ = " << localgg_
            << " local relaxation steps ltot = " << ltot << endl;
  std::clog.flush();
#endif

  return;
}

// not used...
double ThermalStrain::getEleStrain(int i, int j) {
  if ((i >= ns_) || (i < 0) || (j < 0) || (j >= 6)) {
    std::clog << "i should be between 0 and ns_, "
              << "and j should be between 0 and 6." << endl;
    std::cerr << "i should be between 0 and ns_, "
              << "and j should be between 0 and 6." << endl;
    exit(1);
  } else {
    return elestrain_[i][j];
  }
}

double ThermalStrain::getEleStress(int i, int j) {
  if ((i >= ns_) || (i < 0) || (j < 0) || (j >= 6)) {
    std::clog << "i should be between 0 and ns_, "
              << "and j should be between 0 and 6." << endl;
    cerr << "i should be between 0 and ns_, "
         << "and j should be between 0 and 6." << endl;
    exit(1);
  } else {
    return elestress_[i][j];
  }
}

void ThermalStrain::setExpansionCoord(int index, std::vector<int> cor) {
  map<int, vector<int>>::iterator p = exp_.find(index);
  if (p != exp_.end()) {
    p->second = cor;
  } else {
    exp_.insert(make_pair(index, cor));
  }
}
