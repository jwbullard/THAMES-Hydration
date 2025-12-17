/**
@file AppliedStrain.cc
@brief Method definitions for the AppliedStrain derived class.
*/
#include "AppliedStrain.h"

using std::cout;
using std::endl;
using std::string;

AppliedStrain::AppliedStrain(int nx, int ny, int nz, int dim,
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

#ifdef DEBUG
  std::clog << "AppliedStrain::AppliedStrain Constructor" << endl;
  std::clog.flush();
#endif

  exx_ = eyy_ = ezz_ = 0.0;
  exz_ = eyz_ = exy_ = 0.0;

  kmax_ = 80; // 3;
  std::clog
      << endl
      << "AppliedStrain::AppliedStrain - "
         "the number of relaxation steps for elastic computation :  kmax_ = "
      << kmax_ << endl;

  is_[0] = 26;
  is_[1] = 2;
  is_[2] = 1;
  is_[3] = 0;
  is_[4] = 25;
  is_[5] = 18;
  is_[6] = 17;
  is_[7] = 16;
}

void AppliedStrain::femat() {
  // double dndx[8], dndy[8], dndz[8];
  // double g[3][3][3];
  // double es[6][8][3];
  // double delta[8][3];
  // double x, y, z;
  // int nxy = nx * ny;
  // int is[8]; moved in constructor as is_[8]
  double delta[8][3];
  int nx1 = nx_ - 1;
  int ny1 = ny_ - 1;
  int nz1 = nz_ - 1;
  double nx_dbl = nx_;
  double ny_dbl = ny_;
  double nz_dbl = nz_;
  int m;
  double sum = 0;

  /*
  ///
  /// (User) NOTE: complete elastic modulus matrix is used, so an anisotropic
  /// matrix could be directly input at any point, since program is written to
  /// use a general elastic moduli tensor, but is only explicitly implemented
  /// for isotropic materials.
  ///

  ///
  /// Initialize stiffness matrices
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
  /// Set up elastic moduli matrices for each kind of element
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
                  sum = 0.0;
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
  */

  ///
  /// Set up vector for linear term, b, and constant term, C, in the elastic
  /// energy. This is done using the stiffness matrices, and the periodic terms
  /// in the applied strain that come in at the boundary pixels via the periodic
  /// boundary conditions and the condition that an applied macroscopic strain
  /// exists (see Sec. 2.2 in the manual). It is easier to set b up this way
  /// than to analytically write out all the terms involved.
  ///

  ///
  /// Initialize b and C
  ///

  for (int m3 = 0; m3 < 3; m3++) {
    for (int m = 0; m < ns_; m++) {
      b_[m][m3] = 0.0;
    }
  }
  C_ = 0.0;

  ///
  /// For all cases, the correspondence between 0-7 finite element node labels
  /// and 0-26 neighbor labels is (see Table 4 in manual):
  /// 0:ib_[m][26], 1:ib_[m][2],
  /// 2:ib_[m][1], 3:ib_[m][0],
  /// 4:ib_[m][25], 5:ib_[m][18],
  /// 6:ib_[m][17], 7:ib_[m][16].
  ///

  // moved in constructor as is_[8]  //check!
  // is[0] = 26;
  // is[1] = 2;
  // is[2] = 1;
  // is[3] = 0;
  // is[4] = 25;
  // is[5] = 18;
  // is[6] = 17;
  // is[7] = 16;

  ///
  /// x = nx - 1, face
  ///

  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 1) || (i8 == 2) || (i8 == 5) || (i8 == 6)) {
        delta[i8][0] = exx_ * nx_dbl;
        delta[i8][1] = exy_ * nx_dbl;
        delta[i8][2] = exz_ * nx_dbl;
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
              C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                    delta[mm][nn];
            }
          }
          b_[ib_[m][is_[mm]]][nn] += sum;
        }
      }
    }
  }

  ///
  /// y = ny - 1, face
  ///

  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 2) || (i8 == 3) || (i8 == 6) || (i8 == 7)) {
        delta[i8][0] = exy_ * ny_dbl;
        delta[i8][1] = eyy_ * ny_dbl;
        delta[i8][2] = eyz_ * ny_dbl;
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
              C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                    delta[mm][nn];
            }
          }
          b_[ib_[m][is_[mm]]][nn] += sum;
        }
      }
    }
  }

  ///
  /// z = nz - 1, face
  ///

  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 4) || (i8 == 5) || (i8 == 6) || (i8 == 7)) {
        delta[i8][0] = exz_ * nz_dbl;
        delta[i8][1] = eyz_ * nz_dbl;
        delta[i8][2] = ezz_ * nz_dbl;
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
              C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                    delta[mm][nn];
            }
          }
          b_[ib_[m][is_[mm]]][nn] += sum;
        }
      }
    }
  }

  ///
  /// x = nx - 1, y = ny - 1, edge
  ///

  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 1) || (i8 == 5)) {
        delta[i8][0] = exx_ * nx_dbl;
        delta[i8][1] = exy_ * nx_dbl;
        delta[i8][2] = exz_ * nx_dbl;
      }
      if ((i8 == 3) || (i8 == 7)) {
        delta[i8][0] = exy_ * ny_dbl;
        delta[i8][1] = eyy_ * ny_dbl;
        delta[i8][2] = eyz_ * ny_dbl;
      }
      if ((i8 == 2) || (i8 == 6)) {
        delta[i8][0] = exy_ * ny_dbl + exx_ * nx_dbl;
        delta[i8][1] = eyy_ * ny_dbl + exy_ * nx_dbl;
        delta[i8][2] = eyz_ * ny_dbl + exz_ * nx_dbl;
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
            C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                  delta[mm][nn];
          }
        }
        b_[ib_[m][is_[mm]]][nn] += sum;
      }
    }
  }

  ///
  /// x = nx - 1, z = nz - 1, edge
  ///

  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 1) || (i8 == 2)) {
        delta[i8][0] = exx_ * nx_dbl;
        delta[i8][1] = exy_ * nx_dbl;
        delta[i8][2] = exz_ * nx_dbl;
      }
      if ((i8 == 4) || (i8 == 7)) {
        delta[i8][0] = exz_ * nz_dbl;
        delta[i8][1] = eyz_ * nz_dbl;
        delta[i8][2] = ezz_ * nz_dbl;
      }
      if ((i8 == 5) || (i8 == 6)) {
        delta[i8][0] = exz_ * nz_dbl + exx_ * nx_dbl;
        delta[i8][1] = eyz_ * nz_dbl + exy_ * nx_dbl;
        delta[i8][2] = ezz_ * nz_dbl + exz_ * nx_dbl;
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
            C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                  delta[mm][nn];
          }
        }
        b_[ib_[m][is_[mm]]][nn] += sum;
      }
    }
  }

  ///
  /// y = ny - 1, z = nz - 1, edge
  ///

  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if ((i8 == 4) || (i8 == 5)) {
        delta[i8][0] = exz_ * nz_dbl;
        delta[i8][1] = eyz_ * nz_dbl;
        delta[i8][2] = ezz_ * nz_dbl;
      }
      if ((i8 == 2) || (i8 == 3)) {
        delta[i8][0] = exy_ * ny_dbl;
        delta[i8][1] = eyy_ * ny_dbl;
        delta[i8][2] = eyz_ * ny_dbl;
      }
      if ((i8 == 6) || (i8 == 7)) {
        delta[i8][0] = exy_ * ny_dbl + exz_ * nz_dbl;
        delta[i8][1] = eyy_ * ny_dbl + eyz_ * nz_dbl;
        delta[i8][2] = eyz_ * ny_dbl + ezz_ * nz_dbl;
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
            C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                  delta[mm][nn];
          }
        }
        b_[ib_[m][is_[mm]]][nn] += sum;
      }
    }
  }

  ///
  /// x = nx - 1, y = ny - 1, z = nz - 1, corner
  ///

  for (int i3 = 0; i3 < 3; i3++) {
    for (int i8 = 0; i8 < 8; i8++) {
      delta[i8][i3] = 0.0;
      if (i8 == 1) {
        delta[i8][0] = exx_ * nx_dbl;
        delta[i8][1] = exy_ * nx_dbl;
        delta[i8][2] = exz_ * nx_dbl;
      }
      if (i8 == 3) {
        delta[i8][0] = exy_ * ny_dbl;
        delta[i8][1] = eyy_ * ny_dbl;
        delta[i8][2] = eyz_ * ny_dbl;
      }
      if (i8 == 4) {
        delta[i8][0] = exz_ * nz_dbl;
        delta[i8][1] = eyz_ * nz_dbl;
        delta[i8][2] = ezz_ * nz_dbl;
      }
      if (i8 == 7) {
        delta[i8][0] = exy_ * ny_dbl + exz_ * nz_dbl;
        delta[i8][1] = eyy_ * ny_dbl + eyz_ * nz_dbl;
        delta[i8][2] = eyz_ * ny_dbl + ezz_ * nz_dbl;
      }
      if (i8 == 5) {
        delta[i8][0] = exx_ * nx_dbl + exz_ * nz_dbl;
        delta[i8][1] = exy_ * nx_dbl + eyz_ * nz_dbl;
        delta[i8][2] = exz_ * nx_dbl + ezz_ * nz_dbl;
      }
      if (i8 == 2) {
        delta[i8][0] = exx_ * nx_dbl + exy_ * ny_dbl;
        delta[i8][1] = exy_ * nx_dbl + eyy_ * ny_dbl;
        delta[i8][2] = exz_ * nx_dbl + eyz_ * ny_dbl;
      }
      if (i8 == 6) {
        delta[i8][0] = exx_ * nx_dbl + exy_ * ny_dbl + exz_ * nz_dbl;
        delta[i8][1] = exy_ * nx_dbl + eyy_ * ny_dbl + eyz_ * nz_dbl;
        delta[i8][2] = exz_ * nx_dbl + eyz_ * ny_dbl + ezz_ * nz_dbl;
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
          C_ += 0.5 * delta[m8][m3] * dk_[pix_[m]][m8][m3][mm][nn] *
                delta[mm][nn];
        }
      }
      b_[ib_[m][is_[mm]]][nn] += sum;
    }
  }

  return;
}

// double AppliedStrain::energy(int nx, int ny, int nz, int ns) {
double AppliedStrain::energy() {
  double utot = 0.0;
  for (int m3 = 0; m3 < 3; m3++) {
    for (int m = 0; m < ns_; m++) {
      gb_[m][m3] = 0.0;
    }
  }

  ///
  /// Do global matrix multiply via small stiffness matrices, gb_ = A * u_
  /// The long statement below correctly brings in all the terms from the
  /// global matrix A using only the small stiffness matrices.
  ///

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

  utot = C_;
  for (int m3 = 0; m3 < 3; m3++) {
    for (int m = 0; m < ns_; m++) {
      utot += 0.5 * u_[m][m3] * gb_[m][m3] + b_[m][m3] * u_[m][m3];
      gb_[m][m3] += b_[m][m3];
    }
  }

  return utot;
}

// int AppliedStrain::dembx(int ns, double gg, int ldemb, int kkk) {
int AppliedStrain::dembx(int ldemb, int kkk) {
  double lambda, gamma;
  double hAh, gglast;
  int Lstep;

  if (kkk == 0) {
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < ns_; m++) {
        h_[m][m3] = gb_[m][m3];
      }
    }
  }

  ///
  /// Lstep counts the number of conjugate gradient steps taken in each
  /// call to dembx.

  Lstep = 0;

  for (int ijk = 0; ijk < ldemb; ijk++) {
    Lstep += 1;

    // Show CG iteration progress every 10 iterations
    if (Lstep % 10 == 0 || Lstep == 1) {
      std::clog << "    CG iteration " << Lstep << "/" << ldemb
                << " (gg=" << std::scientific << gg_ << ")" << endl;
      std::clog.flush();
    }

    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < ns_; m++) {
        Ah_[m][m3] = 0.0;
      }
    }

    ///
    /// Do global matrix multiply via small stiffness matrices, Ah_ = A_ * h_.
    /// The long statement below correctly brings in all the terms from the
    /// global matrix A_ using only the small stiffness matrices dk_.

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

    hAh = 0.0;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < ns_; m++) {
        hAh += h_[m][m3] * Ah_[m][m3];
      }
    }
    lambda = gg_ / hAh;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < ns_; m++) {
        u_[m][m3] = u_[m][m3] - lambda * h_[m][m3];
        gb_[m][m3] = gb_[m][m3] - lambda * Ah_[m][m3];
      }
    }

    gglast = gg_;
    gg_ = 0.0;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < ns_; m++) {
        gg_ += gb_[m][m3] * gb_[m][m3];
      }
    }
    if (gg_ < gtest_)
      return Lstep;

    gamma = gg_ / gglast;
    for (int m3 = 0; m3 < 3; m3++) {
      for (int m = 0; m < ns_; m++) {
        h_[m][m3] = gb_[m][m3] + gamma * h_[m][m3];
      }
    }
  }
  return Lstep;
}

// void AppliedStrain::stress(int nx, int ny, int nz, int ns) {
void AppliedStrain::stress() {
  double dndx[8], dndy[8], dndz[8];
  double es[6][8][3];
  double uu[8][3];
  double str11, str22, str33, str13, str23, str12;
  double s11, s22, s33, s13, s23, s12;
  int m = 0;
  double nx_dbl = nx_;
  double ny_dbl = ny_;
  double nz_dbl = nz_;
  int nx1 = nx_ - 1;
  int ny1 = ny_ - 1;
  int nz1 = nz_ - 1;

  ///
  /// Set up single element strain matrix
  /// dndx, dndy, and dndz are the components of the average strain matrix
  /// in a pixel.
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
  /// Build averaged strain matrix, follows code in femat, but for average
  /// strain over the pixel, not the strain at a point.
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
  /// Compute components of the average stress and strain tensors in each pixel.
  ///

  strxx_ = stryy_ = strzz_ = strxz_ = stryz_ = strxy_ = 0.0;
  sxx_ = syy_ = szz_ = sxz_ = syz_ = sxy_ = 0.0;

  double exx_nx_dbl = exx_ * nx_dbl;
  double exy_nx_dbl = exy_ * nx_dbl;
  double exz_nx_dbl = exz_ * nx_dbl;

  double exy_ny_dbl = exy_ * ny_dbl;
  double eyy_ny_dbl = eyy_ * ny_dbl;
  double eyz_ny_dbl = eyz_ * ny_dbl;

  double exz_nz_dbl = exz_ * nz_dbl;
  double eyz_nz_dbl = eyz_ * nz_dbl;
  double ezz_nz_dbl = ezz_ * nz_dbl;

  int i, j, k, nxy_k, nx_j;

  double strxx, stryy, strzz, strxz, stryz, strxy;
  double sxx, syy, szz, sxz, syz, sxy;
  for (i = 0; i < nx_; i++) {
    strxx = stryy = strzz = strxz = stryz = strxy = 0.0;
    sxx = syy = szz = sxz = syz = sxy = 0.0;
    for (k = 0; k < nz_; k++) {
      for (j = 0; j < ny_; j++) {
        nxy_k = nxy_ * k;
        nx_j = nxy_k + nx_ * j;
        m = nx_j + i;

        ///
        /// Load in elements of 8-vector using periodic boundary condition.
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
        /// wrong for a pixel on periodic boundary. Since they come from an
        /// opposite face, need to put in applied strain to correct them.
        ///

        // if (i == (nx - 1)) {
        if (i == nx1) {
          uu[1][0] += exx_nx_dbl; // exx_nx_dbl = exx_ * nx_dbl;
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
        // if (j == (ny - 1)) {
        if (j == ny1) {
          uu[2][0] += exy_ny_dbl; // exy_ny_dbl = exy_ * ny_dbl;
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
        // if (k == (nz - 1)) {
        if (k == nz1) {
          uu[4][0] += exz_nz_dbl; // exz_nz_dbl = exz_ * nz_dbl;
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
        /// Local stresses and strains in a pixel
        ///

        str11 = str22 = str33 = str13 = str23 = str12 = 0.0;
        s11 = s22 = s33 = s13 = s23 = s12 = 0.0;

        for (int n3 = 0; n3 < 3; n3++) {
          for (int n8 = 0; n8 < 8; n8++) {
            s11 += es[0][n8][n3] * uu[n8][n3];
            s22 += es[1][n8][n3] * uu[n8][n3];
            s33 += es[2][n8][n3] * uu[n8][n3];
            s13 += es[3][n8][n3] * uu[n8][n3];
            s23 += es[4][n8][n3] * uu[n8][n3];
            s12 += es[5][n8][n3] * uu[n8][n3];
            for (int n = 0; n < 6; n++) {
              str11 += cmod_[pix_[m]][0][n] * es[n][n8][n3] * uu[n8][n3];
              str22 += cmod_[pix_[m]][1][n] * es[n][n8][n3] * uu[n8][n3];
              str33 += cmod_[pix_[m]][2][n] * es[n][n8][n3] * uu[n8][n3];
              str13 += cmod_[pix_[m]][3][n] * es[n][n8][n3] * uu[n8][n3];
              str23 += cmod_[pix_[m]][4][n] * es[n][n8][n3] * uu[n8][n3];
              str12 += cmod_[pix_[m]][5][n] * es[n][n8][n3] * uu[n8][n3];
            }
          }
        }

        strxx += str11;
        stryy += str22;
        strzz += str33;
        strxy += str12;
        strxz += str13;
        stryz += str23;
        sxx += s11;
        syy += s22;
        szz += s33;
        sxy += s12;
        sxz += s13;
        syz += s23;

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
        /// Compute the strain energy for each element
        ///

        for (int ii = 0; ii < 6; ii++) {
          for (int jj = 0; jj < 6; jj++) {
            strainengy_[m] = 0.5 * elestrain_[m][ii] * cmod_[pix_[m]][ii][jj] *
                             elestrain_[m][jj];
          }
        }

        ///
        /// Sum local strains and stresses into global values
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

    /* Averaging depends on whether layer averaging is being done (for ITZ) */
    if (hasAggregateSlab_) {
      /* Layer average of stresses and strains in each layer perpendicular to
       * aggregate slab */

      double nyzd = static_cast<double>(ny_ * nz_);
      strxx /= (nyzd);
      stryy /= (nyzd);
      strzz /= (nyzd);
      strxy /= (nyzd);
      strxz /= (nyzd);
      stryz /= (nyzd);
      sxx /= (nyzd);
      syy /= (nyzd);
      szz /= (nyzd);
      sxy /= (nyzd);
      sxz /= (nyzd);
      syz /= (nyzd);

      K_[i] = (1.0 / 3.0) * (strxx + stryy + strzz) / (sxx + syy + szz);
      G_[i] = (1.0 / 3.0) * ((strxz / sxz) + (stryz / syz) + (strxy / sxy));
      /* GODZILLA */
      /* Ended here and have not provided a way to output this layer-by-layer
       * information */
      /* GODZILLA */
    }
  }

  ///
  /// Volume average of global stresses and strains
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

void AppliedStrain::relax(int kmax) {

  ///
  /// RELAXATION LOOP
  /// (USER) kmax is the maximum number of times dembx will be called, with
  /// ldemb conjugate gradient steps performed during each call. The total
  /// number of conjugate gradient steps allowed for a given elastic computation
  /// is kmax * ldemb.
  ///

#ifdef DEBUG
  std::clog << "AppliedStrain::relax Entering energy" << endl;
  std::clog.flush();
#endif

  // Using ldemb=100 to match VCCTL's elastic.c for equivalent performance
  int ldemb = 100, ltot = 0;
  double utot;
  int Lstep;

  std::clog << "AppliedStrain::relax Starting relaxation (kmax=" << kmax
            << ", ldemb=" << ldemb << ", gtest=" << gtest_ << ")" << endl;
  std::clog.flush();

  ///
  /// Call energy to get initial energy and initial gradient.
  ///

  // utot = energy(nx_, ny_, nz_, ns_);
  utot = energy();

  ///
  /// gg_ is the norm squared of the gradient (gg_ = gb_ * gb_)
  ///

  gg_ = 0.0;
  for (int m3 = 0; m3 < 3; m3++) {
    for (int m = 0; m < ns_; m++) {
      gg_ += (gb_[m][m3] * gb_[m][m3]);
    }
  }

  double gginit = gg_; // Save initial gradient for progress tracking
  std::clog << "  Initial gg = " << gg_ << ", gtest = " << gtest_ << endl;
  std::clog.flush();

  for (int kkk = 0; kkk < kmax; kkk++) {

    ///
    /// Call dembx to go into the conjugate gradient solver
    ///

    // Lstep = dembx(ns_, gg_, ldemb, kkk);
    Lstep = dembx(ldemb, kkk);
    ltot += Lstep;

    ///
    /// Call energy to compute energy after dembx call. If gg_ < gtest_, this
    /// will be the final energy. If gg_ is still larger than gtest_, then this
    /// will give an intermediate energy with which to check how the relaxation
    /// process is coming along.
    ///

    // utot = energy(nx_, ny_, nz_, ns_);
    utot = energy();

    // Always show progress output (not just in verbose mode)
    double percent = 100.0 * ((double)(kkk + 1) / (double)(kmax));
    std::clog << "  Relaxation step " << (kkk + 1) << "/" << kmax << " ("
              << std::fixed << std::setprecision(1) << percent << "%)"
              << " Energy=" << std::scientific << std::setprecision(4) << utot
              << " gg=" << gg_ << endl;
    std::clog.flush();

    if (verbose_) {
      std::clog << "AppliedStrain::relax Number of conjugate steps = " << ltot
                << endl;
      std::clog.flush();
    }

    ///
    /// If relaxation process is finished, jump out of loop.
    ///

    if (gg_ >= gtest_) {

      ///
      /// If relaxation process will continue, compute an d output stresses and
      /// strains as an additional aid to judge how the relaxation procedure is
      /// progressing.
      ///

      // stress(nx_, ny_, nz_, ns_);
      stress();

#ifdef DEBUG
      std::clog
          << "AppliedStrain::relax stresses: xx,    yy,    zz,    xz,    yz,  "
             "  xy"
          << endl;
      std::clog << "AppliedStrain::relax         " << strxx_ << " " << stryy_
                << " " << strzz_ << " " << strxz_ << " " << stryz_ << " "
                << strxy_ << endl;
      std::clog
          << "AppliedStrain::relax strains: xx,    yy,    zz,    xz,    yz,   "
             " xy"
          << endl;
      std::clog << "AppliedStrain::relax         " << sxx_ << " " << syy_ << " "
                << szz_ << " " << sxz_ << " " << syz_ << " " << sxy_ << endl;
      std::clog.flush();
#endif

    } else {
      // Converged - show success message
      std::clog << "  CONVERGED at step " << (kkk + 1) << "/" << kmax
                << " (gg=" << std::scientific << gg_ << " < gtest=" << gtest_
                << ")" << endl;
      std::clog.flush();
      break;
    }
  }

  // stress(nx_, ny_, nz_, ns_);
  stress();

  if (verbose_) {
    std::clog
        << "AppliedStrain::relax stresses: xx,    yy,    zz,    xz,    yz,    "
           "xy"
        << endl;
    std::clog << "AppliedStrain::relax         " << strxx_ << " " << stryy_
              << " " << strzz_ << " " << strxz_ << " " << stryz_ << " "
              << strxy_ << endl;
    std::clog
        << "AppliedStrain::relax strains: xx,    yy,    zz,    xz,    yz,    xy"
        << endl;
    std::clog << "AppliedStrain::relax         " << sxx_ << " " << syy_ << " "
              << szz_ << " " << sxz_ << " " << syz_ << " " << sxy_ << endl;
    std::clog.flush();
  }

  return;
}

void AppliedStrain::calc(std::vector<int> *p_vectPhId, double exx, double eyy,
                         double ezz, double exz, double eyz, double exy) {
  int kmax = kmax_; // 80;
  int m;

  ///
  /// Read in a microstructure in subroutine setMicrostructure, and set up
  /// pix_[m] with the appropriate phase assignments.
  ///

  setMicrostructure(p_vectPhId);

  ///
  /// Count and output the volume fractions of the different phases.
  ///

  // assig(ns_, nphase_);
  assig();

#ifdef DEBUG
  for (int i = 0; i < nphase_; i++) {
    std::clog << "AppliedStrain::calc Phase " << i
              << " bulk = " << phasemod_[i][0] << " shear = " << phasemod_[i][1]
              << endl;
  }

  for (int i = 0; i < nphase_; i++) {
    std::clog << "AppliedStrain::calc Volume fraction of phase " << i << "  is "
              << prob_[i] << endl;
  }
  std::clog.flush();
#endif

  ///
  /// (USER) Set applied strains.
  /// Actual shear strain applied in is exy_, exz_, and eyz_ as given in the
  /// statements below. The engineering shear strain, by which the shear modulus
  /// is usually defined, is twice these values.
  ///

  exx_ = exx;
  eyy_ = eyy;
  ezz_ = ezz;
  exz_ = exz;
  eyz_ = eyz;
  exy_ = exy;

#ifdef DEBUG
  std::clog << "AppliedStrain::calc Applied engineering strains" << endl;
  std::clog
      << "AppliedStrain::calc  exx_,   eyy_,   ezz_,   exz_,   eyz_,   exy_"
      << endl;
  std::clog << "AppliedStrain::calc " << exx_ << "   " << eyy_ << "   " << ezz_
            << "   " << exz_ << "   " << eyz_ << "   " << exy_ << endl;
  std::clog.flush();
#endif

  ///
  /// Set up elastic modulus variables, finite element stiffness matrices,
  /// the constant, C, and vector, b, required for computing the energy.
  /// (USER) If anisotropic elastic moduli tensors are used, these need to be
  /// input in subroutine femat.
  ///

  femat();

  ///
  /// Apply chosen strains as a homogeneous macroscopic strain
  /// as the initial condition.
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
        // u_[m][0] = x * exx + y * exy + z * exz;
        // u_[m][1] = x * exy + y * eyy + z * eyz;
        // u_[m][2] = x * exz + y * eyz + z * ezz;
        u_[m][0] = i * exx + j * exy + k * exz;
        u_[m][1] = i * exy + j * eyy + k * eyz;
        u_[m][2] = i * exz + j * eyz + k * ezz;
      }
    }
  }

  ///
  /// RELAXATION LOOP
  /// (USER) kmax is the maximum number of times dembx will be called, with
  /// ldemb conjugate gradient steps performed during each call. The total
  /// number of conjugate gradient steps allowed for a given elastic computation
  /// is kmax * ldemb.
  ///

  relax(kmax);

  return;
}

double AppliedStrain::getBulkModulus(vector<int> *p_vectPhId) {
  double bulk;
  double Stress, Strain;
  // Stress = Strain = 0.0;

  // std::clog << "AppliedStrain::getBulkModulus - bf-calc" << endl;

  calc(p_vectPhId, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05);

  // std::clog << "AppliedStrain::getBulkModulus - af-calc" << endl;

  // for (int i = 0; i < 3; i++) {
  //   Stress += getStress(i);
  //   Strain += getStrain(i);
  // }
  Stress = getStress(0) + getStress(1) + getStress(2);
  Strain = getStrain(0) + getStrain(1) + getStrain(2);
  bulk = Stress / 3.0 / Strain;

  return bulk;
}
