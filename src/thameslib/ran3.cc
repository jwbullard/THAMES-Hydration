/**
@file ran3.cc
@brief Algorithm for the ran3 portable random number generator.

Taken directly from Press, W.H., Teukolsky, S.A., Vetterling, W.T.,
Flannery, B.P., Numerical Recipes in C, 2nd Edition, University Press,
London, 1997 (Section 7.1). Constants below are Knuth's subtractive
generator seed values from the same reference.

@warning NOT thread-safe: `ma[]`, `inext`, `inextp`, `iff` are `static`
locals holding generator state across calls. If THAMES ever runs kinetic
models in parallel, this must be replaced with a per-thread engine.
Reset with a negative seed to re-initialize the state.
*/
#include "ran3.h"

// Knuth's subtractive generator constants (Numerical Recipes §7.1).
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0 / MBIG)

double ran3(int *idum) {
  static int inext, inextp;
  static int ma[56];
  static int iff = 0;
  int mj, mk;
  int i, ii, k;

  if (*idum < 0 || iff == 0) {
    iff = 1;
    mj = MSEED - (*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for (i = 1; i <= 54; i++) {
      ii = (21 * i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ)
        mk += MBIG;
      mj = ma[ii];
    }
    for (k = 1; k <= 4; k++)
      for (i = 1; i <= 55; i++) {
        ma[i] -= ma[1 + (i + 30) % 55];
        if (ma[i] < MZ)
          ma[i] += MBIG;
      }
    inext = 0;
    inextp = 31;
    *idum = 1;
  }

  if (++inext == 56)
    inext = 1;
  if (++inextp == 56)
    inextp = 1;
  mj = ma[inext] - ma[inextp];
  if (mj < MZ)
    mj += MBIG;
  ma[inext] = mj;

  return mj * FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
