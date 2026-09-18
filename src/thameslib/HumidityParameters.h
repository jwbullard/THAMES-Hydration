/**
@struct HumidityParameters
@brief Relative-humidity dependence of one phase's reaction rate.

The rate of every kinetic model is multiplied by

    f(h) = [ max(0, (h - h0) / (1 - h0)) ]^n,    f = 1 for h >= 1

where h is the Kelvin relative humidity of the pore solution, set once
per cycle by KineticController from the pore-size bin that holds the
meniscus (Lattice::getLargestSaturatedPore). In saturated mode every pore
is full, so h = 1 and f = 1 exactly.

Defaults (h0 = 0.70, n = 1, linear): Patel, Killoh, Parrott & Gutteridge,
Mater. Struct. 21 (1988) 192-197, Table 2 — QXRD degree-of-hydration
gain from 14 to 90 d, relative to saturated curing, is 0.38 / 0.60 / 0.93
at 81 / 91 / 97 % RH and ~0 at <= 69 % RH (linear h0 = 0.70 gives
0.37 / 0.70 / 0.90). Consistent with Killoh, Parrott & Patel, ACI SP-114
(1989) 157-174, where bound water is flat below ~70 % RH and rises
roughly linearly above it. Both are increments over weeks, not
instantaneous rates. The experimental h is an ambient (total) RH; in
Kelvin-only terms h0 shifts to ~h0 / a_w ~ 0.71, which is ignored here.

Documented alternative: Parrott & Killoh, Brit. Ceram. Proc. 35 (1984)
41-53, Fig. 6 — h0 = 0.55, n = 4, a fit to older overall-cement data
applied to all clinker phases. Underestimates the Patel 1988 data by
~3x at 81 % RH.

Per-phase overrides are expected: Killoh 1989 reports that the
pozzolanic reaction of fly ash is severely restricted below 80 % RH, and
Patel 1988 Figs. 1-2 show belite more sensitive, and C3A / C4AF less
sensitive, than alite.
*/

#ifndef SRC_THAMESLIB_HUMIDITYPARAMETERS_H_
#define SRC_THAMESLIB_HUMIDITYPARAMETERS_H_

#include <cmath>

struct HumidityParameters {
  double h0 = 0.70;      /**< RH at and below which the rate is zero */
  double exponent = 1.0; /**< shape of f(h) above h0 */
};

/**
@brief Rate multiplier for relative humidity rh.

@param rh is the Kelvin relative humidity [0, 1]
@param p holds h0 and the exponent
@return f(rh) in [0, 1]
*/
inline double rhRateFactor(double rh, const HumidityParameters &p) {
  if (rh >= 1.0)
    return 1.0;
  double x = (rh - p.h0) / (1.0 - p.h0);
  return (x <= 0.0) ? 0.0 : std::pow(x, p.exponent);
}

#endif // SRC_THAMESLIB_HUMIDITYPARAMETERS_H_
