/**
@struct SaturatingRateParameters
@brief Bullard 2015 / Han 2025 Eq. 7 rate parameters for one direction
(dissolution or precipitation) of a phase.

Populated by `KineticController::parseSaturatingRateBlock` from a
`dissolution` (required) or `precipitation` (optional) sub-block of a
phase's `kinetic_data`. `SaturatingRateModel` holds two of these — one
each for dissolution and precipitation — as `std::optional`, so a
missing precipitation block naturally falls back to the dissolution
values by microscopic reversibility.

Units:
- rateConstant: mol / (m^2 * s)
- B, n: dimensionless
*/

#ifndef SRC_THAMESLIB_SATURATINGRATEPARAMETERS_H_
#define SRC_THAMESLIB_SATURATINGRATEPARAMETERS_H_

struct SaturatingRateParameters {
  double rateConstant;  ///< k in Bullard 2015 / Han 2025 Eq. 7
  double B;             ///< dimensionless exponent scale (Han); = 1/C1 (Bullard)
  double n;             ///< dimensionless exponent (Han n; Bullard r)
};

#endif  // SRC_THAMESLIB_SATURATINGRATEPARAMETERS_H_
