/**
@file PhysicalConstants.h
@brief Universal physical and mathematical constants.

Held in a small dependency-free header (no GEMS, no nlohmann, no C-style headers)
so that pure-math modules like NucleationRate.cc can consume them without
pulling in the full simulation stack. `global.h` re-exports these by including
this header, so any code that already includes `global.h` sees no change.
*/

#ifndef SRC_THAMESLIB_PHYSICALCONSTANTS_H_
#define SRC_THAMESLIB_PHYSICALCONSTANTS_H_

// Ideal gas constant [J/mol/K]
inline const double GASCONSTANT = 8.314;

// Boltzmann constant [J/K] (exact SI-2019 definition)
inline const double BOLTZMANNCONSTANT = 1.380649e-23;

// Avogadro constant [1/mol] (exact SI-2019 definition)
inline const double AVOGADROCONSTANT = 6.02214076e23;

// Pi
inline const double Pi = 3.14159265359;

#endif  // SRC_THAMESLIB_PHYSICALCONSTANTS_H_
