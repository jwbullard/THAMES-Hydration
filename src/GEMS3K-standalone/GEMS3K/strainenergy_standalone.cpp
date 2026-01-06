/**
@file strainenergy_standalone.cpp
@brief Provides a default strainenergy vector for standalone GEMS3K builds.

When GEMS3K is built standalone (without the full THAMES application),
this file provides a zero-initialized strainenergy vector. The actual
strain energy values are only meaningful when running the full THAMES
hydration simulation with elastic calculations enabled.
*/

#include <vector>

// Default strainenergy vector for standalone GEMS3K builds
// This will be zero-initialized, meaning no strain energy contribution
// to the Gibbs free energy calculations
std::vector<double> strainenergy;
