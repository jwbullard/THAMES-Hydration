/**
@file global.h
@brief Declare and assign global variables for THAMES.

*/
#ifndef SRC_THAMESLIB_GLOBAL_H_
#define SRC_THAMESLIB_GLOBAL_H_

// #ifndef DEBUG
// #define DEBUG
// #endif

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <ctime>
#include "../version.h"

#include "../GEMS3K-standalone/GEMS3K/node.h"
#include "../Resources/include/nlohmann/json.hpp"

using json = nlohmann::json;

// String indicating THAMES version in the input microstructure file
inline const char VERSIONSTRING[] = "#THAMES:Version:";

inline const char VERSIONBUGFIX[] = "0";

// String indicating the voxel resolution in micrometers
// for the THAMES input microstructure file
inline const char IMGRESSTRING[] = "#THAMES:Image_Resolution:";

// String indicating the x dimension in voxels
// for the THAMES input microstructure file
inline const char XSIZESTRING[] = "#THAMES:X_Size:";

// String indicating the y dimension in voxels
// for the THAMES input microstructure file
inline const char YSIZESTRING[] = "#THAMES:Y_Size:";

// String indicating the z dimension in voxels
// for the THAMES input microstructure file
inline const char ZSIZESTRING[] = "#THAMES:Z_Size:";

// String indicating the Alite DC name
inline const char AliteDCName[] = "C3S";

// String indicating the Belite DC name
inline const char BeliteDCName[] = "C2S";

// String indicating the Aluminate DC name
inline const char AluminateDCName[] = "C3A";

// String indicating the Ferrite DC name
inline const char FerriteDCName[] = "C4AF";

// String indicating the liquid solution GEM Phase name
// @todo Make this general somehow
inline const char ElectrolyteGEMName[] = "aq_gen";

// String indicating the liquid water DC name
// @todo Make this general somehow
inline const char WaterDCName[] = "H2O@";

// String indicating the CSH GEM Phase name
// @todo Make this general somehow
inline const char CSHGEMName[] = "CSHQ";

// String indicating the AFt GEM Phase name
// @todo Make this general somehow
//  inline const char AFTDCName[] = "ettringite";
inline const char AFTGEMName[] = "ettr";

// String indicating the AFt DC name
// @todo Make this general somehow
inline const char AFTDCName[] = "ettr";

// String indicating the monosulfate GEM Phase name
// @todo Make this general somehow
inline const char MonosulfGEMName[] = "C4AsH14";

// String indicating the monosulfate DC name
// @todo Make this general somehow
//  inline const char MonosulfGEMName[] = "C4AsH14";
inline const char MonosulfDCName[] = "monosulf14";

// String indicating the monocarboaluminate GEM Phase name
// @todo Make this general somehow
inline const char MonocarbGEMName[] = "C4AcH11";

// String indicating the hydrotalcite GEM Phase name
// @todo Make this general somehow
//  inline const char HydrotalcGEMName[] = "OH-hydrotalcite";
inline const char HydrotalcGEMName[] = "hydrotalc-pyro";

// String indicating a generic kinetic model
inline const char GenericType[] = "Generic";

// String indicating the Parrot Killoh Model
inline const char ParrotKillohType[] = "ParrotKilloh";

// String indicating the Standard Dissolution Model
inline const char StandardType[] = "Standard";

// String indicating the C3S Dissolution Model according to
// Bullard et al., Cement and Concrete Research 74(2015)26
inline const char C3SType[] = "C3S_dissKM";

// String indicating the Pozzolanic Reaction Model
inline const char PozzolanicType[] = "Pozzolanic";

// String for the command to convert between image files (Imagemagick)
inline const std::string ConvertCommand = "magick";

// Flag to indicate exiting the program
inline const int QUIT_PROGRAM = 1;

// Flag to indicate simulation of hydration only
inline const int HYDRATION = 2;

// Flag to indicate simulation of leaching only
inline const int LEACHING = 3;

// Flag to indicate simulation of external sulfate only
inline const int SULFATE_ATTACK = 4;

// Flag to indicate kinetic equations handled within code (NOT USED)
inline const int INTERNAL_KINETICS = 0;

// Flag to indicate kinetics is done a priori outside the code (NOT USED)
inline const int EXTERNAL_KINETICS = 1;

// Normal return condition flag (NOT USED)
inline const int RETURN_NORMAL = 0;

// Flag to indicate an element is out of bounds (NOT USED)
inline const int DB_EOB = 1;

// Flag to indicate call to GEM failed (NOT USED)
inline const int GEMRUN_ERROR = 1;

// Flag to indicate that a lattice element is out of bounds (NOT USED)
inline const int MESH_EOB = 1;

// Flag to indicate that an array element is out of bounds (NOT USED)
inline const int ARRAY_EOB = 1;

// Flag to indicate that a file could not be opened (NOT USED)
inline const int FILE_OPEN_ERROR = 2;

// Flag to indicate unexpected end of file (NOT USED)
inline const int PREMATURE_EOF = 3;

// Flag to indicate datum is bad or of wrong kind
inline const int INVALID_INPUT = 4;

// Special phase ids that are important and must always be the same value
inline const int VOIDID = 0;
inline const int ELECTROLYTEID = 1;
inline const int FIRST_SOLID = ELECTROLYTEID + 1;
inline const int NUMCLINKERPHASES = 4;
// inline const double thrPorosity = 0.0355255;

// The number of face, edge, and corner neighbors to a cubic lattice site
inline const int NUM_NEAREST_NEIGHBORS = 6;
inline const int NUM_SECONDNEAREST_NEIGHBORS = 12;
inline const int NUM_THIRDNEAREST_NEIGHBORS = 8;

inline const int NN_NNN = NUM_NEAREST_NEIGHBORS + NUM_SECONDNEAREST_NEIGHBORS;

// Reference temperature for kinetic calculations [K]
inline const double REFTEMP = 298.15;

// Reference lattice resolution [micrometers]
inline const double REFRES = 4.0;

// Ideal gas constant [J/mol/K]
inline const double GASCONSTANT = 8.314;

// Flags for different kinds of probability distributions (NOT USED)
inline const int DELTA = +1;
inline const int UNIFORM = +2;
inline const int GAUSSIAN = +3;
inline const int LOGNORMAL = +4;

//
// Next constant specifies boundary condition configuration
// 0 = periodic boundaries everywhere
// 1 = periodic in y and z, greased in x
// 2 = periodic in x and z, greased in y
// 3 = periodic in x and y, greased in z
//

inline const int BC = 0;

// Growth mode constants
inline const unsigned int DLA = 0;
inline const unsigned int WMC = 1;
inline const unsigned int IWMC = 2;

// Saturation value for colors
inline const int COLORSATVAL = 255;

// Flags for data formats (NOT USED)
inline const int GEMSFORMAT = 0;
inline const int EXTERNALFORMAT = 0;

inline const double ICTHRESH = 1.0e-9;
inline const double limitICTHRESH = 1.0e-3;

inline const double Pi = 3.14159265359;

inline const double S_PER_MINUTE = 60.0;
inline const double S_PER_H = 3600.0;
inline const double S_PER_DAY = 86400.0;
inline const double S_PER_YEAR = 3.15360000e7;
inline const double H_PER_DAY = 24.0;
inline const double DAY_PER_YEAR = 365.0;

#endif // SRC_THAMESLIB_GLOBAL_H_
