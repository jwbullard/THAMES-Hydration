/**
@file Lattice.h
@brief Declaration of the Lattice class for storing the 3D microstructure

THAMES defines a Lattice class that is instantiated to a Lattice
object at the beginning of the program's execution.  The lattice defines the
three-dimensional environment within which a cement paste microstructure
exists, hydrates, and possibly deteriorates.
*/

#ifndef SRC_THAMESLIB_LATTICE_H_
#define SRC_THAMESLIB_LATTICE_H_

#include "global.h"
#include "Exceptions.h"
#include "AppliedStrain.h"
#include "ChemicalSystem.h"
#include "Interface.h"
#include "Isite.h"
#include "RanGen.h"
#include "Site.h"
#include "utils.h"
// #include "../version.h"

#include <climits>

/**
@struct Sitesize
@brief Structure to catalog site domain sizes
*/
struct Sitesize {
  int siteid; /**< ID of the site in the site_ vector */
  int nsize;  /**< Size of the domain of the phase at that site */
};

struct chemElem {
  int z;
  std::string symb;
  double mass;
};

struct structGrowVect {
  int id;
  int posVect;
  int affinityInt;
};

struct structDissVect {
  int id;
  int posVect;
  int wmcInt;
};

/**
@class Lattice
@brief Defines and stores the 3D microstructure as a discrete lattice of voxel
sites.

*/
class Lattice {

private:
  std::string version_; /**< THAMES version for header information */
  std::string thamesVersion_;
  std::string jobRoot_; /**< The root name for output files */
  std::string damageJobRoot_;

  RanGen *rg_; /**< Pointer to random number generator object */
  int latticeRNGseed_;         /**< the seed of the random number
                                    generator */
  long int numRNGcall_0_;      /**< the first number used to keep
                                    the RNG calls track*/
  long int numRNGcallLONGMAX_; /**< the second number used to keep
                                    the RNG calls track */
  double lastRNG_;             /**< last generated random number */

  int xdim_;               /**< Number of sites in the x dimension */
  int ydim_;               /**< Number of sites in the y dimension */
  int zdim_;               /**< Number of sites in the z dimension */
  double resolution_;      /**< Voxel edge length [micrometers] */
  std::vector<Site> site_; /**< 1D list of Site objects (site = voxel) */
  int numSites_;           /**< Total number of sites */

  ChemicalSystem *chemSys_; /**< Pointer to simulation's ChemicalSystem */
  AppliedStrain *FEsolver_; /**< Pointer to simulation's FE elastic solver */
  std::vector<Interface> interface_; /**< List of the different interface objects
                                          in the microstructure */

  double areaPerFace_;            /**< Converts a voxel face to m2 units */
  double volumePerVoxel_;         /**< Converts a voxel to its volume in m3
                                       units */
  double wsRatio_;                /**< Water-to-solids mass ratio */
  std::vector<double> volumeFraction_;      /**< Array of volume fractions of each
                                                 microstructure phase */
  std::vector<double> surfaceArea_;         /**< Array of surface areas of each
                                                 microstructure phase
                                                 (m2 per 100 g of all solid) */
  std::vector<double> specificSurfaceArea_; /**< Array of specific surface areas
                                                 of each microstructure phase
                                                 (m2 per kg of that phase) */
  std::vector<int> count_;                  /**< Number of sites of each different
                                                 type */

  map<int, std::vector<double>> expansion_;      /**< Map of expansion strain of
                                                      each voxel */
  map<int, std::vector<int>> expansion_coordin_; /**< Map of coordinates of sites
                                                      with local expansion strain */
  double waterChange_;                      /**< How much water must be added or
                                                 subtracted due to hydration or
                                                 deterioration */
  double microstructureVolume_;             /**< Microstructure volume in GEM
                                                 volume units */
  double initialMicrostructureVolume_;      /**< Initial microstructure volume in
                                                 GEM volume units */
  double voxelPoreVolume_;                  /**< Total volume of voxel pores */
  double voxelPoreVolumeFraction_;          /**< Total volume fraction of voxel
                                                 pores */
  double voxelPoreVolumeFractionSaturated_; /**< Total volume fraction of
                                                 saturated voxel pores on
                                                 microstructure volume basis*/
  double subvoxelPoreVolume_;               /**< Total volume of subvoxel pores */
  double nonSolidVolume_;                   /**< Total volume not solid */
  double solidVolumeWithPores_;             /**< Total solid volume including their
                                                 internal pore volume */
  double waterVolume_;                      /**< volume of electrolyte in GEM
                                                 volume units */
  double voidVolume_;                       /**< volume of void in GEM volume
                                                 units */
  double voxelWaterVolume_;                 /**< Volume of voxel pore water */
  double voxelVoidVolume_;                  /**< Volume of voxel void space
                                                 (no water) */
  double subvoxelWaterVolume_;              /**< Volume of water in subvoxel
                                                 pores in GEM units */
  double subvoxelPoreVolumeFraction_;       /**< Total volume fraction of subvoxel
                                                 pores */
  double subvoxelPoreVolumeFractionSaturated_; /**< Total volume fraction of
                                                    saturated subvoxel pores on
                                                    microstructure volume basis*/

  std::vector<struct PoreSizeData> masterPoreSizeDist_; /**< Pore size distribution
                                                             and saturation */

  double time_;              /**< The current simulation time [h] */
  double temperature_;       /**< The current simulation temperature [K] */
  double oldtemp_;           /**< The temperature in the previous
                                  time step [K] */
  double sulfateAttackTime_; /**< Simulation time at which to begin
                                  simulation of sulfate attack [h] */
  double leachTime_;         /**< Simulation time at which to begin
                                  simulation of leaching [h] */

  bool depthEffect_; /**< Whether or not PNG images should have
                          depth effect */
  bool verbose_;     /**< Flag to determine verbose output */
  bool warning_;     /**< Flag to determine warning message output */

  std::vector<chemElem> cfgElem_; /**< Holds periodic table information to output
                                       files in cfg format */

  double initSolidMass_;          /**< the initial solid mass of the system*/

  double wcRatio_;                /**< Water-to-cement mass ratio */

  int numMicroPhases_;            /**< Number of microphases */

  double particRadius_;           /**< used for graphical representation */

  std::vector<int> growthInterfaceSize_;      /**< growth interface size of each 
                                                   microphase */
  std::vector<int> dissolutionInterfaceSize_; /**< dissolution interface size of each
                                                   microphase */

  std::vector<int> growingVectSA_;        /**< for SULFATE ATTACK: contains all
                                               microPhaseIds growing due to SA attack */
  int sizeGrowingVectSA_;                 /**< size of growingVectSA_ vector*/

  std::vector<std::vector<int>> shrinking_;    /**< for each microPhaseId in growingVectSA_,
                                                    all the microDhaseIds that can transform
                                                    into this one*/
  std::vector<std::vector<double>> volratios_; /**< for each microPhaseId in growingVectSA_
                                                    and all corresponding microPhaseIds in
                                                    shrinking_, contains the molar volume
                                                    ratios of the corresponding 
                                                    microPhases */

  int waterDCId_;           /**< the DCId coresp to DCName = "H2O@" */
  double waterMolarMass_;   /**< the water molar mass corresp. to waterDCId_ */
  double waterMolarVol_;    /**< the water molar volume corresp. to waterDCId_ */

  int electrolyteIntPorosity_;
  int voidIntPorosity_;
  int convFactDbl2IntPor_;

  std::vector<std::vector<int>> affinityInt_;
  std::vector<std::vector<bool>> growthTemplate_;
  std::vector<int> microPhasePorosityInt_;

  double oneFaceAreaPerHundredGramSolid_; /** surface area of one voxel's face per 100g of
                                              the initial solid mass of the system*/

  // int DAMAGEID_;

public:
  /**
  @brief Constructor without input microstructure file name.

  This constructor simply initializes the dimensions and time to zero, sets
  the temperature to the globally defined reference temperature, and
  sets the lattice resolution to the globally defined reference value.

  @note Not currently used in THAMES

  @param cs is a pointer to the ChemicalSystem object for the simulation
  */
  Lattice(ChemicalSystem *cs);

  /**
  @brief Overloaded constructor with input microstructure file name.

  This constructor initializes the dimensions and time to zero, sets
  the temperature to the globally defined reference temperature, and
  sets the lattice resolution to the globally defined reference value.
  Afterward, the input microstructure file is opened and read, so that
  the voxel phase assignments can be made at each site.

  @param cs is a pointer to the ChemicalSystem object for the simulation
  @param rg is a pointer to the random number generator object
  @param seedRNG is the random number seed
  @param fileName is the name of the file containing the microstructure data
  @param verbose is true if extra messages are to be printed
  @param warning is true if warning messages are to be printed
  */
  Lattice(ChemicalSystem *cs, RanGen *rg, int seedRNG, const std::string &fileName,
          const bool verbose, const bool warning);

  /**
 @brief Destructor.

 This destructor clears out the `interface_` and `site_` vectors, and
 also deletes the allocated memory for the random number generator object,
 since this is the class that allocated the memory for that object.
 */
  ~Lattice();

  void setChemElem(std::string sval);
  chemElem getChemElem(int ival);

  /**
  @brief Set the number of sites in the x dimension.

  @param x is the number of sites in the x dimension
  */
  void setXDim(const int x) {
    xdim_ = x;
    numSites_ = (xdim_ * ydim_ * zdim_);
  }

  /**
  @brief Get the number of sites in the x dimension.

  @return the number of sites in the x dimension
  */
  int getXDim() const { return xdim_; }

  /**
  @brief Set the number of sites in the y dimension.

  @param y is the number of sites in the y dimension
  */
  void setYDim(const int y) {
    ydim_ = y;
    numSites_ = (xdim_ * ydim_ * zdim_);
  }

  /**
  @brief Get the number of sites in the y dimension.

  @return the number of sites in the y dimension
  */
  int getYDim() const { return ydim_; }

  /**
  @brief Set the number of sites in the z dimension.

  @param z is the number of sites in the z dimension
  */
  void setZDim(const int z) {
    zdim_ = z;
    numSites_ = (xdim_ * ydim_ * zdim_);
  }

  /**
  @brief Get the number of sites in the z dimension.

  @return the number of sites in the z dimension
  */
  int getZDim() const { return zdim_; }

  /**
  @brief Get the total number of lattice sites.

  The lattice is rectangular, so the total number of sites is
  `xdim_ * ydim_ * zdim_`, but we store this value as a class member to
  save having to compute it multiple times.

  @return the total number of lattice sites
  */
  int getNumSites() const { return numSites_; }

  /**
  @brief Set the volume fraction of a given microstructure phase.

  @param i is the index of the microstructure phase
  @param vfrac is the volume fraction to assign on a total microstructure basis
  */
  // void setVolumeFraction(const unsigned int i, const double vfrac) {
  //   if (i > -1 && i < volumeFraction_.size()) {
  //     volumeFraction_[i] = vfrac;
  //   } else {
  //     throw EOBException("Lattice", "setVolumeFraction", "volumeFraction_",
  //                        volumeFraction_.size(), i);
  //   }
  // }

  /**
  @brief Set the initial volume fraction of a given microstructure phase.

  @param i is the index of the microstructure phase
  @param vfrac is the volume fraction to assign on a total microstructure basis
  */
  // void setInitVolumeFraction(const int i, const double vfrac) {
  //   if (i > -1 && i < initVolumeFraction_.size()) {
  //     initVolumeFraction_[i] = vfrac;
  //   } else {
  //     throw EOBException("Lattice", "setInitVolumeFraction",
  //                        "initVolumeFraction_", initVolumeFraction_.size(),
  //                        i);
  //   }
  //   initVolumeFraction_[i] = vfrac;
  // }

  /**
  @brief Set the water-solids mass ratio

  @param ws is the water-solids mass ratio
  */
  void setWsRatio(const double ws) {
    wsRatio_ = 0.0;
    if (ws > 0.0) {
      wsRatio_ = ws;
    }
    return;
  }

  /**
  @brief Get the water-solids mass ratio

  @return the water-solids mass ratio
  */
  double getWsRatio(void) const { return wsRatio_; }

  /**
  @brief Get the water-cement mass ratio

  @return the water-cement mass ratio
  */
  double getWcRatio(void) const { return wcRatio_; }

  /**
  @brief Get the volume fraction of a given microstructure phase.

  This is simply the number of sites with a given phase divided by the
  total number of sites.

  @param i is the index of the microstructure phase
  @return the volume fraction of phase i on a total microstructure basis
  */
  double getVolumeFraction(int i) { return (volumeFraction_[i]); }

  /**
  @brief Calculate the subvoxel pore volume

  @param vol is the array of all microstructure phase volumes
  */
  // void calcSubvoxelPoreVolume(std::vector<double> &vol);

  /**
  @brief Calculate the total volume of solids including
  subvoxel pore volume assigned to solids

  @param vol is the array of all microstructure phase volumes
  it
  */
  // void calcSolidVolumeWithPores(std::vector<double> &vol);

  /**
  @brief Get the total volume of solids including
  subvoxel pore volume assigned to solids

  @return the solid volume including subvoxel pore volume
  */
  double getSolidVolumeWithPores(void) const { return solidVolumeWithPores_; }

  /**
  @brief Calculate the non-solid volume

  @param vol is the array of all microstructure phase volumes
  it
  */
  // void calcNonSolidVolume(std::vector<double> &vol);

  /**
  @brief Get or calculate the non-solid volume

  @return the non-solid volume
  */
  double getNonSolidVolume(void) const { return nonSolidVolume_; }

  /**
  @brief Get the number of neighbor sites each site has.

  This is simply the number of sites with a given phase divided by the
  total number of sites.

  @note NOT USED.

  @return the number of neighbor sites each site has
  */
  // unsigned int getSiteNeighbors() const { return siteNeighbors_; }

  /**
  @brief Set the lattice resolution [meters].

  The lattice resolution is the physical length associated with the edge
  length of a site.

  @param res is the lattice resolution [meters]
  */
  void setResolution(const double res);

  /**
  @brief Get the lattice resolution [meters].

  The lattice resolution is the physical length associated with the edge
  length of a site.

  @note NOT USED.

  @return the lattice resolution [meters]
  */
  double getResolution() const { return resolution_; }

  /**
  @brief Get the areaPerFace_ value [m2].

  @return the areaPerFace_ value [m2]
  */
  double getAreaPerFace() const { return areaPerFace_; }

  /**
  @brief Get the volumePerVoxel_ value [m3].

  @return the volumePerVoxel_ value [m3]
  */
  double getVolumePerVoxel() const { return volumePerVoxel_; }

  /**
  @brief Set the simulation time [hours].

  @note NOT USED.

  @param tval is the simulation time [hours]
  */
  void setTime(const double tval) { time_ = tval; }

  /**
  @brief Get the simulation time [hours].

  @note NOT USED.

  @return the simulation time [hours]
  */
  double getTime() const { return time_; }

  /**
  @brief Get the simulation time at which to start sulfate attack simulation
  [hours].

  @note NOT USED.

  @return the simulation time at which to start sulfate attack [hours]
  */
  double getSulfateAttackTime() const { return sulfateAttackTime_; }

  /**
  @brief Set the simulation time at which to start sulfate attack simulation
  [hours].

  @param sattacktime is the simulation time at which to start sulfate attack
  [hours]
  */
  void setSulfateAttackTime(const double sattacktime) {
    sulfateAttackTime_ = sattacktime;
  }

  /**
  @brief Get the simulation time at which to start leaching simulation [hours].

  @note NOT USED.

  @return the simulation time at which to start leaching [hours]
  */
  double getLeachTime() const { return leachTime_; }

  /**
  @brief Set the simulation time at which to start leaching simulation [hours].

  @param leachtime is the simulation time at which to start leaching [hours]
  */
  void setLeachTime(const double leachtime) { leachTime_ = leachtime; }

  /**
  @brief Set the lattice temperature [K].

  @param tmp is the temperature [K]
  */
  void setTemperature(const double tmp) { temperature_ = tmp; }

  /**
  @brief Get the lattice temperature [K].

  @return the temperature [K]
  */
  double getTemperature() const { return temperature_; }

  /**
  @brief Get the version of THAMES

  @note NOT USED.

  @return the version number as a string
  */
  const std::string &getVersion() const { return version_; }

  /**
  @brief Set the root name for simulation output files.

  @param jobname is the root name for simulation output files
  */
  void setJobRoot(std::string jobname) {
    jobRoot_ = jobname;
    // damageJobRoot_ = jobRoot_ + ".damage";
  }

  /**
  @brief Add a site at location (xp,yp,zp) to the lattice.

  The site is checked for valid coordinates.  If valid a new Site object
  is created and pushed back onto the class's `site_` vector.

  @param xp is the x coordinate of the site to add
  @param yp is the y coordinate of the site to add
  @param zp is the z coordinate of the site to add
  */
  void addSite(const int xp, const int yp, const int zp);

  /**
  @brief Get the x coordinate of a site with a given index in the 1D `site_`
  array.

  @param i is the index of the site in the class's `site_` array
  @return the x coordinate
  */
  int getX(const int i) const { return (site_[i].getX()); }

  /**
  @brief Get the y coordinate of a site with a given index in the 1D `site_`
  array.

  @param i is the index of the site in the class's `site_` array
  @return the y coordinate
  */
  int getY(const int i) const { return (site_[i].getY()); }

  /**
  @brief Get the z coordinate of a site with a given index in the 1D `site_`
  array.

  @param i is the index of the site in the class's `site_` array
  @return the x coordinate
  */
  int getZ(const int i) const { return (site_[i].getZ()); }

  /**
  @brief Get a site's index in the 1D `site_` array, given its (x,y,z)
  coordinates.

  @param ix is the x coordinate of the site
  @param iy is the x coordinate of the site
  @param iz is the x coordinate of the site
  @return the index of the site in the `site_` array
  */
  int getIndex(int ix, int iy, int iz) const;

  /**
  @brief Get the collection of site indices neighboring a given site.

  @param sitenum is the index of the site in question
  @param size is the maximum distance defining the neighborhood [sites]
  @return a list of site indices for all neighbors within the maximum distance
  */
  std::vector<int> getNeighborhood(const int sitenum, const int size);

  /**
  @brief Get the collection of site indices (ids), site phases (phs)
  neighboring a given site (sitenum). Calculate also the number of neighbor sites 
  occupied by electrolyte (numW) and the total porosity (totpor). used for
  sulfate attack (SA)

  */
  void getNeighborhood(const int sitenum, std::vector<int> & ids,
                       std::vector<int> & phs, int & numW,
                       int & totPor);

  /**
  @brief Set the collection of site indices neighboring a given site 
  including the site itself. This ids used for sulfate attack (SA)

  */
  void setNeighborhoodSA(void);

  /**
  @brief Get a pointer to a Site object at a given index in the `site_` array.

  @param index is the index of the Site object in the `site_` array
  @return a pointer to the Site object in question
  */
  Site *getSite(int index) { return &site_[index]; }

  /**
  @brief Designate a site as damaged.

  The site to be damaged is specified by its index in the `site_` array.

  @note NOT USED.

  @param index is the index of the Site object in the `site_` array
  */
  void setDamage(int index) { site_[index].setDamage(); }

  /**
  @brief Change the wmc (weighted mean curvature) of a site by a prescribed
  amount.

  @param index is the index of the Site object in the `site_` array
  @param dwmcval is the increment to add to the wmc
  */
  // void dWmc(int index, double dwmcval) {
  //   site_[index].setWmc(site_[index].getWmc() + dwmcval);
  // }

  /**
  @brief set the internal porosity, wmc0_, of a site i.e. its contribution to
  the "weighted mean curvature", wmc_, of the site

  @param index is the index of the Site object in the `site_` array
  @param wmc0val is the value of wmc0_ to assign to the site
  */
  void setWmc0(int index, int wmc0val) { site_[index].setWmc0(wmc0val); }

  /**
  @brief Compute normalized initial microstructure phase masses

  Given the initial masses of all phases in 1 cm3 of microstructure,
  this method scales them to 100 grams of solid instead.  In the process,
  this method also sets the initial moles of water in the
  chemical system definition.

  So when this method is finished, microPhaseMass has units of g per 100 g of
  solid

  @param microPhaseMass is a vector of all the microstructure masses
  @param cementMass is the combined mass of all the cement components
  @param solidMass is the combined mass of all the solids
  */
  void normalizePhaseMasses(std::vector<double> microPhaseMass);

  /**
  @brief Master method to locate the interfaces for each phase in the
  microstructure.

  */
  void findInterfaces(void);

  /**
  @brief Determine if a voxel has at least one neighbor that is a porousorous
  solid

  @param siteId is the id of the site to check
  @param neighborRange is the number of neighbors to check
  @return true if at least one neighbor is a porous solid
  */
  bool hasPorousSolidNeighbor(const int siteID, const int neighborRange);

  /**
  @brief Add (grow i.e. switch from electrolyte) the prescribed number of
  sites of each microphase that has to grow

  @param growPhaseIDVect is the vector of microphase IDs that must grow
  @param numSiteGrowVect is a vector containing the number of voxels to add for
  each microphase ID in growPhaseIDVect
  @param growPhNameVect is a vector containing the name of each microphase in
  growPhaseIDVect
  @param numtoadd_G is the number of sites switched by this call
  @param nucOK is a flag to characterize the nucleation correctness
  @param totalTRC is the total call number of the changeMicrostructure method
  @return the actual number of sites that were changed for each microphase ID
  from the input growPhaseIDVect vector
  */
  std::vector<int> growPhase(std::vector<int> growPhaseIDVect,
                             std::vector<int> numSiteGrowVect,
                             std::vector<std::string> growPhNameVect,
                             int &numadded_G, bool &nucOK, int totalTRC);

  /**
  @brief create a new growth interface by nucleation of numToNucleate sites for a
  given phase (phaseID); this is necessary when the "growth" of all requested
  sites for this phase was not possible because the size of the corresponding
  growth interface was zero. A nucleation site is chosen according to a probability
  computed taking into account the microPhaseId affinity toward the site surrounding.

  @param phaseid is the id of the microstructure phase to nucleate
  @param numToNucleate is the number of sites to nucleate/create for this phase
  @return the number of nucleation events (nuclei) that have been realised
  */
  int nucleatePhaseAff(const int phaseID, const int numToNucleate);

  /**
  @brief create a new growth interface by nucleation of numToNucleate sites for a
  given phase (phaseID); this is necessary when the "growth" of all requested
  sites for this phase was not possible because the size of the corresponding
  growth interface was zero. A nucleation site is chosen with equal probability
  among the electrolyte sites.

  @param phaseid is the id of the microstructure phase to nucleate
  @param numToNucleate is the number of sites to nucleate/create for this phase
  @return the number of nucleation events (nuclei) that have been realised
  */
  int nucleatePhaseRnd(const int phaseID, const int numToNucleate);

  /**
  @brief Remove (dissolve i.e. switch to electrolyte) the prescribed number of
  sites of each microphase that has to dissolve

  @param dissPhaseIDVect is the vector of microphase IDs that must dissolve
  @param numSiteDissVect is a vector containing the number of voxels to dissolve
  for each microphase ID in dissPhaseIDVect
  @param dissPhNameVect is a vector containing the name of each microphase in
  dissPhaseIDVect
  @param numtoadd_D is the number of sites switched by this call
  @param totalTRC is the total call number of the changeMicrostructure method
  @return vector of the number of voxels of each phase that could not dissolve
  */
  std::vector<int> dissolvePhase(std::vector<int> dissPhaseIDVect,
                                 std::vector<int> numSiteDissVect,
                                 std::vector<std::string> dissPhNameVect,
                                 int &numadded_D, int totalTRC);

  /**
  @brief Change the volume fraction of void space

  This is a pass-through function that alters the volume
  fraction of void space in the microstructure, either adding
  or subtracting electrolyte as necessary

  @param voidFracToAdd is the target change in void volume
  fraction, on a total microstructure volume basis
  */
  double changeSaturationState(double voidFracToAdd);

  /**
  @brief Remove a prescribed volume fraction of water from the system

  This master method empties a given volume fraction of electrolyte
  from the system, starting with the largest voxel pores and moving
  to smaller and smaller saturated pores until either (a) the prescribed
  volume is reached or (b) the system is completely dry.

  If the number of solution-filled sites (voxels) is insufficient to balance the
  loss of electrolyte, then we must visit the sub-voxel porosity and draw
  remaining electrolyte from it.

  This is a master method that splits the task into two parts: (a)
  saturated pore voxels and, if necessary (b) saturated sub-voxel pores.

  @param volFracToRemove is the targeted reduction in electrolyte volume
  fraction, on a total microstructure volume basis
  @return the actual volume fraction of electrolyte removed
  */
  double emptyPorosity(double volFracToRemove);

  /**
  @brief Convert a prescribed number of electrolyte voxels to void

  @param numToEmpty is the number of electrolyte voxels to remove, on a total
  microstructure volume basis
  @return the actual number converted
  */
  int emptyVoxelPorosity(int numToEmpty);

  /**
  @brief Empty a prescribed volume fraction of electrolyte from sub-voxel pores

  @param volFracToRemove is the target reduction in electrolyte volume
  fraction, on a total microstructure volume basis
  @return the actual volume fraction of electrolyte removed
  */
  double emptySubVoxelPorosity(double volFracToRemove);

  /**
  @brief Add a prescribed volume fraction of electrolyte to the system

  This master method adds a given volume fraction of electrolyte
  to the system, starting with the smallest subvoxel pores and moving
  to larger and larger unsaturated pores until either (a) the prescribed
  volume is reached or (b) the system is completely full.

  This is a master method that splits the task into two parts: (a)
  unsaturated sub-voxel pores and, if necessary (b) void voxels

  @param volFracToAdd is the target increase in electrolyte volume
  fraction, on a total microstructure volume basis
  @return the actual volume fraction of electrolyte added
  */
  double fillPorosity(double volFracToAdd);

  /**
  @brief Convert a prescribed number of void voxels to electrolyte voxels

  @param numToFill is the number of electrolyte voxels to fill, on a total
  microstructure volume basis
  @return the actual number converted
  */
  int fillVoxelPorosity(int numToFill);

  /**
  @brief Empty a prescribed volume fraction of electrolyte from sub-voxel pores

  @param volFracToAdd is the target increase in electrolyte volume fraction, on
  a total microstructure volume basis
  @return the actual volume fraction of electrolyte removed
  */
  double fillSubVoxelPorosity(double volFracToAdd);

  /**
  @brief Remove the water from a prescribed number of solution-filled sites.

  This method constructs a list of all the <i>potential</i> void sites, based
  on whether there are multiple connected solution-filled sites in a cluster.
  The list is then sorted essentially by the effective pore size.  Only then
  is the list visited and the prescribed number of sites switched to void.

  @param numsites is the number of sites to switch from water to void
  @return the actual number of sites that were changed
  */
  // int emptyPorosity(int numsites, int cyc);

  /**
  @brief Add water to a prescribed number of empty pore sites.

  This method constructs a list of all the void sites, based
  on whether there are multiple connected void sites in a cluster.
  The list is then sorted essentially by the effective pore size.  Only then
  is the list visited and the prescribed number of sites switched to water.

  @param numToFill is the number of sites to switch from void to water
  @return the actual number of sites that were changed
  */
  int fillPorosity(int numToFill);

  /**
  @brief Add water to the all empty pore sites.

  @param cyc (cycle) is the iteration number in main iteration loop
  @return the added number of moles of water
  */
  double fillAllPorosity(const int cyc);

  /**
  @brief Count the number of solution sites within a box centered on a given
  site.

  @param boxsize is the linear dimension of the cubic box neighborhood
  @param siteid is the index of the queried site in the `site_` array
  @return the number of solution-filled sites in the box neighborhood
  */
  int countBox(int boxsize, unsigned int siteid);

  /**
  @brief Check whether a linear coordinate is outside the lattice boundaries.

  If a given coordinate is outside the lattice boundaries, then the additive
  adjustment is returned that will locate the equivalent site within the
  lattice, assuming periodic boundary conditions.

  @param pos is the linear coordinate to check
  @param size is the dimension of the lattice in that dimension (number of
  sites)
  @return the additive adjustment to locate the equivalent coordinate within
  the lattice
  */
  int checkBC(int pos, int size) {
    if (pos >= size)
      return (-size);
    if (pos < 0)
      return (size);
    return (0);
  }

  /**
  @brief Get a pointer to the ChemicalSystem object for the simulation.

  @return a pointer to the ChemicalSystem object for the simulation
  */
  ChemicalSystem *getChemSys() const { return chemSys_; }

  /**
  @brief Set the phase id of a given site, specified by a pointer to the Site
  object.

  @param s is a pointer to the Site object
  @param i is the phase index to set at that site
  */
  void setMicroPhaseId(Site *s, const int i) {
    count_[s->getMicroPhaseId()]--;
    s->setMicroPhaseId(i);
    count_[i]++;
  }

  /**
  @brief Set the phase id of a given site, specified by the site's index number.

  @param sitenum is the index of the site in the `site_` array
  @param i is the phase index to set at that site
  */
  void setMicroPhaseId(const int sitenum, const int i) {
    count_[site_[sitenum].getMicroPhaseId()]--;
    site_[sitenum].setMicroPhaseId(i);
    count_[i]++;
  }

  /**
  @brief Get the phase id of a given site, specified by the site's index number.

  @param sitenum is the index of the site in the `site_` array
  @return the microstructure phase id at the site
  */
  int getMicroPhaseId(const int sitenum) {
    return (site_[sitenum].getMicroPhaseId());
  }

  /**
  @brisef Add a site to the list of sites where dissolution of a given phase can
  occur.

  @param loc is a pointer to the Site object to add to the list of potential
  dissolution sites
  @param pid is the microstructure phase id
  */
  void addDissolutionSite(Site *loc, int pid);

  /**
  @brief Add a site to the list of sites where growth of a given phase can
  occur.

  @param loc is a pointer to the Site object to add to the list of potential
  growth sites
  @param pid is the microstructure phase id
  */
  void addGrowthSite(Site *loc, int pid);

  /**
  @brief Remove a site from the list of sites where dissolution of a given phase
  can occur.

  @param loc is a pointer to the Site object to remove from the list of
  potential dissolution sites
  @param pid is the microstructure phase id
  */
  void removeDissolutionSite(Site *loc, int pid);

  /**
  @brief Remove a site from the list of sites where growth of a given phase can
  occur, during a dissolution process.

  @param loc is a pointer to the Site object to remove from the list of
  potential growth sites
  @param pid is the microstructure phase id
  */
  void removeGrowthSite_diss(Site *loc, int pid);

  /**
  @brief Remove a site from the list of sites where growth of a given phase can
  occur, during a growth process.

  @param loc is a pointer to the Site object to remove from the list of
  potential growth sites
  @param pid is the microstructure phase id
  */
  void removeGrowthSite_grow(Site *ste0, int pid);

  /**
  @brief Remove a site from the all growth interfaces during a nucleation process.

  @param loc is a pointer to the Site object to remove from the list of
  potential growth sites
  */
  void removeGrowthSite_nucleation(Site *loc);

  /**
  @brief Master method to update a microstructure during after a given time
  interval.

  Updating the microstructure includes determining how many sites of each phase
  to add and subtract from the lattice, determining which site locations will be
  used to do that, and then actually causing the switches in phase id to happen
  at those sites. The interfaces and lists of dissolution and growth sites are
  updated accordingly, too.

  @note Water is assumed to be chemically reactive only if it is in voxel
  porosity (microstructure id ELECTROLYTEID).  If the voxel water is
  exhausted then some reaction can still happen with water in nanoporosity, but
  for now we assume that the nanopore water is chemically unreactive and cannot
  be removed.

  @todo Generalize to allow water in nanopores to be chemically reactive

  @param time is is the simulation time [hours]
  @param simtype is the type of simulation (hydration, leaching, etc)
  @param vectPhNumDiff is the vector of maximum number of voxels belonging to
  each microphase, voxels that can be dissolved according to the system
  configuration (lattice)
  @param vectPhIdDiff is the microphase ID for which a the number of voxels that
  can be dissolved is smaller than the number requested by the corresponding
  kinetic model
  @param vectPhNameDiff is the vector of names of microphases
  @param recalls counts the number of changeMicrostructure calls for a given
  cycle (cyc)
  @param stopPrg is a flag used to stop the program when the possible number of
  nucleation events is smaller than the nuclei requested number
  @param cyc (cycle) is the iteration number in main iteration loop in
  Controller::doCycle - each cycle corresponds to a time step

  @return zero if okay or nonzero if not all requested voxels
  for a certain microphase ID (phDiff) can be dissolved
  */
  int changeMicrostructure(double time, const int simtype,
                           std::vector<int> &vectPhNumDiff,
                           std::vector<int> &vectPhIdDiff,
                           std::vector<std::string> &vectPhNameDiff,
                           int recalls, bool &stopPrg, int cyc);

  /**
  @brief Adjust GEMS calculated volumes of microstructure phases

  The volume fractions passed to this function are those coming directly
  from the chemical system.  But the chemical system does not account for
  occluded porosity that may be associated with a solid phase at length
  scales smaller than the lattice spatial resolution.  This method fixes
  those volume fractions, paying special attention to the water distribution.

  @note Water is assumed to be chemically reactive only if it is in voxel
  porosity (microstructure id ELECTROLYTEID).  If the voxel water is
  exhausted then some reaction can still happen with water in nanoporosity, but
  for now we assume that the nanopore water is chemically unreactive and cannot
  be removed.

  @todo Generalize to allow water in nanopores to be chemically reactive

  @param vol is a vector of the pre-adjusted microstructure volumes that come
  from GEMS (not based on voxels)
  @param volSize is the number of elements in the vol vector
  */
  void adjustMicrostructureVolumes(std::vector<double> &vol, int volSize, int cyc);

  /**
  @brief Calculate microstructure volume fractions

  @param names is a vector of the adjusted microstructure volumes
  @param vol is a vector of the pre-adjusted microstructure volumes that come
  from GEMS (not based on voxels)
  @param vfrac will hold the microstructure volume fractions
  @param volSize is the number of elements in the vol vector 
  */
  void adjustMicrostructureVolFracs(std::vector<std::string> &names,
                                    const std::vector<double> vol,
                                    std::vector<double> &vfrac, int volSize);

  /**
  @brief Calculate the pore size distribution data

  */
  void calculatePoreSizeDistribution(void);

  /**
  @brief Get the pore size distribution contribution of each
  phase to the whole microstructure

  This function returns a matrix. Each column corresponds to one of the
  microstructure phases that posssesses an internal porosity and each
  row of a column correponds to a volume fraction for a given pore diameter

  @return A 2D vector of doubles
  */
  std::vector<std::vector<struct PoreSizeData>> getPhasePoreSizeDistributions(void);

  /**
  @brief Get the maximum phase pore diameter from among all the defined
  pore diameters for each microstructure phase

  @param phasePoreSizeDist is a 2D matrix representing the pore size
  distribution of each phase
  @return The maximum pore diameter found (nm) within any phase
  */
  double getMaxPhasePoreDiameter(
      const std::vector<std::vector<struct PoreSizeData>> phasePoreSizeDist);

  /**
  @brief Calculate the overall pore size distribution of the
  microstructure, `masterPoreSizeDist_`, as a histogram with pre-defined
  diameters

  @param histogramDiameters is a vector of diameters (nm) to use for binning
  @param phasePoreSizeDist is a 2D matrix representing the pore size
  distribution of each phase
  */
  void calcMasterPoreSizeDist(
      const std::vector<double> histogramDiameters,
      const std::vector<std::vector<struct PoreSizeData>> phasePoreSizeDist);
  /**
  @brief Set up the binned pore diameters for master pore size distribution

  @return A vector of doubles corresponding holding the binning diameters (nm)
  */
  std::vector<double> setPSDiameters(void);

  /**
  @brief Get the pore volume fractions associated with each phase

  This will be the volume fraction of the microstructure that is occupied
  by pores (no matter the size) of each phase.

  @return A vector of doubles corresponding to the pore volume fractions of each
  phase
  */
  std::vector<double> getPoreVolumeFractions(void);

  /**
  @brief Write the pore size distribution data to a file

  @param curtime is the current time in hours
  @param timeString is the current time resolved into y,d,h,m
  */
  void writePoreSizeDistribution(const double curtime,
                                 const string timeString);

  /**
  @brief Write the microstructure colors to a file

  This is done to save processing the simparams.json file just to get the colors
  and will make post-processing of images easier.

  */
  void writeMicroColors();

  /**
  @brief Write the 3D microstructure to a file.

  The microstructure output file will indicate the phase id at each site.

  @param timeString is the current time resolved into y,d,h,m
  */
  void writeLattice(const std::string timeString);

  /**
  @brief Write the 3D microstructure to a xyz file.

  @param curtime is the current time in hours
  @param timeString is the current time resolved into y,d,h,m
  */
  void writeLatticeXYZ(const double curtime, const std::string timeString);

  /**
  @brief Write one by one the 3D microstructure into the same xyz file (append).

  @param curtime is the current time in hours
  */
  void appendXYZ(double curtime);

  /**
  @brief Write the 3D microstructure to a cfg file that can be read by AtomEye
  program.

  @param timeString is the current time resolved into y,d,h,m
  */
  void writeLatticeCFG(const std::string timeString);

  /**
  @brief Write to a file of a 3D sub-microStructure of the current microStructure.

  The sub-microStructure output file will indicate the microPhaseId at each site.

  @param newZdim is the new z dimension (height starting from  z = 0) of
  the sub-microStructure; x and y dimensions don't change.
  */
  void writeNewLattice(int newZdim);

  /**
  @brief Write the 3D microstructure to a file.

  The damage output file is binary, each site either being damaged or not.

  @param timeString is the current time resolved into y,d,h,m
  */
  void writeDamageLattice(const std::string timeString);

  /**
  @brief Write the 3D microstructure to a png file that can be immediately
  rendered.

  @param timeString is the current time resolved into y,d,h,m
  */
  void writeLatticePNG(const std::string timeString);

  /**
  @brief Write the 3D microstructure to a png file that can be immediately
  rendered.

  The damage output file is binary, each site either being damaged or not.

  @param timeString is the current time resolved into y,d,h,m
  */
  void writeDamageLatticePNG(const std::string timeString);

  /**
  @brief Create files of sequential slices of the microstructure in the x
  direction.

  The slices are individual PPM files of 2D (y,z) microstructure slices,
  written back to back, in the same file.  Once created, the files are each
  converted to GIFs using a system call to Imagemagick, and then the GIFs are
  converted to an animated GIF file using a system call to gifsicle.

  @note NOT USED.

  @warning This method currently depends on system calls
  @warning This method currently depends on having Imagemagick installed
  @warning This method currently depends on having gifsicle installed

  @todo Remove the dependence on system calls, Imagemagick, and gifsicle

  @param root is the root name of the png output file to create
  */
  void makeMovie();

  /**
  @brief Set the expansion strain components of a site specified by its index.

  This function changes the strain components of a site already in the
  list of expansion sites.  If the prescribed site is not already in the
  list of expansion sites, then the site will be added to that list.

  @param index is the index of the site in the `site_` array
  @param val is the vector of expansion strain components to set
  */
  void setExpansion(int index, std::vector<double> val) {
    std::map<int, std::vector<double>>::iterator p = expansion_.find(index);
    if (p != expansion_.end()) {
      p->second = val;
    } else {
      expansion_.insert(make_pair(index, val));
    }
  }

  /**
  @brief Get the expansion strain components of a site specified by its index.

  @param index is the index of the site in the `site_` array
  @return the vector of expansion strain components to set
  */
  std::vector<double> getExpansion(int index) {
    std::map<int, std::vector<double>>::iterator p = expansion_.find(index);
    if (p != expansion_.end()) {
      return p->second;
    } else {
      std::string msg = "Could not find expansion_ match to index provided";
      throw EOBException("Lattice", "getExpansion", msg, expansion_.size(),
                         index);
    }
  }

  /**
  @brief Get the expansion strain components for all strained sites in the
  lattice.

  @return the map of the strain components, keyed to the site index numbers
  */
  std::map<int, std::vector<double>> getExpansion() { return expansion_; }

  /**
  @brief Get the coordinates of local region for calculating expansion stress.

  This gets the coordinates of the center site of a box in the lattice within
  which the expansion strain is calculated in the ThermalStrain model due to
  local crystallization pressure.

  @todo Change the function name to something like getExpansionSiteCoordinates.

  @param index is the index of a site that has crystallization pressure
  @return the (x,y,z) coordinates of the site
  */
  // std::vector<int> getExpansionCoordin(int index) {
  //   std::map<int, std::vector<int>>::iterator p =
  //   expansion_coordin_.find(index); if (p != expansion_coordin_.end()) {
  //     return p->second;
  //   } else {
  //     std::string msg = "Could not find expansion_coordin_ match to index
  //     provided"; throw EOBException("Lattice", "getExpansionCoordin", msg,
  //                        expansion_coordin_.size(), index);
  //   }
  // }

  /**
  @brief Set the coordinates of local site for calculating expansion stress.

  This gets the coordinates of the center site of a box in the lattice within
  which the expansion strain is calculated in the ThermalStrain model due to
  local crystallization pressure.

  @note NOT USED (commented in Controller)

  @todo Change the function name to something like setExpansionSiteCoordinates

  @param index is the index of a site that has crystallization pressure
  @param coordin is the (x,y,z) triple of the site's coordinates
  @return the (x,y,z) coordinates of the site
  */
  // void setExpansionCoordin(int index, std::vector<int> coordin) {
  //   std::map<int, std::vector<int>>::iterator p =
  //   expansion_coordin_.find(index); if (p == expansion_coordin_.end()) {
  //     expansion_coordin_.insert(make_pair(index, coordin));
  //   }
  // }

  /**
  @brief Get the microstructure volume

  @return the microstructure volume (GEMS volume units)
  */
  double getMicrostructureVolume(void) const {
    return (chemSys_->getMicroVolume());
  }

  /**
  @brief Get the initial microstructure volume

  @return the initial microstructure volume (GEMS volume units)
  */
  double getInitialMicrostructureVolume(void) const {
    return (chemSys_->getInitMicroVolume());
  }

  /**
  @brief Get the total voxel pore volume

  @return the volume of voxel pores (GEMS volume units)
  */
  double getVoxelPoreVolume(void) const { return voxelPoreVolume_; }

  /**
  @brief Set the voxel pore volume

  @param voxelporevolume is the voxel pore volume (GEMS volume units)
  */
  void setVoxelPoreVolume(double voxelporevolume) {
    voxelPoreVolume_ = voxelporevolume;
  }

  /**
  @brief Get the total voxel pore volume fraction
  This is calculated on a total system volume basis

  @return the volume fraction of voxel pores (microstructure basis)
  */
  double getVoxelPoreVolumeFraction(void) const {
    return voxelPoreVolumeFraction_;
  }

  /**
  @brief Get the total voxel-scale saturated pore volume fraction
  This is calculated on a total microstructure volume basis

  @return the volume fraction of voxel-scale pores (microstructure basis)
  */
  double getVoxelPoreVolumeFractionSaturated(void) const {
    return voxelPoreVolumeFractionSaturated_;
  }

  /**
  @brief Get the total voxel-scale saturated pore volume fraction
  This is calculated on a total microstructure volume basis

  @return the volume fraction of voxel pores (microstructure basis)
  */
  double getSubvoxelPoreVolumeFractionSaturated(void) const {
    return subvoxelPoreVolumeFractionSaturated_;
  }

  /**
  @brief Set the voxel pore volume fraction
  This is calculated on a total system volume basis

  @param voxelPoreVolumeFraction is the voxel pore volume
  fraction (microstructure basis)
  */
  void setVoxelPoreVolumeFraction(const double voxelPoreVolumeFraction) {
    voxelPoreVolumeFraction_ = voxelPoreVolumeFraction;
  }

  /**
  @brief Set the voxel-scale pore saturated volume fraction
  This is calculated on a total microstructure volume basis

  @param voxelPoreVolumeFractionSaturated is the voxel pore volume
  fraction (microstructure basis)
  */
  void setVoxelPoreVolumeFractionSaturated(
      const double voxelPoreVolumeFractionSaturated) {
    voxelPoreVolumeFractionSaturated_ = voxelPoreVolumeFractionSaturated;
  }

  /**
  @brief Set the subvoxel-scale pore saturated volume fraction
  This is calculated on a total microstructure volume basis

  @param subvoxelPoreVolumeFractionSaturated is the subvoxel pore volume
  fraction (microstructure basis)
  */
  void setSubvoxelPoreVolumeFractionSaturated(
      const double subvoxelPoreVolumeFractionSaturated) {
    subvoxelPoreVolumeFractionSaturated_ = subvoxelPoreVolumeFractionSaturated;
  }

  /**
  @brief Get the total capillary pore volume

  @return the volume of capillary pores (GEMS volume units)
  */
  // double getCapillaryPoreVolume(void) const { return capillaryPoreVolume_; }

  /**
  @brief Set the capillary pore volume

  @param capillaryporevolume is the capillary pore volume (GEMS volume units)
  */
  // void setCapillaryPoreVolume(double capillaryporevolume) {
  //   capillaryPoreVolume_ = capillaryporevolume;
  // }

  /**
  @brief Get the total capillary pore volume fraction
  This is calculated on a total system volume basis

  @return the volume fraction of capillary pores (microstructure basis)
  */
  // double getCapillaryPoreVolumeFraction(void) const {
  //  return capillaryPoreVolumeFraction_;
  // }

  /**
  @brief Set the capillary pore volume fraction
  This is calculated on a total system volume basis

  @param capillaryPoreVolumeFraction is the capillary pore volume
  fraction (microstructure basis)
  */
  // void setCapillaryPoreVolumeFraction(const double capillaryPoreVolumeFraction) {
  //   capillaryPoreVolumeFraction_ = capillaryPoreVolumeFraction;
  // }

  /**
  @brief Get the total subvoxel pore volume

  @return the volume of subvoxel pores (GEMS volume units)
  */
  double getSubvoxelPoreVolume(void) const { return subvoxelPoreVolume_; }

  /**
  @brief Set the subvoxel pore volume

  @param subvoxelporevolume is the subvoxel pore volume (GEMS volume units)
  */
  void setSubvoxelPoreVolume(const double subvoxelporevolume) {
    subvoxelPoreVolume_ = subvoxelporevolume;
  }

  /**
  @brief Set the subvoxel pore volume

  @param subvoxelporevolume is the subvoxel pore volume (GEMS volume units)
  */
  void setNonSolidVolume(const double nonsolidvolume) {
    nonSolidVolume_ = nonsolidvolume;
  }

  /**
  @brief Get the voxel water volume

  @param vol is the volume of each microstructure phase
  */
  // void calcVoxelWaterVolume(std::vector<double> &vol);

  /**
  @brief Get the voxel water volume

  @return the voxel water volume
  */
  double getVoxelWaterVolume(void) const { return voxelWaterVolume_; }

  /**
  @brief Set the voxel water volume

  @param voxelwatervolume is the voxel water volume (GEMS volume units)
  */
  void setVoxelWaterVolume(const double voxelwatervolume) {
    voxelWaterVolume_ = voxelwatervolume;
  }

  /**
  @brief Get the voxel void volume

  @param vol is the volume of each microstructure phase
  @param calc is true only if calculating instead of just returning
  */
  // void calcVoxelVoidVolume(void);

  /**
  @brief Get the voxel void volume

  @return the voxel void volume
  */
  double getVoxelVoidVolume(void) const { return voxelWaterVolume_; }

  /**
  @brief Set the voxel void volume

  @param voxelvoidvolume is the voxel void volume (GEMS volume units)
  */
  void setVoxelVoidVolume(const double voxelvoidvolume) {
    voxelVoidVolume_ = voxelvoidvolume;
  }

  /**
  @brief Get the total subvoxel pore volume fraction
  This is calculated on a total system volume basis

  @return the volume fraction of subvoxel pores (microstructure basis)
  */
  double getSubvoxelPoreVolumeFraction(void) const {
    return subvoxelPoreVolumeFraction_;
  }

  /**
  @brief Set the subvoxel pore volume fraction
  This is calculated on a total system volume basis

  @param subvoxelporevolumefraction is the subvoxel pore volume
  fraction (microstructure basis)
  */
  void setSubvoxelPoreVolumeFraction(const double subvoxelporevolumefraction) {
    subvoxelPoreVolumeFraction_ = subvoxelporevolumefraction;
  }

  /**
  @brief Get the volume fraction saturated  of the idx element of the pore
  volume distribution
  @param idx is the index to get
  @return the volume fraction saturated of that element in the pore size
  distribution
  */
  // double getMasterPoreVolumeVolfrac(const int idx) {
  //   try {
  //     if (idx >= static_cast<int>(masterPoreSizeDist_.size())) {
  //       throw EOBException("Lattice", "getMasterPoreVolumeVolfrac",
  //                          "masterPoreSizeDist_", masterPoreSizeDist_.size(),
  //                          static_cast<int>(idx));
  //     }
  //   } catch (EOBException ex) {
  //     ex.printException();
  //     exit(1);
  //   }
  //   return (masterPoreSizeDist_[idx].volfrac);
  // }

  /**
  @brief Get the largest diameter of pores containing electrolyte
  @return the diameter of the largest pore containing electrolyte
  */
  double getLargestSaturatedPore(void) {
    double capsize = 1000.0; // nm of voxel pores
    int size = masterPoreSizeDist_.size();
    for (int i = 0; i < size; i++) {
      if (masterPoreSizeDist_[i].volfrac < 1.0) {
        return (masterPoreSizeDist_[i].diam);
      }
    }
    return (capsize);
  }

  /**
  @brief Get the number of sites of water that must be added after a time step.

  @note Currently only used in sulfate attack simulations.

  @return the amount of water that must be added [site units]
  */
  double getWaterChange(void) const { return waterChange_; }

  /**
  @brief Set the number of sites of water that must be added after a time step.

  @note NOT USED.

  @param the number of sites of water that must be added [site units]
  */
  void setWaterChange(double waterchangeval) { waterChange_ = waterchangeval; }

  /**
  @brief Increment the number of sites of water that must be added after a time
  step.

  @param the extra number of sites of water that must be added [site units]
  */
  void dWaterChange(double dwaterchangeval) { waterChange_ += dwaterchangeval; }

  /**
  @brief Set a pointer to the AppliedStrain object for the simulation.

  @param elas is a pointer to the AppliedStrain object for the simulation
  */
  void setFEsolver(AppliedStrain *AppliedStrainSolver) {
    FEsolver_ = AppliedStrainSolver;
  }

  /**
  @brief Write a contiguous subvolume of the lattice to a file.

  @param fileName is the file name to which the subvolume will be written
  @param centerste is a pointer to the central site within the subvolume
  @param size is the extent of the subvolume in each direction away from the
  center site
  @return a list of the site indices belonging to the subvolume that was written
  */
  std::vector<int> writeSubVolume(std::string fileName, Site *centerste, int size);

  /**
  @brief Assign isotropic expansion strain at a set of prescribed sites.

  This function changes the strain components of a site already in the
  list of expansion sites.  If the prescribed site is not already in the
  list of expansion sites, then the site will be added to that list.

  @todo Consider changing the name of this method to applyExpansion

  @param alnb is the collection of site indices to which strain will be assigned
  @param exp is the isotropic expansion strain to set
  */
  void applyExpansion(std::vector<int> alnb, double exp);

  /**
  @brief Estimate the surface areas of all solid phases
  with the aqueous solution, in units of m2 per 100 g of
  total solids

  @param solidMass is the mass of all solids in g
  */
  void calcSurfaceAreas(void) {
    double scaledMass;
    surfaceArea_.resize(numMicroPhases_, 0.0);
    specificSurfaceArea_.resize(numMicroPhases_, 0.0);
    for (int i = 0; i < numMicroPhases_; ++i) {
      calcSurfaceArea(i);

      // Calculate specific surface area of this phase by dividing
      // this surface area by the phase mass (g per 100 g of all solid)
      // Units of specific surface are will be m2 per kg of this phase,
      // to make it consistent with legacy Parrot-Killoh model which
      // uses traditional Blaine fineness units
      scaledMass = chemSys_->getMicroPhaseMass(i);
      if (scaledMass > 0.0) {
        // The factor of 1000.0 converts units from m2/g to m2/kg
        specificSurfaceArea_[i] = 1000.0 * surfaceArea_[i] / scaledMass;
      } else {
        specificSurfaceArea_[i] = 0.0;
      }

    }
  }

  /**
  @brief Estimate the surface area of a phase with the aqueous
  solution, in units of m2 per 100 g of total solids

  @param phaseid is the id of the microstructure phase
  */
  void calcSurfaceArea(int phaseid);

  /**
  @brief Return the current surface area of a phase with the aqueous
  solution.

  @param phaseid is the id of the microstructure phase
  @return the estimated surface area [m2 per 100 g of all solid]
  */
  double getSurfaceArea(int phaseid) {
    // if (phaseid > -1 && phaseid < surfaceArea_.size()) {
    return surfaceArea_[phaseid];
    // }
    // return 0.0;
  }

  /**
  @brief Return the current surface area vector, surfaceArea_, of all the phases
  with the aqueous solution.

  @return the estimated surface area [m2 per 100 g of all solid]
  */
  std::vector<double> getSurfaceArea(void) { return surfaceArea_; }

  /**
  @brief Reset the current surface area vector, surfaceArea_, of all the phases
  with the aqueous solution.

  @param vect is the vector containing the surface area of each microPhase in the
  microStructure
  */
  void resetSurfaceArea(std::vector<double> vect) { surfaceArea_ = vect; }

  /**
  @brief Return the current specific surface area of a phase with the aqueous
  solution.

  @param phaseid is the id of the microstructure phase
  @return the estimated specific surface area [m2 per g of this phase]
  */
  // double getSpecificSurfaceArea(int phaseid) {
  //   if (phaseid > -1 && phaseid < specificSurfaceArea_.size()) {
  //     return specificSurfaceArea_[phaseid];
  //   }
  //   return 0.0;
  // }

  /**
  @brief Return the combined specific surface area of cementitious
  components

  @return the estimated specific surface area [m2 per kg of all cement]
  */
  double getCementSpecificSurfaceArea(void) {
    double cemmass = 0.0;
    // double allsurf = 0.0;
    double cemsurf = 0.0;
    // double allsolidmass = 0.0;
    double thismass, thisssa;
    for (int i = 0; i < numMicroPhases_; ++i) {
      if (i != VOIDID && i != ELECTROLYTEID) {
        // allsolidmass += chemSys_->getMicroPhaseMass(i);
        // allsurf += getSurfaceArea(i);
        if (chemSys_->isCementComponent(i)) {
          thisssa = getSurfaceArea(i); // m2 component /(100g solid)
          thismass =
              chemSys_->getMicroPhaseMass(i); // g component/(100 g solid)
          cemsurf += thisssa;
          cemmass += thismass; // g all cement / (100g solid)
        }
      }
    }
    // At this point ssa is the total surface area of cement per 100 g of all
    // solids, and cemmass is the total g of cement per 100 g of all solids.
    // So ssa/cemmass has units of m2 per g of cement
    // Multiply that by 1000.0 to get units of m2/(kg of cement)
    // if (verbose_) {
    //    cout << "URANIUM all solid mass = " << allsolidmass << " g / (100 g
    //    solid)"
    //          << endl;
    //    cout << "URANIUM all surface = " << allsurf << " m2 / (100 g solid)"
    //            << endl;
    //    cout << "URANIUM cement mass = " << cemmass << " g / (100 g solid)" <<
    //    endl; cout << "URANIUM cement surface = " << cemsurf << " m2 / (100 g
    //    solid)"
    //            << endl;
    //    cout.flush();
    // }
    if (cemmass > 0.0) {
      cemsurf *= (1000.0 / cemmass);
    } else {
      cemsurf = 0.0;
    }
    return cemsurf;
  }

  /**
  @brief Get the sorted distribution of domain sizes

  @param phaseid is the id of the phase to query
  @param numsites is the maximum number of sites to store and sort
  @param maxisze is the maxmimum linear size of interest
  @param sortorder is 0 if sorting in descending order, nonzero otherwise
  @return an STL list of the site ids according to the distribution
  */
  std::vector<int> findDomainSizeDistribution(int phaseid, const int numsites,
                                              int maxsize, int sortorder);

  /**
  @brief Estimate the <i>linear size</i> of a domain

  @param siteid is the id of the microstructure phase
  @param maxsize is the maxmimum linear size of interest
  @return the edge length of the maximum cube that contains the same phase
  */
  int findDomainSize(int siteid, int maxsize);

  /**
  @brief Set the verbose flag

  @param isverbose is true if verbose output should be produced
  */
  void setVerbose(const bool isverbose) {
    verbose_ = isverbose;
    return;
  }

  /**
  @brief Get the verbose flag

  @return the verbose flag
  */
  bool getVerbose() const { return verbose_; }

  /**
  @brief Set the warning flag

  @param iswarning is true if warning output should be produced
  */
  void setWarning(const bool iswarning) {
    warning_ = iswarning;
    return;
  }

  /**
  @brief Get the warning flag

  @return the warning flag
  */
  bool getWarning() const { return warning_; }

  /**
  @brief Get the number of sites occupied by each microPhase in the
  microStructure (count_ vector)

  @return the count_ vector
  */
  std::vector<int> getCount(void) { return count_; }

  /**
  @brief Get the number of sites occupied by a given microPhase,
  having a phId microPhaseId, in the microStructure

  @param phId is the microPhaseId
  @return the number of sites occupied by microPhaseId = phId
  */
  int getCount(int phId) { return count_[phId]; }

  /**
  @brief reset the count_ vector, i.e. the number of sites occupied by each
  microPhase in the microStructure

  @param vect is a vector containing the number of sites for each microPhase in
  the microStructure
  */
  void setCount(std::vector<int> vect) { count_ = vect; }

  /**
  @brief Get the dimension of the interface_ vector

  @return the dimension of the interface_ vector
  */
  int getInterfaceSize(void) { return interface_.size(); }

  /**
  @brief Get both, the dissolution and growth interfaces of a given microPhase,
  having the microPhaseId = phId

  @param phId is the microPhaseId
  @return the interfaces (dissolution and growth) of the microPhase with
  microPhaseId = phId
  */
  Interface getInterface(int phId) { return interface_[phId]; }

  /**
  @brief Set the growth interface of the microPhase having the microPhaseId = phId

  @param phId is the microPhaseId
  @param vect is a vector containing all Isite objects belonging to the growth
  interface oh the microPhase with microPhaseId = phId
  */
  void setGrowthSites(int phId, std::vector<Isite> vect) {
    interface_[phId].setGrowthSites(vect);
  }

  /**
  @brief Set the dissolution interface of the microPhase having the microPhaseId = phId

  @param phId is the microPhaseId
  @param vect is a vector containing all Isite objects belonging to the dissolution
  interface oh the microPhase with microPhaseId = phId
  */
  void setDissolutionSites(int phId, std::vector<Isite> vect) {
    interface_[phId].setDissolutionSites(vect);
  }

  // void setInterfaceMicroPhaseId(int i, int mPhId) {
  //   interface_[i].setMicroPhaseId(mPhId);
  // }

  /**
  @brief Get the initial solid mass of the system

  @return the initial solid mass of the system
  */
  double getInitSolidMass(void) { return initSolidMass_; }

  /**
  @brief find the voxels voxels without contact with electrolyte

  @note NOT USED.
  */
  void findIsolatedClusters(void);

  /**
  @brief create a chemical element database (chemical symbol, atomic number and
  atomic mass)

  */
  void populateElementData(void);

  /**
  @brief Get the chemical symbol of a database element.

  @param ind is the index of an chemical element in the chemical element database
  @return the chemical symbol of the chemical element database having
  the index = ind
  */
  std::string getElemSymb(int ind) { return cfgElem_[ind].symb; }

  /**
  @brief call the random number generator; count the number of calls in such a way
  to be able to restore a given state of the random number generator, and update
  the value of lastRNG_.

  @return the generated random number corresponding to this call.
  */
  double callRNG(void) {
    numRNGcall_0_++;
    if (numRNGcall_0_ == LONG_MAX) {
      numRNGcallLONGMAX_++;
      numRNGcall_0_ = 0;
    }
    lastRNG_ = rg_->Ran3();
    return lastRNG_;
  }

  /**
  @brief get the first number used to keep the RNG calls track
  (see callRNG method)

  @return the first number used to keep the RNG calls track
  */
  long int getNumRNGcall_0(void) { return numRNGcall_0_; }

  /**
  @brief get the second number used to keep the RNG calls track
  (see callRNG method)

  @return the second number used to keep the RNG calls track
  */
  long int getNumRNGcallLONGMAX(void) { return numRNGcallLONGMAX_; }

  /**
  @brief get the last generated random number (see callRNG method)

  @return the last generated random number
  */
  double getLastRNG(void) { return lastRNG_; }

  // void setRNGseed(int seed) { rg_->setSeed(seed); }

  // int getRNGseed(void) { return latticeRNGseed_; }

  /**
  @brief reset the random number generator state to the state corresponding to
  the val_0 and valLONGMAX numbers; valRNG is used to check the correctness of
  this reset.

  @param val_0 is the first number used to keep the RNG calls
  @param valLONGMAX is the second number used to keep the RNG calls
  @param valRNG is the last generated number before this reset
  */
  void resetRNG(long int val_0, long int valLONGMAX, double valRNG) {
    // latticeRNGseed_ = seed;
    rg_->setSeed(latticeRNGseed_);
    numRNGcall_0_ = val_0;
    numRNGcallLONGMAX_ = valLONGMAX;
    // long int count_0 = 0, count_1 = 0;
    long int j0, j1, j11;
    double lastRNGreset = 1.e-16;
    for (j1 = 1; j1 <= numRNGcallLONGMAX_; j1++) {
      for (j11 = 1; j11 <= LONG_MAX; j11++) {
        lastRNGreset = rg_->Ran3();
      }
    }
    for (j0 = 1; j0 <= val_0; j0++) {
      lastRNGreset = rg_->Ran3();
    }
    lastRNG_ = lastRNGreset;

    // cout << endl
    //      << "  Lattice::resetRNG cyc/whileCount/latticeRNGseed_: " << cyc
    //      << " / " << whileCount << " / " << latticeRNGseed_ << endl;
    // cout << "  Lattice::resetRNG "
    //         "numRNGcall_0_/numRNGcallLONGMAX_/lastRNGreset/valRNG: "
    //      << numRNGcall_0_ << " / " << numRNGcallLONGMAX_ << " / "
    //      << lastRNGreset << " / " << valRNG << endl;

    if (abs(lastRNGreset - valRNG) > 1.e-16) {
      std::cout << std::endl
                << "Lattice::resetRNG FAILED => exit" << std::endl;
      exit(0);
    }
  }

  // void increaseLatticeVolume(void);

  void checkSite(int stId);

  /**
  @brief Get the growth interface dimensions of all microPhases in the microStructure

  @return vector containing the dimensions of all the growth interfaces for all
  microPhases in the microStructure (each interface corresponding to a given microPhase)
  */
  std::vector<int> getGrowthInterfaceSize(void) { return growthInterfaceSize_; }

  // int getGrowthInterfaceSize(const int phId) {
  //   return growthInterfaceSize_[phId];
  // }

  /**
  @brief Set the growth interface dimension of each microPhase in the microStructure

  @param vect is a vector containing the dimensions of all the growth interfaces,
  each interface corresponding to a given microPhase in the microStructure
  */
  void setGrowthInterfaceSize(std::vector<int> vect) { growthInterfaceSize_ = vect; }

  /**
  @brief Get the dissolution interface dimensions of all microPhases in the microStructure

  @return vector containing the dimensions of all the dissolution interfaces for all
  microPhases in the microStructure (each interface corresponding to a given microPhase)
  */
  std::vector<int> getDissolutionInterfaceSize(void) {
    return dissolutionInterfaceSize_;
  }

  /**
  @brief Get the dissolution interface dimension of a given microPhase in the microStructure

  @param phId is the microPhaseId
  @return the dimension of the dissolution interface for a microPhase having
  the microPhaseId = phId
  */
  int getDissolutionInterfaceSize(int phId) {
    return dissolutionInterfaceSize_[phId];
  }

  /**
  @brief Set the dissolution interface dimension of each microPhase in the microStructure

  @param vect is a vector containing the dimensions of all the dissolution interfaces,
  each interface corresponding to a given microPhase in the microStructure
  */
  void setDissolutionInterfaceSize(std::vector<int> vect) {
    dissolutionInterfaceSize_ = vect;
  }

  /**
  @brief implements a model of etttringite formation (growth) as a result
  of the sulfate attack. Ettringite formation takes place in two steps: a
  conversion of an Al-bearing phase voxel to an ettringite voxel (calling
  the transformSolSol method) and the conversion, if possible, of one ore
  more electrolyte nearest neighbor voxels to ettringite (calling
  transformLiqSol method). It also calculates the volume of free space
  adjacent to a converted site to determine whether crystallization pressure
  will arise. If so, the method calculates the crystallization pressure and
  crystallization stress-free strain. It then applies the expansion strain
  so that the new stress field can be calculated by the ThermalStrain FE model
  object.

  @todo For the moment the model describes only monosulfate to ettringite
  conversion - generalize this for any phase transformation.

  @param ettrid is the microstructure phase id of ettringite
  @param netsitesEttrid is the number of ettringite sites to grow
  @param dissPhaseIDVect is a vector containing all phase Ids that can transform
  in ettringite
  @param numSiteDissVect is a vector containing the number of voxels to dissolve
  (not to transform in ettringite!) for each microphase ID in dissPhaseIDVect
  @param dissPhNameVect is a vector containing the name of each microphase in
  dissPhaseIDVect
  @param volumeRatio is a vector containing the ratios between the molar volume
  of each microphase in dissPhNameVect and the molar volume of ettringite
  @param numtoadd_D is the number of voxels converted from an Al-bearing phase
  to ettringite (first step of the model!)
  @param totalTRC is the total call number of the changeMicrostructure method
  @return a vector having a dimension N equal to the number of phases that can
  transform into ettringite + 1: first N-1 positions contain the number of voxels
  of each phase that must be dissolved further by the "normal" dissolution (using
  dissolvePhase method), while the last position contains the number of ettringite
  voxels added by the current method. If the the number of  added ettringite voxels
  is smaller than netsitesEttrid, this difference will be added calling the
  growPhase method ("normal" growth).
  */
  std::vector<int>
       transformPhase(int ettrid, int netsitesEttrid,
                      std::vector<int> dissPhaseIDVect,
                      std::vector<int> numSiteDissVect,
                      std::vector<std::string> dissPhNameVect,
                      std::vector<double> volumeRatio,
                      int &numadded_D, int totalTRC);

  /**
  @brief implements the first step of the model describing the ettringite growth
  under sulfate attack conditions (see transformPhase method). It converts a solid
  microphase occupying a voxel to another solid microphase.

  @param ste is apointer to the Site object being occupied by a solid microphase
  @param oldPhId is the solid microphase Id occupying the site ste
  @param newPhId is the new solid microphase Id that must occupy the site ste
  @param totalTRC is the total call number of the changeMicrostructure method
  */
  void transformSolSol(Site *ste, int oldPhId, int newPhId, int totalTRC); // sol to sol

  /**
  @brief implements the second step of the model describing the ettringite growth
  under sulfate attack conditions (see transformPhase method). It converts the
  electrolyte Id occupying a voxel to a solid microphase Id.
  @param ste is apointer to the Site object being occupied by electrolyte
  @param growPhID is the new solid microphase Id that must occupy the site ste
  @param totalTRC is the total call number of the changeMicrostructure method
  @return a vector containing the microPhaseIds of the sites that, modifying the ste site
  occupancy, cannot belong to the dissolution interfaces of the microPhases
  occupying these sites
  */
  std::vector<int> transformLiqSol(Site *ste, int growPhID, int totalTRC);  // liq to sol

  /**
  @brief creates the vectors growingVectSA_, shrinking_ and volratios_ vectors:
  growingVectSA_ - contains all microPhaseIds growing due to SA attack (AFt)
  shrinking_     - contains for each microPhaseId in growingVectSA_, all the microDhaseIds
                   that can transform into this one (Monosulfate)
  volratios_     - for each microPhaseId in growingVectSA_ and all corresponding
                   microPhaseIds in shrinking_, contains the molar volume ratios
                   of the corresponding microPhases:
        volratios_[i, j] = molarVolume(growingVectSA_[i])/molarVolume(shrinking_[i, j])

  @todo For the moment the model describes only monosulfate to ettringite
  conversion - generalize this for any phase transformation.
  */
  void createGrowingVectSA(void);

  std::vector<int> getAllSitesPhId(void) { // check!
    std::vector<int> allPhId(numSites_, 0);
    for (int i = 0; i < numSites_; i++) {
      allPhId[i] = site_[i].getMicroPhaseId();
    }
    return allPhId;
  }

  void addSeedCSHQ(bool seedMassCEM, bool seedMassC3S,
                   bool seedMassC2S, double massFraction);

  /**
  @brief Calculate the surface area of one voxel's face per 100g of
  the initial solid mass of the system (in units of m2 per 100 g of 
  total solids)
  */
  void calcOneFaceAreaPerHundredGramSolid(void);

}; // End of Lattice class

#endif // SRC_THAMESLIB_LATTICE_H_

///
/// The functions below are used to aid in comparison of one site to another,
/// by which means lists of the sites can be sorted.
///

#ifndef CMPFUNCS
#define CMPFUNCS

/**
@brief Compare two sites, returning true is the first site is "less than" the
second.

The comparison is made on the basis of what THAMES loosely calls the
<i>weighted mean curvature</i>, (wmc).  A site with high wmc is a site where
dissolution of a phase is likely to occur, and growth of another phase is
unlikely to occur. Conversely, a site with a low wmc is a site where growth of
a phase is likely to occur but dissolution of a phase is unlikely to occur.

@param s1 is a pointer to the first site in the comparison
@param s2 is a pointer to the second site in the comparison
@return true if the first site has lower wmc than the second, false otherwise
*/
bool cmp(const Site *s1, const Site *s2);

/**
@brief Sort two sites based on their affinity for a given phase.

The comparison is made on the basis of what THAMES loosely calls the
<i>affinity</i>.  A site with high affinity is a site where growth of a phase
is more likely to occur because of an affinity between it and the interface.

@param s1 is the first site in the comparison
@param s2 is the second site in the comparison
@return true if the first site has <i>greater</i> affinity than the second,
false otherwise
*/
bool affinitySort(const Isite s1, const Isite s2);

#endif
