/**
@file Site.h
@brief Declaration of the Site class.
*/

#ifndef SRC_THAMESLIB_SITE_H_
#define SRC_THAMESLIB_SITE_H_

#include "ChemicalSystem.h"
#include "Exceptions.h"
#include "global.h"

/**
@class Site
@brief Handle behavior of individual lattice sites.

The Site class manages changes in the phase id of a lattice site,
and whether or not it is a dissolution site or a growth site.  The site
can store a flag that determines whether or not it is damaged by
deterioration, and stores a list of its neighbor sites by their index
number.

*/
class Site {

protected:
  int x_ = 0;            /**< x-coordinate in mesh coordinate frame */
  int y_ = 0;            /**< y-coordinate in mesh coordinate frame */
  int z_ = 0;            /**< y-coordinate in mesh coordinate frame */
  int id_ = 0;           /**< Unique id in the 1D array of all sites */
  int microPhaseId_ = 0; /**< The microstructure phase assignment */
  ChemicalSystem
      *chemSys_; /**< Pointer to simulation's ChemicalSystem object */
  std::vector<int> growth_; /**< Vector of phases that can grow at this site */
  double stressFreeVolume_; /**< Stress-free volume of the site */
  double trueVolume_;       /**< Actual volume of site, accounting for stress */
  bool damage_ = false;     /**< True if site is damaged, false otherwise */
  std::vector<Site *>
      nb_; /**< List of site ids that are neighbors to this site */
  std::vector<int>
      nbSA_; /**< List of site ids that are neighbors to this site */

  /**
  @brief Ranking of potential for dissolution if the site is an interface site.

  The wmc (acronym for weighted mean curvature) is a quantitative ranking of the
  potential for a site's phase to dissolve based on its local environment.  It
  is kind of like mean curvature, but potentially can be weighted to account for
  crystalline anisotropy.

  THAMES is not currently able to quantitatively calculate wmc as defined in
  the metallurgy literature:

      - Taylor, J.E., Cahn, J.W., Handwerker, C.E., Geometric models of crystal
          growth, Acta metall. mater. 40 [7] (1992) 1443--1474.
      - Taylor, J.E., II--Mean curvature and weighted mean curvature,
          Acta metall. mater. 40 [7] (1992) 1475--1485.

  Instead, THAMES uses the digital template method to calculate a quantity that
  is roughly linearly proportional to mean curvature:

      - Bullard, J.W., Garboczi, E.J., Carter, W.C., Fuller, E.J., Numerical
          methods for computing interfacial mean curvature, Comput. Mater. Sci.
          4 (1995) 103--116.

  This provides a ranking of dissolution potential only.
  */
  int wmc_;  /**< total porosity ("surface curvature") at this site */
  int wmc0_; /**< this site internal porosity
                  (its own contribution at wmc_ value) */

  double expstrain_; /**< Assigned expansion strain by phase
                          constrained transformation or an applied load */

  bool verbose_; /**< Flag to determine verbose output */

  std::vector<int>
      inGrowInterfacePos_; /**< vector of the site position in each
                                growth interface (-1 if the site doesn't
                                belong to a growth interface) */
  int inDissInterfacePos_; /**< site position in the corresponding
                                dissolution interface (-1 if the site
                                doesn't belong to the dissolution
                                interface) */

  std::vector<int>
      inGrowthVectorPos_;      /**< vector of the site position in the
                                    growth vector (-1 if the site doesn't
                                    belong to the growth vector) */
  int inDissolutionVectorPos_; /**< vector of the site position in the
                                    dissolution vector (-1 if the site doesn't
                                    belong to the growth vector) */

  int visit_; /**< flag used to avoid acting twice or more
                   on the site (or its neighborhood) during
                   a current action. Set to the current action's ID
                   when the site is first touched in an iteration;
                   compared against the same ID on subsequent visits
                   to short-circuit re-processing. Caller responsible
                   for clearing (or bumping the action ID) between
                   distinct iterations. */

public:
  /**
  @brief Default constructor.

  @note NOT USED.
  */
  Site();

  /**
  @brief Overloaded constructor.

  This constructor takes arguments for all the member variables so that they can
  be assigned during construction.

  The (x,y,z) coordinates and the (x,y,z) dimensions of the lattice are used
  to calculate and assign a unique index number for the site in a 1D array
  for the lattice.

  @param xp is the x-coordinate of the site in the mesh coordinate frame
  @param yp is the y-coordinate of the site in the mesh coordinate frame
  @param zp is the z-coordinate of the site in the mesh coordinate frame
  @param xs is the number of sites in the x dimension of the lattice
  @param ys is the number of sites in the y dimension of the lattice
  @param zs is the number of sites in the z dimension of the lattice
  @param neigh is the number of adjacent sites to be considered as neighbors
  @param csys is a pointer to the simulation's ChemicalSystem object
  @param verbose is true if verbose output should be produced
  */
  Site(int xp, int yp, int zp, int xs, int ys, int zs, int neigh,
       ChemicalSystem *csys, const bool verbose = false);

  /**
  @brief Set the visit flag to the sv value

  @param sv is 1 if the site has been already checked and 0 if not
  */
  void setVisit(int sv) { visit_ = sv; }

  /**
  @brief Get the visit flag

  @return the visit flag
  */
  int getVisit(void) { return visit_; }

  /**
  @brief Set the site position (pos) on a growth interface corresponding to a
  given microPhase (having microPhaseId = phId)

  @param phId is the microPhaseId
  @param pos is the site position on the growth interface of the given
  microPhase, having microPhaseId = phId
  */
  void setInGrowInterfacePos(int phId, int pos) {
    inGrowInterfacePos_[phId] = pos;
  }

  /**
  @brief Get the site position on a growth interface corresponding to a given
  microPhase (having microPhaseId = phId)

  @param phId is the microPhaseId
  @return the site position on the growth interface corresponding to the
  microPhase having microPhaseId = phId
  */
  int getInGrowInterfacePos(int phId) { return inGrowInterfacePos_[phId]; }

  /**
  @brief Set the site positions on the growth interfaces at some previous values
  given by the vector vect

  @param vect is the vector containing some previous site positions on the
  corresponding growth interfaces (-1 if the site doesn't belong to a given
  microPhase growth interface)
  */
  void setInGrowInterfacePosVector(std::vector<int> vect) {
    inGrowInterfacePos_ = vect;
  }

  /**
  @brief Get the vector containing the site positions on the growth interfaces
  correponding to a given microStructure configuration (a given time)

  @return the vector containing the site positions on the growth interfaces
  correponding to a given microStructure configuration (a given time)
  */
  std::vector<int> getInGrowInterfacePosVector() { return inGrowInterfacePos_; }

  /**
  @brief Set the site position (pos) on the dissolution interface corresponding
  to the microPhase occupying this site

  @param pos is the site position on the dissolution interface of the microPhase
  occupying this site
  */
  void setInDissInterfacePos(int pos) { inDissInterfacePos_ = pos; }

  /**
  @brief Get the site position (inDissInterfacePos_) on the dissolution
  interface corresponding to the microPhase occupying this site

  @return the site position (inDissInterfacePos_) on the dissolution interface
  corresponding to the microPhase occupying this site
  */
  int getInDissInterfacePos(void) { return inDissInterfacePos_; }

  /**
  @brief Set the site position (pos) corresponding to a given microPhase
  (having microPhaseId = phId) on the current growth vector (-1 if the site
  doesn't belong, for the given phId, to the current growth vector)

  @param phId is the microPhaseId
  @param pos is the site position corresponding to the given microPhase, having
  microPhaseId = phId, on the current growth vector (-1 if the site doesn't
  belong, for the given phId, to the current growth vector)
  */
  void setInGrowthVectorPos(int phId, int pos) {
    inGrowthVectorPos_[phId] = pos;
  }

  /**
  @brief Get the site position corresponding to a given microPhase
  (having microPhaseId = phId) on the current growth vector (-1 if the site
  doesn't belong, for the given phId, to the current growth vector)

  @param phId is the microPhaseId
  @return the site position corresponding to the given microPhase, having
  microPhaseId = phId, on the current growth vector (-1 if the site doesn't
  belong, for the given phId, to the current growth vector)
  */
  int getInGrowthVectorPos(int phId) { return inGrowthVectorPos_[phId]; }

  /**
  @brief Set the site position (pos) on the dissolution vector

  @param pos is the site position on the dissolution vector
  */
  void setInDissolutionVectorPos(int pos) { inDissolutionVectorPos_ = pos; }

  /**
  @brief Get the site position on the dissolution vector

  @return the site position on the dissolution vector
  */
  int getInDissolutionVectorPos(void) { return inDissolutionVectorPos_; }

  /**
  @brief Get a pointer to a given site in the site's neighborhood.

  The neighbors are stored in a 1D vector of site index values.  The
  construction of the neighbor table happens in the Lattice class.

  @param pos is the index of the neighbor in the neighbor table
  @return a pointer to the neighboring site
  */
  Site *nb(const int pos) const {
    if (pos >= static_cast<int>(nb_.size())) {
      throw EOBException("Site", "nb", "nb_", nb_.size(), pos);
    }
    return nb_[pos];
  }

  /**
  @brief Get the size of the neighbor table (number of neighbors).

  @param dist is the distance (site dimensions) away from the site to consider
  @return the number of neighbors within this distance
  */
  int nbSize(int dist = 3) const {
    switch (dist) {
    case 0:
      return 1;
      break;
    case 1:
      return NUM_NEAREST_NEIGHBORS;
      break;
    case 2:
      return (NN_NNN);
      break;
    default:
      return nb_.size();
      break;
    }
    return nb_.size();
  }

  /**
  @brief Set a neighbor position with a particular site pointer.

  This method is not used to push neighbors onto the neighbor vector.  The
  creation of the neighbor table happens elsewhere in the Lattice class.

  @param i is the index in the allocated neighbor table.
  @param neigh is a pointer to the site that is to be assigned to index i of the
  neighbor table
  */
  void setNb(int i, Site *neigh) {
    if (i >= static_cast<int>(nb_.size()))
      throw EOBException("Site", "setNb", "nb_", nb_.size(), i);
    nb_[i] = neigh;
    return;
  }

  /**
  @brief Get the neighbor vector of this site.

  @return the neighbor vector (pointers) of this site
  */
  std::vector<Site *> getNb() { return nb_; }

  /**
  @brief Set a neighbor position with a particular site.

  This method is not used to push neighbors onto the neighbor vector.  The
  creation of the neighbor table happens elsewhere in the Lattice class.

  @param ind is the index in the allocated neighbor table.
  */
  void setNbSA(int ind) { nbSA_.push_back(ind); }

  /**
  @brief Get the neighbor vector of this site.

  @return the neighbor vector of this site
  */
  std::vector<int> getNbSA() { return nbSA_; }

  /**
  @brief Get the index number of the site (position in the 1D Lattice vector).

  @return the index number of the site
  */
  int getId() const { return id_; }

  /**
  @brief Get the microstructure phase id number assigned to the site.

  @return the microstructure phase id number
  */
  int getMicroPhaseId() const { return microPhaseId_; }

  /**
  @brief Get the microstructure phase name assigned to the site.

  @return the microstructure phase name
  */
  // string getMicroPhaseName() const {
  //   return chemSys_->getMicroPhaseName(microPhaseId_);
  // }

  /**
  @brief Set the microstructure phase id number assigned to the site.

  @param pid is the microstructure phase id number to assign
  */
  void setMicroPhaseId(const int pid) {
    microPhaseId_ = pid;
    stressFreeVolume_ = 1.0;
    return;
  }

  /**
  @brief Get the x-coordinate of the site in the mesh coordinate frame.

  @return the x-coordinate of the site in the mesh coordinate frame
  */
  int getX() const { return x_; }

  /**
  @brief Get the y-coordinate of the site in the mesh coordinate frame.

  @return the y-coordinate of the site in the mesh coordinate frame
  */
  int getY() const { return y_; }

  /**
  @brief Get the z-coordinate of the site in the mesh coordinate frame.

  @return the z-coordinate of the site in the mesh coordinate frame
  */
  int getZ() const { return z_; }

  /**
  @brief Get the vector of x, y and z coordinates of this site in the
  mesh coordinate frame.

  @return the vector of x, y and z coordinates of this site in the
  mesh coordinate frame.
  */
  std::vector<int> getXYZ();

  /**
  @brief Get the "weighted mean curvature" of the site.

  The wmc (acronym for weighted mean curvature) is a quantitative ranking of the
  potential for a site's phase to dissolve based on its local environment.  It
  is kind of like mean curvature, but potentially can be weighted to account for
  crystalline anisotropy.

  THAMES is not currently able to quantitatively calculate wmc as defined in
  the metallurgy literature:

      - Taylor, J.E., Cahn, J.W., Handwerker, C.E., Geometric models of crystal
          growth, Acta metall. mater. 40 [7] (1992) 1443--1474.
      - Taylor, J.E., II--Mean curvature and weighted mean curvature,
          Acta metall. mater. 40 [7] (1992) 1475--1485.

  Instead, THAMES uses the digital template method to calculate a quantity that
  is roughly linearly proportional to mean curvature:

      - Bullard, J.W., Garboczi, E.J., Carter, W.C., Fuller, E.J., Numerical
          methods for computing interfacial mean curvature, Comput. Mater. Sci.
          4 (1995) 103--116.

  This provides a ranking of dissolution potential only.

  @return the relative potential for dissolution at this site
  */
  int getWmc(void) const { return wmc_; }

  /**
  @brief get the internal porosity, wmc0_, of this site i.e. its contribution to
  the "weighted mean curvature", wmc_, of the site

  @return the value of wmc0_ corresponding to the microPhase occupying this site
  */
  int getWmc0(void) const { return wmc0_; }

  /**
  @brief Determine if a site is occupied by a porous solid

  A porous solid is any phase that is not pure electrolyte but which
  has some internal porosity (<i>i.e.<i>, 0 < porosity < 1 )

  @return true if the phase occupying this site is a porous solid
  */
  bool isPorousSolid(void) {
    if ((wmc0_ > 0) && (microPhaseId_ > ELECTROLYTEID))
      return (true);
    return (false);
  }

  /**
  @brief Set the "weighted mean curvature" of the site.

  @param wmcval is the value of wmc_ to assign to the site
  */
  void setWmc(int wmcval) { wmc_ = wmcval; }

  /**
  @brief set the internal porosity, wmc0_, of this site i.e. its contribution to
  the "weighted mean curvature", wmc_, of the site

  @param wmcval is the value of wmc0_ to assign to the site
  */
  void setWmc0(int wmcval) { wmc0_ = wmcval; }

  /**
  @brief Increment the "weighted mean curvature" of the site.

  @param dwmcval is the increment to make to the existing value of wmc
  */
  void dWmc(int dwmcval) {
    wmc_ += dwmcval;
    if (wmc_ < 0) {
      std::clog << std::endl
                << "Site::dWmc error :  wmc_ < 0! i.e. wmc_ = " << wmc_
                << std::endl;
      std::clog << std::endl << "exit" << std::endl;
      exit(0);
    }
  }

  /**
  @brief Calculate the "weighted mean curvature" of the site.

  The wmc (acronym for weighted mean curvature) is a quantitative ranking of the
  potential for a site's phase to dissolve based on its local environment.  It
  is kind of like mean curvature, but potentially can be weighted to account for
  crystalline anisotropy.

  THAMES is not currently able to quantitatively calculate wmc as defined in
  the metallurgy literature:

      - Taylor, J.E., Cahn, J.W., Handwerker, C.E., Geometric models of crystal
          growth, Acta metall. mater. 40 [7] (1992) 1443--1474.
      - Taylor, J.E., II--Mean curvature and weighted mean curvature,
          Acta metall. mater. 40 [7] (1992) 1475--1485.

  Instead, THAMES uses the digital template method to calculate a quantity that
  is roughly linearly proportional to mean curvature:

      - Bullard, J.W., Garboczi, E.J., Carter, W.C., Fuller, E.J., Numerical
          methods for computing interfacial mean curvature, Comput. Mater. Sci.
          4 (1995) 103--116.

  This provides a ranking of dissolution potential only.
  */
  // void calcWmc(void);

  /**
  @brief Designate the site as a potential dissolution site for a particular
  phase.

  @param pid is the microstructure phase id of the phase that can dissolve at
  the site
  */
  void clearGrowth(void) { growth_.clear(); }

  /**
  @brief Designate the site as a potential growth site for a particular phase.

  @param pid is the microstructure phase id of the phase that can grow at the
  site
  */
  void setGrowthSite(int pid);

  /**
  @brief add a microPhaseId to the growth_ vector belonging to this site;
  growth_ vector contains the microPhaseIds of all microPhases that can
  grow at this site.

  @param pid is the microstructure phase id of the phase that can grow at the
  site
  */
  void addGrowthPhaseId(int pid) { growth_.push_back(pid); }

  /**
  @brief Remove a phase from the list of phases that can grow at the site.

  If a site is currently occupied by aqueous solution, then any of several solid
  phases may be able to grow there.  Therefore, the list of phases in the
  `growth_` vector can have multiple members.

  @param pid is the id of the microstructure phase to remove
  */
  void removeGrowthSite(int pid);

  /**
  @brief Get the entire list of all microPhases that can grow at the site.

  @return the list of ids of all microPhases that can grow at the site
  */
  std::vector<int> getGrowthPhases() const { return growth_; }

  /**
  @brief Set the entire list of all microPhases that can grow at the site.

  @param vect is the vector of ids of all microPhases that can grow at the site
  */
  void setGrowthPhases(std::vector<int> vect) {
    growth_.clear();
    growth_ = vect;
  }

  /**
  @brief Find out if the site is designated as damaged by some kind of
  deterioration.

  @todo Change the method name to isDamaged.

  @return true if the site is damaged (damage_ flag set), or false otherwise
  */
  bool IsDamage() { return damage_; }

  /**
  @brief Designate the site as damaged by some kind of deterioration.

  */
  void setDamage() { damage_ = true; }

  /**
  @brief Set the stress-free volume of the site.

  @note NOT USED.
  @todo Change the method name to something like setStressfreevolume.

  @param vol is the stress-free volume of the site to assign, normalized by
  reference site volume
  */
  void setStressFreeVolume(double vol);

  /**
  @brief Get the stress-free volume of the site.

  @note NOT USED.
  @todo Change the method name to something like getStressfreevolume.

  @return the stress-free volume of the site, normalized by reference site
  volume
  */
  double getStressFreeVolume() { return stressFreeVolume_; }

  /**
  @brief Get the true volume of the site.

  @return the actual volume of the site, normalized by the strain-free site
  volume
  */
  double getTrueVolume() { return trueVolume_; }

  /**
  @brief Set the true volume of the site.

  @param vol is the actual volume of the site, normalized by the strain-free
  site volume
  */
  void setTrueVolume(double vol);

  /**
  @brief Set the expansion strain at the site.

  @param val is the isotropic (<i>i.e.</i>, volumetric) strain at the site
  */
  void setExpansionStrain(double val) {
    if (val > expstrain_) {
      expstrain_ = val;
    }
    return;
  }

  /**
  @brief Get the expansion strain at the site.

  @return the isotropic (<i>i.e.</i>, volumetric) strain at the site
  */
  double getExpansionStrain() { return expstrain_; }

  /**
  @brief One site is equal to another iff their wmc values are equal.

  This kind of function is used to provide a comparison for sorting a list of
  Site objects
  */
  friend bool operator==(const Site &s1, const Site &s2) {
    return (s1.getWmc() == s2.getWmc());
  }

  /**
  @brief One site is less than another iff its wmc value is less.

  This kind of function is used to provide a comparison for sorting a list of
  Site objects
  */
  friend bool operator<(const Site &s1, const Site &s2) {
    return (s1.getWmc() < s2.getWmc());
  }

  /**
  @brief Set the verbose flag

  @param isverbose is true if verbose output should be produced
  */
  void setVerbose(const bool isverbose) { verbose_ = isverbose; }

  /**
  @brief Get the verbose flag

  @return the verbose flag
  */
  bool getVerbose() const { return verbose_; }

}; // End of the Site class

#endif // SRC_THAMESLIB_SITE_H_
