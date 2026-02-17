/**
@file Controller.h
@brief Declaration of the Controller class.
*/

#ifndef SRC_THAMESLIB_CONTROLLER_H_
#define SRC_THAMESLIB_CONTROLLER_H_

#include "AdaptiveTimeController.h"
#include "Exceptions.h"
#include "KineticController.h"
#include "Lattice.h"
#include "Site.h"
#include "ThermalStrain.h"
#include "global.h"
#include <memory>

struct RestoreSite {
  // for each site in site_:
  int microPhaseId;        /**< The microstructure phase assignment */
  std::vector<int> growth; /**< vector of phases that can grow at this
                                site */
  std::vector<int> inGrowInterfacePos; /**< vector of the site position in each
                                            growth interface */
  int inDissInterfacePos; /**< site position in the corresponding dissolution
                               interface */
  int wmc;                /**< total porosity ("surface curvature") at this
                               site */
  int wmc0;               /**< this site internal porosity (its own contribution
                               at wmc_ value) */
  int visit;              /**< reset to 0 */
};

struct RestoreInterface {
  //  from Interface
  unsigned int microPhaseId;           /**< The phase id of the voxels at this
                                            interface */
  std::vector<Isite> growthSites;      /**< The list of all sites eligible for
                                            adjacent growth */
  std::vector<Isite> dissolutionSites; /**< The list of sites eligible for
                                            self-dissolution */
};

struct RestoreSystem {
  // from ChemicalSystem (in fact from KineticController):
  // std::vector<double> ICMoles;
  std::vector<double> DCMoles;
  // from Lattice:
  std::vector<int> count;
  std::vector<int> growthInterfaceSize;
  std::vector<int> dissolutionInterfaceSize;
  std::vector<RestoreSite> site; /**< 1D list of Site objects (site = voxel) */
  // from Interface
  std::vector<RestoreInterface> interface;

  long int numRNGcall_0;
  long int numRNGcallLONGMAX;
  double lastRNG;
};

// from ChemicalSystem:
//   double *ICMoles_;     /**< List of number of moles of each IC in system */
//   double *DCMoles_;             /**< List of moles of each DC */
//   double *prevGEMPhaseMoles_; /**< List of moles of each phase in the system
//   in the previous time step */
//  double *prevGEMPhaseMass_;  /**< List of mass of each phase in the system in
//  the previous time step */ double *prevGEMPhaseVolume_; /**< List of volume
//  of each phase in the system in the previous time step */
//
// from Lattice:
//   std::vector<Site> site_;     /**< 1D list of Site objects (site = voxel) */
//   for each site in site_:
//     unsigned int microPhaseId_;   // The microstructure phase assignment
//     std::vector<unsigned int> growth_; // Vector of phases that can grow at
//     this site
//    double wmc_;                  // total porosity ("surface curvature") at
//    this site double wmc0_;                 // this site internal porosity
//    (its own contribution at wmc_ value)
//    >>int visit_;<<               // reset to 0
// std::vector<Interface> interface_;     //
//   from Interface
//     microPhaseId_; /**< The phase id of the voxels at this interface */
//     std::vector<Isite> growthSites_; /**< The list of all sites eligible
//     foradjacent growth */ std::vector<Isite> dissolutionSites_; /**< The list
//     of sites eligible for self-dissolution */ for each Isite:
//       unsigned int id_; /**< The id of the corresponding Site */
//       int affinity_;    /**< The affinity for growth of a phase at the site
//       */
//       */ bool verbose_;    /**< Flag for whether to produce verbose output */
//       double prob_;     /**< The growth probability of a phase at this site
//       (computed according the affinity) */ double probIni_;
// std::vector<int> count_;               // recreate or restored

/**
@class Controller
@brief Controls the running of simulation iterations.

The `Controller` class is the hub for THAMES simulations.  It has pointers
to all the other major objects that are instantiated by THAMES, including
the KineticController, the ThermalStrain, the ChemicalSystem, and the Lattice
that are associated with the system.

The `Controller` object is also responsible for running each iteration
of the simulation and deciding which modules to run based on whether
hydration, leaching, or sulfate attack is desired.

In THAMES, the `Lattice` class can be thought of as the system.
By itself, it can identify what materials it contains, the properties of
the materials, the temperature, current age, degree of hydration, etc.
However, the `Lattice` class does not possess a <i>driver</i> to
control the microstructure development.  This latter
functionality is contained in the `Controller` class, which operates
directly on the `Lattice` to determine how to modify the lattice at
each time step specified in the phase input file.

THAMES keeps track of physical and chemical data
about individual material phases, including
specific gravity, internal porosity, composition, molar volume, etc.  All the
phases are stored in a phase database.

The ultimate objective of the `Controller` class is to cycle through
the phase input file data, and to modify the lattice accordingly at each
time step.
*/

class Controller {

protected:
  std::string jobRoot_; /**< Root name for all output files */
  Lattice *lattice_;    /**< Pointer to microstructure lattice object */
  KineticController *kineticController_; /**< Pointer to kinetic controller
                                              object */
  ThermalStrain *thermalstr_;            /**< Pointer to the finite element
                                              model object */
  ChemicalSystem *chemSys_;              /**< Pointer to `ChemicalSystem`
                                              object */
  std::vector<double> time_;             /**< List of simulation times for
                                              each iteration */
  std::vector<double> timeInitial_;      /**< List of simulation times for
                                              each iteration */
  vector<double> outputImageTime_;       /**< List of times to output image */
  double outputImageTimeInterval_;       /**< Frequency to output images */

  int simType_; /**< Hydration, leaching, or sulfate attack for now */

  bool attack_;                     /**< for sulfate attack */
  double beginAttackTime_;          /**< Simulation time at which to begin the
                                         attack (leach/sulfate attack) */
  double endAttackTime_;            /**< Simulation time at which to stop the
                                         attack (leach/sulfate attack) */
  double attackTimeInterval_;       /**< Simulation time interval to do the
                                         attack (leach/sulfate attack) */
  std::vector<int> isParrotKilloh_; /**< all microPhaseIds for microPhases
                                         controlled by Parrot-Killoh model */
  int sizePK_;                      /**< size of isParrotKilloh_ vector */
  bool notPKPhase = true;           /**< flag saying if a microPhases is or not
                                         controlled by Parrot-Killoh model */

private:
  double sulfateAttackTime_; /**< Simulation time at which to begin sulfate
                                attack, in hours */
  double leachTime_;         /**< Simulation time at which to begin leaching,
                                      in hours */
  int oldDamageCount_; /**< Number of pixels in the system that were already
                            damaged before the current cycle */
  int newDamageCount_; /**< Number of pixels in the system damaged during the
                            current cycle */
  int allDamageCount_; /**< Total number of pixels in the system that are
                          damaged at the end of the current cycle */

  bool verbose_; /**< Flag for verbose output */
  bool warning_; /**< Flag for warning output */
  // bool xyz_;   /**< Flag for 3D movie data output */
  bool xyzMovie_; /**< Flag for 3D movie data output */
  bool xyzFiles_; /**< Flag for 3D files (.xyz) data output */

  int numMicroPhases_;     /**< Number of microPhases */
  int numGEMPhases_;       /**< Number of GEM phases in the CSD */
  int numICs_;             /**< Number of independent components (IC) */
  int numDCs_;             /**< Number of dependent components (DC) */
  double temperature_;     /**< Temperature [K]*/
  int waterDCId_;          /**< the DCId coresp to DCName = "H2O@" */
  double waterMolarMass_;  /**< the water molar mass corresp. to waterDCId_ */
  int numSites_;           /**< Total number of microStructure voxels */
  double initMicroVolume_; /**< Initial absolute volume of the microStructure */

  double finalHydrationTime_;

  double deltaTime_;
  double lastGoodTime_;
  double stepTimeTHR_;

  /**
  @brief Adaptive time stepping controller.

  Manages timestep selection based on GEMS solver feedback.
  Replaces fixed linear time stepping with intelligent adaptation.
  Only used when useAdaptiveTimeStepping_ is true.
  */
  std::unique_ptr<AdaptiveTimeController> adaptiveTimeController_;

  /**
  @brief Flag to enable/disable adaptive time stepping.

  When true, uses AdaptiveTimeController for timestep decisions.
  When false, uses legacy pre-generated time array.
  */
  bool useAdaptiveTimeStepping_;

  /**
  @brief Maximum relative change in DC moles per timestep (kinetics constraint).

  Used with computeKineticsBasedMaxTimestep() to limit timesteps when
  dissolution/precipitation rates are high. Read from simparams.json
  adaptive_stepping.max_relative_change. Default: 0.05 (5%).
  */
  double maxRelativeChange_;

public:
  /**
  @brief The constructor.

  This is the only Controller constructor provided.  It requires that all the
  auxiliary objects be defined already, including

      - The lattice object
      - The kinetic model object
      - The chemical system object (interface between GEM and THAMES
      - The finite element model for tracking strain and stress

  @param msh is a pointer to the already-instantiated `Lattice` object
  @param kc is a pointer to the already-instantiated `KineticController`
  @param cs is a pointer to the already-instantiated `ChemicalSystem` object
  @param thmstr is a pointer to the already-instantiated `ThermalStrain` object
  @param simtype is the type of simulation to run
  @param jsonFileName is the name of the input parameter file
  @param jobname is the root name to give to all output files
  @param verbose is true if verbose output should be produced
  @param warning is true if warning output should be produced
  @param xyz is true if 3D visualization data should be produced
  */
  Controller(Lattice *msh, KineticController *kc, ChemicalSystem *cs,
             ThermalStrain *thmstr, const int simtype,
             const std::string &jsonFileName, const std::string &jobname,
             const bool verbose, const bool warning, const bool xyz);

  /**
  @brief Run a computational iteration.

  This method launches one computational iteration, which includes

      - Consulting the kinetic controller to determine the number of moles of
          each independent component (IC) to add or subtract from the system
      - Running the GEM thermodynamic calculation
      - Running the finite element code (optionally) to update stress and strain
  states
      - Updating the lattice to reflect the new microstructure

  @param elemTimeInterval is used by the time advance algorithm in case of GEMS
  failure
  */
  // void doCycle(const std::string &statfilename, int choice, double
  // elemTimeInterval);
  void doCycle(double elemTimeInterval);

  /**
  @brief Write exit status to a JSON file for UI consumption.

  This method writes an exit_status.json file to the Result directory
  indicating how the simulation terminated. The UI reads this file to
  display appropriate messages to the user.

  @param exit_code 0 for normal completion, non-zero for abnormal exit
  @param exit_reason Human-readable description of why simulation ended
  @param diagnostics Additional diagnostic information (optional)
  */
  void writeExitStatus(int exit_code, const std::string &exit_reason,
                       const std::string &diagnostics = "");

  /**
  @brief Calculate the state of the system (called by doCycle).

  This method calculates the change in state of the system during a cycle,
  including

      - Consulting the kinetic model to get IC moles dissolved or added
      - (Optionally) determining IC moles to add from an external sulfate
  solution
      - Launching a thermodynamic calculation
      - (Optionally) determining the AFt saturation index for crystallization
  pressure calculations
      - Updating the microstructure
      - Outputting the microstructure phase volume fractions and other data

  @param time is the simulation time [hours]
  @param dt is the change in simulation time used by the kinetic model [hours]
  @param isFirst is true if this is the first state calculation
  (initialization)
  @param cyc is the cycle number for the main controller loop (iteration over
  time)
  @return the node status handle (from ChemicalSystem::calculateState)
  */
  // void calculateState(double time, double dt, bool isFirst, int cyc);
  int calculateState(double &time, double dt, bool isFirst, int cyc);

  /**
  @brief Parse the input JSON file specifying Controller parameters to use.

  Controller parameters that need to be input are

      - Length of time to calculate [hours]
      - Frequency to output microstructure images

  @todo Need some error checking for what to do if a required field is not found

  @param docname is the name of the JSON input file containing the Controller
  parameters
  */
  void parseDoc(const std::string &docname);

  /**
  @brief Set the simulation time at which to begin sulfate attack simulation.

  @param sattacktime is the simulation time to begin sulfate attack [hours]
  */
  void setSulfateAttackTime(const double sattacktime) {
    sulfateAttackTime_ = sattacktime;
  }

  /**
  @brief Get the simulation time at which to begin sulfate attack simulation.

  @return the simulation time to begin sulfate attack [hours]
  */
  double getSulfateAttackTime(void) const { return sulfateAttackTime_; }

  /**
  @brief Set the simulation type

  @param simtype is the simulation type
  */
  // void setSimType(const double simtype) { simType_ = simtype; }

  /**
  @brief Get the simulation type.

  @return the simulation type
  */
  double getSimType(void) const { return simType_; }

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

  /**
  @brief Set the warning flag

  @param iswarning is true if warning output should be produced
  */
  void setWarning(const bool iswarning) { warning_ = iswarning; }

  /**
  @brief Get the warning flag

  @return the warning flag
  */
  bool getWarning() const { return warning_; }

  /**
  @brief Enable or disable adaptive time stepping.

  When enabled, the controller uses GEMS solver feedback to
  intelligently adjust timestep size. When disabled, uses
  the legacy pre-generated time array.

  @param enable true to enable adaptive stepping
  */
  void setUseAdaptiveTimeStepping(bool enable) {
    useAdaptiveTimeStepping_ = enable;
  }

  /**
  @brief Check if adaptive time stepping is enabled.

  @return true if adaptive stepping is enabled
  */
  bool getUseAdaptiveTimeStepping() const { return useAdaptiveTimeStepping_; }

  /**
  @brief Get pointer to the adaptive time controller.

  @return pointer to AdaptiveTimeController, or nullptr if not initialized
  */
  AdaptiveTimeController *getAdaptiveTimeController() {
    return adaptiveTimeController_.get();
  }

  /**
  @brief Set the xyz flag

  @param isxyz is true if 3D visualization data output should be produced
  */
  // void setXyz(const bool isxyz) { xyz_ = isxyz; }

  /**
  @brief Get the xyz flag

  @return the xyz flag
  */
  // bool getXyz() const { return xyz_; }

  /**
  @brief Master function for writing ascii text files

  @param time is the simulation time
  */
  void writeTxtOutputFiles(double time);

  /**
  @brief Master function for writing ascii text files of ICs and DCs

  @param time is the simulation time
  */
  void writeTxtOutputFiles_onlyICsDCs(double time);

  /**
  @brief Convert time in hours to y,d,h,m (string format ***y***d**h**m)

  @param curtime is the simulation time in h
  @return string "***y***d**h**m"
  */
  std::string getTimeString(const double curtime);

}; // End of Controller class

#endif // SRC_THAMESLIB_CONTROLLER_H_
