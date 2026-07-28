/**
@file KineticController.h
@brief Declaration of the KineticController class.

@section Introduction
The `KineticController` class keeps track of all the
different kinetic models that govern the rate of hydration.

*/

#ifndef SRC_THAMESLIB_KINETICCONTROLLER_H_
#define SRC_THAMESLIB_KINETICCONTROLLER_H_

#include "global.h"
// #include "../Resources/include/nlohmann/json.hpp"
#include "ChemicalSystem.h"
#include "JMAKGrowth.h"
#include "JMAKParameters.h"
#include "KineticData.h"
#include "KineticModel.h"
#include "Lattice.h"
#include "ParrotKillohModel.h"
#include "PozzolanicModel.h"
#include "StandardKineticModel.h"
#include "Exceptions.h"

// using json = nlohmann::json;

/**
@class KineticController
@brief Manages kinetic models

THAMES allows some flexibility in defining different types of kinetic models.
*/

class KineticController {

private:
  int numPhases_;            /**< Total number of phases in the kinetic model */
  ChemicalSystem *chemSys_;  /**< Pointer to the ChemicalSystem object for this
                                  simulation */
  Lattice *lattice_;         /**< Pointer to the lattice object holding the
                                  microstructure */
  std::vector<KineticModel *> phaseKineticModel_; /***< Kinetic model for each phase */
  double temperature_;       /**< Temperature [K] */
  double refT_;              /**< Reference temperature for PK model [K] */
  double sulfateAttackTime_; /**< Time at which sulfate attack simulation starts
                                  [hours] */
  double leachTime_;         /**< Time at which leaching simulation starts [hours] */

  std::vector<std::string> name_;              /**< List of names of phases in the
                                                    kinetic model */
  std::vector<int> microPhaseId_;              /**< List of microstructure ids that
                                                    are in kinetic model */
  std::vector<double> initScaledMass_;         /**< List of initial scaled masses */
  std::vector<double> scaledMass_;             /**< List of scaled masses */
  std::vector<double> specificSurfaceArea_;    /**< List of specific surface areas */
  std::vector<double> refSpecificSurfaceArea_; /**< List of reference specific surface
                                                    areas */
  std::vector<bool> isKinetic_;    /**< vector setting the isKinetic property of each
                                        microPhase in the system; true for a kinetic
                                        controlled microPhase */
  int DCNum_;                      /**< Number of DCs in chemical system */
  int GEMPhaseNum_;                /**< Number of GEM phases in chemical system */
  bool verbose_;                   /**< Flag for verbose output */
  bool warning_;                   /**< Flag for warnining output */
  bool useNucleationKinetics_ = false;
      /**< CNT global switch — set by Controller after parseDoc */

  // ---------- JMAK per-voxel growth state (2026-07-28) ------------------
  //
  // Per-phase JMAK bookkeeping. Populated per-phase during
  // parseKineticData / makeModel and updated cycle-by-cycle inside
  // calculateKineticStep. Only active for phases with BOTH a nucleation
  // sub-block AND a jmak sub-block in the phase's kinetic_data JSON.
  //
  // Non-JMAK CNT phases (nucleation block but no jmak block) continue
  // to use the pre-existing "1 voxel = 1 nucleation event" placement
  // path via accumulateNucleation / drainNucleationInteger — the two
  // paths coexist and are selected per-phase via jmakEnabled_[midx].
  //
  // See docs/POST_ALPHA_TODOS.md "CNT growth model needs JMAK-per-voxel"
  // for the physical rationale and mathematical model. Free functions
  // that do the per-cycle math live in JMAKGrowth.{h,cc}.

  /**
  @brief Per-phase flag: true if this phase uses JMAK growth (nucleation
  block + jmak block both present in kinetic_data).
  */
  std::vector<bool> jmakEnabled_;

  /**
  @brief Per-phase JMAK parameters (Avrami n and morphology alpha).
  */
  std::vector<JMAKParameters> jmakParams_;

  /**
  @brief Per-phase JMAK global moment accumulator (M_0..M_3, G_acc).
  Updated once per cycle in calculateKineticStep by adding J(t)·dt and
  G(t)·dt contributions.
  */
  std::vector<jmak::GlobalMoments> jmakGlobals_;

  /**
  @brief Per-generation snapshot record.

  When a new generation of voxels is seeded (Poisson-triggered by the
  per-phase seed accumulator crossing 1.0), we record how many voxels
  it contains plus the global-moment snapshot at that moment. From then
  on, the generation's extended-volume-per-voxel is a difference-of-
  moments computation against the current global moments (see
  jmak::extendedVolumePerVoxel).
  */
  struct JMAKGeneration {
    int Nvoxels;                             /**< voxel count in this generation */
    jmak::GenerationMomentsAtSeed seed;      /**< moment snapshot at seed time */
  };

  /**
  @brief Per-phase list of generations. Each generation carries an
  independent JMAK growth trajectory driven by the same global moments.
  */
  std::vector<std::vector<JMAKGeneration>> jmakGenerations_;

  /**
  @brief Per-phase fractional-voxel seeding accumulator.

  Each cycle we compute the expected number of newly-seeded voxels as
  `J(t) * V_unseeded_electrolyte * dt`. This is usually less than 1 per
  cycle, so we accumulate the fraction and only create a new generation
  when the accumulator crosses 1.0.
  */
  std::vector<double> jmakSeedAccumulator_;

  /**
  @brief Per-phase count of voxels already placed in the lattice via the
  JMAK path.

  Lattice sync compares this to `floor(sum_g N_g * X_g)` each cycle and
  places any positive delta via `nucleatePhaseRnd`. Kept as an integer
  so the delta-based placement is clean.
  */
  std::vector<int> jmakVoxelsInLattice_;

  std::vector<double> DCMoles_;    /**< vector of all DC moles - after the dissolution
                                        corresponding to the current time step - to be
                                        sent to GEMS */
  std::vector<double> DCMolesIni_; /**< vector of all DC moles - before to start the
                                        dissolution corresponding to the current time
                                        step*/
  std::vector<double> scaledMassIni_; /**< List of scaled masses before a given time
                                           step*/

  std::vector<int> impurityDCID_;     /**< vector of the DCIds of all impurities
                                           contained and able to be eliberated by
                                           dissolution of each kinetic controlled
                                           microPhases (in order, DCIds of: K2O, Na2O,
                                           Per, SO3) */
  // Per-cycle impurity release vectors [mol]. Indexed by kinetic-model
  // index (midx) — one entry per kinetic-controlled microphase. Populated
  // inside calculateKineticStep from the phase's K2O/Na2O/MgO/SO3 mass
  // fractions and its dissolved mass this cycle, then added to the
  // corresponding IC (impurityDCID_) columns of DCMoles_.
  std::vector<double> impurity_K2O_;
  std::vector<double> impurity_Na2O_;
  std::vector<double> impurity_Per_;
  std::vector<double> impurity_SO3_;

  int pKMsize_;                  /**< dimension of the phaseKineticModel_ vector */

  double initScaledCementMass_;  /**< initial scaled cement mass i.e. the sum of all
                                      scalled masses corresponding to the microPhases
                                      controlled by the Parrot-Killoh model */
  double hydTimeIni_;            /**< the time elapsed before the current time step */
  int waterDCId_;                /**< the DCId coresp to DCName = "H2O@" */
  double beginAttackTime_;       /**< Simulation time at which to begin the attack
                                      (sulfate attack for now); hydration stops when
                                      the current time equqls beginAttackTime_ */

  std::vector<double> surfaceAreaIni_; /**< vector of surface areas of each microPhase
                                       before to start the dissolution for a given
                                       time step */

  std::vector<double> initialMicroPhaseSI_; /**< Saturation indices for each microPhase
                                                 from pre-equilibration GEMS calculation.
                                                 Used for SI-aware initial timestep
                                                 estimation. */
  bool hasSIData_; /**< Flag indicating if SI data from pre-equilibration is available */

public:
  /**
  @brief Default constructor.

  This constructor is not used in THAMES.  It just establishes default values
  for all the member variables.

  @note NOT USED.
  */
  KineticController();

  /**
  @brief Overloaded constructor.

  This constructor is used in THAMES.

  @param cs is a pointer to the ChemicalSystem object for the simulation
  @param lattice is a pointer to the Lattice object holding the microstructure
  @param jsonFileName is the name of the JSON file with the input for the
  kinetic model
  @param verbose is true if verbose output should be produced
  @param warning is false if suppressing warning output
  */
  KineticController(ChemicalSystem *cs, Lattice *lattice,
                    const std::string &jsonFileName, const bool verbose,
                    const bool warning);

  /**
  @brief Destructor does nothing.
  */
  virtual ~KineticController();

  /**
  @brief Initialize the kinetic data structure
  */
  void initKineticData(struct KineticData &kineticData) {
    kineticData.name = "";
    kineticData.microPhaseId = kineticData.GEMPhaseId = kineticData.DCId = 0;
    kineticData.type = "thermo";
    kineticData.scaledMass = 0.0;
    kineticData.surfaceAreaMultiplier = 1.0;
    kineticData.temperature = kineticData.reftemperature = 293.15;
    kineticData.k1 = kineticData.k2 = kineticData.k3 = 1.0;
    kineticData.n1 = kineticData.n3 = 1.0;
    kineticData.critDOR = 0.0;
    kineticData.dissolutionRateConst = 0.0;
    kineticData.diffusionRateConstEarly = 0.0;
    kineticData.diffusionRateConstLate = 0.0;
    kineticData.siexp = kineticData.dfexp = 0.0;
    kineticData.dorexp = kineticData.ohexp = 0.0;
    kineticData.dissolvedUnits = 1.0;
    kineticData.activationEnergy = 0.0;
    kineticData.loi = kineticData.sio2 = kineticData.al2o3 = kineticData.cao =
        0.0;
    // Reset the CNT block so a per-phase parse of one phase does not leak
    // into the next phase's KineticData (which is reused across the parse loop).
    kineticData.nucleation.reset();
    // Same reason: reset the JMAK block (2026-07-28). JMAK is a growth-
    // side extension of CNT and inherits the same per-phase leak risk.
    kineticData.jmak.reset();
    // Same reason: reset the SaturatingRate blocks. Otherwise a preceding
    // SaturatingRate phase's parameters would silently apply to a
    // Standard/Pozzolanic/ParrotKilloh phase parsed after it.
    kineticData.saturatingDissolution.reset();
    kineticData.saturatingPrecipitation.reset();
  }

  /**
  @brief Master method controlling the parsing of JSON input to the kinetic
  model.

  @param docName is the name of the (purported) JSON input file
  */
  void parseDoc(const std::string &docName);

  /**
  @brief Parse the input data for one phase in the JSON input file.

  @param cdi is an iterator over the JSON data
  @param numEntry is the number of solid entries in the JSON file, will be
  incremented
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseMicroPhases(const json::iterator cdi, int &numEntry,
                        struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for one phase in the JSON input file.

  @todo Need error checking for what to do if a required entry is not present

  @param p is an iterator over the JSON data
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticData(const json::iterator p,
                        struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for the Parrot-Killoh kinetic model.

  @todo Need error checking for what to do if a required entry is not present

  @param pp is an iterator over the JSON data
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticDataForParrotKilloh(const json::iterator pp,
                                       struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for the standard kinetic model.

  @todo Need error checking for what to do if a required entry is not present

  @param pp is an iterator over the JSON data
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticDataForStandard(const json::iterator pp,
                                   struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for the pozzolanic kinetic model.

  @todo Need error checking for what to do if a required entry is not present

  @param pp is an iterator over the JSON data
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticDataForPozzolanic(const json::iterator pp,
                                     struct KineticData &kineticData);

  /**
  @brief Parse the SaturatingRate-model kinetic data block.

  Fills `kineticData.saturatingDissolution` (required) and
  `kineticData.saturatingPrecipitation` (optional) from a phase's
  `kinetic_data` JSON block plus the scalar fields shared with other
  model types (`surfaceAreaMultiplier`, `dissolvedUnits`,
  `activationEnergy`). Also invokes `parseNucleationBlock` at the tail
  so SaturatingRate phases can carry a `nucleation` block just like
  Standard and Pozzolanic phases.

  @param pp is an iterator into the phase's `kinetic_data` JSON block
  @param kineticData is the KineticData struct being populated
  */
  void parseKineticDataForSaturatingRate(const json::iterator pp,
                                         struct KineticData &kineticData);

  /**
  @brief Parses a single `dissolution` or `precipitation` sub-block
  inside a SaturatingRate phase's `kinetic_data`.

  Expected JSON (all three fields required):
      { "rateConstant": <double>, "B": <double>, "n": <double> }

  If the sub-block is absent, the target `std::optional` is left
  unchanged (callers decide whether the absence is an error).

  @param p is the iterator into the phase's `kinetic_data` JSON block
  @param blockName is either "dissolution" or "precipitation"
  @param target is the optional to populate on success
  */
  void parseSaturatingRateSubBlock(
      const json::iterator p, const std::string &blockName,
      std::optional<SaturatingRateParameters> &target);

  /**
  @brief Parses the optional `nucleation` sub-block inside a phase's `kinetic_data`.

  Reads only the `.value` field of each parameter (`gamma`, `theta`, `A0`);
  the accompanying `range` and `provenance` fields are user-facing metadata
  and are ignored at runtime. If no `nucleation` block is present, the
  `kineticData.nucleation` optional is left empty (CNT disabled for the phase).

  @param pp is the iterator into the phase's `kinetic_data` JSON block
  @param kineticData is the KineticData struct being populated
  */
  void parseNucleationBlock(const json::iterator pp,
                            struct KineticData &kineticData);

  /**
  @brief Per-cycle JMAK update for one JMAK-enabled phase.

  Advances the phase's global moment accumulator, seeds new generations
  based on Poisson-triggered accumulator crossings, evaluates all
  generations' transformed fractions X_g, and syncs the lattice to
  match `floor(sum_g N_g * X_g)` — placing new voxels via
  `Lattice::nucleatePhaseRnd` and updating DCMoles / microPhaseMass
  through the same "correct per-100g scaling" path as the classical
  CNT placement block does. See docs/POST_ALPHA_TODOS.md
  "CNT growth model needs JMAK-per-voxel" for the mathematical model.

  Called from the CNT dispatch inside `calculateKineticStep` when
  `jmakEnabled_[midx]` is true, in place of the classical
  computeNucleationVoxels/accumulateNucleation/drainNucleationInteger
  path.

  @param midx        kinetic-model index (per-phase array offset)
  @param timestep    cycle timestep [hours]
  @param cyc         cycle number (for logging)
  */
  void updateJMAKPhase(int midx, double timestep, int cyc);

  /**
  @brief Get the scaled mass of the phase in the kinetic model.

  The scaled mass of a phase is its mass percent on a total solids basis.

  @return the vector of scaled masses [percent solids]
  */
  std::vector<double> getScaledMass() const { return scaledMass_; }

  /**
  @brief Get the scaled mass of one phase

  The scaled mass of a phase is its mass percent on a total solids basis.

  @note NOT USED.

  @param midx is the microstructure id of the phase to query
  @return the vector of scaled masses [percent solids]
  */

  double getScaledMass(const int midx) { return scaledMass_[midx]; }

  /**
  @brief Get the <i>initial</i> mass of the microstructure phases

  The scaled mass of a phase is its mass percent on a total solids basis.

  @return the vector of initial scaled masses [percent solids]
  */
  std::vector<double> getInitScaledMass() const { return initScaledMass_; }

  /**
  @brief Get the <i>initial</i> scaled mass of one microstructure phase

  The scaled mass of a phase is its mass percent on a total solids basis.

  @note NOT USED.

  @param midx is the microstructure id of the phase to query
  @return the initial scaled mass [percent solids]
  */
  double getInitScaledMass(const int midx) { return initScaledMass_[midx]; }

  /**
  @brief Compute normalized initial microstructure phase masses

  @note NOT USED

  Given the initial masses of all phases in the microstructure,
  this method scales them to 100 grams of solid.  In the process,
  this method also sets the initial moles of water in the
  chemical system definition.
  */
  void calcPhaseMasses(void);

  /**
  @brief Get sum of phase masses

  */
  double getSolidMass(void);

  /**
  @brief Get the list of specific surface areas

  @return the vector of specific surface areas [m2/kg]
  */
  std::vector<double> getSpecificSurfaceArea() const { return specificSurfaceArea_; }

  /**
  @brief Get the specific surface area of one phase


  @param midx is the microstructure id of the phase to query
  @return the specific surface area [m2/kg]
  */
  double getSpecificSurfaceArea(const int midx) { return specificSurfaceArea_[midx]; }

  /**
  @brief Get the list of reference specific surface areas

  @return the vector of reference specific surface areas [m2/kg]
  */
  std::vector<double> getRefSpecificSurfaceArea() const {
    return refSpecificSurfaceArea_;
  }

  /**
  @brief Get the reference specific surface area of one phase


  @param midx is the microstructure id of the phase to query
  @return the reference specific surface area [m2/kg]
  */
  double getRefSpecificSurfaceArea(const int midx) {
    return refSpecificSurfaceArea_[midx];
  }

  /**
  @brief Make a kinetic model for a given phase

  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void makeModel(struct KineticData &kineticData);

  /**
  @brief Get the ChemicalSystem object for the simulation used by the kinetic
  model.

  @note NOT USED.

  @return a pointer to the ChemicalSystem object
  */
  ChemicalSystem *getChemSys() const { return chemSys_; }

  /**
  @brief Set the simulation time at which to begin external sulfate attack.

  @param sattacktime is the simulation time to begin sulfate attack [hours]
  */
  // void setSulfateAttackTime(double sattacktime) {
  //   sulfateAttackTime_ = sattacktime;
  // }

  /**
  @brief Get the simulation time at which to begin external sulfate attack.

  @note NOT USED.

  @return the simulation time to begin sulfate attack [hours]
  */
  // double getSulfateAttackTime(void) const { return sulfateAttackTime_; }

  /**
  @brief Set the simulation time at which to begin leaching.

  @param leachtime is the simulation time to begin leaching [hours]
  */
  // void setLeachTime(double leachtime) { leachTime_ = leachtime; }

  /**
  @brief Get the simulation time at which to begin leaching.

  @note NOT USED.

  @return the simulation time to begin leaching [hours]
  */
  // double getLeachTime(void) const { return leachTime_; }

  /**
  @brief Get the list of phase names used by the kinetic model.

  @note NOT USED.

  @return the vector of names of phases in the kinetic model
  */
  // std::vector<std::string> getName() const { return name_; }

  /**
  @brief Get the name of phase with a given index in the kinetic model.

  @param i is the index of the phase in the kinetic model
  @return the name of the phase with index i
  */
  // std::string getName(const unsigned int i) const { return name_[i]; }

  /**
  @brief Set kinetic model DC moles

  */
  // void setKineticDCMoles() {
  //   int size = phaseKineticModel_.size();
  //   for (int i = 0; i < size; ++i) {
  //     phaseKineticModel_[i]->setKineticDCMoles();
  //   }
  //   return;
  // }

  /**
  @brief Zero kinetic model DC moles

  */
  // void zeroKineticDCMoles() {
  //   int size = phaseKineticModel_.size();
  //   for (int i = 0; i < size; ++i) {
  //     phaseKineticModel_[i]->zeroKineticDCMoles();
  //   }
  //   return;
  // }

  /**
  @brief Set the effect of pozzolans on Parrot-Killoh kinetics

  */
  void setPozzEffectOnPK(void);

  /**
  @brief Master method for implementing phase kinetics

  In a given time step, a certain number of moles of each phase may
  dissolve or precipitate.  This function determines the number of moles of each
  phase to change, based on the time interval being simulated.

  @remark This method is very long and several parts are hard-coded when they
  should be made more general.

  @todo Split this method into more convenient chunks
  @todo Make the methods more general, less hardwiring of parameters
  @todo Make the local variable names more descriptive

  @param time is the current time of the simulation [hours]
  @param timestep is the time interval to simulate [hours]
  @param cyc is the iteration number in main iteration loop in
  Controller::doCycle
  */
  void calculateKineticStep(double time, const double timestep, int cyc);

  /**
  @brief reset the dissolved number of moles for a kinetic controlled microPhase,
  having microPhaseId = pId, when the lattice configuration cannot be changed
  according to the kinetic models/GEMS previsions.

  @param cyc is the iteration number in main iteration loop in
  Controller::doCycle
  @param pId is the microPhaseId of a kinetic controlled microPhase
  @param scaledMass is the amount from the initial mass of this microPhase, mass
  computed by the previuos call of calculateKineticStep, amount that cannot be
  dissolved because of the current lattice configuration
  @param timestep is the time interval to simulate [hours]
  */
  void updateKineticStep(int cyc, int pId, double scaledMass, double timestep);

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

  @param iswarning is true if verbose output should be produced
  */
  void setWarning(const bool iswarning) { warning_ = iswarning; }

  /**
  @brief Get the warning flag

  @return the warning flag
  */
  bool getWarning() const { return warning_; }

  /**
  @brief Get the number of moles of every dependent component (DC) in the
  system.

  @return a vector containing the number of moles of every dependent component
  (DC) in the system.
  */
  std::vector<double> getDCMoles(void) { return DCMoles_; }

  // vector<bool> getIsKinetic(void) { return isKinetic_; }

  /**
  @brief Set the initial hydration time (hydTimeIni_) at its previous value
  when the lattice configuration cannot be changed according to the
  kinetic models/GEMS previsions.

  @param val is the previos value of the initial hydration time (hydTimeIni_).
  */
  void setHydTimeIni(double val) { hydTimeIni_ = val; }

  /**
  @brief Set the time at which to begin sulfate attack, in hours.

  @param val is the time at which to begin sulfate attack, in hours.
  */
  void setIniAttackTime(const double val) { beginAttackTime_ = val; }

  /**
  @brief Set the initial saturation indices from pre-equilibration GEMS calculation.

  This method stores the saturation indices for each microPhase, which are used
  to provide more accurate initial dissolution rate estimates. The SI values
  allow the kinetic models to compute realistic driving forces rather than
  assuming SI ≈ 0 (maximum driving force).

  @param microPhaseSI vector of saturation indices indexed by microPhaseId
  */
  void setInitialSaturationIndices(const std::vector<double> &microPhaseSI);

  /**
  @brief Estimate the maximum initial dissolution rate across all kinetic models.

  This method queries each kinetic model for its estimated initial dissolution
  rate and returns the maximum rate found. This is used for adaptive timestep
  initialization - the fastest-reacting phase determines the initial timestep.

  If SI data from pre-equilibration is available (via setInitialSaturationIndices),
  the SI-aware rate estimation is used for more accurate results.

  @return the maximum estimated dissolution rate [1/hour], or 0 if no kinetic
          models are defined
  */
  double getMaxInitialDissolutionRate() const;

  /**
  @brief Compute the maximum timestep based on kinetics to limit DC mole changes.

  This method queries each kinetic model for its current molar rate, then computes
  the timestep that would limit the relative change in DC moles to a specified
  maximum (default 5%). This provides a physics-based upper bound on the timestep
  to prevent overshooting equilibrium positions.

  The algorithm:
  1. For each kinetically-controlled DC with non-zero moles:
     - Get current molar rate from kinetic model [mol/100g/h]
     - Get current DC moles from chemical system
     - Compute dt that would cause maxRelativeChange: dt = maxChange * moles / |rate|
  2. Return the minimum dt across all DCs (most restrictive constraint)

  This complements the GEMS-feedback-based adaptive time stepping by preventing
  large chemical changes that could cause GEMS solver difficulties.

  @param maxRelativeChange Maximum allowed relative change per timestep (default 5%)
  @return Maximum timestep in hours, or a large value (1e6) if no kinetics constraint
  */
  double computeKineticsBasedMaxTimestep(double maxRelativeChange = 0.05) const;

  /**
  @brief CNT-driven per-phase cap on dt.

  Iterates over kinetic models. For each phase carrying a nucleation
  block, applies two caps and takes the tighter:

    (1) Electrolyte-fraction cap: `capFraction * count_[ELECTROLYTEID]`.
        Prevents CNT placement from overshooting available saturated
        pore volume in one cycle.

    (2) Mass-availability cap: for each IC that the phase's DC has
        non-zero stoichiometry, `aqICMoles[ic] / (vVoxel/vMolar_DC *
        stoich)` — the number of voxels the currently-aqueous IC moles
        could support. The min across ICs. Prevents CNT from placing
        voxels that `Lattice::changeMicrostructure` then rolls back
        because the mass balance can't feed them. Added 2026-07-24 after
        S4 validation of `SaturatingRateModel` found the CNT/GEMS
        placement/roll-back oscillation; see `docs/SATURATING_RATE.md`
        §6 and the corresponding POST_ALPHA_TODOS entry.

  If `computeNucleationVoxels(dtProposed)` exceeds the effective cap,
  shrinks dt so N lands exactly at the cap; returns the min of the
  phase-level shrunk dts and dtProposed. Fixed-J-within-cycle scaling
  is exactly what the cap exists to protect, so predictive rescaling
  is safe; no retry loop needed.

  Returns dtProposed unchanged when `useNucleationKinetics_` is false
  or no phase carries a nucleation block.

  @param dtProposedHours the proposed cycle timestep [hours]
  @param capFraction     N_cap_electrolyte / N_electrolyte, from
                         simparams.json (`nucleationCapFraction`)
  @return dt to actually use [hours]; always <= dtProposedHours
  */
  double computeNucleationBasedMaxTimestep(double dtProposedHours,
                                            double capFraction) const;

  /**
  @brief Enable/disable Classical Nucleation Theory placement at cycle time.

  When false, the CNT placement block in calculateKineticStep is a no-op
  even when phases carry nucleation blocks in their kinetic_data, and
  computeNucleationBasedMaxTimestep returns dtProposed unchanged.
  Set by Controller after parsing simparams.json.
  */
  void setUseNucleationKinetics(bool enable) {
    useNucleationKinetics_ = enable;
  }

  /**
  @brief Query the CNT global switch.
  */
  bool getUseNucleationKinetics() const { return useNucleationKinetics_; }

}; // End of KineticController class

#endif // SRC_THAMESLIB_KINETICCONTROLLER_H_
