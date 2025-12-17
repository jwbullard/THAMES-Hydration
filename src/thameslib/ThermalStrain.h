/**
@file ThermalStrain.h
@brief Declare the ThermalStrain class.

This class implements a finite element model for soliving linear elastic
problems in which the deformation is due to thermal mismatch or some other
form of stress-free strain, rather than an externally applied strain as
implemented in the AppliedStrain class.
*/

#ifndef SRC_THAMESLIB_THERMALSTRAIN_H_
#define SRC_THAMESLIB_THERMALSTRAIN_H_

#include "ElasticModel.h"
#include "global.h"

/**
@class ThermalStrain
@brief Solves the linear elastic state for systems that are thermally strained.

This class is dervied from the ElasticModel base class, just as the
AppliedStrain class is.

@section intro_sec Background
Program adjusts dimensions of unit cell,
[(1 + macrostrain) times dimension],
in response to phases that have a non-zero eigenstrain and arbitrary
elastic moduli tensors.  All six macrostrains can adjust their values
(3D model), and are stored in the last two positions in the displacement
vector `u_`, as described below.  Periodic boundaries are maintained.  In
the comments below, notes are given for sections of code
that the user might have to change for a particular problem.

@subsection p_anv_v Problem and Variable Definitions
The problem being solved is the minimization of the elastic energy

@f[
    \frac{1}{2} u \cdot A \cdot u + b \cdot u + C + T \cdot u + Y
@f]
where <i>A</i> is the Hessian matrix composed of the stiffness matrices
(`dk`) for each voxel/element, <i>b</i> is a constant vector and
<i>C</i> is a scalar constant, both the latter being functions of the strain
and the periodic boundary conditions.  Up to this point, the energy
is exactly the same as for any applied strain problem (see the `AppliedStrain`
class documentation).  In addition, <i>T</i> is the thermal energy term
that is linear in the displacement, <i>u</i> is the displacement vector
field.

Some of the class variables are used to build the components of this equation.
    - `zcon` is a small array that computes the thermal strain energy
            associated with macrostrains (<i>C</i> term),
    - `u` is the displacement field
    - `gb` is the energy gradient vector
    - `h` and `Ah` are auxiliary vectors
    - `dk` is the single pixel stiffness matrix
    - `pix` is the phase identification vector
    - `ib` is the integer matrix for mapping labels from the 1--27
         nearest neighbor labelling to the 1D system labelling.
    - The array `prob` contains the volume fractions of the <i>i</i>-th phase
    - `strxx`, etc. are the six independent (Voigt notation) volume
        averaged stresses
    - `sxx`, etc. are the six independent (Voigt notation) volume averaged
        strains (not counting the thermal strains)
    - `cmod[i][6][6]` gives the elastic moduli tensor of the <i>i</i>-th phase,
    - `eigen[i][6] gives the six independent elements of the eigenstrain
        tensor for the <i>i</i>-th phase (Voigt notation)
    - `dk([i][8][3][8][3]` is the stiffness matrix of the
        <i>i</i>th phase.
    - `nphase` is the number of phases being considered in the problem,
        and is set by the user.

@subsection dim Dimensions
The main arrays of the problem, `u`, `gb`, `h`, `Ah`, `b`, and `T` are
dimensioned  as `(nx*ny*nz)+2`, which is the number of nodal displacements
plus two for the macrostrains.  Currently the program assumes the number
of different phases is 100, since `phasemod` and `eigen` (the moduli and
eigenstrains for each phase) and `dk` are dimensioned to have at most 100
different kinds. This is easily changed, by changing the dimension of
these three variables manually throughout the program. The parameter `nphase`
gives the number of phases.

Many of the variables and methods needed for this class are already defined
in the base class `ElasticModel`, so this header file only adds a relatively
small number of members needed to handle thermal strain.

@note Program is set up to allow the macrostrains,
which control the overall size of the system, to be dynamic
variables, which are adjusted in order to minimize the overall
energy. That means that if there are no eigenstrains specified
for any of the phases, the overall strain will always relax to
zero. If it is desired to simply apply a fixed strain, with no
eigenstrains, then in subroutines Energy and Dembx, one must
zero out the elements of gb (in Energy and in Dembx) that
correspond to the macrostrains. This is easily done.
This will fix the gradients of the macrostrains to always to be
zero, so that they will not change, so the applied strain (initial
values of the macrostrains) will remain fixed.

@note Manual available at NISTIR 6269 from NTIS or at:
http://ciks.cbt.nist.gov/~garbocz/manual/man.html

@warning Read the manual before using this class.
*/
class ThermalStrain : public ElasticModel {

protected:
  /**
  @brief Should the displacement vector be initialized?

  If the energy is being relaxed in the simulation for the very first time,
  then we initialize the displacement vector to zero.  Otherwise, if this
  is just an incremental recalculation, we keep the displacement vector equal to
  its relaxed state the previous time, which should shorten the number of
  conjugate gradient steps needed to relax the incremented system.
  */
  bool isFirst_;

  double Y_;    /**< Constant energy term needed
                    to correct for the apparent
                    strain caused by periodic BCs */
  int boxsize_; /**< Defines neighborhood within
                    which to perform relaxation */
  int boxnum_;  /**< Number of elements in the local
                    neighborhood */
  /**
  @brief The stopping criterion for local relaxation within a subvolume.

  This is the maximum allowed square of the energy gradient fo which the
  whole system can be considered to have converged to the energy minimum.  This
  is the value to which the quantity `loalgg_` (the computed square of the
  local subvolume energy gradient, `gb_`) is compared.  In most cases,
  `localgtest_` will be set equal to a value <i>abc</i> times the total
  number of elements in the subvolume, so that when `localgg_` is less than
  `localgtest_`, the RMS gradient per element is less than
  @f$\sqrt{abc}@f$.
  */
  double localgtest_;
  double localgg_;                         /**< The square of the gradient
                                               `gb_`<sup>2</sup> */
  std::vector<double> tstrength_;          /**< Assumed tensile strength of
                                        various phases [MPa] */
  std::vector<std::vector<double>> b0_;    /**< xx component of linear b vector
                                     three terms for each phase, so
                                     b0_[nphase][3] */
  std::vector<std::vector<double>> b1_;    /**< yy component of linear b vector
                                     three terms for each phase, so
                                     b0_[nphase][3] */
  std::vector<std::vector<double>> b2_;    /**< zz component of linear b vector
                                     three terms for each phase, so
                                     b0_[nphase][3] */
  std::vector<std::vector<double>> b3_;    /**< xz component of linear b vector
                                     three terms for each phase, so
                                     b0_[nphase][3] */
  std::vector<std::vector<double>> b4_;    /**< yz component of linear b vector
                                     three terms for each phase, so
                                     b0_[nphase][3] */
  std::vector<std::vector<double>> b5_;    /**< xy component of linear b vector
                                     three terms for each phase, so
                                     b0_[nphase][3] */
  std::map<int, std::vector<int>> exp_;    /**< @todo Find out what this is :
                                      list of sites ids and their
                                     coordinates where local expansion
                                     occurs */
  std::vector<std::vector<double>> eigen_; /**< Six components of the
                                   eigenstrain tensor for each mesh element, so
                                   eigen_[ns][6] */
  std::vector<std::vector<double>> T_;     /**< Linear (in displacement) thermal
                                 energy term, one for each phase,
                                 so T_[nphase][3] */
  std::vector<std::vector<std::vector<std::vector<double>>>>
      zcon_; /**< Stores thermal energy associated
                with macrostrains */
  std::vector<std::vector<std::vector<double>>> ss_; /**< Shear stress tensor
                                         components, one at every element, so
                                         the dimensions are ss_[ns][3][3] */
  int kmax_; /**< the number of relaxation steps for a given elastic computation
              */

public:
  /**
  @brief Constructor.

  The constructor initializes the finite element mesh to dimensions specified
  by the user, and sets the maximum number of phases that can be handled.

  @param nx is the number of elements along the x dimension
  @param ny is the number of elements along the y dimension
  @param nz is the number of elements along the z dimension
  @param dim is the total number of mesh elements
  @param cs is a pointer to the governing ChemicalSystem object
  @param hasAggregateSlabe is true if the microstructure has an aggregate slab
  @param verbose is true if verbose output should be produced
  @param warning is false if warning messages should be suppressed
  */
  ThermalStrain(int nx, int ny, int nz, int dim, ChemicalSystem *cs,
                const bool hasAggregateSlab, const bool verbose,
                const bool warning);

  /**
  @brief Destructor.

  */
  ~ThermalStrain() { tstrength_.clear(); }

  /**
  @brief Set up the stiffness matrices and parameters for computing elastic
  energy.

  @todo Change the nskip variable to be boolean.

  @param nx is the number of elements in the x dimension
  @param ny is the number of elements in the y dimension
  @param nz is the number of elements in the z dimension
  @param ns is the total number of mesh elements
  @param nphase is the maximum number of phases in the system
  @param iskip is 0 if this is the first cycle, nonzero otherwise
  */
  // void femat(int nx, int ny, int nz, int ns, int nphase, int iskip);
  void femat(int iskip);

  /**
  @brief Compute gradient of the b vector with respect to the macrostrains.

  Since `b_` is linear with respect to the macrostrains, the derivative with
  respect to any one of them can be computed simply by setting that macrostrain,
  within the method, be equal to one, and setting all the other macrostrains to
  zero.  The structure of the calculation is similar to the loop in `femat` for
  b.

  @param nx is the number of elements in the x dimension
  @param ny is the number of elements in the y dimension
  @param nz is the number of elements in the z dimension
  @param ns is the total number of mesh elements
  @param exx is the xx macrostrain
  @param eyy is the yy macrostrain
  @param ezz is the zz macrostrain
  @param exz is the xz macrostrain
  @param eyz is the yz macrostrain
  @param exy is the xy macrostrain
  */
  // void bgrad(int nx, int ny, int nz, int ns, double exx, double eyy, double
  // ezz,
  //            double exz, double eyz, double exy);
  void bgrad(double exx, double eyy, double ezz, double exz, double eyz,
             double exy);

  /**
  @brief Compute the quadratic term in the macrostrains.

  This corrects for the quadratic contribution to the energy that comes from
  the periodic boundary conditions.  It sets it up as a (2,3) x (2,3) matrix
  that couples to the six macrostrains.

  @todo Change the name of this method to setZcon or something similar.

  @param ns is the total number of mesh elements
  @param nx is the number of elements in the x dimension
  @param ny is the number of elements in the y dimension
  @param nz is the number of elements in the z dimension
  */
  // void constfunc(int ns, int nx, int ny, int nz);
  void constfunc(void);

  /**
  @brief Conjugate gradient relaxation of the elastic energy of the whole mesh.

  @param ns is the total number of mesh elements
  @param gg is the energy gradient
  @param ldemb is the maximum number of conjugate gradient steps to use during
  this call
  @param kkk is the number of times the method has been called
  @return the number of conjugate gradient steps used in this call
  */
  // int dembx(int ns, double gg, int ldemb, int kkk);
  int dembx(int ldemb, int kkk);

  /**
  @brief Conjugate gradient relaxation of the elastic energy of a mesh
  subvolume.

  @param boxsize is the total number of elements in the subvolume
  @param x is the number of elements in the x dimension in this subvolume
  @param y is the number of elements in the y dimension in this subvolume
  @param z is the number of elements in the z dimension in this subvolume
  @param localldemb is the maximum number of conjugate gradient steps to use
  during this call
  @param kkk is zero if this is the first call, causing `h` to be built up.
  @return the number of conjugate gradient steps used in this call
  */
  // int localDembx(int boxsize, int x, int y, int z, int localldemb, int kkk);
  int localDembx(int xlo, int xhi, int ylo, int yhi, int zlo, int zhi,
                 int localldemb, int kkk);

  /**
  @brief Compute the total elastic energy and the gradient in the energy.

  The energy, `utot`, and the gradient, `gb_` are computed for the regular
  displacements and for the macrostrains.

  @param nx is the number of elements in the x dimension
  @param ny is the number of elements in the y dimension
  @param nz is the number of elements in the z dimension
  @param ns is the total number of mesh elements
  @return the computed elastic energy [units, perhaps J]
  */
  // double energy(int nx, int ny, int nz, int ns);
  double energy(void);

  /**
  @brief Compute the stress components for the relaxed system.

  @param nx is the number of elements in the x dimension
  @param ny is the number of elements in the y dimension
  @param nz is the number of elements in the z dimension
  @param ns is the total number of mesh elements
  */
  // void stress(int nx, int ny, int nz, int ns);
  void stress();

  /**
  @brief Control function for finding the minimum energy state.

  This function calculates the energy and the energy gradient by calling
  the energy function, and compares the calculated gradient to the stopping
  criterion.  If the absolute value of the gradient is still larger than
  the stopping criterion, then the conjugate gradient is called again to relax
  the system further, and the process is repeated, stopping only if the total
  number of conjugate gradient iterations exceeds a prescribed maximum.

  @note Argument time is NOT USED.

  @param cyc is the iteration number in main iteration loop in
  Controller::doCycle
  @param time is the simulation time [hours] (not used currently)
  @param kmax is the maximum number of total conjugate gradient iterations
  allowed
  */
  // void relax(double time, int kmax);
  void relax(int cyc, int kmax);

  /**
  @brief Control function for elastic relaxation within a subvolume.

  This function controls the process of finding a minimum elastic energy
  configuration within a cubic subvolume of the microstructure.  It operates
  in the same way as the `relax` function, but is passed the size and
  coordinates of the subvolume.

  @note Argument index is NOT USED.

  @param boxsize is the edge length of the cubic subvolume
  @param x is the x-coordinate of the center element of the subvolume
  @param y is the y-coordinate of the center element of the subvolume
  @param z is the z-coordinate of the center element of the subvolume
  @param index is another variable that is not used
  */
  void localRelax(int boxsize, int x, int y, int z, int index);

  void localRelax(int xlo, int xhi, int ylo, int yhi, int zlo, int zhi, int x,
                  int y, int z, int index);

  /**
  @brief Master function for executing the finite element calculation.

  This is the main controlling function of the entire calculation.  It
  sets up the displacement field, deciding whether to initialize at some
  default or--- if the system has already been solved prior to the current
  perturbation---preconditions the displacement field to be equal to the
  previous relaxed value for computational efficiency.  Once the displacement
  field is initialized, the relaxation control functions are called.
  First, the local relaxations are done in those areas where the greatest
  perturbation to the displacement field is expected, and then the
  global relaxation is executed.

  @note Argument index is NOT USED.

  @param cyc is the iteration number in main iteration loop in
  Controller::doCycle
  @param time is the simulation time [hours]
  @param p_vectPhId is a pointer to a vector containing the prior equilibrium
  displacement field if needed
  @param exx is the xx component of the macrostrain
  @param eyy is the yy component of the macrostrain
  @param ezz is the zz component of the macrostrain
  @param exz is the xz component of the macrostrain
  @param eyz is the yz component of the macrostrain
  @param exy is the xy component of the macrostrain
  */
  // void Calc(double time, std::string fileName, double exx, double eyy, double
  // ezz,
  //           double exz, double eyz, double exy);
  void calc(int cyc, double time, std::vector<int> *p_vectPhId, double exx,
            double eyy, double ezz, double exz, double eyz, double exy);

  /**
  @brief Set the six eigenstrain components of a phase.

  @todo Check the validity of the siteid and throw an exception of out of
  bounds.

  @param siteid is the index of a site in the FE mesh (1D ordering)
  @param xx is the xx component of the eigenstrain (component 1 in Voigt
  notation)
  @param yy is the yy component of the eigenstrain (component 2 in Voigt
  notation)
  @param zz is the zz component of the eigenstrain (component 3 in Voigt
  notation)
  @param xz is the xz component of the eigenstrain (component 4 in Voigt
  notation)
  @param yz is the yz component of the eigenstrain (component 5 in Voigt
  notation)
  @param xy is the xy component of the eigenstrain (component 6 in Voigt
  notation)
  */
  void setEigen(int siteid, double xx, double yy, double zz, double xz,
                double yz, double xy) {
    eigen_[siteid][0] = xx;
    eigen_[siteid][1] = yy;
    eigen_[siteid][2] = zz;
    eigen_[siteid][3] = xz;
    eigen_[siteid][4] = yz;
    eigen_[siteid][5] = xy;
  }

  /**
  @brief Set the six eigenstrain components of every site to zero.

  */
  void setEigen(void) {
    eigen_.clear();
    eigen_.resize(ns_);
    for (int i = 0; i < ns_; i++) {
      eigen_[i].resize(6, 0.0);
    }
  }

  /**
  @brief Set the coordinates for a local expansion strain site.

  The coordinates for an expansion strain site are stored in the `Lattice`
  object, as a private variable called `expansion_coordin_`.  The coordinates
  are originally set by this object in the `Lattice::applyExp` method.  The
  `Controller` object gets the values of these coordinates from the `Lattice`
  object and then calls this function to set the coordinates in the `exp_`
  member of the `ThermalStrain` class.

  @note This seems like an unnecessarily obscure way to associate (x,y,z)
  coordinates with a lattice site, by calling it "expansion coordinates".
  Why not just make a `ThermalStrain` member to store the coordinates for
  each site, calling it something like `elementCoordinates`, and then
  copy it straight from the `Lattice` object to the `ThermalStrain` object?

  @note The `exp_` map class member stores a list of all the sites ids and
  their associated (x,y,z) coordinates where local expansion will occur
  due to, for example, a phase transformation during sulfate attack.

  @todo Fix the way that site coordinates are stored and accessed throughout
  THAMES, so that it is more transparent.  This will involve reworking or
  renaming some of the `Lattice` class and `Controller` class members.  We
  can do this by turning the `exp_` map into a simple vector of site ids, and
  then using a separate vector of vectors to store the (x,y,z) coordinates.

  @param index is the id of the element in the 1D ordering of elements
  @param cor is the (x,y,z) coordinates of that site
  */
  void setExpansionCoord(int index, std::vector<int> cor);

  /**
  @brief Clear the map of expansion sites.

  @note NOT USED.
  */
  void setExpansionCoord(void) { exp_.clear(); }

  /**
  @brief Get a component of the elastic stress at an element.

  @todo Determine what units are returned for stress.

  @todo Change the name of this method to something like
  getElementStressComponent.

  @param i is the element id in the 1D ordering of elements
  @param j is the component of the stress to retrieve (Voigt notation from 0 to
  5)
  @return the stress component [GPa]
  */
  double getEleStress(int i, int j);

  /**
  @brief Get the vector of the elastic stress at an element.

  @todo Determine what units are returned for stress.

  @todo Change the name of this method to something like
  getElementStressComponent.

  @param i is the element id in the 1D ordering of elements (a vector in Voigt
  notation having 6)
  @return the stress vector [GPa]
  */
  std::vector<double> getEleStress(int i) { return elestress_[i]; }

  /**
  @brief Get a component of the the elastic strain at an element.

  @todo Change the name of this method to something like
  getElementStrainComponent.

  @param i is the element id in the 1D ordering of elements
  @param j is the component of the strain to retrieve (Voigt notation from 0 to
  5)
  @return the strain component
  */
  double getEleStrain(int i, int j);

  /**
  @brief Get the tensile strength at a site in the 3D element mesh.

  This method is passed the index of an element in the 1D ordering of mesh
  elements.  It then looks up the phase residing at that element using the
  object's `pix_` array, and returns the value stored in the `tstrength_` array
  for that phase.

  @todo Determine what units are returned for tensile strength.

  @todo Change the name of thie method to getElementTensileStrength.

  @param index is the index of the element in the 1D ordering of elements
  @return the tensile strength of the phase at that element
  */
  double getTstrength(int index) {
    int phaseid = pix_[index];
    return tstrength_[phaseid];
  }

}; // End of ThermalStrain class

#endif // SRC_THAMESLIB_THERMALSTRAIN_H_
