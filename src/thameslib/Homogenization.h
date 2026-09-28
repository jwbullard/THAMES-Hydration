/**
@file Homogenization.h
@brief Effective elastic moduli of a multiphase microstructure by the
       self-consistent scheme.

Supplies the drained bulk modulus a poromechanical estimate of autogenous
shrinkage needs at every output time. The finite element solver in
ElasticModel gives the exact answer for a given microstructure but costs
minutes per call, which rules it out on a per-output basis; this costs
microseconds.

**Why self-consistent.** Every phase is treated on the same footing, embedded
in the effective medium itself. No phase has to be nominated as the matrix,
which matters because a paste at the age autogenous shrinkage is largest has
no matrix: it is a percolating granular network at forty to fifty percent
porosity. Mori-Tanaka would need a matrix chosen, and the natural choice
changes from pore fluid to C-S-H as hydration proceeds, putting a step in the
middle of the curve. The self-consistent estimate also goes to zero as the
solid stops percolating, which is the right limit at the setting point where
ASTM C1698 starts measuring.

**What it costs in exchange.** The estimate is not a bound, where Mori-Tanaka
coincides with a Hashin-Shtrikman bound for well-ordered two-phase
composites. Its stiffness vanishes at a solid fraction of one half for
spherical phases, and that threshold is a property of the approximation, not
of cement: THAMES computes the real rigidity percolation threshold of the
same microstructure, so the two can be compared. At late ages, where the
paste genuinely is C-S-H with clinker and portlandite inclusions,
Mori-Tanaka is the better motivated choice, which is why the multiscale
models in the literature use the self-consistent scheme at the gel level and
Mori-Tanaka at the paste level.

Phases are taken as spherical. Reference: Hill 1965; Berryman 1980; for the
application to cement, Bernard, Ulm & Lemarchand, Cem. Concr. Res. 33 (2003)
1293.

No THAMES dependencies. Testable in isolation; see
`src/unit_tests/test_homogenization.cc`.
*/

#ifndef SRC_THAMESLIB_HOMOGENIZATION_H_
#define SRC_THAMESLIB_HOMOGENIZATION_H_

#include <vector>

namespace homog {

/**
@struct Phase
@brief One constituent of the microstructure.
*/
struct Phase {
  double volumeFraction = 0.0; /**< need not be normalized; the caller's
                                    fractions are rescaled to sum to one */
  double bulkModulus = 0.0;    /**< K, any consistent unit */
  double shearModulus = 0.0;   /**< G, same unit as K */
};

/**
@struct Result
@brief Effective moduli and how the solve went.
*/
struct Result {
  double bulkModulus = 0.0;  /**< effective K, in the phases' unit */
  double shearModulus = 0.0; /**< effective G, in the phases' unit */
  bool converged = false;    /**< false if the iteration ran out first */
  int iterations = 0;        /**< how many passes were taken */
};

/**
@brief Get the effective moduli of a phase assemblage.

Solves the self-consistent equations for spherical phases,

    Sum_r f_r (K_r - Keff) / (K_r + 4 Geff / 3) = 0
    Sum_r f_r (G_r - Geff) / (G_r + Feff)       = 0,
    Feff = (Geff / 6) (9 Keff + 8 Geff) / (Keff + 2 Geff)

(starless names on purpose: the usual K-star notation puts a comment
terminator in the middle of a Doxygen block)

by fixed point iteration from the Voigt average. An assemblage with too
little solid to percolate converges to zero stiffness, which is the answer,
not a failure.

@param phases is the assemblage; entries with zero volume fraction are
       ignored, and phases with zero moduli (pores) are perfectly acceptable
@param tolerance is the relative change at which to stop
@param maxIterations is the cap on passes
@return the effective moduli, with converged false if the cap was reached
*/
Result selfConsistent(const std::vector<Phase> &phases,
                      double tolerance = 1.0e-8, int maxIterations = 500);

} // namespace homog

#endif // SRC_THAMESLIB_HOMOGENIZATION_H_
