/**
@file Homogenization.cc
@brief Method definitions for self-consistent effective moduli.

See Homogenization.h for why the self-consistent scheme is used here and what
is given up by choosing it.
*/

#include "Homogenization.h"

#include <cmath>

namespace homog {

Result selfConsistent(const std::vector<Phase> &phases, double tolerance,
                      int maxIterations) {
  Result result;

  ///
  /// Normalize, and take the Voigt average as the starting guess. Voigt is
  /// the stiffest possible answer, so the iteration approaches the solution
  /// from above and cannot wander into negative stiffness on the way.
  ///

  double totalFraction = 0.0;
  for (const Phase &p : phases) {
    if (p.volumeFraction > 0.0)
      totalFraction += p.volumeFraction;
  }
  if (totalFraction <= 0.0)
    return result;

  double kEff = 0.0;
  double gEff = 0.0;
  for (const Phase &p : phases) {
    if (p.volumeFraction <= 0.0)
      continue;
    const double f = p.volumeFraction / totalFraction;
    kEff += f * p.bulkModulus;
    gEff += f * p.shearModulus;
  }

  ///
  /// Keep the Voigt average as a fixed yardstick for convergence. Measuring
  /// the change against the CURRENT estimate fails exactly where it matters:
  /// an assemblage below the percolation limit walks down towards zero
  /// stiffness, and a change of 1e-34 against an estimate of 1e-34 is a
  /// hundred percent no matter how converged the answer is.
  ///

  const double reference = (kEff > gEff) ? kEff : gEff;

  if (kEff <= 0.0 && gEff <= 0.0) {
    /// Nothing but pore space: no stiffness, and nothing to iterate on.
    result.converged = true;
    return result;
  }

  for (int iter = 1; iter <= maxIterations; ++iter) {
    result.iterations = iter;

    ///
    /// F* is the shear analogue of 4 G*/3 in the bulk equation. It vanishes
    /// with G*, which is exactly what happens once the solid stops
    /// percolating, so it is guarded rather than assumed positive.
    ///

    const double denomFStar = kEff + 2.0 * gEff;
    const double fStar =
        (denomFStar > 0.0)
            ? (gEff / 6.0) * (9.0 * kEff + 8.0 * gEff) / denomFStar
            : 0.0;

    double kNumerator = 0.0, kDenominator = 0.0;
    double gNumerator = 0.0, gDenominator = 0.0;

    for (const Phase &p : phases) {
      if (p.volumeFraction <= 0.0)
        continue;
      const double f = p.volumeFraction / totalFraction;

      const double kWeight = p.bulkModulus + (4.0 / 3.0) * gEff;
      if (kWeight > 0.0) {
        kNumerator += f * p.bulkModulus / kWeight;
        kDenominator += f / kWeight;
      }

      const double gWeight = p.shearModulus + fStar;
      if (gWeight > 0.0) {
        gNumerator += f * p.shearModulus / gWeight;
        gDenominator += f / gWeight;
      }
    }

    ///
    /// A zero denominator means every phase that could contribute has both
    /// zero modulus and zero weight, so the effective medium has no
    /// stiffness of that kind left. Converging to zero is the answer.
    ///

    const double kNew = (kDenominator > 0.0) ? kNumerator / kDenominator : 0.0;
    const double gNew = (gDenominator > 0.0) ? gNumerator / gDenominator : 0.0;

    const double change =
        (reference > 0.0)
            ? (std::fabs(kNew - kEff) + std::fabs(gNew - gEff)) / reference
            : 0.0;

    ///
    /// Under-relax. The plain fixed point oscillates for an assemblage of
    /// many phases with widely different moduli, which is every real cement
    /// paste, and near the scheme's own percolation limit it can oscillate
    /// without ever settling. Taking a partial step damps that out at the
    /// cost of more iterations, which are cheap here.
    ///

    const double relaxation = 0.5;
    kEff += relaxation * (kNew - kEff);
    gEff += relaxation * (gNew - gEff);

    if (change < tolerance) {
      result.converged = true;
      break;
    }
  }

  result.bulkModulus = kEff;
  result.shearModulus = gEff;
  return result;
}

} // namespace homog
