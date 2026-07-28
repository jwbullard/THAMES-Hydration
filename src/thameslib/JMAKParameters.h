/**
@struct JMAKParameters
@brief Johnson-Mehl-Avrami-Kolmogorov parameters for one phase.

Populated by KineticController from a `jmak` sub-block inside a phase's
`kinetic_data`. Read at simulation time by the JMAK growth machinery in
KineticController that evolves each cohort's transformed volume fraction.

The Avrami exponent n governs the transformed-fraction curve
X(t) = 1 - exp(-K t^n) in the constant-J, constant-G limit. Canonical
values:

  - n = 4: 3D isotropic growth, continuous nucleation
  - n = 3: 3D isotropic growth, site-saturated nucleation (all nuclei
           formed instantaneously at t = 0)
  - n = 3: 2D disk growth, continuous nucleation
  - n = 2: 2D disk growth, site-saturated
  - n = 1 + d * beta, general form (d = growth dimensionality,
           beta = 0 site-saturated, beta = 1 continuous)

For CNT-nucleated phases inside a 1 um^3 lattice voxel with non-ideal
(possibly fractal or non-spherical) morphology, values in [2.5, 4] are
physically defensible. The parser validates n in this range and rejects
values outside; users who need to test outside-range behavior should
edit the source or extend the validated range.

The morphology coefficient alpha appears in the extended-volume
integrand for the fully-3D case:

  Y_c(t) = alpha * V_voxel * integral_{t_c}^{t} J(tau) * r(tau, t)^3 dtau

where r(tau, t) = integral_{tau}^{t} G(s) ds is the growth radius. For
3D isotropic spheres alpha = 4*pi/3. Users would rarely need to change
alpha; it is included for completeness and to allow overriding for
non-spherical crystal habits (e.g., plate-like Ca(OH)2 might use a
smaller effective alpha).

@note For the initial implementation only n = 4 (integer exponent) is
supported cleanly via the moment-decomposition method in
JMAKGrowth.{h,cc}. Non-integer n values require per-cohort numerical
integration of the growth kernel, which is planned but not landed.
*/

#ifndef SRC_THAMESLIB_JMAKPARAMETERS_H_
#define SRC_THAMESLIB_JMAKPARAMETERS_H_

struct JMAKParameters {
  double n;      /**< Avrami exponent; defensible in [2.5, 4] */
  double alpha;  /**< Morphology coefficient; 3D isotropic = 4*pi/3 */
};

#endif  // SRC_THAMESLIB_JMAKPARAMETERS_H_
