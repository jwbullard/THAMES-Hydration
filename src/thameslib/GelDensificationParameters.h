/**
@struct GelDensificationParameters
@brief How the gel porosity of C-S-H follows the space available to it.

THAMES assigns each CSHQ voxel a sub-voxel gel porosity phi, and converts
the GEMS CSHQ solid volume to voxels as V_solid / (1 - phi). Rather than
fixing phi from the solid-solution composition, phi follows the saturated
gel density that Konigsberger, Hellmich and Pichler, Cem. Concr. Res. 88
(2016) 170-183, found to be a single function of the specific precipitation
space

    psi = V_w / (V_sCSH + V_w)

for w/c 0.32, 0.40 and 0.48 (their Fig. 3, fitted to the 1H NMR data of
Muller et al., J. Phys. Chem. C 117 (2013) 403). V_sCSH is the solid C-S-H
on the basis of Allen, Thomas and Jennings, Nat. Mater. 6 (2007) 311
(C1.7-S-H1.8, 2.604 g/cm3) and V_w is all liquid water, capillary plus gel.
psi is a property of the paste, not of the C-S-H: only once capillary water
is gone does all liquid water sit in gel pores, and then psi IS the gel
porosity on the Allen basis.

Target saturated gel density:

    regime II  (psiII_III <= psi):  rho = (a - b psi) rhoSolid      (Eq. 42)
    regime III (psi < psiII_III):   rho = rhoSolid (1 - psi) + rho_w psi  (Eq. 13)

Above psiI_II the Eq. 42 value at psiI_II is held. Konigsberger's regime I
has zero gel porosity there, but Muller finds the earliest C-S-H to be the
least dense, so the onset value is the better-supported choice.

Saturated gel density means the same physical thing whatever "solid" is
taken to be, so phi on the GEMS-solid basis follows exactly as

    phi = (rho_s,GEMS - rho) / (rho_s,GEMS - rho_w),   clamped to [0, 1)

The GEMS CSHQ end-members hold ~2.9 H2O per Si (rho_s,GEMS ~ 2.25 g/cm3),
wetter than the Allen solid, so phi reaches 0 while psi is still ~0.3.
*/

#ifndef SRC_THAMESLIB_GELDENSIFICATIONPARAMETERS_H_
#define SRC_THAMESLIB_GELDENSIFICATIONPARAMETERS_H_

struct GelDensificationParameters {
  bool enabled = true;     /**< false restores the composition-based phi */
  double a = 0.901;        /**< Eq. 42 intercept (Konigsberger 2016) */
  double b = 0.411;        /**< Eq. 42 slope (Konigsberger 2016) */
  double rhoSolid = 2.604; /**< Allen et al. 2007 solid C-S-H [g/cm3] */
  double psiI_II = 0.942;  /**< regime I/II boundary (Table 2) */
  double psiII_III = 0.426; /**< regime II/III boundary (Table 2) */
};

/**
@brief Target saturated C-S-H gel density at specific precipitation space psi.

@param psi is Konigsberger's specific precipitation space [0, 1]
@param p holds the fitted coefficients and regime boundaries
@param rhoWater is the density of liquid water [g/cm3]
@return the saturated gel density [g/cm3]
*/
inline double gelDensityTarget(double psi,
                               const GelDensificationParameters &p,
                               double rhoWater) {
  if (psi > p.psiI_II)
    psi = p.psiI_II;
  if (psi >= p.psiII_III)
    return (p.a - p.b * psi) * p.rhoSolid;
  return p.rhoSolid * (1.0 - psi) + rhoWater * psi;
}

/// Upper bound on gel porosity, so the solid always occupies some volume.
constexpr double GEL_POROSITY_MAX = 0.99;

/**
@brief Gel porosity of C-S-H for one step, including the envelope rule.

The target density is an average over all the gel, so on its own it would
re-densify existing C-S-H at once and pull the gel envelope (solid plus gel
pores, the volume that becomes voxels) inward. Precipitation into gel pores
cannot do that: in Konigsberger's model the gel volume grows in regime II and
is constant in regime III. That holds there by an exact space balance that an
approximate psi does not reproduce, so it is enforced directly:

    solid growing   -> the envelope may not shrink (it may grow: new gel
                       also precipitates into open space)
    solid shrinking -> the envelope may neither grow nor shrink faster than
                       the solid. Its two limits are gel dissolving at its
                       current porosity (envelope shrinks in proportion) and
                       solid leaving from inside the gel, as in
                       decalcification (envelope unchanged, porosity rises).
                       Dissolving gel cannot push outward into its
                       neighbours. (carbonation, leaching)

These reduce to bounds on phi. There is no adjustable parameter.

All volumes share one frame (any units); only their ratios matter.

@param psi is Konigsberger's specific precipitation space
@param solidDensity is the density of the GEMS C-S-H solid [g/cm3]
@param waterDensity is the density of liquid water [g/cm3]
@param solidVolume is the C-S-H solid volume now
@param solidCommitted is the solid volume at the last accepted step (0: none)
@param envelopeCommitted is the envelope volume at the last accepted step
@param p holds the fitted coefficients and regime boundaries
@return phi in [0, GEL_POROSITY_MAX]
*/
inline double densifiedGelPorosity(double psi, double solidDensity,
                                   double waterDensity, double solidVolume,
                                   double solidCommitted,
                                   double envelopeCommitted,
                                   const GelDensificationParameters &p) {
  const double target = gelDensityTarget(psi, p, waterDensity);
  double phi = (solidDensity - target) / (solidDensity - waterDensity);
  if (phi < 0.0)
    phi = 0.0;

  if (envelopeCommitted > 0.0 && solidCommitted > 0.0) {
    // Holding the envelope fixed gives this phi.
    const double phiSameEnvelope = 1.0 - solidVolume / envelopeCommitted;
    if (solidVolume >= solidCommitted) {
      if (phi < phiSameEnvelope)
        phi = phiSameEnvelope;
    } else {
      // Shrinking in proportion keeps the committed porosity.
      const double phiProportional = 1.0 - solidCommitted / envelopeCommitted;
      if (phi < phiProportional)
        phi = phiProportional;
      if (phi > phiSameEnvelope)
        phi = phiSameEnvelope;
    }
  }
  if (phi > GEL_POROSITY_MAX)
    phi = GEL_POROSITY_MAX;
  return phi;
}

#endif // SRC_THAMESLIB_GELDENSIFICATIONPARAMETERS_H_
