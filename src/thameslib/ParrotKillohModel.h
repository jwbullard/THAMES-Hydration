/**
@file ParrotKillohModel.h
@brief Declaration of the ParrotKillohModel class.

@section Introduction
This class implements the Parrot and Killoh (PK) model of
1984 [1]---in the same way as described by Lothenbach and
Winnefeld [2]---for cement clinker phases.

The PK model accounts only for dissolution of the four major clinker minerals:
alite (impure C<sub>3</sub>S), belite (impure C<sub>2</sub>S), aluminate
(C<sub>3</sub>A), and ferrite (C<sub>4</sub>AF).  The model provides
mathematical equations for the rates of three broad rate-controlling phenomena:

    -# nucleation and growth,
    -# early-age diffusion, and
    -# late-age diffusion.

For each clinker mineral <i>i</i>, the three rate equations are

@f{eqnarray*}
R_{ng}(i) &=& \frac{A K_1(i)}{N_1(i)} \left( 1 - \alpha \right) \left( - \ln
\left( 1 - \alpha \right) \right)^{1 - N_1(i)} \\
R_{de}(i) &=& \frac{K_2(i) \left( 1 - \alpha \right)^{2/3}}{1 - \left( 1 -
\alpha \right)^{1/3}} \\ R_{dl}(i) &=& K_3(i) \left( 1 - \alpha \right)^{N_3(i)}
@f}

where <i>A</i> is the overall surface area of the cement powder
(cm<sup>2</sup>/kg) and \f$\alpha\f$ is the overall degree of hydration on a
mass basis. <i>K</i><sub>1</sub>, <i>K</i><sub>2</sub>, <i>K</i><sub>3</sub>,
<i>N</i><sub>1</sub>, and <i>N</i><sub>3</sub> are constants defined for each
clinker mineral.  The values of these constants used by Parrot and Killoh in
Ref. [1] are shown in the table.  In any particular time interval, the predicted
rate of dissolution of a clinker mineral

@f[
R(i) = \min (R_{ng}(i),R_{de}(i),R_{dl}(i)) \cdot f(\text{RH}) \cdot g(w/c)
@f]

where \f$f(\text{RH})\f$ and \f$g(w/c)\f$ account empirically for the influences
of relative humidity and water-cement mass ratio (w/c), respectively, according
to

@f{eqnarray*}
f(\text{RH}) &=& \left( \frac{ \text{RH} - 0.55}{0.45} \right)^4 \\
g(w/c) &=&
\begin{cases}
1 & \text{if}\ \alpha \le 1.333\, w/c \\
(1 + 4.444 w/c - 3.333 \alpha)^4 & \text{otherwise}
\end{cases}
@f}

The new degree of hydration at the end of the time interval is calculated
according to the difference equation

@f[
\alpha(t+\Delta t) = \alpha(t) + R(t) \Delta t
@f]

<table>
<caption id="multi_row">Empirical constants used by Parrot and Killoh</caption>
<tr><th>Parameter               <th>Alite           <th>Belite <th>Aluminate
<th>Ferrite <tr><td><i>K</i><sub>1</sub>    <td>1.5             <td>0.5 <td>1.0
<td>0.37 <tr><td><i>N</i><sub>1</sub>    <td>0.7             <td>1.0 <td>0.85
<td>0.7 <tr><td><i>K</i><sub>2</sub>    <td>0.05            <td>0.006 <td>0.05
<td>0.015 <tr><td><i>K</i><sub>3</sub>    <td>1.1             <td>0.2 <td>1.0
<td>0.4 <tr><td><i>N</i><sub>3</sub>    <td>3.3             <td>5.0 <td>3.2
<td>3.7 \end{center}
</table>

@section References

    -# Parrot, L.J., Killoh, D.C., Prediction of cement hydration, British
Ceramic Proceedings 35 (1984) 41-53.
    -# Lothenbach, B., Winnefeld, F., Thermodynamic modelling of the hydration
of portland cement, Cement and Concrete Research 36 (2006) 209--226.

*/

#ifndef SRC_THAMESLIB_PARROTKILLOHMODEL_H_
#define SRC_THAMESLIB_PARROTKILLOHMODEL_H_

#include "global.h"
#include "Exceptions.h"
#include "ChemicalSystem.h"
#include "KineticController.h"
#include "KineticData.h"
#include "KineticModel.h"
#include "Lattice.h"

// Ref specific surface area adjusted downward from
// 385.0 m2/kg (the published value) to 372.0 m2/kg to
// better agree with kinetics when calculating actual
// specific surface area from the microstructure.

const double RefSpecificSurfaceArea = 372.0;

/**
@class ParrotKillohModel
@brief Handles the Parrot and Killoh(1984) kinetic model of clinker ractions

The Parrot and Killoh model [1] is used to empirically estimate the
mass fraction of each clinker phase that dissolves in a unit time.  Eventually
this can be expanded to handle other kinetically controlled phases outside the
Parrot and Killoh model, such as the growth of portlandite or C--S--H.
*/

class ParrotKillohModel : public KineticModel {

protected:
  double wsRatio_; /**< water-solid mass ratio */
  double wcRatio_; /**< water-cement mass ratio */

  double k1_; /**< List of Parrot and Killoh <i>K</i><sub>1</sub> values */
  double k2_; /**< List of Parrot and Killoh <i>K</i><sub>2</sub> values */
  double k3_; /**< List of Parrot and Killoh <i>K</i><sub>3</sub> values */
  double n1_; /**< List of Parrot and Killoh <i>N</i><sub>1</sub> values */
  double n3_; /**< List of Parrot and Killoh <i>N</i><sub>3</sub> values */
  double dorHcoeff_;
  double critDOR_;   /**< List of critical degrees of hydration for w/c
                          effect in the Parrot and Killoh model */
  double pfk_;       /**< Multiplicative factor for k's to account for
                          effects of pozzolanic additions */
  double rh_;        /**< relative humidity */
  double rhFactor_;  /**< relative humidity factor, i.e. the correction
                          of the hydration rate taking into account the
                          ambient relative humidity */
  double arrhenius_; /**< arrhenius factor */

public:
  /**
  @brief Default constructor.

  This constructor is not used in THAMES.  It just establishes default values
  for all the member variables.

  @note NOT USED.
  */
  ParrotKillohModel();

  /**
  @brief Overloaded constructor.

  This constructor is the one invoked by THAMES.  It can only be called once the
  various other objects for the simulation are allocated and constructed.

  @param cs is a pointer to the ChemicalSystem object for the simulation
  @param lattice is a pointer to the Lattice object holding the microstructure
  @param kineticData is the collection of kinetic parameters already stored
  @param verbose is true if verbose output should be produced
  @param warning is false if suppressing warning output
  */
  ParrotKillohModel(ChemicalSystem *cs, Lattice *lattice,
                    struct KineticData &kineticData, const bool verbose,
                    const bool warning);

  /**
  @brief Get the type of kinetic model

  @return a string indicating the model type
  */
  std::string getType() const { return (ParrotKillohType); }

  /**
  @brief Set the w/s mass ratio of the system for the kinetic model equations.

  @note NOT USED.

  @param wsr is the w/s mass ratio to set
  */
  void setWsRatio(double wsr) { wsRatio_ = wsr; }

  /**
  @brief Get the w/s mass ratio of the system used by the kinetic model
  equations.

  @return the w/s mass ratio
  */
  double getWsRatio() const { return wsRatio_; }

  /**
  @brief Get the list of <i>K</i><sub>1</sub> values for clinker phases in the
  PK model.

  @note NOT USED.

  @return the <i>K</i><sub>1</sub> value for the phase in the PK model
  */
  double getK1() const { return k1_; }

  /**
  @brief Get the list of <i>K</i><sub>2</sub> values for clinker phases in the
  PK model.

  @note NOT USED.

  @return the <i>K</i><sub>2</sub> value for the phase in the PK model
  */
  double getK2() const { return k2_; }

  /**
  @brief Get the list of <i>K</i><sub>3</sub> values for clinker phases in the
  PK model.

  @note NOT USED.

  @return the <i>K</i><sub>2</sub> value for the phase in the PK model
  */
  double getK3() const { return k3_; }

  /**
  @brief Set the multiplicative factor for k's due to pozzolanic effects

  @param pfk is the multiplicative factor to set.
  */
  void setPfk(const double pfk) {
    pfk_ = pfk;
    if (pfk_ < 1.0e-5)
      pfk_ = 1.0e-5;
  }

  /**
  @brief Get the multiplicative factor for k's due to pozzolanic effects

  @note NOT USED.

  @return the pfk_ value for the phase in the PK model
  */
  double getPfk() const { return pfk_; }

  /**
  @brief Get the list of <i>N</i><sub>1</sub> values for clinker phases in the
  PK model.

  @note NOT USED.

  @return the <i>N</i><sub>2</sub> value for the phase in the PK model
  */
  double getN1() const { return n1_; }

  /**
  @brief Get the list of <i>N</i><sub>3</sub> values for clinker phases in the
  PK model.

  @note NOT USED.

  @return the <i>N</i><sub>3</sub> value for the phase in the PK model
  */
  double getN3() const { return n3_; }

  /**
  @brief Get the critical degrees of reaction for w/c effects in the kinetic
  model.

  @note NOT USED.

  @return the critical degree of reaction for this phase
  */
  double getCritDOR() const { return critDOR_; }

  /**
  @brief Master method for implementing one kinetic time step.

  Overloaded from base class to handle Parrot and Killoh model.

  @todo Split this method into more convenient chunks
  @todo Make the methods more general, less hardwiring of parameters
  @todo Make the local variable names more descriptive

  @param timestep is the time interval to simulate [h]
  @param temperature is the absolute temperature during this step [K]
  @param rh is the internal relative humidity
  @param scaledMass is C-style array of the normalized mass of each
  microstructure phase [g/100 g]
  @param massDissolved is the C-style array of dissolved mass of each
  microstructure phase [g/100g]
  @param cyc is the cycle number (iteration of main loop)
  @param totalDOR is the total degree of reaction [dimensionless]
  */
  virtual void calculateKineticStep(const double timestep, double &scaledMass,
                                    double &massDissolved, int cyc,
                                    double totalDOR, bool doTweak,
                                    bool &doNotModif);

  /**
  @brief Estimate the initial dissolution rate for this clinker phase.

  Uses the Parrot-Killoh rate equations at DOR ≈ 0 to estimate the
  initial dissolution rate. At very low DOR, the nucleation and growth
  rate typically dominates.

  @return the estimated dissolution rate [1/hour]
  */
  double estimateInitialDissolutionRate() const override;

  /**
  @brief Get the current molar rate of dissolution.

  Calculates the current dissolution rate based on the Parrot-Killoh model
  equations. The rate depends on the degree of reaction (DOR), not on
  saturation index.

  Note: The PK model only handles dissolution (positive rates). Precipitation
  is handled by GEMS thermodynamics for the hydration products.

  @param scaledMass the current scaled mass of this phase [g/100g solids]
  @return the current molar rate [mol/hour per 100g solids]
          Always positive (dissolution only) for clinker phases
  */
  double getCurrentMolarRate(double scaledMass) const override;

}; // End of ParrotKillohModel class

#endif // SRC_THAMESLIB_PARROTKILLOHMODEL_H_
