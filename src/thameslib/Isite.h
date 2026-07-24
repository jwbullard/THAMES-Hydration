/**
@file Isite.h
@brief Declare a class that keeps track of the chemical potential of sites

This class instantiates to objects that have

    - an id that defines a corresponding `Site` object,
    - an <i>affinity</i> value, which is a qualitative chemical potential
variable, similar in some ways to weighted mean curvature of an interface
*/

#ifndef SRC_THAMESLIB_ISITE_H_
#define SRC_THAMESLIB_ISITE_H_

#include "global.h"
#include "Exceptions.h"

/**
@class Declaration of the Isite class.

*/
class Isite {

private:
  int id_;          /**< The id of the corresponding Site */
  int affinityInt_; /**< The affinity for growth of a phase at the site */
  bool verbose_;    /**< Flag for whether to produce verbose output */
  double prob_;     /**< Growth probability of a phase at this site (computed
                         from affinity). @note DECLARED BUT UNUSED — no
                         getter, no setter, initialized via the constructor's
                         default `prb=0` and then never referenced. Candidate
                         for deletion; kept for now to preserve constructor
                         ABI. */

public:
  /**
  @brief Default constructor initializes members to zero.

  @note NOT USED.
  */
  Isite();

  /**
  @brief Overloaded constructor sets the members to prescribed values at
  construction time.

  @param idval is the id of the corresponding Site object
  @param aftyval is the prescribed value of the affinity to set
  @param verbose is the flag for verbose output
  */
  Isite(int idval, int aftyvalInt, const bool verbose = false, double prb = 0);

  /**
  @brief Copy constructor.

  @param The Isite object to copy
  */
  Isite(const Isite &obj);

  Isite &operator=(const Isite &obj); // copy assignment operator

  /**
  @brief Get the id number of the this Isite object
  (the id_ of the corresponding Site object).

  @return the id number of this Isite object
  (the id_ of the corresponding Site object).
  */
  int getId(void) const { return id_; }

  /**
  @brief Set the id number of this Isite object
  (the id_ of the corresponding Site object).

  @todo Maybe the argument should be declared const

  @param idval is the id number of this Isite object
  (the id_ of the corresponding Site object).
  */
  void setId(int idval) { id_ = idval; }

  /**
  @brief Get the growth affinity of this Isite object.

  @return the growth affinity of this Isite object
  */
  int getAffinityInt(void) const { return affinityInt_; }

  /**
  @brief Set the growth affinity of this Isite object.

  @note NOT USED.

  @param num is the growth affinity of this Isite object
  */
  void setAffinityInt(int num) { affinityInt_ = num; }

  /**
  @brief Update the growth affinity of this Isite object.

  @param afty is the value that must be added to the already growth
  affinity value of this Isite object
  */
  void updateAffinityInt(int afty) { affinityInt_ += afty; }

  /**
  @brief Set the verbose flag

  @param isverbose is true if verbose output should be produced
  */
  void setVerbose(const bool isverbose) { verbose_ = isverbose; }

  /**
  @brief Get the verbose flag

  @return the verbose flag
  */
  bool getVerbose(void) const { return verbose_; }

}; // End of the Isite class
#endif // SRC_THAMESLIB_ISITE_H_
