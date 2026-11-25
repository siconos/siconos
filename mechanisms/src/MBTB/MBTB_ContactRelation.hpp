// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef MBTB_CONTACTRELATION
#define MBTB_CONTACTRELATION

#include "NewtonEulerNormalVelocityR.hpp"  // Base class

namespace siconos::mechanisms {

class MBTB_Contact;

//! It is a relation dedicated for the simple unilateral (ie: without  Coulomb
//! friction).
/*!
  Aggregation to the class MBTB_Contact, the member _pContact contains the CAD
  information. It derivates from Siconos::NewtonEulerNormalVelocityR. This class does the
  link between CAD and Siconos.
 */
class MBTB_ContactRelation : public siconos::modeling::NewtonEulerNormalVelocityR {
 protected:
  std::shared_ptr<MBTB_Contact> _pContact{nullptr};

 public:
  /** Constructor
   *  \param pC [in] a pointer to the MBTB_Contact. Must be allocated/free by
   * the caller.
   */
  MBTB_ContactRelation(std::shared_ptr<MBTB_Contact> pC);

  /**
     to compute the output y = h(q) of the Relation

      \param[in] q1 generalized coordinates vector of the fist dynamical system involved
      in the relation
      \param[in] q2 generalized coordinates vector of the second dynamical system
      involved in the relation
      \param[in,out] y the resulting vector
  */
  virtual void computeh(
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  virtual ~MBTB_ContactRelation() noexcept = default;
};

}  // namespace siconos::mechanisms
#endif
