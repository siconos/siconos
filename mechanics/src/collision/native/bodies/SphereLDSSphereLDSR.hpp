// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file SphereLDSSphereLDSR.hpp
  \brief Two spheres relation - Inherits from LagrangianScleronomousR
*/

#ifndef SphereLDSSphereLDSR_h
#define SphereLDSSphereLDSR_h

#include "LagrangianScleronomousR.hpp"

namespace siconos::collision::native::bodies {

class SphereLDSSphereLDSR : public siconos::modeling::LagrangianScleronomousR,
                            public std::enable_shared_from_this<SphereLDSSphereLDSR> {
 private:
  ACCEPT_SERIALIZATION(SphereLDSSphereLDSR);

  double r1, r2, r1pr2;

 public:
  /** Constructor
   *
   *  \param r1 disk1 radius
   *  \param r2 disk2 radius
   */
  SphereLDSSphereLDSR(double r1, double r2);

  double distance(double, double, double, double, double, double, double, double);

  /**
    to compute the output y = h(q) of the Relation

    \param q coordinates of the dynamical systems involved in the relation
    \param y the resulting vector
  */
  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /** Computes \f$ \nabla^\top_q h(q) \f$
   * \param q coordinates of the dynamical systems involved in the relation
   */
  void computeContactMatrix(const siconos::algebra::BlockVector& q) override;

  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::collision::native::bodies
#endif /* SphereLDSSphereLDSR_h */
