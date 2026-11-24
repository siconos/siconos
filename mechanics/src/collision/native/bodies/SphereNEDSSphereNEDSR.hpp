// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file SphereNEDSSphereNEDSR.hpp
  \brief Two spheres relation with the Newton Euler formalism and quaternions.
*/

#ifndef SphereNEDSSphereNEDSR_h
#define SphereNEDSSphereNEDSR_h

#include "NewtonEulerVelocityR.hpp"

namespace siconos::collision::native::bodies {
class SphereNEDSSphereNEDSR : public siconos::modeling::NewtonEulerVelocityR,
                              public std::enable_shared_from_this<SphereNEDSSphereNEDSR> {
 private:
  ACCEPT_SERIALIZATION(SphereNEDSSphereNEDSR);

  double radius1_{0.};
  double radius2_{0.};
  double radius_sum{0.};

 public:
  ~SphereNEDSSphereNEDSR() noexcept = default;

  /** Constructor
   *
   *  \param r1 disk1 radius
   *  \param r2 disk2 radius
   */
  SphereNEDSSphereNEDSR(double r1, double r2);

  double distance(double, double, double, double, double, double, double, double);

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

  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::collision::native::bodies
#endif /* SphereNEDSSphereNEDSR_h */
