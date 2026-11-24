// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file SphereNEDSPlanR.hpp

  \brief Newton Euler sphere relation with a plan Ax+By+Cz+D=0

*/

#ifndef SphereNEDSPlanR_h
#define SphereNEDSPlanR_h

#include "NewtonEulerVelocityR.hpp"
#include "SiconosVector.hpp"

namespace siconos::collision::native::bodies {
class SphereNEDSPlanR : public siconos::modeling::NewtonEulerVelocityR,
                        public std::enable_shared_from_this<SphereNEDSPlanR> {
 private:
  ACCEPT_SERIALIZATION(SphereNEDSPlanR);

  double sphere_radius_{0.};
  /* ax + by + cz + d = 0 */
  siconos::algebra::SiconosVector3 plane_normal_ = siconos::algebra::SiconosVector3::Zero();
  double plane_d_{0};

 public:
  /** Constructor
   *
   *  \param r disk radius
   *  \param A
   *  \param B
   *  \param C
   *  \param D
   */
  SphereNEDSPlanR(double r, double A, double B, double C, double D);

  ~SphereNEDSPlanR() noexcept = default;

  double distance(double, double, double, double);

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

  bool equal(double A, double B, double C, double D, double r) const {
    return (A == plane_normal_.x() && B == plane_normal_.y() && C == plane_normal_.z() &&
            D == plane_d_ && r == sphere_radius_);
  }
  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::collision::native::bodies
#endif /* SphereNEDSPlanR_h */
