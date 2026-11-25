/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "SphereNEDSPlanR.hpp"

#include <cmath>

#include "SiconosVector.hpp"

siconos::collision::native::bodies::SphereNEDSPlanR::SphereNEDSPlanR::SphereNEDSPlanR(
    double r, double A, double B, double C, double D)
    : siconos::modeling::NewtonEulerVelocityR{}, sphere_radius_{r}, plane_d_{D} {
  plane_normal_ << A, B, C;
  double norm = plane_normal_.norm();
  if (norm < 1e-12) {
    THROW_EXCEPTION("SphereNEDSPlanR: normal to the plan has 0 norm");
  } else {
    plane_normal_.normalize();
    plane_d_ /= norm;
  }
}

double siconos::collision::native::bodies::SphereNEDSPlanR::SphereNEDSPlanR::distance(
    double x, double y, double z, double rad) {
  return (
      fabs(plane_normal_.x() * x + plane_normal_.y() * y + plane_normal_.z() * z + plane_d_) -
      rad);
}

void siconos::collision::native::bodies::SphereNEDSPlanR::SphereNEDSPlanR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  assert(!q2);
  double q_0 = q1(0);
  double q_1 = q1(1);
  double q_2 = q1(2);

  y(0) = distance(q_0, q_1, q_2, sphere_radius_);
  contactPoint1_ = q1.head<3>() - sphere_radius_ * plane_normal_;
  contactPoint2_ = contactPoint1_;
  nc_ = plane_normal_;
}
