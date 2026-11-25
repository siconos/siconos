// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SphereNEDSSphereNEDSR.hpp"

#include "SiconosVector.hpp"

siconos::collision::native::bodies::SphereNEDSSphereNEDSR::SphereNEDSSphereNEDSR::
    SphereNEDSSphereNEDSR(double r, double rr)
    : siconos::modeling::NewtonEulerVelocityR{},
      radius1_{r},
      radius2_{rr},
      radius_sum{r + rr} {}

double siconos::collision::native::bodies::SphereNEDSSphereNEDSR::distance(
    double x1, double y1, double z1, double r1, double x2, double y2, double z2, double r2) {
  double dx = x1 - x2;
  double dy = y1 - y2;
  double dz = z1 - z2;

  return (sqrt(dx * dx + dy * dy + dz * dz) - radius_sum);
}

void siconos::collision::native::bodies::SphereNEDSSphereNEDSR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  assert(q2);
  y(0) = distance(q1(0), q1(1), q1(2), radius1_, (*q2)(0), (*q2)(1), (*q2)(2), radius2_);
  // Approximation _Pc1=_Pc2
  contactPoint1_ = (radius1_ * q1.head<3>() + radius2_ * q2->head<3>()) / radius_sum;
  contactPoint2_ = contactPoint1_;
  nc_ = (q1.head<3>() - q2->head<3>()) / (y(0) + radius_sum);
}
