/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "SphereNEDSSphereNEDSR.hpp"

#include "BlockVector.hpp"
#include "SiconosVector.hpp"

siconos::collision::native::bodies::SphereNEDSSphereNEDSR::SphereNEDSSphereNEDSR::
    SphereNEDSSphereNEDSR(double r, double rr)
    : siconos::modeling::NewtonEuler3DR{} {
  r1 = r;
  r2 = rr;
  r1pr2 = r1 + r2;
}

double siconos::collision::native::bodies::SphereNEDSSphereNEDSR::distance(
    double x1, double y1, double z1, double r1, double x2, double y2, double z2, double r2) {
  double dx = x1 - x2;
  double dy = y1 - y2;
  double dz = z1 - z2;

  return (sqrt(dx * dx + dy * dy + dz * dz) - r1pr2);
}

void siconos::collision::native::bodies::SphereNEDSSphereNEDSR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  assert(q2);
  double q_0 = q1(0);
  double q_1 = q1(1);
  double q_2 = q1(2);
  double q_7 = 0;
  double q_8 = 0;
  double q_9 = 0;
  if (q2) {
    q_7 = (*q2)(0);
    q_8 = (*q2)(1);
    q_9 = (*q2)(2);
  }
  y(0) = distance(q_0, q_1, q_2, r1, q_7, q_8, q_9, r2);
  // Approximation _Pc1=_Pc2
  (*_Pc1)(0) = (r1 * q_0 + r2 * q_7) / (r1 + r2);
  (*_Pc1)(1) = (r1 * q_1 + r2 * q_8) / (r1 + r2);
  (*_Pc1)(2) = (r1 * q_2 + r2 * q_9) / (r1 + r2);
  (*_Pc2)(0) = (r1 * q_0 + r2 * q_7) / (r1 + r2);
  (*_Pc2)(1) = (r1 * q_1 + r2 * q_8) / (r1 + r2);
  (*_Pc2)(2) = (r1 * q_2 + r2 * q_9) / (r1 + r2);
  (*_Nc)(0) = (q_0 - q_7) / (y(0) + r1pr2);
  (*_Nc)(1) = (q_1 - q_8) / (y(0) + r1pr2);
  (*_Nc)(2) = (q_2 - q_9) / (y(0) + r1pr2);
}
