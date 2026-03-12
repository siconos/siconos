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

#include "SphereNEDSPlanR.hpp"

#include <cmath>

#include "BlockVector.hpp"
#include "SiconosVector.hpp"

siconos::collision::native::bodies::SphereNEDSPlanR::SphereNEDSPlanR::SphereNEDSPlanR(
    double r, double A, double B, double C, double D)
    : siconos::modeling::NewtonEuler3DR{}, r{r}, A{A}, B{B}, C{C}, D{D} {
  nN = sqrt(A * A + B * B + C * C);

  n1 = A / nN;
  n2 = B / nN;
  n3 = C / nN;
}

double siconos::collision::native::bodies::SphereNEDSPlanR::SphereNEDSPlanR::distance(
    double x, double y, double z, double rad) {
  return (fabs(A * x + B * y + C * z + D) / nN - rad);
}

void siconos::collision::native::bodies::SphereNEDSPlanR::SphereNEDSPlanR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  double q_0 = q1(0);
  double q_1 = q1(1);
  double q_2 = q1(2);

  y(0) = distance(q_0, q_1, q_2, r);
  (*_Pc1)(0) = q_0 - r * n1;
  (*_Pc1)(1) = q_1 - r * n2;
  (*_Pc1)(2) = q_2 - r * n3;
  (*_Pc2)(0) = q_0 - r * n1;
  (*_Pc2)(1) = q_1 - r * n2;
  (*_Pc2)(2) = q_2 - r * n3;
  (*_Nc)(0) = n1;
  (*_Nc)(1) = n2;
  (*_Nc)(2) = n3;
}
