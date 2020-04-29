/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include <cmath>
#include "SphereNEDSPlanR.hpp"
#include <Interaction.hpp>
#include <BlockVector.hpp>

SphereNEDSPlanR::SphereNEDSPlanR(double r, double A, double B, double C, double D)
  : NewtonEuler3DR(), r(r), A(A), B(B), C(C), D(D)
{
  nN = sqrt(A * A + B * B + C * C);

  n1 = A / nN;
  n2 = B / nN;
  n3 = C / nN;



}

double SphereNEDSPlanR::distance(double x, double y, double z, double rad)
{

  return (fabs(A * x + B * y + C * z + D) / nN - rad);
}


void SphereNEDSPlanR::computeh(double time, const BlockVector& q0, SiconosVector& y)
{

  double q_0 = q0(0);
  double q_1 = q0(1);
  double q_2 = q0(2);

  y.setValue(0, distance(q_0, q_1, q_2, r));
  _Pc1->setValue(0, q_0 - r * n1);
  _Pc1->setValue(1, q_1 - r * n2);
  _Pc1->setValue(2, q_2 - r * n3);
  _Pc2->setValue(0, q_0 - r * n1);
  _Pc2->setValue(1, q_1 - r * n2);
  _Pc2->setValue(2, q_2 - r * n3);
  _Nc->setValue(0, n1);
  _Nc->setValue(1, n2);
  _Nc->setValue(2, n3);
}
