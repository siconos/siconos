/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "SphereNEDSSphereNEDSR.hpp"
#include <Interaction.hpp>
#include <BlockVector.hpp>

SphereNEDSSphereNEDSR::SphereNEDSSphereNEDSR(double r,
    double rr)
  : NewtonEulerFrom3DLocalFrameR()
{
  r1 = r;
  r2 = rr;
  r1pr2 = r1 + r2;
}

double SphereNEDSSphereNEDSR::distance(double x1, double y1, double z1, double r1,
                                       double x2, double y2, double z2, double r2)
{
  double dx = x1 - x2;
  double dy = y1 - y2;
  double dz = z1 - z2;

  return (sqrt(dx * dx + dy * dy + dz * dz) - r1pr2);
}


void SphereNEDSSphereNEDSR::computeh(double time, BlockVector& q0, SiconosVector& y)
{


  double q_0 = q0(0);
  double q_1 = q0(1);
  double q_2 = q0(2);
  double q_7 = q0(7);
  double q_8 = q0(8);
  double q_9 = q0(9);

  y.setValue(0, distance(q_0, q_1, q_2, r1, q_7, q_8, q_9, r2));
  //Approximation _Pc1=_Pc2
  _Pc1->setValue(0, (r1 * q_0 + r2 * q_7) / (r1 + r2));
  _Pc1->setValue(1, (r1 * q_1 + r2 * q_8) / (r1 + r2));
  _Pc1->setValue(2, (r1 * q_2 + r2 * q_9) / (r1 + r2));
  _Pc2->setValue(0, (r1 * q_0 + r2 * q_7) / (r1 + r2));
  _Pc2->setValue(1, (r1 * q_1 + r2 * q_8) / (r1 + r2));
  _Pc2->setValue(2, (r1 * q_2 + r2 * q_9) / (r1 + r2));
  _Nc->setValue(0, (q_0 - q_7) / (y.getValue(0) + r1pr2));
  _Nc->setValue(1, (q_1 - q_8) / (y.getValue(0) + r1pr2));
  _Nc->setValue(2, (q_2 - q_9) / (y.getValue(0) + r1pr2));
  //std::cout<<" SphereNEDSSphereNEDSR::computeh dist:"<<y->getValue(0)<<"\n";
  //std::cout<<"_Pc:\n";
  //_Pc->display();
  //std::cout<<"_Nc:\n";
  //_Nc->display();
}
