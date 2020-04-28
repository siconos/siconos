/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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


#include "Lagrangian2d1DR.hpp"
#include <boost/math/quaternion.hpp>
#include "NewtonEulerDS.hpp"
#include "Interaction.hpp"
#include "BlockVector.hpp"
#include <boost/math/quaternion.hpp>

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"


void Lagrangian2d1DR::initialize(Interaction& inter)
{
  LagrangianR::initialize(inter);
  //proj_with_q  _jachqProj.reset(new SimpleMatrix(_jachq->size(0),_jachq->size(1)));

  if((inter.getSizeOfDS() !=3) and (inter.getSizeOfDS() !=6))
  {
    RuntimeException::selfThrow("Lagrangian2d1DR::initialize(Interaction& inter). The size of ds must of size 3");
  }
  unsigned int qSize = 3 * (inter.getSizeOfDS() / 3);
  _jachq.reset(new SimpleMatrix(1, qSize));
}

void Lagrangian2d1DR::computeJachq(const BlockVector& q, BlockVector& z)
{
  DEBUG_BEGIN("Lagrangian2d1DR::computeJachq(Interaction& inter, SP::BlockVector q0 \n");

  double Nx = _Nc->getValue(0);
  double Ny = _Nc->getValue(1);
  double Px = _Pc1->getValue(0);
  double Py = _Pc1->getValue(1);
  double G1x = q.getValue(0);
  double G1y = q.getValue(1);

  _jachq->setValue(0,0,Nx);
  _jachq->setValue(0,1,Ny);
  _jachq->setValue(0,2,(G1y-Py)*Nx - (G1x-Px)*Ny);


  if(q.size() ==6)
  {
    DEBUG_PRINT("take into account second ds\n");
    double G2x = q.getValue(3);
    double G2y = q.getValue(4);

    _jachq->setValue(0,3,-Nx);
    _jachq->setValue(0,4,-Ny);
    _jachq->setValue(0,5,- ((G2y-Py)*Nx - (G2x-Px)*Ny));
  }
  DEBUG_EXPR(_jachq->display(););
  DEBUG_END("Lagrangian2d1DR::computeJachq(Interaction& inter, SP::BlockVector q0) \n");

}

double Lagrangian2d1DR::distance() const
{
  DEBUG_BEGIN("Lagrangian2d1DR::distance(...)\n")
  SiconosVector dpc(*_Pc2 - *_Pc1);
  DEBUG_END("Lagrangian2d1DR::distance(...)\n")
  return dpc.norm2() * (inner_prod(*_Nc, dpc) >= 0 ? -1 : 1);

}

void Lagrangian2d1DR::computeh(const BlockVector& q, BlockVector& z, SiconosVector& y)
{
  DEBUG_BEGIN("Lagrangian2d1DR::computeh(...)\n");
  DEBUG_EXPR(q.display());
  // Contact points and normal are stored as relative to q1 and q2, if
  // no q2 then pc2 and normal are absolute.

  // Update pc1 based on q and relPc1

  double angle= q(2);
  (*_Pc1)(0) = q(0) + cos(angle) * (*_relPc1)(0)+ sin(angle) * (*_relPc1)(1);
  (*_Pc1)(1) = q(1) + sin(angle) * (*_relPc1)(0)- cos(angle) * (*_relPc1)(1);
  if(q.size() == 6)
  {
    // To be checked
    DEBUG_PRINT("take into account second ds\n");
    angle = q(5);
    (*_Pc2)(0) = q(3) + cos(angle) * (*_relPc2)(0)+ sin(angle) * (*_relPc2)(1);
    (*_Pc2)(1) = q(4) + sin(angle) * (*_relPc2)(0)- cos(angle) * (*_relPc2)(1);
    (*_Nc)(0) =  cos(angle) * (*_relNc)(0)+ sin(angle) * (*_relNc)(1);
    (*_Nc)(1) =  sin(angle) * (*_relNc)(0)+ cos(angle) * (*_relNc)(1);
  }
  else
  {
    *_Pc2 = *_relPc2;
    *_Nc = *_relNc;
  }


  // SP::SiconosVector q1 = (q0.getAllVect())[0];
  // ::boost::math::quaternion<double> qq1((*q1)(3), (*q1)(4), (*q1)(5), (*q1)(6));
  // ::boost::math::quaternion<double> qpc1(0,(*_relPc1)(0),(*_relPc1)(1),(*_relPc1)(2));

  // // apply q1 rotation and add
  // qpc1 = qq1 * qpc1 / qq1;
  // (*_Pc1)(0) = qpc1.R_component_2() + (*q1)(0);
  // (*_Pc1)(1) = qpc1.R_component_3() + (*q1)(1);
  // (*_Pc1)(2) = qpc1.R_component_4() + (*q1)(2);

  // if (q0.numberOfBlocks() > 1)
  // {
  //   // Update pc2 based on q0 and relPc2
  //   SP::SiconosVector q2 = (q0.getAllVect())[1];
  //   ::boost::math::quaternion<double> qq2((*q2)(3), (*q2)(4), (*q2)(5), (*q2)(6));
  //   ::boost::math::quaternion<double> qpc2(0,(*_relPc2)(0),(*_relPc2)(1),(*_relPc2)(2));

  //   // apply q2 rotation and add
  //   qpc2 = qq2 * qpc2 / qq2;
  //   (*_Pc2)(0) = qpc2.R_component_2() + (*q2)(0);
  //   (*_Pc2)(1) = qpc2.R_component_3() + (*q2)(1);
  //   (*_Pc2)(2) = qpc2.R_component_4() + (*q2)(2);

  //   // same for normal
  //   ::boost::math::quaternion<double> qnc(0, (*_relNc)(0), (*_relNc)(1), (*_relNc)(2));
  //   qnc = qq2 * qnc / qq2;
  //   (*_Nc)(0) = qnc.R_component_2();
  //   (*_Nc)(1) = qnc.R_component_3();
  //   (*_Nc)(2) = qnc.R_component_4();
  // }
  // else
  // {
  //   *_Pc2 = *_relPc2;
  //   *_Nc = *_relNc;
  // }

  LagrangianScleronomousR::computeh(q, z, y);
  y.setValue(0, distance());
  DEBUG_EXPR(y.display(););
  DEBUG_EXPR(display(););
  DEBUG_END("Lagrangian2d1DR::computeh(...)\n")
}
void Lagrangian2d1DR::display() const
{
  LagrangianR::display();

  std::cout << " _Pc1 :" << std::endl;
  if(_Pc1)
    _Pc1->display();
  else
    std::cout << " nullptr :" << std::endl;

  std::cout << " _Pc2 :" << std::endl;
  if(_Pc2)
    _Pc2->display();
  else
    std::cout << " nullptr :" << std::endl;

  std::cout << " _relPc1 :" << std::endl;
  if(_relPc1)
    _relPc1->display();
  else
    std::cout << " nullptr :" << std::endl;

  std::cout << " _relPc2 :" << std::endl;
  if(_relPc2)
    _relPc2->display();
  else
    std::cout << " nullptr :" << std::endl;

  std::cout << " _Nc :" << std::endl;
  if(_Nc)
    _Nc->display();
  else
    std::cout << " nullptr :" << std::endl;
  std::cout << " _relNc :" << std::endl;
  if(_relNc)
    _relNc->display();
  else
    std::cout << " nullptr :" << std::endl;

}
