/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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


#include "Lagrangian2d3DR.hpp"
#include "Interaction.hpp"
#include "BlockVector.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"



void Lagrangian2d3DR::initialize(Interaction& inter)
{
  LagrangianR::initialize(inter);
  //proj_with_q  _jachqProj.reset(new SimpleMatrix(_jachq->size(0),_jachq->size(1)));

  if((inter.getSizeOfDS() !=3) and (inter.getSizeOfDS() !=6))
  {
    THROW_EXCEPTION("Lagrangian2d3DR::initialize(Interaction& inter). The size of ds must of size 3");
  }
  unsigned int qSize = 3 * (inter.getSizeOfDS() / 3);
  _jachq.reset(new SimpleMatrix(3, qSize));
}

void Lagrangian2d3DR::computeJachq(const BlockVector& q, BlockVector& z)
{
  DEBUG_BEGIN("Lagrangian2d3DR::computeJachq(Interaction& inter, SP::BlockVector q0 \n");

  double Nx = _Nc->getValue(0);
  double Ny = _Nc->getValue(1);
  double Px = _Pc1->getValue(0);
  double Py = _Pc1->getValue(1);
  double G1x = q.getValue(0);
  double G1y = q.getValue(1);


  /* construct tangent vector */
  double Tx = Ny;
  double Ty = - Nx;


  double lever_arm_x = G1x-Px ;
  double lever_arm_y = G1y-Py ;
  DEBUG_PRINTF("lever_arm_x = %4.2e,\t lever_arm_ y = %4.2e\n", lever_arm_x, lever_arm_y);

  _jachq->setValue(0,0,Nx);
  _jachq->setValue(0,1,Ny);
  _jachq->setValue(0,2,lever_arm_y*Nx - lever_arm_x*Ny);

  _jachq->setValue(1,0,Tx);
  _jachq->setValue(1,1,Ty);
  _jachq->setValue(1,2,lever_arm_y*Tx - lever_arm_x*Ty);

  _jachq->setValue(2,0,0.0);
  _jachq->setValue(2,1,0.0);
  _jachq->setValue(2,2,1.0);


  if(q.size() ==6)
  {
    DEBUG_PRINT("take into account second ds\n");
    double G2x = q.getValue(3);
    double G2y = q.getValue(4);

    _jachq->setValue(0,3,-Nx);
    _jachq->setValue(0,4,-Ny);
    _jachq->setValue(0,5,- ((G2y-Py)*Nx - (G2x-Px)*Ny));

    _jachq->setValue(1,3,-Tx);
    _jachq->setValue(1,4,-Ty);
    _jachq->setValue(1,5,-((G2y-Py)*Tx - (G2x-Px)*Ty));

    _jachq->setValue(2,3,0.0);
    _jachq->setValue(2,4,0.0);
    _jachq->setValue(2,5,-1.0);

  }
  DEBUG_EXPR(_jachq->display(););
  DEBUG_END("Lagrangian2d3DR::computeJachq(Interaction& inter, SP::BlockVector q0) \n");

}

double Lagrangian2d3DR::distance() const
{
  DEBUG_BEGIN("Lagrangian2d3DR::distance(...)\n")
  SiconosVector dpc(*_Pc2 - *_Pc1);
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_EXPR(_Pc2->display(););
  DEBUG_EXPR(dpc.display(););
  DEBUG_END("Lagrangian2d3DR::distance(...)\n")
  return dpc.norm2() * (inner_prod(*_Nc, dpc) >= 0 ? -1 : 1);

}

void Lagrangian2d3DR::computeh(const BlockVector& q, BlockVector& z, SiconosVector& y)
{
  DEBUG_BEGIN("Lagrangian2d3DR::computeh(...)\n");
  DEBUG_EXPR(q.display());
  // Contact points and normal are stored as relative to q1 and q2, if
  // no q2 then pc2 and normal are absolute.

  // Update pc1 based on q and relPc1

  double angle= q(2);
  DEBUG_PRINTF("angle (ds1)= %e\n", angle);
  (*_Pc1)(0) = q(0) + cos(angle) * (*_relPc1)(0)- sin(angle) * (*_relPc1)(1);
  (*_Pc1)(1) = q(1) + sin(angle) * (*_relPc1)(0)+ cos(angle) * (*_relPc1)(1);
  if(q.size() == 6)
  {
    // To be checked
    DEBUG_PRINT("take into account second ds\n");
    angle = q(5);
    DEBUG_PRINTF("angle (ds2) = %e\n", angle);
    (*_Pc2)(0) = q(3) + cos(angle) * (*_relPc2)(0)- sin(angle) * (*_relPc2)(1);
    (*_Pc2)(1) = q(4) + sin(angle) * (*_relPc2)(0)+ cos(angle) * (*_relPc2)(1);
    (*_Nc)(0) =  cos(angle) * (*_relNc)(0)- sin(angle) * (*_relNc)(1);
    (*_Nc)(1) =  sin(angle) * (*_relNc)(0)+ cos(angle) * (*_relNc)(1);
  }
  else
  {
    *_Pc2 = *_relPc2;
    *_Nc = *_relNc;
  }
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_EXPR(_Pc2->display(););
  DEBUG_EXPR(_Nc->display(););

  LagrangianScleronomousR::computeh(q, z, y);
  y.setValue(0, distance());
  DEBUG_EXPR(y.display(););
  DEBUG_EXPR(display(););
  DEBUG_END("Lagrangian2d3DR::computeh(...)\n")
}
void Lagrangian2d3DR::display() const
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
