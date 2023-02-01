/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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


#include "Cable2d3DR.hpp"
#include "Interaction.hpp"
#include "BlockVector.hpp"
#include "SiconosException.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"



void siconos::mechanics::fem::Cable2d3DR::initialize(Interaction& inter)
{
  unsigned int qSize = inter.getSizeOfDS();
  _jachq.reset(new SimpleMatrix(2, qSize));
}

void siconos::mechanics::fem::Cable2d3DR::computeJachq(const BlockVector& q, BlockVector& z)
{
  DEBUG_BEGIN("Cable2d3DR::computeJachq(const BlockVector& q, BlockVector& z \n");

  double Nx = _Normal->getValue(0);
  double Ny = _Normal->getValue(1);
  double Nz = _Normal->getValue(2);

  double Tx = _Tangent->getValue(0);
  double Ty = _Tangent->getValue(1);
  double Tz = _Tangent->getValue(2);
  
  DEBUG_PRINTF("N_x = %4.2e,\t N_y = %4.2e,\t N_z= %4.2e\n", Nx, Ny, Nz);
  DEBUG_PRINTF("T_x = %4.2e,\t T_y = %4.2e,\t T_z= %4.2e\n", Tx, Ty, Tz);

  _jachq->setValue(0,_node_dof_index,  Nx);
  _jachq->setValue(0,_node_dof_index+1,Ny);
  _jachq->setValue(0,_node_dof_index+2,Nz);
  
  _jachq->setValue(1,_node_dof_index,  Tx);
  _jachq->setValue(1,_node_dof_index+1,Ty);
  _jachq->setValue(1,_node_dof_index+2,Tz);

  if(q.size() ==6)
  {
    DEBUG_PRINT("take into account second ds\n");
    THROW_EXCEPTION("Cable2d3DR is not implemented for cable/cable contact");
  }
  DEBUG_EXPR(_jachq->display(););
  DEBUG_END("Cable2d3DR::computeJachq(const BlockVector& q, BlockVector& z) \n");

}

double siconos::mechanics::fem::Cable2d3DR::distance() const
{
  DEBUG_BEGIN("Cable2d3DR::distance(...)\n")
  SiconosVector dpc(*_Pc2 - *_Pc1);
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_EXPR(_Pc2->display(););
  DEBUG_EXPR(dpc.display(););
  DEBUG_END("Cable2d3DR::distance(...)\n")
  return dpc.norm2() * (inner_prod(*_Normal, dpc) >= 0 ? -1 : 1);

}

void siconos::mechanics::fem::Cable2d3DR::computeh(const BlockVector& q, BlockVector& z, SiconosVector& y)
{
  DEBUG_BEGIN("Cable2d3DR::computeh(...)\n");

  // LagrangianScleronomousR::computeh(q, z, y);
  SiconosVector & position = *((q.getAllVect())[0]);
  _Pc1->setValue(0, position(_node_dof_index));
  _Pc1->setValue(1, position(_node_dof_index+1));
  _Pc1->setValue(2, position(_node_dof_index+2));
  y.setValue(0, distance());
  DEBUG_EXPR(y.display(););
  DEBUG_EXPR(display(););
  DEBUG_END("Cable2d3DR::computeh(...)\n")
}

void siconos::mechanics::fem::Cable2d3DR::display() const
{
  LagrangianR::display();

  std::cout << " _node_dof_index :" << _node_dof_index<< std::endl;
  
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

  std::cout << " _Normal :" << std::endl;
  if(_Normal)
    _Normal->display();
  else
    std::cout << " nullptr :" << std::endl;
  
  std::cout << " _Tangent :" << std::endl;
  if(_Tangent)
    _Tangent->display();
  else
    std::cout << " nullptr :" << std::endl;

}
