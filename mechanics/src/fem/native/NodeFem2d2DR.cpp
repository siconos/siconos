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

#include "NodeFem2d2DR.hpp"
#include "Interaction.hpp"
#include "BlockVector.hpp"
#include "SiconosException.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::mechanics::fem::NodeFem2d2DR::initialize(modeling::Interaction& inter)
{
  unsigned int qSize = inter.getSizeOfDS();
  _jachq.reset(new siconos::algebra::SimpleMatrix(2, qSize));
}

void siconos::mechanics::fem::NodeFem2d2DR::computeJachq(const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z)
{
  DEBUG_BEGIN("NodeFem2d2DR::computeJachq(const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z \n");

  double Nx = _Normal->getValue(0);
  double Ny = _Normal->getValue(1);

  double Tx = _Tangent->getValue(0);
  double Ty = _Tangent->getValue(1);
  
  // double Px = _Pc1->getValue(0);
  // double Py = _Pc1->getValue(1);

  DEBUG_PRINTF("N_x = %4.2e,\t N_y = %4.2e\n", Nx, Ny);
  DEBUG_PRINTF("T_x = %4.2e,\t T_y = %4.2e\n", Tx, Ty);
 
  _jachq->setValue(0,(*_node->dofIndex())[0],Nx);
  _jachq->setValue(0,(*_node->dofIndex())[1],Ny);
  
  _jachq->setValue(1,(*_node->dofIndex())[0],Tx);
  _jachq->setValue(1,(*_node->dofIndex())[1],Ty);

  if(q.size() ==6)
  {
    DEBUG_PRINT("take into account second ds\n");
    THROW_EXCEPTION("NodeFem2d2DR is not implemented for cable/cable contact");
  }
  //DEBUG_EXPR(_jachq->display(););
  DEBUG_END("NodeFem2d2DR::computeJachq(const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z) \n");

}

double siconos::mechanics::fem::NodeFem2d2DR::distance() const
{
  DEBUG_BEGIN("NodeFem2d2DR::distance(...)\n")
  siconos::algebra::SiconosVector dpc(*_Pc2 - *_Pc1);
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_EXPR(_Pc2->display(););
  DEBUG_EXPR(dpc.display(););
  DEBUG_END("NodeFem2d2DR::distance(...)\n")
  return dpc.norm2() * (inner_prod(*_Normal, dpc) >= 0 ? -1 : 1);

}

void siconos::mechanics::fem::NodeFem2d2DR::computeh(const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z, siconos::algebra::SiconosVector& y)
{
  DEBUG_BEGIN("NodeFem2d2DR::computeh(...)\n");

  LagrangianScleronomousR::computeh(q, z, y);
  siconos::algebra::SiconosVector & displacement = *((q.getAllVect())[0]);
  _Pc1->setValue(0, displacement((*_node->dofIndex())[0])+_node->x());
  _Pc1->setValue(1, displacement((*_node->dofIndex())[1])+_node->y());
  y.setValue(0, distance());
  DEBUG_PRINTF("distance = %e\n",distance());
  DEBUG_EXPR(y.display(););
  DEBUG_EXPR(display(););
  DEBUG_END("NodeFem2d2DR::computeh(...)\n")
}

void siconos::mechanics::fem::NodeFem2d2DR::display() const
{
  LagrangianR::display();


  std::cout << " _node :" << std::endl;
  if(_node)
    _node->display();
  else
    std::cout << " nullptr :" << std::endl;
  
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
