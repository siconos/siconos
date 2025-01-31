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

#include "NodeFem2d3DR.hpp"
#include "Interaction.hpp"
#include "BlockVector.hpp"
#include "SiconosException.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"



void siconos::mechanics::fem::NodeFem2d3DR::initialize(siconos::modeling::Interaction& inter)
{
  //proj_with_q  _jachqProj.reset(new SimpleMatrix(_jachq->size(0),_jachq->size(1)));
  auto sizeDS = inter.getSizeOfDS();

  if (!jacobianhOver_q_internal_storage_) {
    jacobianhOver_q_internal_storage_ = std::make_unique<std::vector<double>>(2 * sizeDS);
  }
  
  if((sizeDS!=3) and (sizeDS !=6))
  {
    THROW_EXCEPTION("NodeFem2d3DR::initialize(Interaction& inter). The size of ds must of size 3");
  }
  unsigned int qSize = 3 * (inter.getSizeOfDS() / 3);
  jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianhOver_q_internal_storage_->data(), 2, sizeDS);
}

void siconos::mechanics::fem::NodeFem2d3DR::computeJacobianhOver_q(const siconos::algebra::BlockVector &q)
{
  DEBUG_BEGIN("NodeFem2d3DR::computeJachq(const BlockVector& q, BlockVector& z \n");

  double Nx = _Normal->getValue(0);
  double Ny = _Normal->getValue(1);
  double Nz = _Normal->getValue(2);

  double Tx = _Normal->getValue(0);
  double Ty = _Normal->getValue(1);
  double Tz = _Normal->getValue(2);
  
  double Px = _Pc1->getValue(0);
  double Py = _Pc1->getValue(1);

  DEBUG_PRINTF("N_x = %4.2e,\t N_y = %4.2e,\t N_z= %4.2e\n", Nx, Ny, Nz);
  DEBUG_PRINTF("T_x = %4.2e,\t T_y = %4.2e,\t T_z= %4.2e\n", Tx, Ty, Tz);

  jacobianhOver_q_view_->setValue(0,_node_index*3,  Nx);
  jacobianhOver_q_view_->setValue(0,_node_index*3+1,Ny);
  jacobianhOver_q_view_->setValue(0,_node_index*3+2,Nz);
  
  jacobianhOver_q_view_->setValue(1,_node_index*3,  Tx);
  jacobianhOver_q_view_->setValue(1,_node_index*3+1,Ty);
  jacobianhOver_q_view_->setValue(1,_node_index*3+2,Tz);

  if(q.size() ==6)
  {
    DEBUG_PRINT("take into account second ds\n");
    THROW_EXCEPTION("NodeFem2d3DR is not implemented for cable/cable contact");
  }
  DEBUG_EXPR(jacobianhOver_q_view_->display(););
  DEBUG_END("NodeFem2d3DR::computeJachq(const BlockVector& q, BlockVector& z) \n");

}

double siconos::mechanics::fem::NodeFem2d3DR::distance() const
{
  DEBUG_BEGIN("NodeFem2d3DR::distance(...)\n")
  siconos::algebra::SiconosVector dpc(*_Pc2 - *_Pc1);
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_EXPR(_Pc2->display(););
  DEBUG_EXPR(dpc.display(););
  DEBUG_END("NodeFem2d3DR::distance(...)\n")
    return dpc.norm2() * (_Normal->dot(dpc) >= 0 ? -1 : 1);

}

void siconos::mechanics::fem::NodeFem2d3DR::computeh(const siconos::algebra::BlockVector &q,
                Eigen::Ref<siconos::algebra::SiconosVector> y)
{
  DEBUG_BEGIN("NodeFem2d3DR::computeh(...)\n");

  siconos::modeling::LagrangianScleronomousR::computeh(q, y);
  y.setValue(0, distance());
  DEBUG_EXPR(y.display(););
  DEBUG_EXPR(display(););
  DEBUG_END("NodeFem2d3DR::computeh(...)\n")
}

void siconos::mechanics::fem::NodeFem2d3DR::display() const
{
  LagrangianR::display();

  std::cout << " _node_index :" << _node_index<< std::endl;
  
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
  if(_Normal)
    _Normal->display();
  else
    std::cout << " nullptr :" << std::endl;

}
