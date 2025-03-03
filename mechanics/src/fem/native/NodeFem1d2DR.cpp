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

#include "NodeFem1d2DR.hpp"

#include "BlockVector.hpp"
#include "FENode.hpp"
#include "Interaction.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::mechanics::fem::NodeFem1d2DR::initialize(modeling::Interaction& inter) {
  auto sizeDS = inter.getSizeOfDS();

  if (!jacobianhOver_q_internal_storage_) {
    jacobianhOver_q_internal_storage_ = std::make_unique<std::vector<double>>(1 * sizeDS);
  }
  jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianhOver_q_internal_storage_->data(), 1, sizeDS);
}

void siconos::mechanics::fem::NodeFem1d2DR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& q) {
  DEBUG_BEGIN(
      "NodeFem1d2DR::computeJacobianhOver_q(const siconos::algebra::BlockVector& q, "
      "siconos::algebra::BlockVector& z \n");

  double Nx = (*_Normal)(0);
  double Ny = (*_Normal)(1);

  DEBUG_PRINTF("N_x = %4.2e,\t N_y = %4.2e\n", Nx, Ny);

  jacobianhOver_q_view_->setValue(0, (*_node->dofIndex())[0], Nx);
  jacobianhOver_q_view_->setValue(0, (*_node->dofIndex())[1], Ny);

  if (q.size() == 6) {
    DEBUG_PRINT("take into account second ds\n");
    THROW_EXCEPTION("NodeFem1d2DR is not implemented for node/node contact");
  }
  DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_););
  DEBUG_END(
      "NodeFem1d2DR::computeJacobianhOver_q(const siconos::algebra::BlockVector& q, "
      "siconos::algebra::BlockVector& z) \n");
}

double siconos::mechanics::fem::NodeFem1d2DR::distance() const {
  DEBUG_BEGIN("NodeFem1d2DR::distance(...)\n")
  siconos::algebra::SiconosVector dpc(*_Pc2 - *_Pc1);
  DEBUG_EXPR(siconos::algebra::print(*_Pc1););
  DEBUG_EXPR(siconos::algebra::print(*_Pc2););
  DEBUG_EXPR(siconos::algebra::print(dpc););
  DEBUG_END("NodeFem1d2DR::distance(...)\n")
  return dpc.norm() * (_Normal->dot(dpc) >= 0 ? -1 : 1);
}

void siconos::mechanics::fem::NodeFem1d2DR::computeh(
    const siconos::algebra::BlockVector& q, Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_BEGIN("NodeFem1d2DR::computeh(...)\n");

  LagrangianScleronomousR::computeh(q, y);
  siconos::algebra::SiconosVector& displacement = *((q.getAllVect())[0]);
  (*_Pc1)(0) = displacement((*_node->dofIndex())[0]) + _node->x();
  (*_Pc1)(1) = displacement((*_node->dofIndex())[1]) + _node->y();
  y(0) = distance();

  DEBUG_EXPR(siconos::algebra::print(y););
  DEBUG_EXPR(display(););
  DEBUG_END("NodeFem1d2DR::computeh(...)\n")
}

void siconos::mechanics::fem::NodeFem1d2DR::updateContactPoint(
    siconos::algebra::SiconosVector& pc2, siconos::algebra::SiconosVector& normal) {
  *_Pc2 = pc2;
  *_Normal = normal;
};

void siconos::mechanics::fem::NodeFem1d2DR::updateContactPoint(double pc2[2],
                                                               double normal[2]) {
  (*_Pc2)(0) = pc2[0];
  (*_Pc2)(1) = pc2[1];
  (*_Normal)(0) = normal[0];
  (*_Normal)(1) = normal[1];
};

void siconos::mechanics::fem::NodeFem1d2DR::display() const {
  LagrangianR::display();

  std::cout << " _node :\n";
  if (_node)
    _node->display();
  else
    std::cout << " nullptr :\n";

  std::cout << " _Pc1 :\n";
  if (_Pc1)
    siconos::algebra::print(*_Pc1);
  else
    std::cout << " nullptr :\n";

  std::cout << " _Pc2 :\n";
  if (_Pc2)
    siconos::algebra::print(*_Pc2);
  else
    std::cout << " nullptr :\n";

  std::cout << " _Normal :\n";
  if (_Normal)
    siconos::algebra::print(*_Normal);
  else
    std::cout << " nullptr\n";
}
