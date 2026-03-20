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

#include "BlockVector.hpp"
#include "FENode.hpp"
#include "Interaction.hpp"
#include "SiconosException.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::mechanics::fem::NodeFem2d3DR::initialize(siconos::modeling::Interaction &inter) {
  auto sizeDS = inter.getSizeOfDS();

  if (!jacobianhOver_q_internal_storage_) {
    jacobianhOver_q_internal_storage_ = std::make_unique<std::vector<double>>(2 * sizeDS);
  }

  if ((sizeDS != 3) and (sizeDS != 6)) {
    THROW_EXCEPTION(
        "NodeFem2d3DR::initialize(Interaction& inter). The size of ds must of size 3");
  }
  jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianhOver_q_internal_storage_->data(), 2, sizeDS);
  jacobianhOver_q_view_->setZero();
}

void siconos::mechanics::fem::NodeFem2d3DR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector &q) {
  DEBUG_BEGIN("NodeFem2d3DR::computeJachq(const BlockVector& q, BlockVector& z \n");

  jacobianhOver_q_view_->setValue(0, _node_index * 3, nc_.x());
  jacobianhOver_q_view_->setValue(0, _node_index * 3 + 1, nc_.y());
  jacobianhOver_q_view_->setValue(0, _node_index * 3 + 2, nc_.z());

  jacobianhOver_q_view_->setValue(1, _node_index * 3, tangent_.x());
  jacobianhOver_q_view_->setValue(1, _node_index * 3 + 1, tangent_.y());
  jacobianhOver_q_view_->setValue(1, _node_index * 3 + 2, tangent_.z());

  if (q.size() == 6) {
    DEBUG_PRINT("take into account second ds\n");
    THROW_EXCEPTION("NodeFem2d3DR is not implemented for cable/cable contact");
  }
  DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_view_););
  DEBUG_END("NodeFem2d3DR::computeJachq(const BlockVector& q, BlockVector& z) \n");
}

double siconos::mechanics::fem::NodeFem2d3DR::distance() const {
  DEBUG_BEGIN("NodeFem2d3DR::distance(...)\n")
  siconos::algebra::SiconosVector dpc = contactPoint2_ - contactPoint1_;
  DEBUG_EXPR(siconos::algebra::print(contactPoint1_););
  DEBUG_EXPR(siconos::algebra::print(contactPoint2_););
  DEBUG_EXPR(siconos::algebra::print(dpc););
  DEBUG_END("NodeFem2d3DR::distance(...)\n")
  return dpc.norm() * (nc_.dot(dpc) >= 0 ? -1 : 1);
}

void siconos::mechanics::fem::NodeFem2d3DR::computeh(
    const siconos::algebra::BlockVector &q, Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_BEGIN("NodeFem2d3DR::computeh(...)\n");

  siconos::modeling::LagrangianScleronomousR::computeh(q, y);
  y(0) = distance();
  DEBUG_EXPR(siconos::algebra::print(y););
  DEBUG_EXPR(display(););
  DEBUG_END("NodeFem2d3DR::computeh(...)\n")
}

void siconos::mechanics::fem::NodeFem2d3DR::display() const {
  LagrangianR::display();

  std::cout << " _node_index :" << _node_index << std::endl;

  std::cout << " _Pc1 :" << std::endl;
  siconos::algebra::print(contactPoint1_);

  std::cout << " _Pc2 :" << std::endl;
  siconos::algebra::print(contactPoint2_);

  std::cout << " _Normal :" << std::endl;
  siconos::algebra::print(nc_);

  std::cout << " _Tangent :" << std::endl;
  siconos::algebra::print(tangent_);
}
