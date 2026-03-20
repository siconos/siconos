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

#include "NodeFem2d2DR.hpp"

#include "BlockVector.hpp"
#include "FENode.hpp"  // global_dof_index
#include "Interaction.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::mechanics::fem::NodeFem2d2DR::initialize(modeling::Interaction& inter) {
  setComputeJacobianhOver_qFunction([this](const siconos::algebra::BlockVector& q,
                                           Eigen::Ref<siconos::algebra::MapType> result) {
    DEBUG_BEGIN(
        "NodeFem2d2DR::computeJacobianhOver_q(const siconos::algebra::BlockVector& q, "
        "siconos::algebra::BlockVector& z \n");

    result.setValue(0, node_->global_dof_index()[0], nc_.x());
    result.setValue(0, node_->global_dof_index()[1], nc_.y());

    result.setValue(1, node_->global_dof_index()[0], tangent_.x());
    result.setValue(1, node_->global_dof_index()[1], tangent_.y());

    if (q.size() == 6) {
      DEBUG_PRINT("take into account second ds\n");
      THROW_EXCEPTION("NodeFem2d2DR is not implemented for cable/cable contact");
    }
    // DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_););
    DEBUG_END(
        "NodeFem2d2DR::computeJacobianhOver_q(const siconos::algebra::BlockVector& q, "
        "siconos::algebra::BlockVector& z) \n");
  });
  LagrangianScleronomousR::initialize(inter);
}

double siconos::mechanics::fem::NodeFem2d2DR::distance() const {
  DEBUG_BEGIN("NodeFem2d2DR::distance(...)\n")
  siconos::algebra::SiconosVector dpc = contactPoint2_ - contactPoint1_;
  DEBUG_EXPR(siconos::algebra::print(contactPoint1););
  DEBUG_EXPR(siconos::algebra::print(contactPoint2););
  DEBUG_EXPR(siconos::algebra::print(dpc););
  DEBUG_END("NodeFem2d2DR::distance(...)\n")
  return dpc.norm() * (nc_.dot(dpc) >= 0 ? -1 : 1);
}

void siconos::mechanics::fem::NodeFem2d2DR::computeh(
    const siconos::algebra::BlockVector& q, Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_BEGIN("NodeFem2d2DR::computeh(...)\n");

  LagrangianScleronomousR::computeh(q, y);
  auto& displacement = *((q.getAllVect())[0]);
  contactPoint1_(0) = displacement(node_->global_dof_index()[0]) + node_->x();
  contactPoint1_(1) = displacement(node_->global_dof_index()[1]) + node_->y();
  y(0) = distance();
  DEBUG_PRINTF("distance = %e\n", distance());
  DEBUG_EXPR(siconos::algebra::print(y););
  DEBUG_EXPR(display(););
  DEBUG_END("NodeFem2d2DR::computeh(...)\n")
}

void siconos::mechanics::fem::NodeFem2d2DR::updateContactPoint(
    const siconos::algebra::SiconosVector2& pc2,
    const siconos::algebra::SiconosVector2& normal,
    const siconos::algebra::SiconosVector2& tangent) {
  contactPoint2_ = pc2;
  nc_ = normal;
  tangent_ = tangent;
};

/** update the contact points from array
 */
void siconos::mechanics::fem::NodeFem2d2DR::updateContactPoint(double pc2[2], double normal[2],
                                                               double tangent[2]) {
  contactPoint2_(0) = pc2[0];
  contactPoint2_(1) = pc2[1];
  nc_(0) = normal[0];
  nc_(1) = normal[1];
  tangent_(0) = tangent[0];
  tangent_(1) = tangent[1];
};
void siconos::mechanics::fem::NodeFem2d2DR::display() const {
  LagrangianR::display();

  std::cout << " Node:\n";
  if (node_)
    node_->display();
  else
    std::cout << " nullptr :\n";

  std::cout << " _Pc1 :\n";
  siconos::algebra::print(contactPoint1_);

  std::cout << " _Pc2 :\n";
  siconos::algebra::print(contactPoint2_);
  std::cout << " _Normal :\n";
  siconos::algebra::print(nc_);

  std::cout << " _Tangent :\n";
  siconos::algebra::print(tangent_);
}
