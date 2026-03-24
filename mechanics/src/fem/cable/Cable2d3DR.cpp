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

#include "Cable2d3DR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::fem::cable::Cable2d3DR::initialize(siconos::modeling::Interaction& inter) {
  auto qSize = inter.getSizeOfDS();

  if (!jacobianhOver_q_internal_storage_) {
    jacobianhOver_q_internal_storage_ = std::make_unique<std::vector<double>>(2 * qSize);
  }
  jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianhOver_q_internal_storage_->data(), 2, qSize);
  jacobianhOver_q_view_->setZero();
}

void siconos::fem::cable::Cable2d3DR::computeh(const siconos::algebra::BlockVector& q,
                                               Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_BEGIN("Cable2d3DR::computeh(...)\n");

  // LagrangianScleronomousR::computeh(q, z, y);
  auto& position = *((q.getAllVect())[0]);
  contactPoint1_(0) = position(_node_dof_index);
  contactPoint1_(1) = position(_node_dof_index + 1);
  contactPoint1_(2) = position(_node_dof_index + 2);
  y(0) = distance();
  DEBUG_EXPR(siconos::algebra::print(y););
  DEBUG_EXPR(display(););
  DEBUG_END("Cable2d3DR::computeh(...)\n")
}

void siconos::fem::cable::Cable2d3DR::computeJacobianhOver_q(
    const siconos::algebra::BlockVector& q) {
  DEBUG_BEGIN(
      "Cable2d3DR::computeJacobianhOver_q(const siconos::algebra::BlockVector& q, "
      "siconos::algebra::BlockVector& z \n");

  jacobianhOver_q_view_->row(0).segment(_node_dof_index, 3) = nc_;
  jacobianhOver_q_view_->row(1).segment(_node_dof_index, 3) = tangent_;

  if (q.size() == 6) {
    DEBUG_PRINT("take into account second ds\n");
    THROW_EXCEPTION("Cable2d3DR is not implemented for cable/cable contact");
  }

  DEBUG_END(
      "Cable2d3DR::computeJacobianhOver_q(const siconos::algebra::BlockVector& q, "
      "siconos::algebra::BlockVector& z) \n");
}

double siconos::fem::cable::Cable2d3DR::distance() const {
  DEBUG_BEGIN("Cable2d3DR::distance(...)\n")
  siconos::algebra::SiconosVector dpc = contactPoint2_ - contactPoint1_;
  DEBUG_EXPR(siconos::algebra::print(contactPoint1_));
  DEBUG_EXPR(siconos::algebra::print(contactPoint2_););
  DEBUG_EXPR(siconos::algebra::print(dpc););
  DEBUG_END("Cable2d3DR::distance(...)\n")
  return dpc.norm() * (nc_.dot(dpc) >= 0 ? -1 : 1);
}

/** update the contact points from references
 */
void siconos::fem::cable::Cable2d3DR::updateContactPoints(
    const siconos::algebra::SiconosVector3& pc1, const siconos::algebra::SiconosVector3& pc2,
    const siconos::algebra::SiconosVector3& normal,
    const siconos::algebra::SiconosVector3& tangent) {
  contactPoint1_ = pc1;
  contactPoint2_ = pc2;
  nc_ = normal;
  tangent_ = tangent;
};

void siconos::fem::cable::Cable2d3DR::updateContactPoints(double pc1[3], double pc2[3],
                                                         double normal[3], double tangent[3]) {
  contactPoint1_(0) = pc1[0];
  contactPoint1_(1) = pc1[1];
  contactPoint1_(2) = pc1[2];
  contactPoint2_(0) = pc2[0];
  contactPoint2_(1) = pc2[1];
  contactPoint2_(2) = pc2[2];
  nc_(0) = normal[0];
  nc_(1) = normal[1];
  nc_(2) = normal[2];
  tangent_(0) = tangent[0];
  tangent_(1) = tangent[1];
  tangent_(2) = normal[2];
};

void siconos::fem::cable::Cable2d3DR::display() const {
  LagrangianR::display();

  std::cout << " _node_dof_index :" << _node_dof_index << "\n";

  std::cout << " _Pc1: \n";
  siconos::algebra::print(contactPoint1_);
  std::cout << " _Pc2 :\n";
  siconos::algebra::print(contactPoint2_);
  std::cout << " _Normal:\n";
  siconos::algebra::print(nc_);

  std::cout << " _Tangent:\n";
  siconos::algebra::print(tangent_);
}
