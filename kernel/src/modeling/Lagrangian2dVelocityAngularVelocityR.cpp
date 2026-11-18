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

#include "Lagrangian2dVelocityAngularVelocityR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "SiconosException.hpp"
#include "siconos_debug.h"

void siconos::modeling::Lagrangian2dVelocityAngularVelocityR::initialize(Interaction& inter) {
  // proj_with_q  jacobianhOver_q_Proj =
  // std::make_shared<siconos::algebra::SiconosMatrix>(jacobianhOver_q_->rows(),jacobianhOver_q_->cols()));

  if ((inter.getSizeOfDS() != 3) and (inter.getSizeOfDS() != 6)) {
    THROW_EXCEPTION(
        "siconos::modeling::Lagrangian2dVelocityAngularVelocityR::initialize(Interaction& "
        "inter). The size of ds "
        "must of size 3 or 6");
  }
  auto qSize = 3 * (inter.getSizeOfDS() / 3);
  if (!H_storage_) {
    H_storage_ = std::make_unique<std::vector<double>>(3 * qSize);
  }
  H_view_ = std::make_shared<siconos::algebra::MapType>(H_storage_->data(), 3, qSize);
  H_view_->setZero();
}

void siconos::modeling::Lagrangian2dVelocityAngularVelocityR::computeH(
    const siconos::algebra::BlockVector& q) {
  DEBUG_BEGIN(
      "siconos::modeling::Lagrangian2dVelocityAngularVelocityR::computeH(Interaction& inter, "
      "siconos::algebra::BlockVector q0 \n");
  siconos::algebra::SiconosVector2 tangent;
  tangent << -nc_.y(), nc_.x();
  siconos::algebra::SiconosVector2 lever_arm;
  lever_arm << contactPoint1_.x() - q(0), contactPoint1_.y() - q(1);

  DEBUG_PRINTF("lever_arm_x = %4.2e,\t lever_arm_ y = %4.2e\n", lever_arm.x(), lever_arm.y());
  H_view_->row(0).segment(0, 2) = nc_;
  (*H_view_)(0, 2) = lever_arm.x() * nc_.y() - lever_arm.y() * nc_.x();
  H_view_->row(1).segment(0, 2) << tangent.x(), tangent.y();
  (*H_view_)(1, 2) = lever_arm.x() * tangent.y() - lever_arm.y() * tangent.x();

  H_view_->row(2).setZero();
  (*H_view_)(2, 2) = 1.0;

  if (q.size() == 6) {
    (*H_view_)(2, 5) = -1.0;
    lever_arm << contactPoint1_.x() - q(3), contactPoint1_.y() - q(4);

    H_view_->row(0).segment(3, 2) = -nc_;
    (*H_view_)(0, 5) = -lever_arm.x() * nc_.y() + lever_arm.y() * nc_.x();
    H_view_->row(1).segment(3, 2) << -tangent.x(), -tangent.y();
    (*H_view_)(1, 5) = -lever_arm.x() * tangent.y() + lever_arm.y() * tangent.x();
  }
  DEBUG_EXPR(std::cout << H_view_ << "\n";);
  DEBUG_END(
      "siconos::modeling::Lagrangian2dVelocityAngularVelocityR::computeH(Interaction& inter, "
      "siconos::algebra::BlockVector q0) \n");
}

double siconos::modeling::Lagrangian2dVelocityAngularVelocityR::distance() const {
  DEBUG_BEGIN("siconos::modeling::Lagrangian2dVelocityAngularVelocityR::distance(...)\n")
  siconos::algebra::SiconosVector dpc = contactPoint2_ - contactPoint1_;
  DEBUG_EXPR(siconos::algebra::print(contactPoint1_););
  DEBUG_EXPR(siconos::algebra::print(contactPoint2_););
  DEBUG_EXPR(siconos::algebra::print(dpc););
  DEBUG_END("siconos::modeling::Lagrangian2dVelocityAngularVelocityR::distance(...)\n")
  return dpc.norm() * (nc_.dot(dpc) >= 0 ? -1 : 1);
}

void siconos::modeling::Lagrangian2dVelocityAngularVelocityR::computeh(
    const siconos::algebra::BlockVector& q, Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_BEGIN("siconos::modeling::Lagrangian2dVelocityAngularVelocityR::computeh(...)\n");
  DEBUG_EXPR(siconos::algebra::print(q));

  DEBUG_EXPR(siconos::algebra::print(contactPoint1_););
  DEBUG_EXPR(siconos::algebra::print(contactPoint2_););
  DEBUG_EXPR(siconos::algebra::print(nc_););

  LagrangianScleronomousR::computeh(q, y);
  y(0) = distance();
  DEBUG_EXPR(siconos::algebra::print(y););
  DEBUG_EXPR(display(););
  DEBUG_END("siconos::modeling::Lagrangian2dVelocityAngularVelocityR::computeh(...)\n")
}

void siconos::modeling::Lagrangian2dVelocityAngularVelocityR::display() const {
  LagrangianR::display();

  std::cout << " _Pc1 :" << std::endl;
  siconos::algebra::print(contactPoint1_);
  std::cout << " _Pc2 :" << std::endl;
  siconos::algebra::print(contactPoint2_);
  std::cout << " _Nc :" << std::endl;
  siconos::algebra::print(nc_);
}
