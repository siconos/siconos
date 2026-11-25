// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
  if (!contactMatrix_storage_) {
    contactMatrix_storage_ = std::make_unique<std::vector<double>>(3 * qSize);
  }
  contactMatrix_view_ =
      std::make_shared<siconos::algebra::MapType>(contactMatrix_storage_->data(), 3, qSize);
  contactMatrix_view_->setZero();
}

void siconos::modeling::Lagrangian2dVelocityAngularVelocityR::computeContactMatrix(
    const siconos::algebra::BlockVector& q) {
  DEBUG_BEGIN(
      "siconos::modeling::Lagrangian2dVelocityAngularVelocityR::computeContactMatrix("
      "Interaction& inter, "
      "siconos::algebra::BlockVector q0 \n");
  siconos::algebra::SiconosVector2 tangent;
  tangent << -nc_.y(), nc_.x();
  siconos::algebra::SiconosVector2 lever_arm;
  lever_arm << contactPoint1_.x() - q(0), contactPoint1_.y() - q(1);

  DEBUG_PRINTF("lever_arm_x = %4.2e,\t lever_arm_ y = %4.2e\n", lever_arm.x(), lever_arm.y());
  contactMatrix_view_->row(0).segment(0, 2) = nc_;
  (*contactMatrix_view_)(0, 2) = lever_arm.x() * nc_.y() - lever_arm.y() * nc_.x();
  contactMatrix_view_->row(1).segment(0, 2) << tangent.x(), tangent.y();
  (*contactMatrix_view_)(1, 2) = lever_arm.x() * tangent.y() - lever_arm.y() * tangent.x();

  contactMatrix_view_->row(2).setZero();
  (*contactMatrix_view_)(2, 2) = 1.0;

  if (q.size() == 6) {
    (*contactMatrix_view_)(2, 5) = -1.0;
    lever_arm << contactPoint1_.x() - q(3), contactPoint1_.y() - q(4);

    contactMatrix_view_->row(0).segment(3, 2) = -nc_;
    (*contactMatrix_view_)(0, 5) = -lever_arm.x() * nc_.y() + lever_arm.y() * nc_.x();
    contactMatrix_view_->row(1).segment(3, 2) << -tangent.x(), -tangent.y();
    (*contactMatrix_view_)(1, 5) = -lever_arm.x() * tangent.y() + lever_arm.y() * tangent.x();
  }
  DEBUG_EXPR(std::cout << contactMatrix_view_ << "\n";);
  DEBUG_END(
      "siconos::modeling::Lagrangian2dVelocityAngularVelocityR::computeContactMatrix("
      "Interaction& inter, "
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
