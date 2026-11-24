// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "Lagrangian2dNormalVelocityR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::modeling::Lagrangian2dNormalVelocityR::initialize(Interaction& inter) {
  // proj_with_q  jacobianhOver_q_Proj =
  // std::make_shared<siconos::algebra::SiconosMatrix>(jacobianhOver_q_->rows(),jacobianhOver_q_->cols()));

  if ((inter.getSizeOfDS() != 3) and (inter.getSizeOfDS() != 6)) {
    THROW_EXCEPTION(
        "siconos::modeling::Lagrangian2dNormalVelocityR::initialize(Interaction& inter). The "
        "size of ds "
        "must of size 3 or 6");
  }

  auto qSize = 3 * (inter.getSizeOfDS() / 3);
  if (!contactMatrix_storage_) {
    contactMatrix_storage_ = std::make_unique<std::vector<double>>(qSize);
  }
  contactMatrix_view_ =
      std::make_shared<siconos::algebra::MapType>(contactMatrix_storage_->data(), 1, qSize);
  contactMatrix_view_->setZero();
}

void siconos::modeling::Lagrangian2dNormalVelocityR::computeContactMatrix(
    const siconos::algebra::BlockVector& q) {
  DEBUG_BEGIN(
      "siconos::modeling::Lagrangian2dNormalVelocityR::computeContactMatrix(Interaction& "
      "inter, "
      "siconos::algebra::BlockVector q0 \n");

  double G1x = q(0);
  double G1y = q(1);

  contactMatrix_view_->row(0).segment(0, 2) = nc_;

  (*contactMatrix_view_)(0, 2) =
      (G1y - contactPoint1_.y()) * nc_.x() - (G1x - contactPoint1_.x()) * nc_.y();

  if (q.size() == 6) {
    DEBUG_PRINT("take into account second ds\n");
    double G2x = q(3);
    double G2y = q(4);

    contactMatrix_view_->row(0).segment(3, 2) = -nc_;
    (*contactMatrix_view_)(0, 5) =
        (G2y - contactPoint1_.y()) * nc_.x() - (G2x - contactPoint1_.x()) * nc_.y();
  }
  DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_););
  DEBUG_END(
      "siconos::modeling::Lagrangian2dNormalVelocityR::computeContactMatrix(Interaction& "
      "inter, "
      "siconos::algebra::BlockVector q0) \n");
}

double siconos::modeling::Lagrangian2dNormalVelocityR::distance() const {
  DEBUG_BEGIN("siconos::modeling::Lagrangian2dNormalVelocityR::distance(...)\n")
  siconos::algebra::SiconosVector2 dpc = contactPoint2_ - contactPoint1_;
  DEBUG_END("siconos::modeling::Lagrangian2dNormalVelocityR::distance(...)\n")
  return dpc.norm() * (nc_.dot(dpc) >= 0 ? -1 : 1);
}

void siconos::modeling::Lagrangian2dNormalVelocityR::computeh(
    const siconos::algebra::BlockVector& q, Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_BEGIN("siconos::modeling::Lagrangian2dNormalVelocityR::computeh(...)\n");
  DEBUG_EXPR(siconos::algebra::print(q));

  LagrangianScleronomousR::computeh(q, y);
  y(0) = distance();
  DEBUG_EXPR(siconos::algebra::print(y););
  DEBUG_EXPR(display(););
  DEBUG_END("siconos::modeling::Lagrangian2dNormalVelocityR::computeh(...)\n")
}
void siconos::modeling::Lagrangian2dNormalVelocityR::display() const {
  LagrangianR::display();

  std::cout << " _Pc1 :" << std::endl;
  siconos::algebra::print(contactPoint1_);
  std::cout << " _Pc2 :" << std::endl;

  siconos::algebra::print(contactPoint2_);

  std::cout << " _Nc :" << std::endl;
  siconos::algebra::print(nc_);
}
