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
#include "LagrangianLinearDiagonalDS.hpp"

#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <iostream>

#include "siconos_debug.h"

// --- Constructor for the undamped system with identity mass matrix
siconos::modeling::LagrangianLinearDiagonalDS::LagrangianLinearDiagonalDS(
    Eigen::Ref<siconos::algebra::SiconosVector> q0,
    Eigen::Ref<siconos::algebra::SiconosVector> v0,
    Eigen::Ref<siconos::algebra::SiconosVector> stiffness_diag)
    : LagrangianDS(q0, v0) {
  stiffnessMatrix_view_ = std::make_shared<siconos::algebra::MapVectorType>(
      stiffness_diag.data(), stiffness_diag.size());
}

void siconos::modeling::LagrangianLinearDiagonalDS::initRhs(double time) {
  THROW_EXCEPTION(
      "siconos::modeling::LagrangianLinearDiagonalDS::initRhs - not yet implemented for "
      "LagrangianLinearDiagonalDS.");
}

void siconos::modeling::LagrangianLinearDiagonalDS::setDampingMatrix(
    Eigen::Ref<siconos::algebra::SiconosVector> damping) {
  dampingMatrix_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(damping.data(), damping.size());
}

void siconos::modeling::LagrangianLinearDiagonalDS::setMassMatrix(
    Eigen::Ref<siconos::algebra::SiconosVector> mass) {
  massMatrix_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(mass.data(), mass.size());
}

void siconos::modeling::LagrangianLinearDiagonalDS::computeTotalForces(
    const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time) {
  DEBUG_PRINT("LagrangianLinearTIDS::computeTotalForces(v,q,t) \n");

  if (!totalForces_) {
    totalForces_ = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  }

  *totalForces_ = -1. * stiffnessMatrix_view_->cwiseProduct(q);
  if (dampingMatrix_view_) *totalForces_ -= dampingMatrix_view_->cwiseProduct(velocity);

  if (fext_view_) {
    computeFext(time);
    *totalForces_ += *fext_view_;
  }
}

void siconos::modeling::LagrangianLinearDiagonalDS::display(bool brief) const {
  LagrangianDS::display();
  std::cout << "===== Lagrangian Linear Diagonal System display ===== \n ";
  std::cout << "- Mass Matrix M : \n";
  if (massMatrix_view_)
    std::cout << *massMatrix_view_ << "\n";
  else
    std::cout << "-> nullptr \n";
  std::cout << "- Stiffness Matrix K : " << *stiffnessMatrix_view_ << "\n";

  std::cout << "- Damping Matrix C : \n";
  if (dampingMatrix_view_)
    std::cout << *dampingMatrix_view_ << "\n";
  else
    std::cout << "-> nullptr\n";
  std::cout << "=========================================================== \n";
}
