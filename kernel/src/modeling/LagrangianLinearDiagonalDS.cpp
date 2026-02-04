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

#include "StorageTools.hpp"
#include "siconos_debug.h"

// --- Constructor for the undamped system with identity mass matrix
siconos::modeling::LagrangianLinearDiagonalDS::LagrangianLinearDiagonalDS(
    Eigen::Ref<siconos::algebra::SiconosVector> q0,
    Eigen::Ref<siconos::algebra::SiconosVector> v0,
    Eigen::Ref<siconos::algebra::SiconosVector> stiffness_diag, siconos::algebra::AliasTag)
    : LagrangianDS(q0, v0, siconos::algebra::alias_t) {
  setStiffnessMatrix(stiffness_diag, siconos::algebra::alias_t);
}

siconos::modeling::LagrangianLinearDiagonalDS::LagrangianLinearDiagonalDS(
    const siconos::algebra::SiconosVector& q0, const siconos::algebra::SiconosVector& v0,
    const siconos::algebra::SiconosVector& stiffness_diag, siconos::algebra::CopyTag tag)
    : LagrangianDS(q0, v0, siconos::algebra::copy_t) {
  setStiffnessMatrix(stiffness_diag, siconos::algebra::copy_t);
}

void siconos::modeling::LagrangianLinearDiagonalDS::initRhs(double time) {
  THROW_EXCEPTION(
      "siconos::modeling::LagrangianLinearDiagonalDS::initRhs - not yet implemented for "
      "LagrangianLinearDiagonalDS.");
}

void siconos::modeling::LagrangianLinearDiagonalDS::setStiffnessMatrix(
    const siconos::algebra::SiconosVector& newValue, siconos::algebra::CopyTag tag) {
  assert(newValue.size() == ndof_);

  stiffnessMatrix_storage_ = std::make_unique<siconos::algebra::SiconosVector>(newValue);
  hasFint_ = true;
}

void siconos::modeling::LagrangianLinearDiagonalDS::setStiffnessMatrix(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue, siconos::algebra::AliasTag tag) {
  assert(newValue.size() == ndof_);

  stiffnessMatrix_storage_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), ndof_);
  hasFint_ = true;
}

void siconos::modeling::LagrangianLinearDiagonalDS::setDampingMatrix(
    const siconos::algebra::SiconosVector& newValue, siconos::algebra::CopyTag tag) {
  assert(newValue.size() == ndof_);

  dampingMatrix_storage_ = std::make_unique<siconos::algebra::SiconosVector>(newValue);
  hasFint_ = true;
}

void siconos::modeling::LagrangianLinearDiagonalDS::setDampingMatrix(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue, siconos::algebra::AliasTag tag) {
  assert(newValue.size() == ndof_);

  dampingMatrix_storage_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), ndof_);
  hasFint_ = true;
}

void siconos::modeling::LagrangianLinearDiagonalDS::setMassMatrix(
    const siconos::algebra::SiconosVector& newValue, siconos::algebra::CopyTag tag) {
  assert(newValue.size() == ndof_);

  massMatrix_storage_ = std::make_unique<siconos::algebra::SiconosVector>(newValue);
  hasMass_ = true;
}

void siconos::modeling::LagrangianLinearDiagonalDS::setMassMatrix(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue, siconos::algebra::AliasTag tag) {
  assert(newValue.size() == ndof_);

  massMatrix_storage_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), ndof_);
  hasMass_ = true;
}

void siconos::modeling::LagrangianLinearDiagonalDS::computeTotalForces(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time) {
  DEBUG_PRINT("LagrangianLinearTIDS::computeTotalForces(v,q,t) \n");

  if (!totalForces_) {
    totalForces_ = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  }

  use_stiffnessMatrix([&](auto const& K) { *totalForces_ -= K.cwiseProduct(q); });
  if (hasDampingMatrix())
    use_dampingMatrix([&](auto const& C) { *totalForces_ -= C.cwiseProduct(velocity); });

  if (hasFext_) {
    computeFext(time);
    use_fext([&](auto const& fext) { *totalForces_ += fext; });
  }
}

void siconos::modeling::LagrangianLinearDiagonalDS::display(bool brief) const {
  LagrangianDS::display();
  std::cout << "===== Lagrangian Linear Diagonal System display ===== \n ";
  std::cout << "- Mass Matrix M : \n";
  if (hasStiffnessMatrix()) {
    std::cout << "- Stiffness matrix\n ";
    use_stiffnessMatrix([&](const auto& M) { siconos::algebra::print(M); });
    std::cout << "\n";
  }
  if (hasDampingMatrix()) {
    std::cout << "- Damping matrix\n ";
    use_dampingMatrix([&](const auto& M) { siconos::algebra::print(M); });
    std::cout << "\n";
  }
  if (hasMassMatrix()) {
    std::cout << "- Mass matrix\n ";
    use_massMatrix([&](const auto& M) { siconos::algebra::print(M); });
    std::cout << "\n";
  }

  std::cout << "=========================================================== \n";
}
