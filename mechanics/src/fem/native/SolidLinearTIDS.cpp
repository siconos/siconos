/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include "SolidLinearTIDS.hpp"

#include <stdio.h>

#include <iostream>

#include "FiniteElementLinearTIDS.hpp"
#include "FiniteElementModel.hpp"
#include "Mesh.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES

#include "siconos_debug.h"

siconos::mechanics::fem::SolidLinearTIDS::SolidLinearTIDS(
    std::shared_ptr<Mesh> mesh, const std::map<int, const Material>& materials)
    : FiniteElementLinearTIDS(mesh, materials) {
  DEBUG_BEGIN("SolidLinearTIDS::SolidLinearTIDS(...)\n");

  auto nElements = FEModel_->elements().size();
  auto dim = FEModel_->mesh()->dim();

  // Computes the size of the stress tensor
  size_t dimStressInElem = (dim == 2) ? 3 : 6;  // Works also for Beam cases
  dimStress_ = siconos::algebra::to_index(dimStressInElem * nElements);

  x_size_ = ndof_ + dimStress_;
  // aka the true unknowns are v and sigma
  //    _q0.reset(new siconos::algebra::SiconosVector(ndof_,0.0));
  //    _velocity0.reset(new siconos::algebra::SiconosVector(ndof_,0.0));

  sigma0_ = std::make_unique<siconos::algebra::SiconosVector>(dimStress_);
  sigma0_->setZero();
  epsilon_p0_ = std::make_unique<siconos::algebra::SiconosVector>(dimStress_);
  epsilon_p0_->setZero();

  // -- Memory allocation for stress --
  //    _sigma = std::make_shared<siconos::algebra::SiconosVector>(*sigma0_);
  stress_[0] = std::make_shared<siconos::algebra::SiconosVector>(*sigma0_);
  stress_[1] = std::make_shared<siconos::algebra::SiconosVector>(*sigma0_);
  stress_[2] = std::make_shared<siconos::algebra::SiconosVector>(*sigma0_);
  epsilon_p_[0] = std::make_shared<siconos::algebra::SiconosVector>(dimStress_);
  epsilon_p_[1] = std::make_shared<siconos::algebra::SiconosVector>(*epsilon_p0_);

  // -- Initialize the elasticity matrix --
  elasticity_matrix_storage_ =
      std::make_unique<siconos::algebra::SiconosSparseMatrix>(dimStress_, dimStress_);
  siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
      elasticity_matrix_storage_,
      [&](auto& matrix) { FEModel_->computeElasticityMatrix(matrix, materials_); },
      "elasticity_matrix_storage_");

  // -- Initialize the B matrix --
  B_matrix_storage_ =
      std::make_unique<siconos::algebra::SiconosSparseMatrix>(dimStress_, ndof_);
  siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
      B_matrix_storage_, [&](auto& matrix) { FEModel_->computeBMatrix(matrix, materials_); },
      "B_matrix_storage_");

  DEBUG_END(
      "SolidLinearTIDS::SolidLinearTIDS(std::shared_ptr<Mesh> mesh, std::shared_ptr<Material> "
      "material\n");
}

// --- Functions for memory handling ---
void siconos::mechanics::fem::SolidLinearTIDS::initMemory(
    siconos::algebra::blocks::size_type steps) {
  LagrangianSparseDS::initMemory(steps);
  stressMemory_.setMemorySize(steps, dimStress_);
  // plasticDeformationMemory_.setMemorySize(steps, dimStress_);
  // plasticRate_Memory_.setMemorySize(steps, dimStress_);
}

void siconos::mechanics::fem::SolidLinearTIDS::swapInMemory() {
  LagrangianSparseDS::swapInMemory();
  stressMemory_.swap(*stress_[0]);
  // stressRateMemory_.swap(*stress_[1]);
  // plasticDeformationMemory_.swap(*epsilon_p_[0]);
  // plasticRate_Memory_.swap(*epsilon_p_[1]);
}

void siconos::mechanics::fem::SolidLinearTIDS::display(bool brief) const {
  std::cout << "===== SolidLinearTIDS display ===== " << std::endl;
  FiniteElementLinearTIDS::display();
}