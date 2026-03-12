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
#include "FiniteElementLinearTIDS.hpp"

#include "BoundaryCondition.hpp"
#include "FiniteElementModel.hpp"
#include "Material.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #include "siconos_debug.h"

siconos::mechanics::fem::FiniteElementLinearTIDS::FiniteElementLinearTIDS(
    std::shared_ptr<Mesh> mesh, const std::map<int, const Material>& materials)
    : LagrangianSparseLinearTIDS(), mesh_(mesh), materials_(materials) {
  // Warning FP: the DS is built from default empty constructor in an unusual
  // way. Care must be taken to properly set all attributes in DS, SecondOrder,
  // Lagrangian ... It may be better to :
  // - build FEModel from mesh
  // - compute ndof
  // - build ds from ndof or q0, v0

  FEModel_ = std::make_shared<FiniteElementModel>(mesh);
  ndof_ = FEModel_->init();

  q0_storage_ = std::make_unique<siconos::algebra::SiconosVector>(ndof_);
  velocity0_storage_ = std::make_unique<siconos::algebra::SiconosVector>(ndof_);
  use_q0([&](auto& v) { v.setZero(); });
  use_velocity0([&](auto& v) { v.setZero(); });

  // -- Memory allocation for vector and matrix members --
  state_q_[0] = std::make_shared<siconos::algebra::SiconosVector>(q0());
  state_q_[1] = std::make_shared<siconos::algebra::SiconosVector>(velocity0());

  p_[1] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  p_[1]->setZero();

  x_size_ = 2 * ndof_;

  // Mass ...
  // Deal with 'plugged' mass later
  mass_storage_ = std::make_unique<siconos::algebra::SiconosSparseMatrix>(ndof_, ndof_);
  hasConstantMass_ = true;
  hasMass_ = true;
  computemass_ = nullptr;
  hasLUMass_ = false;
  siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
      mass_storage_, [&](auto& matrix) { FEModel_->computeMassMatrix(matrix, materials_); },
      "mass_storage_");

  stiffnessMatrix_storage =
      std::make_unique<siconos::algebra::SiconosSparseMatrix>(ndof_, ndof_);
  siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
      stiffnessMatrix_storage,
      [&](auto& matrix) { FEModel_->computeStiffnessMatrix(matrix, materials_); },
      "stiffnessMatrix_storage");

  // dampingMatrix_ = std::make_shared<siconos::algebra::SiconosSparseMatrix>(ndof_, ndof_);
}

void siconos::mechanics::fem::FiniteElementLinearTIDS::applyDirichletBoundaryConditions(
    int physical_entity_tag, const std::vector<int>& node_dof_index, double imposedVelocity) {
  if (!boundaryConditions_)
    boundaryConditions_ = std::make_shared<siconos::modeling::BoundaryCondition>(
        siconos::modeling::BoundaryCondition::Indices{});

  FEModel_->applyDirichletBoundaryConditions(physical_entity_tag, node_dof_index,
                                             boundaryConditions_, imposedVelocity);

  reactionToBoundaryConditions_ = std::make_shared<siconos::algebra::SiconosVector>(
      boundaryConditions_->velocityIndices().size());
};

void siconos::mechanics::fem::FiniteElementLinearTIDS::applyNodalForces(
    int physical_entity_tag, const siconos::algebra::SiconosVector& nodal_forces) {
  if (!std::holds_alternative<siconos::algebra::OwnedDenseVector>(fext_storage_)) {
    fext_storage_ = std::make_unique<siconos::algebra::SiconosVector>(ndof_);
  }
  hasFext_ = true;
  hasConstantFext_ = true;  // set later by applyNodalForces
  //  We should probably use something like setComputeFext(applyNodalForce)
  computefext_ = nullptr;

  use_fext([&](auto& fext) {
    fext.setZero();
    FEModel_->applyNodalForces(physical_entity_tag, nodal_forces, fext);
  });
};

double siconos::mechanics::fem::FiniteElementLinearTIDS::elasticPotentialEnergy() const {
  siconos::algebra::SiconosVector tmp{ndof_};
  useStiffness([&](const auto& K) { tmp = K * *state_q_[0]; });
  return 0.5 * state_q_[0]->dot(tmp);
}

void siconos::mechanics::fem::FiniteElementLinearTIDS::display(bool brief) const {
  std::cout << "===== FiniteElementLinearTIDS display ===== " << std::endl;
  LagrangianSparseLinearTIDS::display();
  FEModel_->display(brief);
}
