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
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::mechanics::fem::FiniteElementLinearTIDS::FiniteElementLinearTIDS(
    std::shared_ptr<Mesh> mesh, std::map<unsigned int, std::shared_ptr<Material>> materials,
    siconos::algebra::StorageType storageType)
    : LagrangianLinearTIDS::LagrangianLinearTIDS(),
      _mesh(mesh),
      _materials(materials),
      _storageType(storageType) {
  DEBUG_BEGIN(
      "FiniteElementLinearTIDS::FiniteElementLinearTIDS(std::shared_ptr<Mesh> "
      "mesh, "
      "std::shared_ptr<Material> material\n");

  // Warning FP: the DS is built from default empty constructor in an unusual
  // way. Care must be taken to properly set all attributes in DS, SecondOrder,
  // Lagrangian ... It may be better to :
  // - build FEModel from mesh
  // - compute ndof
  // - build ds from ndof or q0, v0

  _FEModel = std::make_shared<FiniteElementModel>(mesh);
  ndof_ = _FEModel->init();

  q0_internal_storage_ = std::make_unique<std::vector<double>>(ndof_);
  q0_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(q0_internal_storage_->data(), ndof_);
  q0_view_->setZero();
  velocity0_internal_storage = std::make_unique<std::vector<double>>(ndof_);
  velocity0_view_ = std::make_shared<siconos::algebra::MapVectorType>(
      velocity0_internal_storage->data(), ndof_);
  velocity0_view_->setZero();

  // -- Memory allocation for vector and matrix members --
  state_q_[0] = std::make_shared<siconos::algebra::SiconosVector>(*q0_view_);
  state_q_[1] = std::make_shared<siconos::algebra::SiconosVector>(*velocity0_view_);

  p_[1] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  p_[1]->setZero();
  //   _zeroPlugin();
  x_size_ = 2 * ndof_;

  // Mass ...
  // Deal with 'plugged' mass later
  hasConstantMass_ = true;
  hasMass_ = true;
  computemass_ = nullptr;
  mass_internal_storage_ = std::make_unique<std::vector<double>>(ndof_ * ndof_);
  mass_view_ = std::make_shared<siconos::algebra::MapType>(mass_internal_storage_->data(),
                                                           ndof_, ndof_);
  // _mass->setIsSymmetric(true);
  // _mass->setIsPositiveDefinite(true);

  _FEModel->computeMassMatrix(mass_view_, _materials);
  //

  // Stiffness must be set with setStiffnessMatrix ...
  // stiffnessMatrix_ = std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_);
  // stiffnessMatrix_->setIsSymmetric(true);
  // stiffnessMatrix_->setIsPositiveDefinite(true);

  _FEModel->computeStiffnessMatrix(*stiffnessMatrix_view_, _materials);

  // if(!_C)
  // {
  //   _C = std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_,
  //   _storageType);
  // }
  // _C->setZero();

  DEBUG_END(
      "FiniteElementLinearTIDS::FiniteElementLinearTIDS(std::shared_ptr<Mesh> "
      "mesh, "
      "std::shared_ptr<Material> material\n");
}

void siconos::mechanics::fem::FiniteElementLinearTIDS::applyDirichletBoundaryConditions(
    int physical_entity_tag, std::shared_ptr<std::vector<int>> node_dof_index) {
  if (!boundaryConditions_)
    boundaryConditions_ = std::make_shared<siconos::modeling::BoundaryCondition>(
        siconos::modeling::BoundaryCondition::Indices{});

  _FEModel->applyDirichletBoundaryConditions(physical_entity_tag, node_dof_index,
                                             boundaryConditions_);

  reactionToBoundaryConditions_ = std::make_shared<siconos::algebra::SiconosVector>(
      boundaryConditions_->velocityIndices().size());
};

void siconos::mechanics::fem::FiniteElementLinearTIDS::applyNodalForces(
    int physical_entity_tag, std::shared_ptr<siconos::algebra::SiconosVector> nodal_forces) {
  if (!fext_view_) {
    if (!fext_internal_storage_) {
      fext_internal_storage_ = std::make_unique<std::vector<double>>(ndof_);
    }
    fext_view_ = std::make_shared<siconos::algebra::MapVectorType>(
        fext_internal_storage_->data(), ndof_);  // TODOSAM : what to do here ?
  }
  _FEModel->applyNodalForces(physical_entity_tag, nodal_forces, fext_view_);
};

double siconos::mechanics::fem::FiniteElementLinearTIDS::elasticPotentialEnergy() const {
  siconos::algebra::SiconosVector tmp{ndof_};
  tmp = *stiffnessMatrix_view_ * *state_q_[0];
  return 0.5 * state_q_[0]->dot(tmp);
}

void siconos::mechanics::fem::FiniteElementLinearTIDS::display(bool brief) const {
  std::cout << "===== FiniteElementLinearTIDS display ===== " << std::endl;
  LagrangianLinearTIDS::display();
  _FEModel->display(brief);
}
