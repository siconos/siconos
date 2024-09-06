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
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"
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
  _ndof = _FEModel->init();

  _q0 = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  _q0->zero();
  _velocity0 = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  _velocity0->zero();

  // -- Memory allocation for vector and matrix members --
  _q[0] = std::make_shared<siconos::algebra::SiconosVector>(*_q0);
  _q[1] = std::make_shared<siconos::algebra::SiconosVector>(*_velocity0);

  _p[1] = std::make_shared<siconos::algebra::SiconosVector>(_ndof);

  _zeroPlugin();
  _n = 2 * _ndof;

  if (!_mass) {
    _mass = std::make_shared<MapType>(nullptr, _ndof, _ndof); // TODOSAM : handle mass creation
    // _mass->setIsSymmetric(true);
    // _mass->setIsPositiveDefinite(true);
  }
  _FEModel->computeMassMatrix(_mass, _materials);

  if (!_K) {
    _K = std::make_shared<Matrix>(_ndof, _ndof);
    // _K->setIsSymmetric(true);
    // _K->setIsPositiveDefinite(true);
  }
  _FEModel->computeStiffnessMatrix(_K, _materials);

  // if(!_C)
  // {
  //   _C = std::make_shared<siconos::algebra::SiconosMatrix>(_ndof, _ndof,
  //   _storageType);
  // }
  // _C->zero();

  DEBUG_END(
      "FiniteElementLinearTIDS::FiniteElementLinearTIDS(std::shared_ptr<Mesh> "
      "mesh, "
      "std::shared_ptr<Material> material\n");
}

void siconos::mechanics::fem::FiniteElementLinearTIDS::applyDirichletBoundaryConditions(
    int physical_entity_tag, std::shared_ptr<std::vector<int>> node_dof_index) {
  if (!_boundaryConditions)
    _boundaryConditions = std::make_shared<siconos::modeling::BoundaryCondition>(
        siconos::modeling::BoundaryCondition::Indices{});

  _FEModel->applyDirichletBoundaryConditions(physical_entity_tag, node_dof_index,
                                             _boundaryConditions);

  _reactionToBoundaryConditions = std::make_shared<siconos::algebra::SiconosVector>(
      _boundaryConditions->velocityIndices().size());
};

void siconos::mechanics::fem::FiniteElementLinearTIDS::applyNodalForces(
    int physical_entity_tag, std::shared_ptr<siconos::algebra::SiconosVector> nodal_forces) {
  if (!_fExt) {
    _fExt = std::make_shared<siconos::algebra::SiconosVector>(dimension());
  }

  _FEModel->applyNodalForces(physical_entity_tag, nodal_forces, _fExt);
};

double siconos::mechanics::fem::FiniteElementLinearTIDS::elasticPotentialEnergy() const {
  auto tmp = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  siconos::algebra::prod(*_K, *q(), *tmp, true);
  return 0.5 * siconos::algebra::inner_prod(*q(), *tmp);
}

void siconos::mechanics::fem::FiniteElementLinearTIDS::display(bool brief) const {
  std::cout << "===== FiniteElementLinearTIDS display ===== " << std::endl;
  LagrangianLinearTIDS::display();
  _FEModel->display(brief);
}
