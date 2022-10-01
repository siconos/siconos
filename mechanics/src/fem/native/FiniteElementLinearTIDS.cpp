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

#include "FiniteElementLinearTIDS.hpp"
#include "FiniteElementModel.hpp"
#include "Material.hpp"
#include "SiconosMatrix.hpp"
#include "BoundaryCondition.hpp"



// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

#include <stdio.h>
#include <iostream>

siconos::mechanics::fem::native::FiniteElementLinearTIDS::FiniteElementLinearTIDS(std::shared_ptr<Mesh> mesh,
    std::map<unsigned int, std::shared_ptr<Material> > materials,
    Siconos::UBLAS_TYPE storageType):
  LagrangianLinearTIDS::LagrangianLinearTIDS(), _mesh(mesh), _materials(materials), _storageType(storageType)
{
  DEBUG_BEGIN("FiniteElementLinearTIDS::FiniteElementLinearTIDS(std::shared_ptr<Mesh> mesh, std::shared_ptr<Material> material\n");

  _FEModel.reset(new FiniteElementModel(mesh));
  _ndof = _FEModel->init();

  _q0.reset(new SiconosVector(_ndof,0.0));
  _velocity0.reset(new SiconosVector(_ndof,0.0));

  LagrangianDS::_init(_q0,_velocity0);

  if(!_mass)
  {
    _mass.reset(new SimpleMatrix(_ndof, _ndof, _storageType));
    _mass->setIsSymmetric(true);
    _mass->setIsPositiveDefinite(true);
  }
  _FEModel->computeMassMatrix(_mass, _materials);

  if(!_K)
  {
    _K.reset(new SimpleMatrix(_ndof, _ndof, _storageType));
    _K->setIsSymmetric(true);
    _K->setIsPositiveDefinite(true);
  }
  _FEModel->computeStiffnessMatrix(_K, _materials);

  // if(!_C)
  // {
  //   _C.reset(new SimpleMatrix(_ndof, _ndof, _storageType));
  // }
  // _C->zero();

  DEBUG_END("FiniteElementLinearTIDS::FiniteElementLinearTIDS(std::shared_ptr<Mesh> mesh, std::shared_ptr<Material> material\n");

}
void siconos::mechanics::fem::native::FiniteElementLinearTIDS::applyDirichletBoundaryConditions(int physical_entity_tag, std::shared_ptr<IndexInt> node_dof_index)
{

  if(!_boundaryConditions)
  {
    std::shared_ptr<IndexInt> bdIndex = std::make_shared<IndexInt>(0);
    std::shared_ptr<SiconosVector> bdPrescribedVelocity = std::make_shared<SiconosVector>(0);
    _boundaryConditions = std::make_shared<BoundaryCondition>(bdIndex,bdPrescribedVelocity);
  }

  _FEModel->applyDirichletBoundaryConditions(physical_entity_tag, node_dof_index, _boundaryConditions);
  _reactionToBoundaryConditions = std::make_shared<SiconosVector>(_boundaryConditions->velocityIndices()->size());

};

void siconos::mechanics::fem::native::FiniteElementLinearTIDS::applyNodalForces(int physical_entity_tag, std::shared_ptr<SiconosVector> nodal_forces)
{

  if(!_fExt)
  {
    _fExt =  std::make_shared<SiconosVector>(dimension());
    _fExt->zero();
  }

  _FEModel->applyNodalForces(physical_entity_tag, nodal_forces, _fExt);

};

void siconos::mechanics::fem::native::FiniteElementLinearTIDS::display(bool brief) const
{
  std::cout << "===== FiniteElementLinearTIDS display ===== " <<std::endl;
  LagrangianLinearTIDS::display();
  _FEModel->display(brief);
}
