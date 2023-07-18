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

#include "FiniteElementModel.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "Mesh.hpp"
#include "SiconosVectorOp.hpp"
#include "SimpleMatrix.hpp"
//#include "SiconosAlgebraProd.hpp"
//#include "SimpleMatrixFriends.hpp"


// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "siconos_debug.h"

#include <stdio.h>
#include <iostream>

//siconos::mechanics::fem::SolidLinearTIDS::SolidLinearTIDS(std::shared_ptr<Mesh> mesh,
//    std::map<unsigned int, std::shared_ptr<Material> > materials,
//    siconos::algebra::UblasType storageType):
//  FiniteElementLinearTIDS::FiniteElementLinearTIDS(mesh, materials, storageType)
//{
//  DEBUG_BEGIN("SolidLinearTIDS::SolidLinearTIDS(std::shared_ptr<Mesh> mesh, std::shared_ptr<Material> material\n");

//  int dim = _FEModel->mesh()->dim();
//  int dimStress = dim*(dim+1)/2;
//  int nElements = _FEModel->elements().size();
////  _q0.reset(new siconos::algebra::SiconosVector(_ndof,0.0));
////  _velocity0.reset(new siconos::algebra::SiconosVector(_ndof,0.0));

//////  LagrangianDS::_init(_q0,_velocity0);

//  if(!_S)
//  {
//    _S.reset(new siconos::algebra::SimpleMatrix(nElements*dimStress, nElements*dimStress, _storageType));
//  }
//  _FEModel->computeBMatrix(_B, _materials);


//  if(!_B)
//  {
//    _B.reset(new siconos::algebra::SimpleMatrix(nElements*dimStress, _ndof, _storageType));
//  }
//  _FEModel->computeBMatrix(_B, _materials);

//  DEBUG_END("SolidLinearTIDS::SolidLinearTIDS(std::shared_ptr<Mesh> mesh, std::shared_ptr<Material> material\n");

//}

siconos::mechanics::fem::SolidLinearTIDS::SolidLinearTIDS(std::shared_ptr<Mesh> mesh,
    std::map<unsigned int, std::shared_ptr<Material> > materials,
    siconos::algebra::UblasType storageType) : FiniteElementLinearTIDS::FiniteElementLinearTIDS(mesh,materials,storageType)
{
    DEBUG_BEGIN("FiniteElementLinearTIDS::FiniteElementLinearTIDS(std::shared_ptr<Mesh> mesh, std::shared_ptr<Material> material\n");
//    _mesh = mesh;
//    _materials = materials;
//    _storageType = storageType;
    std::cout << "In SolidLinearTIDS constructor !" << std::endl;
    _FEModel.reset(new FiniteElementModel(mesh));
    _ndof = _FEModel->init();
    int nElements = _FEModel->elements().size();
    int dim = _FEModel->mesh()->dim();
    int dimStress = dim*(dim+1)/2;
    int dimState = _ndof*2 + dimStress*nElements;
    std::cout << "Fe model inintialized !" << std::endl;
    _q0.reset(new siconos::algebra::SiconosVector(_ndof,0.0));
    _velocity0.reset(new siconos::algebra::SiconosVector(_ndof,0.0));

//  LagrangianDS::_init(_q0,_velocity0);


//   // Mass matrix:
//    if(!_M)
//    {
//      _M.reset(new siconos::algebra::SimpleMatrix(_ndof, _ndof));
//    }
//  std::cout << "M resized !" << std::endl;
//  _FEModel->computeMassMatrix(_M, _materials);
//  std::cout << "M computed !" << std::endl;
  if(!_S)
  {
      _S = std::make_shared<siconos::algebra::SimpleMatrix>(dimStress*nElements,dimStress*nElements,
                                                               _storageType);
  }
  _FEModel->computeSMatrix(_S,_materials);
  std::cout << "S computed in constructor!" << std::endl;
  _S->display();
//  if(!_mass)
//  {
//      _mass.reset(new siconos::algebra::SimpleMatrix(dimState, dimState, _storageType));
//      _mass->setIsSymmetric(true);
//      _mass->setIsPositiveDefinite(true);
//  }
//  std::cout << "Setting mass matrix..." << std::endl;
//  for (int i=0;i<_ndof;i++){
//      for (int j=0;j<_ndof;j++){
//          _mass->setValue(i, j, _M->getValue(i,j));
//      }
//  }
//  for (int i=0;i<_ndof;i++){
//      _mass->setValue(i+_ndof, i+_ndof, 1.0);
//  }
//  for (int i=0;i<dimStress*nElements;i++){
//        for (int j=0;j<dimStress*nElements;j++){
//      _mass->setValue(2*_ndof + i, 2*_ndof + j, _S->getValue(i,j));
//        }
//  }
//  std::cout << "mass matrix set !" << std::endl;

    if(!_B)
    {
        _B = std::make_shared<siconos::algebra::SimpleMatrix>(nElements*dimStress, _ndof,
                                                                 _storageType);
    }
    _FEModel->computeBMatrix(_B, _materials);

//    if(!_C)
//    {
//        _C.reset(new siconos::algebra::SimpleMatrix(dimState, dimState, _storageType));
//    }
//    std::cout << "C matrix computation..." << std::endl;
//  for (int i=0;i<_ndof;i++){
//      for (int j=0;j<3*nElements;j++){
//          _C->setValue(i, j+2*_ndof, _B->getValue(j,i));
//      }
//  }
//  for (int i=0;i<3*nElements;i++){
//      for (int j=0;j<_ndof;j++){
//          _C->setValue(i+2*_ndof, j, -_B->getValue(i,j));
//      }
//  }
//  for (int i=0;i<_ndof;i++){
//      _C->setValue(i+_ndof, i, 1.0);
//  }
//  std::cout << "C matrix set !" << std::endl;

//  std::cout << "display of K in constructor:" << std::endl;


//  if(!_K)
//  {
//      _K.reset(new siconos::algebra::SimpleMatrix(dimState, dimState, _storageType));
//  }
//  _K->zero();
//  _K->display();
  // if(!_C)
  // {
  //   _C.reset(new SimpleMatrix(_ndof, _ndof, _storageType));
  // }
  // _C->zero();

  DEBUG_END("FiniteElementLinearTIDS::FiniteElementLinearTIDS(std::shared_ptr<Mesh> mesh, std::shared_ptr<Material> material\n");

}



double siconos::mechanics::fem::SolidLinearTIDS::kineticEnergy() const
{

  siconos::algebra::SiconosVector * tmp = new siconos::algebra::SiconosVector(_ndof);
  siconos::algebra::SiconosVector& velocity = *(_q[1]);
  siconos::algebra::prod(*_mass, velocity, *tmp, true);
  double kineticEnergy = 0.5*inner_prod(velocity,*tmp);
  return kineticEnergy;
}

double siconos::mechanics::fem::SolidLinearTIDS::elasticPotentialEnergy() const
{

  siconos::algebra::SiconosVector * tmp = new siconos::algebra::SiconosVector(_ndof);
  siconos::algebra::SiconosVector& displacement = *(_q[0]);
  siconos::algebra::prod(*_K, displacement, *tmp, true);
  double potentialEnergy = 0.5*inner_prod(displacement,   *tmp);
  return potentialEnergy;
  
}
void siconos::mechanics::fem::SolidLinearTIDS::display(bool brief) const
{
  std::cout << "===== SolidLinearTIDS display ===== " <<std::endl;
  LagrangianLinearTIDS::display();
  _FEModel->display(brief);
}
