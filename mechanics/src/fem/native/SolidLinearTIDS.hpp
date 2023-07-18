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

/*! \file SolidLinearTIDS.hpp

 */
#ifndef SOLIDLAGRANGIANTIDS_H
#define SOLIDLAGRANGIANTIDS_H

#include "FiniteElementLinearTIDS.hpp"

/** Finite Element discretization of an elastic solids that inherits from FEM with time invariant coefficients and exhibits sigma
 *  to enable easy plasticty implementation using a stress relation.
 * - \f$M\dot v + K u +  B^T \sigma = F_{ext}(t,z) + p
 * v = \dot u
 * \dot \sigma = B v \f$
 */
namespace siconos::mechanics::fem::native
{
class SolidLinearTIDS : public FiniteElementLinearTIDS
{

protected:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(LagrangianLinearTIDS);
  using siconos::mechanics::fem::native::FiniteElementLinearTIDS::_mesh;
  using siconos::mechanics::fem::native::FiniteElementLinearTIDS::_materials;
  using siconos::mechanics::fem::native::FiniteElementLinearTIDS::_storageType;
  /** M Matrix */
  std::shared_ptr<siconos::algebra::SimpleMatrix> _M;
  /** Elasticity Matrix */
  std::shared_ptr<siconos::algebra::SimpleMatrix> _S;
  /** B matrix from FEM */
  std::shared_ptr<siconos::algebra::SimpleMatrix> _B;
//std::shared_ptr<Mesh> _mesh;
  /** default constructor */
  SolidLinearTIDS():FiniteElementLinearTIDS() {};


public:

  /** constructor from initial state and all matrix operators.
   * \param mesh the mesh that defined the spatial discretization
   * \param material
   */
  SolidLinearTIDS(std::shared_ptr<Mesh> mesh,
                          std::map<unsigned int, std::shared_ptr<Material> > materials,
                          siconos::algebra::UblasType storageType=siconos::algebra::UblasType::DENSE);


  /** destructor */
  ~SolidLinearTIDS() {};

  void setStorageType(siconos::algebra::UblasType type)
  {
    _storageType=type;
  }


  std::shared_ptr<FiniteElementModel> FEModel()
  {
    return _FEModel;
  };


  void applyDirichletBoundaryConditions(int physical_entity_tag, std::shared_ptr<IndexInt> node_dof_index);

  void applyNodalForces(int physical_entity_tag, std::shared_ptr<siconos::algebra::SiconosVector> nodal_forces);

  /** Compute kinetic energy
   */
  double kineticEnergy() const;
  double elasticPotentialEnergy() const;

  void display(bool brief) const override;




};

} // namespace siconos::mechanics::fem::native
#endif // SOLIDLAGRANGIANTIDS_H
