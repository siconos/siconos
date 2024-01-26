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

/*! \file FiniteElementLinearTIDS.hpp */
#ifndef FINITEELEMENTLAGRANGIANTIDS_H
#define FINITEELEMENTLAGRANGIANTIDS_H

#include "LagrangianLinearTIDS.hpp"
#include "SiconosAlgebraTypes.hpp"  // For UblasType
/** Finite Element discretization of elastic solids that inherits from
 * Lagrangian Linear Systems with time invariant coefficients
 * - \f$M\dot v + Cv + Kq = F_{ext}(t,z) + p \f$
 */
namespace siconos::mechanics::fem {

class Mesh;
class FiniteElementModel;
class Material;

class FiniteElementLinearTIDS : public modeling::LagrangianLinearTIDS {
 protected:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(LagrangianLinearTIDS);

  /** a mesh */
  std::shared_ptr<Mesh> _mesh{nullptr};

  /** a material */
  std::map<unsigned int, std::shared_ptr<Material>> _materials;

  /** a finite element model */
  std::shared_ptr<FiniteElementModel> _FEModel{nullptr};

  /* Storage type for the matrices */
  siconos::algebra::UblasType _storageType{siconos::algebra::UblasType::DENSE};

 public:
  /** constructor from initial state and all matrix operators.
   * \param mesh the mesh that defined the spatial discretization
   * \param material
   */
  FiniteElementLinearTIDS(
      std::shared_ptr<Mesh> mesh, std::map<unsigned int, std::shared_ptr<Material>> materials,
      siconos::algebra::UblasType storageType = siconos::algebra::UblasType::DENSE);

  /** destructor */
  ~FiniteElementLinearTIDS() noexcept = default;

  void setStorageType(siconos::algebra::UblasType type) { _storageType = type; }

  std::shared_ptr<FiniteElementModel> FEModel() { return _FEModel; };

  void applyDirichletBoundaryConditions(int physical_entity_tag,
                                        std::shared_ptr<std::vector<int>> node_dof_index);

  void applyNodalForces(int physical_entity_tag,
                        std::shared_ptr<siconos::algebra::SiconosVector> nodal_forces);

  double elasticPotentialEnergy() const;

  void display(bool brief = true) const override;
};

}  // namespace siconos::mechanics::fem
#endif  // FINITEELEMENTLAGRANGIANTIDS_H
