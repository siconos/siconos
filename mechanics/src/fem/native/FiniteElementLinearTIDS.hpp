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

#include "LagrangianSparseLinearTIDS.hpp"

/** Finite Element discretization of elastic solids that inherits from
 * Lagrangian Linear Systems with time invariant coefficients
 * - \f$M\dot v + Cv + Kq = F_{ext}(t,z) + p \f$
 */
namespace siconos::mechanics::fem {

class Mesh;
class FiniteElementModel;
class Material;

class FiniteElementLinearTIDS : public modeling::LagrangianSparseLinearTIDS {
 protected:
  /** a mesh */
  std::shared_ptr<Mesh> mesh_{nullptr};

  /** materials */
  std::map<unsigned int, const Material> materials_;

  /** a finite element model */
  std::shared_ptr<FiniteElementModel> FEModel_{nullptr};

 public:
  /** constructor from initial state and all matrix operators.
   * \param mesh the mesh that defined the spatial discretization
   * \param material
   */
  FiniteElementLinearTIDS(std::shared_ptr<Mesh> mesh,
                          const std::map<unsigned int, const Material>& materials);

  /** destructor */
  ~FiniteElementLinearTIDS() noexcept = default;

  std::shared_ptr<FiniteElementModel> FEModel() { return FEModel_; };

  /** @brief Apply Dirichlet Boundary conditions for all nodes of all elements matching a tag
   *
   * @param[in] physical_entity_tag the tag to be matched
   * @param[in] node_dof_index list of dof to be constrained (local index in the node)
   **/
  void applyDirichletBoundaryConditions(int physical_entity_tag,
                                        const std::vector<int>& node_dof_index);

  /** @brief Apply Neuman Boundary conditions (nodal forces) for all elements matching a tag
   *
   * @param[in] physical_entity_tag the tag to be matched
   * @param[in] nodal forces input vector (contains only the values to be set to tagged elem)
   **/
  void applyNodalForces(int physical_entity_tag,
                        const siconos::algebra::SiconosVector& nodal_forces);

  double elasticPotentialEnergy() const;

  void display(bool brief = true) const override;
};

}  // namespace siconos::mechanics::fem
#endif  // FINITEELEMENTLAGRANGIANTIDS_H
