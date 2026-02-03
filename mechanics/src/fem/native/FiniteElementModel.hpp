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

/*! \file FiniteElementLinearTIDS.hpp

 */
#ifndef FINITEELEMENTMODEL_H
#define FINITEELEMENTMODEL_H

#include <map>
#include <memory>
#include <span>
#include <unordered_set>
#include <vector>

#include "FETypes.hpp"
#include "FiniteElement.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::modeling {
class BoundaryCondition;
}

namespace siconos::mechanics::fem {

class FENode;
class FiniteElement;
class Mesh;
class MeshVertex;
class MeshElement;
class Material;

/** A finite element node */

constexpr double I_TET2_X[] = {0.13819660112501052, 0.13819660112501052, 0.13819660112501052,
                               0.58541019662496840};
constexpr double I_TET2_Y[] = {0.13819660112501052, 0.13819660112501052, 0.58541019662496840,
                               0.13819660112501052};
constexpr double I_TET2_Z[] = {0.13819660112501052, 0.58541019662496840, 0.13819660112501052,
                               0.13819660112501052};
constexpr double I_TET2_W[] = {0.04166666666666666, 0.04166666666666666, 0.04166666666666666,
                               0.04166666666666666};

/** a finite element model */
class FiniteElementModel {
 protected:
  /** a mesh */
  std::shared_ptr<Mesh> mesh_{nullptr};

  /** nodes */
  std::vector<std::shared_ptr<FENode>> nodes_ = {};

  /** elements */
  std::vector<std::shared_ptr<FiniteElement>> elements_ = {};

  /** vertex to node map **/
  std::map<std::shared_ptr<MeshVertex>, std::shared_ptr<FENode>> vertexToNode_;

  /** MeshElement to FiniteElement map **/
  std::map<std::shared_ptr<MeshElement>, std::shared_ptr<FiniteElement>> mElementTOFiniteElement_;

  /** Rule of five */
  FiniteElementModel() = delete;
  FiniteElementModel(FiniteElementModel&) = delete;
  FiniteElementModel& operator=(const FiniteElementModel&) = delete;
  FiniteElementModel(FiniteElementModel&&) = delete;
  FiniteElementModel& operator=(FiniteElementModel&&) = delete;

  // Helper function to check if the element type is valid (i.e. properly implemented in
  // siconos in the current FiniteElementModel
  inline bool is_valid_element(FiniteElementType value) {
    static const std::unordered_set<FiniteElementType> valid_values = {
        FiniteElementType::T3, FiniteElementType::TH4, FiniteElementType::B2,
        FiniteElementType::B3};
    return valid_values.count(value) > 0;
  }

 public:
  /** Constructor
      \param m the mesh
      \param ndof dof number
      \param e mesh element
   */
  FiniteElementModel(std::shared_ptr<Mesh> m) : mesh_(m) {};

  ~FiniteElementModel() noexcept = default;

  auto mesh() { return mesh_; }

  const std::span<const std::shared_ptr<FiniteElement>> elements() const { return elements_; }

  const std::span<const std::shared_ptr<FENode>> nodes() { return nodes_; }

  std::shared_ptr<FENode> vertexToNode(std::shared_ptr<MeshVertex> v) const;

  /** @brief creates the FEM model from the mesh and the element type
   *
   * \return the number of dof */
  siconos::algebra::Index init();

  /* Assembly method for elementary matrix */
  void assembleElementaryMatrix(siconos::algebra::SiconosSparseMatrix& M,
                                const siconos::algebra::SiconosDenseMatrix& Me,
                                const FiniteElement& fe);

  // should be computeMass of LagrangianDS ?

  /* Assembly method specific for elementary fem B matrix */
  void assembleElementary_B_Matrix(siconos::algebra::SiconosSparseMatrix& M,
                                   const siconos::algebra::SiconosDenseMatrix& Be,
                                   const FiniteElement& fe, size_t elemCnt);

  void assembleElementary_S_Matrix(siconos::algebra::SiconosSparseMatrix& S,
                                   const siconos::algebra::SiconosDenseMatrix& Se,
                                   const FiniteElement& fe, size_t elem_cnt);

  /** compute Mass Matrix
   **/
  void computeMassMatrix(siconos::algebra::SiconosSparseMatrix& mass,
                         const std::map<int, const Material>& mat);

  /** compute elementary Mass Matrix
   **/
  void computeElementaryMassMatrix(siconos::algebra::SiconosDenseMatrix&, const FiniteElement& fe,
                                   double massDensity);

  void computeBeamElementaryMassMatrix_direct(
      siconos::algebra::SiconosDenseMatrix& Me, FiniteElement& fe,
      const std::map<int, const Material>& materials);

  /** compute Stiffness Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeStiffnessMatrix(siconos::algebra::SiconosSparseMatrix&,
                              const std::map<int, const Material>& mat);

  /** compute elementary Stiffness Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeElementaryStiffnessMatrix(siconos::algebra::SiconosDenseMatrix& Me, FiniteElement& fe,
                                        const siconos::algebra::SiconosDenseMatrix& D,
                                        double thickness);

  /** compute elementary Stiffness Matrix with a direct method
   * for linear element
   **/
  void computeElementaryStiffnessMatrix_direct(siconos::algebra::SiconosDenseMatrix& Me,
                                               FiniteElement& fe,
                                               const siconos::algebra::SiconosDenseMatrix& D,
                                               double thickness);

  void computeBeamElementaryStiffnessMatrix_direct(
      siconos::algebra::SiconosDenseMatrix& Ke, FiniteElement& fe,
      const std::map<int, const Material>& materials);

  void computeElementary_B_Matrix(FiniteElement& fe, siconos::algebra::SiconosDenseMatrix& B,
                                  double length);

  void computeBeamElementaryBMatrix_direct(
      FiniteElement& fe, siconos::algebra::SiconosDenseMatrix& Be,
      const std::map<int, const Material>& materials);

  void computeElementaryBMatrix_direct(FiniteElement& fe, siconos::algebra::SiconosDenseMatrix& Be,
                                       double thickness);

  void computeBMatrix(siconos::algebra::SiconosSparseMatrix& B,
                      const std::map<int, const Material>& mat);

  void computeElasticityMatrix(siconos::algebra::SiconosSparseMatrix& S,
                               const std::map<int, const Material>& mat);

  /** @brief Apply Dirichlet Boundary conditions for all nodes of all elements matching a tag
   *
   * @param[in] physical_entity_tag the tag to be matched
   * @param[in] node_dof_index list of dof to be constrained (local index in the node)
   * @param[in,out] bc boundary conditions to be updated (global indices will be appended if
   * @param[in, optional] imposedVelocity, default = 0. Value chosen for the constrained dofs
   *not already present)
   **/
  void applyDirichletBoundaryConditions(
      int physical_entity_tag, const std::vector<int>& node_dof_index,
      std::shared_ptr<siconos::modeling::BoundaryCondition> bc, double imposedVelocity = 0);

  /** @brief Apply Neuman Boundary conditions (nodal forces) for all elements matching a tag
   *
   * @param[in] physical_entity_tag the tag to be matched
   * @param[in] nodal forces input vector (contains only the values to be set to tagged elem)
   * @param[in,out] forces output vector (all elements/nodes)
   **/
  void applyNodalForces(int physical_entity_tag,
                        const siconos::algebra::SiconosVector& nodal_forces,
                        Eigen::Ref<siconos::algebra::SiconosVector> forces);

  /** @brief get the list of possible contacting nodes for a given tag on element.
   *  @param[in] contact_entity_tag tag to be matched
   *  @return a list of FENodes (shared ptr)
   */
  std::vector<std::shared_ptr<FENode>> contactingNodes(int contact_entity_tag);

  void display(bool brief) const;
};

}  // namespace siconos::mechanics::fem

#endif  // FINITEELEMENTMODEL_H
