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

#include <list>
#include <map>
#include <memory>
#include <vector>

namespace siconos::algebra {
class SimpleMatrix;
class SiconosVector;
class SiconosMatrix;
}  // namespace siconos::algebra

namespace siconos::modeling {
class BoundaryCondition;
}

namespace siconos::mechanics::fem {

class FENode;
class FElement;
class Mesh;
class MVertex;
class MElement;
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
  std::shared_ptr<Mesh> _mesh{nullptr};

  /** nodes */
  std::vector<std::shared_ptr<FENode>> _nodes = {};

  /** elements */
  std::vector<std::shared_ptr<FElement>> _elements = {};

  /** vertex to node map **/
  std::map<std::shared_ptr<MVertex>, std::shared_ptr<FENode>> _vertexToNode;

  /** MElement to FElement map **/
  std::map<std::shared_ptr<MElement>, std::shared_ptr<FElement>> _mElementTOFElement;

  /** Rule of five */
  FiniteElementModel() = delete;
  FiniteElementModel(FiniteElementModel&) = delete;
  FiniteElementModel& operator=(const FiniteElementModel&) = delete;
  FiniteElementModel(FiniteElementModel&&) = delete;
  FiniteElementModel& operator=(FiniteElementModel&&) = delete;

 public:
  /** Constructor
      \param m the mesh
      \param ndof dof number
      \param e mesh element
   */
  FiniteElementModel(std::shared_ptr<Mesh> m) : _mesh(m) {};

  ~FiniteElementModel() noexcept = default;

  auto mesh() { return _mesh; }

  auto& elements() { return _elements; }

  auto& nodes() { return _nodes; }

  std::shared_ptr<FENode> vertexToNode(std::shared_ptr<MVertex> v);

  /* create the FEM model from the mesh and the element type
   * \return the number of dof */
  unsigned int init();

  /* Assembly method for elemetary matrix */
  void AssembleElementaryMatrix(std::shared_ptr<siconos::algebra::SiconosMatrix> M,
                                siconos::algebra::SimpleMatrix& Me, FElement& fe);

  /* Assembly method specific for elementary fem B matrix */
  void AssembleElementary_B_Matrix(std::shared_ptr<siconos::algebra::SiconosMatrix> M,
                                     siconos::algebra::SimpleMatrix& Be, FElement& fe, int elemCnt);

  void AssembleElementary_S_Matrix(std::shared_ptr<siconos::algebra::SiconosMatrix> S,
                              siconos::algebra::SimpleMatrix& Se, FElement& fe, int elem_cnt);

  /** compute Mass Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeMassMatrix(std::shared_ptr<siconos::algebra::SiconosMatrix>,
                         std::map<unsigned int, std::shared_ptr<Material>>& mat);

  /** compute elementary Mass Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeElementaryMassMatrix(siconos::algebra::SimpleMatrix& Me, FElement& fe,
                                   double massDensity);

  void computeBeamElementaryMassMatrix_direct(siconos::algebra::SimpleMatrix &Me,
                                              FElement &fe, std::map<unsigned int, std::shared_ptr<Material>> &materials);

  /** compute Stiffness Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeStiffnessMatrix(std::shared_ptr<siconos::algebra::SiconosMatrix>,
                              std::map<unsigned int, std::shared_ptr<Material>>& mat);

  /** compute elementary Stiffness Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeElementaryStiffnessMatrix(siconos::algebra::SimpleMatrix& Me, FElement& fe,
                                        std::shared_ptr<siconos::algebra::SimpleMatrix> D,
                                        double thickness);

  /** compute elementary Stiffness Matrix with a direct method
   * for linear element
   **/
  void computeElementaryStiffnessMatrix_direct(
      siconos::algebra::SimpleMatrix& Me, FElement& fe,
      std::shared_ptr<siconos::algebra::SimpleMatrix> D, double thickness);

  void computeBeamElementaryStiffnessMatrix_direct(siconos::algebra::SimpleMatrix &Ke, FElement &fe, std::map<unsigned int, std::shared_ptr<Material>> &materials);

  void computeElementary_B_Matrix(FElement& fe, siconos::algebra::SimpleMatrix& B, double length);

  void computeBeamElementaryBMatrix_direct(FElement& fe, siconos::algebra::SimpleMatrix& Be,  std::map<unsigned int, std::shared_ptr<Material>> &materials);

  void computeElementaryBMatrix_direct(FElement& fe, siconos::algebra::SimpleMatrix& Be, double thickness);

  void computeBMatrix(std::shared_ptr<siconos::algebra::SiconosMatrix> B,
    std::map<unsigned int, 	std::shared_ptr<Material> > & mat);

  void computeSMatrix(std::shared_ptr<siconos::algebra::SiconosMatrix> S,
                      std::map<unsigned int, 	std::shared_ptr<Material> > & mat);
  void computeBeam_S_Matrix(std::shared_ptr<siconos::algebra::SiconosMatrix> S,
                            std::map<unsigned int, 	std::shared_ptr<Material> > & mat);

  /** apply Dirichlet Boundary conditions for a given tag on element.
   **/
  void applyDirichletBoundaryConditions(
      int physical_entity_tag, std::shared_ptr<std::vector<int>> node_dof_index,
      std::shared_ptr<siconos::modeling::BoundaryCondition> _boundaryCondition);
  void applyUniformDirichletBoundaryConditions(
      int physical_entity_tag, std::shared_ptr<std::vector<int>> node_dof_index,
      std::shared_ptr<siconos::modeling::BoundaryCondition> _boundaryCondition, double imposedVelocity);

  /** apply Neuman Boundary conditions (nodal forces) for a given tag on
   *element.
   **/
  void applyNodalForces(int physical_entity_tag,
                        std::shared_ptr<siconos::algebra::SiconosVector> nodal_forces,
                        std::shared_ptr<siconos::algebra::SiconosVector> forces);
  /** get the list of possible contacting nodes for a given tag on element.
   **/
  std::shared_ptr<std::list<std::shared_ptr<FENode>>> contactingNodes(int contact_entity_tag);

  void display(bool brief) const;
};

}  // namespace siconos::mechanics::fem

#endif  // FINITEELEMENTMODEL_H
