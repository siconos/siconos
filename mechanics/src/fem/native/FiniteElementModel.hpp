// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
/** 2D elasticity tensor D (sigma = D epsilon, Voigt notation (xx, yy, xy) with the
 * engineering shear strain) of an isotropic material, for its analysis type (plane strain
 * or plane stress). Throws for other analysis types.
 */
siconos::algebra::SiconosMatrix33 elasticityTensor2D(const Material& mat);

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
  std::map<std::shared_ptr<MeshElement>, std::shared_ptr<FiniteElement>>
      mElementTOFiniteElement_;

  /** Rule of five */
  FiniteElementModel() = delete;
  FiniteElementModel(FiniteElementModel&) = delete;
  FiniteElementModel& operator=(const FiniteElementModel&) = delete;
  FiniteElementModel(FiniteElementModel&&) = delete;
  FiniteElementModel& operator=(FiniteElementModel&&) = delete;

  // Helper function to check if the element type is valid (i.e. properly implemented in
  // siconos in the current FiniteElementModel and for the given mesh dimension
  inline bool is_valid_element(FiniteElementType value, int mesh_dimension) {
    static const std::unordered_set<FiniteElementType> valid_2d_values = {
        FiniteElementType::T3, FiniteElementType::B2};
    static const std::unordered_set<FiniteElementType> valid_3d_values = {
        FiniteElementType::TH4, FiniteElementType::B3};
    if (mesh_dimension == 2)
      return valid_2d_values.count(value) > 0;
    else if (mesh_dimension == 3)
      return valid_3d_values.count(value) > 0;
    else
      return 0;
  }

 public:
  /** Constructor
      \param m the mesh
      \param ndof dof number
      \param e mesh element
   */
  FiniteElementModel(std::shared_ptr<Mesh> m) : mesh_(m) {};

  ~FiniteElementModel() noexcept = default;

  auto mesh() const { return mesh_; }

  const std::span<const std::shared_ptr<FiniteElement>> elements() const { return elements_; }

  const std::span<const std::shared_ptr<FENode>> nodes() const { return nodes_; }

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
   * In 2D, the mass density is multiplied by the thickness of the material.
   * @param[out] mass the assembled mass matrix
   * @param[in] mat materials (one per element tag)
   * @param[in] massType consistent (default) or lumped (row sums of the consistent
   *            elementary matrices; throws for beam elements or a non-positive diagonal)
   **/
  void computeMassMatrix(siconos::algebra::SiconosSparseMatrix& mass,
                         const std::map<int, const Material>& mat,
                         MassType massType = MassType::consistent);

  /** compute elementary Mass Matrix
   **/
  void computeElementaryMassMatrix(siconos::algebra::SiconosDenseMatrix&,
                                   const FiniteElement& fe, double massDensity);

  void computeBeamElementaryMassMatrix_direct(siconos::algebra::SiconosDenseMatrix& Me,
                                              FiniteElement& fe,
                                              const std::map<int, const Material>& materials);

  /** compute Stiffness Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeStiffnessMatrix(siconos::algebra::SiconosSparseMatrix&,
                              const std::map<int, const Material>& mat);

  /** compute elementary Stiffness Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeElementaryStiffnessMatrix(siconos::algebra::SiconosDenseMatrix& Me,
                                        FiniteElement& fe,
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

  void computeBeamElementaryBMatrix_direct(FiniteElement& fe,
                                           siconos::algebra::SiconosDenseMatrix& Be,
                                           const std::map<int, const Material>& materials);

  /** @brief strain operator of a T3 element for the stress unknowns of SolidLinearTIDS:
   * scaled by sqrt(Area * thickness), in Mandel notation (exx, eyy, sqrt2 exy)
   **/
  void computeElementaryBMatrix_direct(FiniteElement& fe,
                                       siconos::algebra::SiconosDenseMatrix& Be,
                                       double thickness);

  /** @brief Compute the plain (unscaled) strain-displacement matrix for a T3
   * element, i.e. epsilon = B * u with no additional area/thickness scaling.
   *
   * Unlike computeElementaryBMatrix_direct(), which scales its result by
   * sqrt(Area * thickness) and uses the Mandel notation for the global sparse B/S
   * assembly (K = B^T * S^{-1} * B), this returns the physical strain-displacement matrix
   * in Voigt notation (exx, eyy, gamma_xy), suitable for evaluating element strain/stress
   * from a nodal displacement vector.
   *
   * @param[in] fe the (T3) finite element
   * @param[out] Be the 3x6 strain-displacement matrix
   **/
  void computeElementaryStrainDisplacementMatrix_direct(
      FiniteElement& fe, siconos::algebra::SiconosDenseMatrix& Be);

  void computeBMatrix(siconos::algebra::SiconosSparseMatrix& B,
                      const std::map<int, const Material>& mat);

  /** compute the compliance matrix S (block diagonal, one block per element): D_M^{-1}
   * for the continuum elements, with D_M = T D T, T = diag(1, 1, sqrt2), the elasticity
   * tensor in Mandel notation (D = elasticityTensor2D, Voigt notation, plane strain or plane
   * stress according to the material), the beam compliance for B2/B3. With B scaled by
   * sqrt(area * thickness) (Mandel notation for the continuum elements), K = B^T S^{-1} B.
   **/
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

  /** @brief Get contact segments on a given physical boundary.
   *
   *  @param contact_tag Physical tag identifier (e.g., 4 for "Contact Nodes" in gmsh)
   *  @return Vector of segments, each an array of 4 (2D) or 6 (3D) global DOF indices
   *  @throws If no elements with the given tag are found or element type is unsupported
   */
  std::vector<std::array<siconos::algebra::Index, 4>> getContactSegments(
      int contact_tag) const;

  /** @brief Convert contact vertex coordinates to global DOF indices.
   *
   *  @param contact_coords Nx2 matrix (2D) or Nx3 matrix (3D) of vertex coordinates, ordered
   * along the boundary (CCW)
   *  @return Flattened vector of global DOF indices for the contact boundary
   *  @throws If a coordinate does not match any mesh vertex
   */
  std::vector<siconos::algebra::Index> getContactGlobalIndices(
      const Eigen::Ref<const Eigen::MatrixXd>& contact_coords) const;

  void display(bool brief) const;
};

}  // namespace siconos::mechanics::fem

#endif  // FINITEELEMENTMODEL_H
