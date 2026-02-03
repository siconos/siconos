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

#include "FiniteElementModel.hpp"

#include <iostream>
#include <memory>

#include "BoundaryCondition.hpp"
#include "FENode.hpp"
#include "FETypes.hpp"
#include "FiniteElement.hpp"
#include "Material.hpp"
#include "Mesh.hpp"  // MeshVertex, MeshElement ...
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Tools.hpp"
#include "op3x3.h"  // det3x3
#define DEBUG_STDOUT
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#include "siconos_debug.h"

std::shared_ptr<siconos::mechanics::fem::FENode>
siconos::mechanics::fem::FiniteElementModel::vertexToNode(
    std::shared_ptr<MeshVertex> v) const {
  if (vertexToNode_.find(v) == vertexToNode_.end()) {
    return nullptr;  // the vertex is not in the map
  }
  return vertexToNode_.at(v);
}

siconos::algebra::Index siconos::mechanics::fem::FiniteElementModel::init() {
  DEBUG_BEGIN("FiniteElement::init()\n");

  // nodes_.resize(mesh_->vertices().size());

  siconos::algebra::Index dofIdx = 0;
  siconos::algebra::Index ndof = 0;
  std::size_t num_node = 0;

  std::vector<std::shared_ptr<MeshElement>> ignored_elements;

  // Parses the elements (MeshElement) of the mesh to create FiniteElements
  for (auto elem : mesh_->elements()) {
    DEBUG_PRINTF("MeshElement->num() : %zu\n", elem->num());

    if (is_valid_element(elem->type(), mesh_->dim())) {
      // Add nodes of the current element into the global vector of nodes
      auto ndofPerNode = number_of_dof_per_node(elem->type());
      std::vector<std::shared_ptr<FENode>> fenodes;
      for (auto v : elem->vertices()) {
        // only if the node has not already been registered
        if (vertexToNode_.find(v) == vertexToNode_.end()) {
          // Total number of dofs
          ndof += ndofPerNode;
          // Global index of the nodes of the current element
          std::vector<std::size_t> dofIndex(ndofPerNode);
          for (int d = 0; d < ndofPerNode; d++) dofIndex[d] = dofIdx++;
          DEBUG_PRINTF(
              "  create node num_mode: %lu with dofIndex.size() = %lu, for vertex "
              "num = %lu\n",
              num_node, dofIndex.size(), v->num());
          nodes_.push_back(std::make_shared<FENode>(num_node++, v, dofIndex));
          vertexToNode_[v] = nodes_.back();
        } else {
          DEBUG_PRINTF("  node already exists for vertex : %zu \n", v->num());
        }
        // Add node into the list of the current finite element
        fenodes.push_back(vertexToNode_.at(v));
      }
      // Creates the finite element
      elements_.push_back(std::make_shared<FiniteElement>(elem, fenodes));
      mElementTOFiniteElement_[elem] = elements_.back();
    } else {
      ignored_elements.push_back(elem);
      continue;
    }
  }

  std::cout << "Element type not recognised or ignored : [";
  for (auto e : ignored_elements) {
    std::cout << " " << e->num() << "(type: " << tools::enum_to_string(e->type()) << "),";
  }
  std::cout << "]" << std::endl;

  DEBUG_PRINTF("number of nodes : %lu\n", nodes_.size());
  DEBUG_PRINTF("number of elements : %lu\n", elements_.size());

  DEBUG_END("FiniteElement::init()\n");
  return ndof;
}

void siconos::mechanics::fem::FiniteElementModel::assembleElementaryMatrix(
    siconos::algebra::SiconosSparseMatrix& M, const siconos::algebra::SiconosDenseMatrix& Me,
    const FiniteElement& fe) {
  // FP: We use coeffRef to insert element in M. Insert should be used at first call.
  // Reminder:
  // - use insert if the component does not exist
  // - use coeffRef to update the value of an existing component (will insert if it does not
  // exist but less perf.)
  size_t node1_cnt = 0;
  for (const auto& node1 : fe.nodes()) {
    auto dofIndex1 = node1->global_dof_index();
    size_t node2_cnt = 0;
    for (const auto& node2 : fe.nodes()) {
      auto dofIndex2 = node2->global_dof_index();
      for (std::size_t i = 0; i < dofIndex1.size(); i++) {
        for (std::size_t j = 0; j < dofIndex2.size(); j++) {
          M.coeffRef(dofIndex1[i], dofIndex2[j]) +=
              Me(i + node1_cnt * dofIndex1.size(), j + node2_cnt * dofIndex2.size());
        }
      }
      node2_cnt++;
    }
    node1_cnt++;
  }
  // DEBUG_EXPR(M->display());
}

void siconos::mechanics::fem::FiniteElementModel::assembleElementary_B_Matrix(
    siconos::algebra::SiconosSparseMatrix& Bmatrix,
    const siconos::algebra::SiconosDenseMatrix& Be, const FiniteElement& fe, size_t elem_cnt) {
  size_t node1_cnt = 0;

  auto dim_stress = fe.dimStress();
  for (const auto& node1 : fe.nodes()) {
    auto dofIndex1 = node1->global_dof_index();
    for (size_t i = 0; i < dofIndex1.size(); i++)
      for (size_t j = 0; j < dofIndex1.size(); j++) {
        Bmatrix.coeffRef(dim_stress * elem_cnt + i, dofIndex1[j]) +=
            Be(i, j + node1_cnt * dofIndex1.size());
      }
    node1_cnt++;
  }
}

void siconos::mechanics::fem::FiniteElementModel::assembleElementary_S_Matrix(
    siconos::algebra::SiconosSparseMatrix& S, const siconos::algebra::SiconosDenseMatrix& Se,
    const FiniteElement& fe, size_t elem_cnt) {
  auto dim_stress = fe.dimStress();
  for (siconos::algebra::Index i = 0; i < dim_stress; i++)
    for (siconos::algebra::Index j = 0; j < dim_stress; j++) {
      S.coeffRef(dim_stress * elem_cnt + i, dim_stress * elem_cnt + j) += Se(i, j);
    }
}

void siconos::mechanics::fem::FiniteElementModel::computeElementaryMassMatrix(
    siconos::algebra::SiconosDenseMatrix& Me, const FiniteElement& fe, double massDensity) {
  Me.setZero();

  auto nnodes = fe.nodes().size();

  int dim = mesh_->dim();

  std::vector<double> N(nnodes);
  std::vector<double> Nksi(nnodes);  // only in 2D
  std::vector<double> Neta(nnodes);
  std::vector<double> Nzeta(nnodes);

  std::vector<double> J(dim * dim);

  /** We perform integration by summing over the gauss points
   * this could be simplified by explicit formulae
   */

  const auto& nodes = fe.nodes();
  int integrationOrder = 2;
  for (const auto& gp : fe.GaussPoints(integrationOrder)) {
    if (dim == 2 and fe.family() == FiniteElementFamily::isoparametric) {
      // Compute shape function and derivatives of shape function
      double gp_eta = gp[0];
      double gp_ksi = gp[1];
      double gp_w = gp[2];
      fe.shapeFunctionIso2D(gp_eta, gp_ksi, N, Nksi, Neta);
      // Compute element determinant
      for (int i = 0; i < 4; i++) J[i] = 0.0;
      for (decltype(nnodes) n = 0; n < nnodes; n++) {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n,
        // Neta[n]); DEBUG_PRINTF(" x = %e\t y = %e\n", nodes[n]->_mVertex->x(),
        // nodes[n]->_mVertex->y());
        J[0] += Nksi[n] * nodes[n]->x();
        J[1] += Nksi[n] * nodes[n]->y();
        J[2] += Neta[n] * nodes[n]->x();
        J[3] += Neta[n] * nodes[n]->y();
      }
      double detJ = J[0] * J[3] - J[1] * J[2];
      DEBUG_PRINTF("detJ = %e\n", detJ);

      if (detJ <= 0) THROW_EXCEPTION("computeElementaryMassMatrix. detJ <=0");

      // DEBUG_EXPR(std::cout << "Gauss points : "<< gp_eta << " "  << gp_ksi <<
      // " "  << gp_w
      // << " "   << std::endl;);

      double coeff = gp_w * massDensity * detJ;

      /* M += (coeff * Nt N)*/
      for (decltype(nnodes) i = 0; i < nnodes; i++) {
        for (decltype(nnodes) j = 0; j < nnodes; j++) {
          // DEBUG_PRINTF(" N[%i] = %e\t N[%i] = %e\t entry = %e \n", i, N[i],
          // j, N[j], coeff* N[i]*N[j]);
          Me(i * 2, j * 2) += coeff * N[i] * N[j];
          Me(i * 2 + 1, j * 2 + 1) += coeff * N[i] * N[j];
        }
      }
    } else if (dim == 3 and fe.family() == FiniteElementFamily::isoparametric)  // Ugly
    {
      // Compute shape function and derivatives of shape function
      double gp_eta = gp[0];
      double gp_ksi = gp[1];
      double gp_zeta = gp[2];
      double gp_w = gp[3];
      fe.shapeFunctionIso3D(gp_eta, gp_ksi, gp_zeta, N, Nksi, Neta, Nzeta);
      // Compute element determinant
      for (int i = 0; i < 9; i++) J[i] = 0.0;
      for (decltype(nnodes) n = 0; n < nnodes; n++) {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\t Nzeta[%i] = %e\n", n,
        // Nksi[n], n, Neta[n], n, Nzeta[n]); DEBUG_PRINTF(" x = %e\t y = %e, z
        // = %e\n", nodes[n]->x(), nodes[n]->y(), nodes[n]->z());
        J[0] = J[0] + Nksi[n] * nodes[n]->x();
        J[1] = J[1] + Nksi[n] * nodes[n]->y();
        J[2] = J[2] + Nksi[n] * nodes[n]->z();
        J[3] = J[3] + Neta[n] * nodes[n]->x();
        J[4] = J[4] + Neta[n] * nodes[n]->y();
        J[5] = J[5] + Neta[n] * nodes[n]->z();
        J[6] = J[6] + Nzeta[n] * nodes[n]->x();
        J[7] = J[7] + Nzeta[n] * nodes[n]->y();
        J[8] = J[8] + Nzeta[n] * nodes[n]->z();
      }
      double detJ = det3x3(J.data());
      DEBUG_PRINTF("detJ = %e\n", detJ);
      // DEBUG_EXPR(std::cout << "Gauss points : "<< gp_eta << " "  << gp_ksi <<
      // " "  << gp_w
      // << " "   << std::endl;);

      double coeff =
          gp_w * massDensity * detJ / 6.0;  // we divide again by 6.0 since the reference
                                            // element has volume equal to 1/6.0

      DEBUG_EXPR(std::cout << "coeff: " << coeff << std::endl;);
      /* M += (coeff * Nt N)*/
      for (size_t i = 0; i < nnodes; i++) {
        for (size_t j = 0; j < nnodes; j++) {
          // DEBUG_PRINTF(" N[%i] = %e\t N[%i] = %e\t entry = %e \n", i, N[i],
          // j, N[j], coeff* N[i]*N[j]);
          Me(i * 3, j * 3) += coeff * N[i] * N[j];
          Me(i * 3 + 1, j * 3 + 1) += coeff * N[i] * N[j];
          Me(i * 3 + 2, j * 3 + 2) += coeff * N[i] * N[j];
        }
      }
    }
  }
}

void siconos::mechanics::fem::FiniteElementModel::computeBeamElementaryMassMatrix_direct(
    siconos::algebra::SiconosDenseMatrix& Me, FiniteElement& fe,
    const std::map<int, const Material>& materials) {
  int dim = mesh_->dim();
  double length = fe.length();
  double length2 = length * length;
  double length3 = length2 * length;
  const Material& mat = materials.at(fe.mElement()->tags(0));
  double massDensity = mat.massDensity();
  double A = mat.crossSectionArea();
  double J = mat.secondMomentOfArea();
  auto Te = fe.TeMatrix();
  auto Te_transpose = fe.TeMatrix_transpose();

  if (dim == 2) {
    siconos::algebra::SiconosMatrix66 Me_loc;
    Me_loc.setZero();
    Me_loc(0, 0) = massDensity * A * length * 2 / 6;
    Me_loc(0, 3) = massDensity * A * length * 1 / 6;
    Me_loc(3, 0) = Me_loc(0, 3);

    Me_loc(1, 1) = massDensity * A * length * 156 / 420;
    Me_loc(1, 2) = massDensity * A * length2 * 21 / 420;
    Me_loc(1, 4) = massDensity * A * length * 54 / 420;
    Me_loc(1, 5) = massDensity * A * length2 * (-13) / 420;
    Me_loc(2, 1) = Me_loc(1, 2);
    Me_loc(4, 1) = Me_loc(1, 4);
    Me_loc(5, 1) = Me_loc(1, 5);

    Me_loc(2, 2) = massDensity * A * length3 * 4 / 420;
    Me_loc(2, 4) = massDensity * A * length2 * 13 / 420;
    Me_loc(2, 5) = massDensity * A * length3 * (-3) / 420;
    Me_loc(4, 2) = Me_loc(2, 4);
    Me_loc(5, 2) = Me_loc(2, 5);

    Me_loc(3, 3) = massDensity * A * length * 2 / 6;

    Me_loc(4, 4) = massDensity * A * length * 156 / 420;
    Me_loc(4, 5) = massDensity * A * length2 * (-21) / 420;
    Me_loc(5, 4) = Me_loc(4, 5);

    Me_loc(5, 5) = massDensity * A * length3 * 4 / 420;
    Me.noalias() = Te_transpose * Me_loc * Te;

  } else if (dim == 3) {
    siconos::algebra::SiconosMatrix Me_loc{12, 12};
    // Axial
    Me_loc(0, 0) = massDensity * A * length * 2 / 6;
    Me_loc(0, 6) = massDensity * A * length * 1 / 6;
    Me_loc(6, 0) = Me_loc(0, 3);
    Me_loc(6, 6) = Me_loc(0, 0);

    // Torsional
    Me_loc(3, 3) = massDensity * A * length * J / (3 * A);
    Me_loc(9, 9) = Me_loc(3, 3);
    Me_loc(3, 9) = massDensity * A * length * J / (6 * A);
    Me_loc(9, 3) = Me_loc(3, 9);

    Me_loc(1, 1) = massDensity * A * length * 156 / 420;
    Me_loc(2, 2) = Me_loc(1, 1);

    Me_loc(1, 5) = massDensity * A * length2 * 22 / 420;
    Me_loc(1, 7) = massDensity * A * length * 54 / 420;
    Me_loc(1, 11) = massDensity * A * length2 * (-13) / 420;
    Me_loc(5, 1) = Me_loc(1, 5);
    Me_loc(7, 1) = Me_loc(1, 7);
    Me_loc(11, 1) = Me_loc(1, 11);
    Me_loc(2, 4) = -Me_loc(1, 5);
    Me_loc(2, 8) = Me_loc(1, 7);
    Me_loc(2, 10) = -Me_loc(1, 11);
    Me_loc(4, 2) = Me_loc(2, 4);
    Me_loc(8, 2) = Me_loc(2, 8);
    Me_loc(10, 2) = Me_loc(2, 10);

    Me_loc(4, 4) = massDensity * A * length3 * 4 / 420;
    Me_loc(5, 5) = Me_loc(4, 4);

    Me_loc(4, 8) = -massDensity * A * length2 * 13 / 420;
    Me_loc(4, 10) = massDensity * A * length3 * (-3) / 420;
    Me_loc(8, 4) = Me_loc(4, 8);
    Me_loc(10, 4) = Me_loc(4, 10);
    Me_loc(5, 7) = -Me_loc(4, 8);
    Me_loc(5, 11) = Me_loc(4, 10);
    Me_loc(7, 5) = Me_loc(5, 7);
    Me_loc(11, 5) = Me_loc(5, 11);

    Me_loc(7, 7) = Me_loc(1, 1);
    Me_loc(8, 8) = Me_loc(1, 1);

    Me_loc(7, 11) = -massDensity * A * length2 * 22 / 420;
    Me_loc(11, 7) = Me_loc(7, 11);
    Me_loc(8, 10) = -Me_loc(7, 11);
    Me_loc(10, 8) = Me_loc(8, 10);

    Me_loc(10, 10) = Me_loc(4, 4);
    Me_loc(11, 11) = Me_loc(4, 4);
    Me.noalias() = Te_transpose * Me_loc * Te;
  }
}

void siconos::mechanics::fem::FiniteElementModel::computeMassMatrix(
    siconos::algebra::SiconosSparseMatrix& M, const std::map<int, const Material>& mat) {
  /* loop over the elements */
  for (auto fe : elements()) {
    // to be optimized if all the element are similar
    siconos::algebra::Index ndofElement = fe->ndof();
    siconos::algebra::SiconosDenseMatrix Me{ndofElement, ndofElement};
    // Note FP: consider the case where all elements have the same ndof
    // --> use a buffer for Me to avoid allocation at each call

    if (fe->mElement()->type() == FiniteElementType::B2 ||
        fe->mElement()->type() == FiniteElementType::B3) {
      computeBeamElementaryMassMatrix_direct(Me, *fe, mat);
    } else {
      std::cout << "TATATATATA GGS " << fe->mElement()->tags(0) << "\n";
      double massDensity = mat.at(fe->mElement()->tags(0)).massDensity();
      computeElementaryMassMatrix(Me, *fe, massDensity);
    }
    assembleElementaryMatrix(M, Me, *fe);
  }
  M.makeCompressed();
}

void siconos::mechanics::fem::FiniteElementModel::computeElementaryStiffnessMatrix_direct(
    siconos::algebra::SiconosDenseMatrix& Ke, FiniteElement& fe,
    const siconos::algebra::SiconosDenseMatrix& D, double thickness) {
  DEBUG_BEGIN(
      "siconos::mechanics::fem::FiniteElementModel::computeElementaryStiffnessMatrix_direct("
      ")"
      "\n");

  Ke.setZero();

  // Compute element determinant
  auto ndof = fe.ndof();

  const auto& nodes = fe.nodes();

  /** We perform integration by summing over the gauss points
   * this could be simplified by explicit formulae
   */

  if (mesh_->dim() == 2 and fe.family() == FiniteElementFamily::isoparametric)  // Ugly
  {
    // Direct computation without Gauss Integration
    double x1 = nodes[0]->x();
    double x2 = nodes[1]->x();
    double x3 = nodes[2]->x();

    double y1 = nodes[0]->y();
    double y2 = nodes[1]->y();
    double y3 = nodes[2]->y();

    double x21 = x2 - x1;
    double x31 = x3 - x1;
    double x32 = x3 - x2;

    double y21 = y2 - y1;
    double y31 = y3 - y1;
    double y32 = y3 - y2;

    double twoA = x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3 + x1 * y2 - x2 * y1;
    DEBUG_PRINTF("twoA = %e\n", twoA);
    siconos::algebra::SiconosDenseMatrix B{3, ndof};
    B.setZero();
    B(0, 0) = -y32;
    B(0, 2) = y31;
    B(0, 4) = -y21;

    B(1, 1) = x32;
    B(1, 3) = -x31;
    B(1, 5) = x21;

    B(2, 0) = x32;
    B(2, 1) = -y32;
    B(2, 2) = -x31;
    B(2, 3) = y31;
    B(2, 4) = x21;
    B(2, 5) = -y21;

    B *= 1.0 / twoA;

    DEBUG_EXPR(siconos::algebra::print(B));

    // Compute BT D B
    siconos::algebra::SiconosDenseMatrix DB{3, ndof};
    DB.noalias() = D * B;
    B.transposeInPlace();
    siconos::algebra::SiconosDenseMatrix BTDB{ndof, ndof};
    BTDB.noalias() = B * DB;
    Ke.noalias() = twoA / 2.0 * thickness * BTDB;
  } else if (mesh_->dim() == 3 and fe.family() == FiniteElementFamily::isoparametric)  // Ugly
  {
    /* Construct the B matrix (its form is consistent with the choice
     * of the representation of strain) */
    siconos::algebra::SiconosDenseMatrix B{6, ndof};
    B.setZero();
    // Direct computation without Gauss Integration
    double x1 = nodes[0]->x();
    double x2 = nodes[1]->x();
    double x3 = nodes[2]->x();
    double x4 = nodes[3]->x();

    double y1 = nodes[0]->y();
    double y2 = nodes[1]->y();
    double y3 = nodes[2]->y();
    double y4 = nodes[3]->y();

    double z1 = nodes[0]->z();
    double z2 = nodes[1]->z();
    double z3 = nodes[2]->z();
    double z4 = nodes[3]->z();

    double x21 = x2 - x1;

    double x31 = x3 - x1;
    //    double x32 = x3 - x2;

    double x41 = x4 - x1;
    // double x42 = x4 - x2;
    // double x43 = x4 - x3;

    double y21 = y2 - y1;

    double y31 = y3 - y1;
    double y32 = y3 - y2;

    double y41 = y4 - y1;
    double y42 = y4 - y2;
    double y43 = y4 - y3;

    double z21 = z2 - z1;

    double z31 = z3 - z1;
    double z32 = z3 - z2;

    double z41 = z4 - z1;
    double z42 = z4 - z2;
    double z43 = z4 - z3;

    double a1 = y2 * z43 - y3 * z42 + y4 * z32;
    double a2 = -y1 * z43 + y3 * z41 - y4 * z31;
    double a3 = y1 * z42 - y2 * z41 + y4 * z21;
    double a4 = -y1 * z32 + y2 * z31 - y3 * z21;

    double b1 = -x2 * z43 + x3 * z42 - x4 * z32;
    double b2 = x1 * z43 - x3 * z41 + x4 * z31;
    double b3 = -x1 * z42 + x2 * z41 - x4 * z21;
    double b4 = x1 * z32 - x2 * z31 + x3 * z21;

    double c1 = x2 * y43 - x3 * y42 + x4 * y32;
    double c2 = -x1 * y43 + x3 * y41 - x4 * y31;
    double c3 = x1 * y42 - x2 * y41 + x4 * y21;
    double c4 = -x1 * y32 + x2 * y31 - x3 * y21;

    double sixV = x21 * (y31 * z41 - y41 * z31) + y21 * (x41 * z31 - x31 * z41) +
                  z21 * (x31 * y41 - x41 * y31);

    DEBUG_PRINTF("V = %e\n", sixV / 6.0);

    B(0, 0) = a1;
    B(0, 3) = a2;
    B(0, 6) = a3;
    B(0, 9) = a4;

    B(1, 1) = b1;
    B(1, 4) = b2;
    B(1, 7) = b3;
    B(1, 10) = b4;

    B(2, 2) = c1;
    B(2, 5) = c2;
    B(2, 8) = c3;
    B(2, 11) = c4;

    B(3, 0) = b1;
    B(3, 1) = a1;
    B(3, 3) = b2;
    B(3, 4) = a2;
    B(3, 6) = b3;
    B(3, 7) = a3;
    B(3, 9) = b4;
    B(3, 10) = a4;

    B(4, 1) = c1;
    B(4, 2) = b1;
    B(4, 4) = c2;
    B(4, 5) = b2;
    B(4, 7) = c3;
    B(4, 8) = b3;
    B(4, 10) = c4;
    B(4, 11) = b4;

    B(5, 0) = c1;
    B(5, 2) = a1;
    B(5, 3) = c2;
    B(5, 5) = a2;
    B(5, 6) = c3;
    B(5, 8) = a3;
    B(5, 9) = c4;
    B(5, 11) = a4;

    DEBUG_EXPR(siconos::algebra::print(B););
    B *= 1. / sixV;

    // Compte BT D B
    siconos::algebra::SiconosDenseMatrix DB{6, ndof};
    DB.noalias() = D * B;
    B.transposeInPlace();
    siconos::algebra::SiconosDenseMatrix BTDB{ndof, ndof};
    BTDB.noalias() = B * DB;
    DEBUG_EXPR(siconos::algebra::print(BTDB););
    Ke.noalias() = sixV / 6.0 * BTDB;
  }

  DEBUG_EXPR(siconos::algebra::print(Ke););

  DEBUG_END(
      "siconos::mechanics::fem::FiniteElementModel::"
      "computeElementaryStiffnessMatrix_direct()\n");
}

void siconos::mechanics::fem::FiniteElementModel::computeBeamElementaryStiffnessMatrix_direct(
    siconos::algebra::SiconosDenseMatrix& Ke, FiniteElement& fe,
    const std::map<int, const Material>& materials) {
  int dim = mesh_->dim();
  const Material& mat = materials.at(fe.mElement()->tags(0));
  double E = mat.elasticYoungModulus();
  double I = mat.momentOfInertia();
  double A = mat.crossSectionArea();
  double G = mat.shearModulus();
  double J = mat.secondMomentOfArea();

  double length = fe.length();
  double length2 = length * length;
  double length3 = length2 * length;

  // views, no copy, read-only
  auto Te = fe.TeMatrix();
  auto Te_transpose = fe.TeMatrix_transpose();

  double axial, shearing, bending, torsion;
  double shearOnBend, shearOnshear, bendOnBend;
  axial = E * A / length;
  shearing = E * I * 12 / length3;
  bending = E * I * 4 / length;
  torsion = G * J / length;
  shearOnBend = E * I * 6 / length2;
  shearOnshear = shearing;
  bendOnBend = E * I * 2 / length;

  if (dim == 2) {
    siconos::algebra::SiconosMatrix66 Ke_loc;

    Ke_loc(0, 0) = axial;
    Ke_loc(0, 3) = -axial;
    Ke_loc(3, 0) = Ke_loc(0, 3);

    Ke_loc(1, 1) = shearing;
    Ke_loc(1, 2) = shearOnBend;
    Ke_loc(1, 4) = -shearOnshear;
    Ke_loc(1, 5) = shearOnBend;
    Ke_loc(2, 1) = Ke_loc(1, 2);
    Ke_loc(4, 1) = Ke_loc(1, 4);
    Ke_loc(5, 1) = Ke_loc(1, 5);

    Ke_loc(2, 2) = bending;
    Ke_loc(2, 4) = -shearOnBend;
    Ke_loc(2, 5) = bendOnBend;
    Ke_loc(4, 2) = Ke_loc(2, 4);
    Ke_loc(5, 2) = Ke_loc(2, 5);

    Ke_loc(3, 3) = axial;

    Ke_loc(4, 4) = shearing;
    Ke_loc(4, 5) = -shearOnBend;
    Ke_loc(5, 4) = Ke_loc(4, 5);

    Ke_loc(5, 5) = bending;
    Ke.noalias() = Te_transpose * Ke_loc * Te;

  } else if (dim == 3) {
    siconos::algebra::SiconosMatrix Ke_loc{12, 12};

    Ke_loc(0, 0) = axial;
    Ke_loc(0, 6) = -axial;
    Ke_loc(6, 0) = Ke_loc(0, 6);

    Ke_loc(1, 1) = shearing;
    Ke_loc(1, 5) = shearOnBend;
    Ke_loc(1, 7) = -shearOnshear;
    Ke_loc(1, 11) = shearOnBend;
    Ke_loc(2, 2) = shearing;
    Ke_loc(2, 4) = -shearOnBend;
    Ke_loc(2, 8) = -shearOnshear;
    Ke_loc(2, 10) = -shearOnBend;
    Ke_loc(5, 1) = Ke_loc(1, 5);
    Ke_loc(7, 1) = Ke_loc(1, 7);
    Ke_loc(11, 1) = Ke_loc(1, 11);
    Ke_loc(4, 2) = Ke_loc(2, 4);
    Ke_loc(8, 2) = Ke_loc(2, 8);
    Ke_loc(10, 2) = Ke_loc(2, 10);

    Ke_loc(3, 3) = torsion;
    Ke_loc(3, 9) = -torsion;
    Ke_loc(9, 3) = Ke_loc(3, 9);

    Ke_loc(4, 4) = bending;
    Ke_loc(4, 8) = shearOnBend;
    Ke_loc(4, 10) = bendOnBend;
    Ke_loc(8, 4) = Ke_loc(4, 8);
    Ke_loc(10, 4) = Ke_loc(4, 10);
    Ke_loc(5, 5) = bending;
    Ke_loc(5, 7) = -shearOnBend;
    Ke_loc(5, 11) = bendOnBend;
    Ke_loc(7, 5) = Ke_loc(5, 7);
    Ke_loc(11, 5) = Ke_loc(5, 11);

    Ke_loc(6, 6) = axial;

    Ke_loc(7, 7) = shearing;
    Ke_loc(7, 11) = -shearOnBend;
    Ke_loc(11, 7) = Ke_loc(7, 11);
    Ke_loc(8, 8) = shearing;
    Ke_loc(8, 10) = shearOnBend;
    Ke_loc(10, 8) = Ke_loc(8, 10);

    Ke_loc(9, 9) = torsion;

    Ke_loc(10, 10) = bending;
    Ke_loc(11, 11) = bending;
    Ke.noalias() = Te_transpose * Ke_loc * Te;
  }
}

void siconos::mechanics::fem::FiniteElementModel::computeElementaryStiffnessMatrix(
    siconos::algebra::SiconosDenseMatrix& Ke, FiniteElement& fe,
    const siconos::algebra::SiconosDenseMatrix& D, double thickness) {
  DEBUG_BEGIN(
      "siconos::mechanics::fem::FiniteElementModel::computeElementaryStiffnessMatrix()\n");

  Ke.setZero();
  // Compute element determinant
  auto ndof = fe.ndof();
  const auto& nodes = fe.nodes();
  auto nnodes = fe.nodes().size();

  int dim = mesh_->dim();
  std::vector<double> N(nnodes);
  std::vector<double> Nksi(nnodes);
  std::vector<double> Neta(nnodes);
  std::vector<double> Nzeta(nnodes);
  std::vector<double> Nx(nnodes);
  std::vector<double> Ny(nnodes);
  std::vector<double> Nz(nnodes);
  std::vector<double> J(dim * dim);
  std::vector<double> Jinv(dim * dim);

  // double b[3];

  /** We perform integration by summing over the gauss points
   * this could be simplified by explicit formulae
   */

  int integrationOrder = 1;
  for (const auto& gp : fe.GaussPoints(integrationOrder)) {
    if (mesh_->dim() == 2 and fe.family() == FiniteElementFamily::isoparametric)  // Ugly
    {
      // Compute shape function and derivatives of shape function
      double gp_eta = gp[0];
      double gp_ksi = gp[1];
      double gp_w = gp[2];
      fe.shapeFunctionIso2D(gp_eta, gp_ksi, N, Nksi, Neta);
      // Compute element determinant
      for (size_t i = 0; i < 4; i++) J[i] = 0.0;
      for (size_t n = 0; n < nnodes; n++) {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n,
        // Neta[n]); DEBUG_PRINTF(" x = %e\t y = %e\n", nodes[n]->x(),
        // nodes[n]->y());
        J[0] = J[0] + Nksi[n] * nodes[n]->x();
        J[1] = J[1] + Nksi[n] * nodes[n]->y();
        J[2] = J[2] + Neta[n] * nodes[n]->x();
        J[3] = J[3] + Neta[n] * nodes[n]->y();
      }
      double detJ = J[0] * J[3] - J[1] * J[2];
      DEBUG_PRINTF("detJ = %e\n", detJ);

      // compute inverse of the Jacobian
      Jinv[0] = J[3] / detJ;
      Jinv[1] = -J[1] / detJ;
      Jinv[2] = -J[2] / detJ;
      Jinv[3] = J[0] / detJ;

      // Compute the derivative w.r.t x and y of the shape function
      for (size_t n = 0; n < nnodes; n++) {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n,
        // Neta[n]);
        Nx[n] = Jinv[0] * Nksi[n] + Jinv[1] * Neta[n];
        Ny[n] = Jinv[2] * Nksi[n] + Jinv[3] * Neta[n];
        // DEBUG_PRINTF(" Nx[%i] = %e\t Ny[%i] = %e\n", n, Nx[n], n, Ny[n]);
      }

      // Construct the B matrix (its form is consistent with the choice of the
      // representation of strain)
      siconos::algebra::SiconosDenseMatrix B{3, ndof};
      B.setZero();
      for (size_t n = 0; n < nnodes; n++) {
        B(0, 2 * n) = Nx[n];
        B(1, 2 * n) = 0.0;
        B(2, 2 * n) = Ny[n];
        B(0, 2 * n + 1) = 0.0;
        B(1, 2 * n + 1) = Ny[n];
        B(2, 2 * n + 1) = Nx[n];
      }

      // Compte BT D B
      siconos::algebra::SiconosDenseMatrix DB{3, ndof};
      DB.noalias() = D * B;

      B.transposeInPlace();

      siconos::algebra::SiconosDenseMatrix BTDB{ndof, ndof};
      BTDB.noalias() = B * DB;
      double coeff = gp_w * detJ * thickness;
      Ke += (coeff * BTDB);

      // // check with direct computation (see IFEM Chap 15 Felippa)
      // std::shared_ptr<siconos::algebra::SiconosDenseMatrix> Ke_direct =
      // std::make_shared<siconos::algebra::SiconosDenseMatrix>(ndof,ndof);
      // computeElementaryStiffnessMatrix_direct(*Ke_direct, fe, D, thickness );
      // std::cout << "diff " <<   siconos::algebra::normInf(*Ke_direct- Ke) << std::endl;
    } else if (mesh_->dim() == 3 and
               fe.family() == FiniteElementFamily::isoparametric)  // Ugly
    {
      // Compute shape function and derivatives of shape function
      double gp_eta = gp[0];
      double gp_ksi = gp[1];
      double gp_zeta = gp[2];
      double gp_w = gp[3];
      fe.shapeFunctionIso3D(gp_eta, gp_ksi, gp_zeta, N, Nksi, Neta, Nzeta);
      // Compute element determinant
      for (int i = 0; i < 9; i++) J[i] = 0.0;
      for (size_t n = 0; n < nnodes; n++) {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n,
        // Neta[n]); DEBUG_PRINTF(" x = %e\t y = %e\n", nodes[n]->x(),
        // nodes[n]->y());
        J[0] = J[0] + Nksi[n] * nodes[n]->x();
        J[1] = J[1] + Nksi[n] * nodes[n]->y();
        J[2] = J[2] + Nksi[n] * nodes[n]->z();
        J[3] = J[3] + Neta[n] * nodes[n]->x();
        J[4] = J[4] + Neta[n] * nodes[n]->y();
        J[5] = J[5] + Neta[n] * nodes[n]->z();
        J[6] = J[6] + Nzeta[n] * nodes[n]->x();
        J[7] = J[7] + Nzeta[n] * nodes[n]->y();
        J[8] = J[8] + Nzeta[n] * nodes[n]->z();
      }
      double detJ = det3x3(J.data());
      DEBUG_PRINTF("detJ = %e\n", detJ);
      // for (int j = 0; j < 3; j++) {
      //   for (int i = 0; i < 3; i++) b[i] = 0.0;
      //   b[j] = 1.0;
      //   // int info = solv3x3(J, &Jinv[j * 3], b);
      // }

      // Compute the derivative w.r.t x and y of the shape function
      for (size_t n = 0; n < nnodes; n++) {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n,
        // Neta[n]);
        Nx[n] = Jinv[0] * Nksi[n] + Jinv[1] * Neta[n] + Jinv[2] * Nzeta[n];
        Ny[n] = Jinv[3] * Nksi[n] + Jinv[4] * Neta[n] + Jinv[5] * Nzeta[n];
        Nz[n] = Jinv[6] * Nksi[n] + Jinv[7] * Neta[n] + Jinv[8] * Nzeta[n];
        // DEBUG_PRINTF(" Nx[%i] = %e\t Ny[%i] = %e\n", n, Nx[n], n, Ny[n]);
      }

      /* Construct the B matrix (its form is consistent with the choice
       * of the representation of strain) */
      siconos::algebra::SiconosDenseMatrix B{6, ndof};
      B.setZero();
      for (size_t n = 0; n < nnodes; n++) {
        B(0, 3 * n) = Nx[n];
        B(1, 3 * n) = 0.0;
        B(2, 3 * n) = 0.0;
        B(3, 3 * n) = Ny[n];
        B(4, 3 * n) = 0.0;
        B(5, 3 * n) = Nz[n];

        B(0, 3 * n + 1) = 0.0;
        B(1, 3 * n + 1) = Ny[n];
        B(2, 3 * n + 1) = 0.0;
        B(3, 3 * n + 1) = Nx[n];
        B(4, 3 * n + 1) = Nz[n];
        B(5, 3 * n + 1) = 0.0;

        B(0, 3 * n + 2) = 0.0;
        B(1, 3 * n + 2) = 0.0;
        B(2, 3 * n + 2) = Nz[n];
        B(3, 3 * n + 2) = 0.0;
        B(4, 3 * n + 2) = Ny[n];
        B(5, 3 * n + 2) = Nx[n];
      }
      DEBUG_EXPR(siconos::algebra::print(B););
      // Compte BT D B
      siconos::algebra::SiconosDenseMatrix DB = D * B;
      B.transposeInPlace();
      siconos::algebra::SiconosDenseMatrix BTDB = B * DB;
      DEBUG_EXPR(siconos::algebra::print(BTDB););

      double coeff = 0.0;
      coeff = gp_w * detJ / 6.0;  // we divide again by 6.0 since the reference
                                  // element has volume equal to 1/6.0
      Ke += (coeff * BTDB);

      // // check with direct computation (see AFEM Chap 16 Felippa)
      // auto Ke_direct =
      // std::make_shared<siconos::algebra::SiconosDenseMatrix>(ndof,ndof);
      // computeElementaryStiffnessMatrix_direct(*Ke_direct, fe, D, thickness );
      // std::cout << "diff " <<   siconos::algebra::normInf(*Ke_direct- Ke) << std::endl;
    }
  }
  DEBUG_EXPR(siconos::algebra::print(Ke););

  DEBUG_END(
      "siconos::mechanics::fem::FiniteElementModel::computeElementaryStiffnessMatrix()\n");
}

void siconos::mechanics::fem::FiniteElementModel::computeStiffnessMatrix(
    siconos::algebra::SiconosSparseMatrix& K, const std::map<int, const Material>& materials) {
  DEBUG_BEGIN("siconos::mechanics::fem::FiniteElementModel::computeStiffnessMatrix()\n");

  // We start to compute the D matrix. Warning:  to be adpated if several
  // materials.

  /* loop over the elements */
  for (auto fe : elements()) {
    if (fe->type() == FiniteElementType::B2 || fe->type() == FiniteElementType::B3) {
      auto ndofElement = fe->ndof();
      auto Ke =
          std::make_shared<siconos::algebra::SiconosDenseMatrix>(ndofElement, ndofElement);
      computeBeamElementaryStiffnessMatrix_direct(*Ke, *fe, materials);
      assembleElementaryMatrix(K, *Ke, *fe);
    } else {
      const Material& mat = materials.at(fe->mElement()->tags(0));
      auto ndofElement = fe->ndof();
      siconos::algebra::SiconosDenseMatrix Ke{ndofElement, ndofElement};
      // to be optimized if all the element are similar
      double E = mat.elasticYoungModulus();
      double nu = mat.Poisson_s_ratio();

      if (mesh_->dim() == 2) {
        siconos::algebra::SiconosMatrix33 Dmat;
        Dmat.setZero();

        if (mat.analysisType2D() == AnalysisType2D::plane_strain) {
          double coef = E / ((1 + nu) * (1 - 2. * nu));
          Dmat(0, 0) = coef * (1. - nu);
          Dmat(0, 1) = coef * nu;

          Dmat(1, 0) = Dmat(0, 1);
          Dmat(1, 1) = Dmat(0, 0);

          Dmat(2, 2) = 0.5 * coef * (1.0 - 2 * nu);
        } else if (mat.analysisType2D() == AnalysisType2D::plane_stress) {
          double coef = E / (1 - nu * nu);
          Dmat(0, 0) = coef;
          Dmat(0, 1) = coef * nu;

          Dmat(1, 0) = Dmat(0, 1);
          Dmat(1, 1) = Dmat(0, 0);
          Dmat(2, 2) = 0.5 * coef * (1.0 - nu);
        } else
          THROW_EXCEPTION(
              "siconos::mechanics::fem::FiniteElementModel::"
              "computeStiffnessMatrix. Type of analysis not yet implemented");
        computeElementaryStiffnessMatrix(Ke, *fe, Dmat, mat.thickness());
      } else if (mesh_->dim() == 3) {
        // Compute 3D elastic tensor.
        siconos::algebra::SiconosDenseMatrix Dmat{6, 6};
        Dmat.setZero();
        double coef = E / ((1 + nu) * (1 - 2. * nu));

        Dmat(0, 0) = coef * (1. - nu);
        Dmat(0, 1) = coef * nu;
        Dmat(0, 2) = coef * nu;

        Dmat(1, 1) = Dmat(0, 0);

        Dmat(1, 0) = Dmat(0, 1);
        Dmat(1, 2) = coef * nu;

        Dmat(2, 2) = Dmat(0, 0);

        Dmat(2, 0) = Dmat(0, 2);
        Dmat(2, 1) = Dmat(1, 2);

        Dmat(3, 3) = coef * (1. - 2. * nu) / 2.;
        Dmat(4, 4) = coef * (1. - 2. * nu) / 2.;
        Dmat(5, 5) = coef * (1. - 2. * nu) / 2.;
        computeElementaryStiffnessMatrix(Ke, *fe, Dmat, mat.thickness());
      }

      assembleElementaryMatrix(K, Ke, *fe);
    }
  }
  // DEBUG_EXPR(siconos::algebra::print(*K););
  DEBUG_END("siconos::mechanics::fem::FiniteElementModel::computeStiffnessMatrix()\n");
}

void siconos::mechanics::fem::FiniteElementModel::computeElementary_B_Matrix(
    FiniteElement& fe, siconos::algebra::SiconosDenseMatrix& B, double length) {
  auto nnodes = fe.nodes().size();
  B.setZero();
  size_t cpt = 0;
  for (const auto& gp : fe.GaussPoints(3)) {
    if (mesh_->dim() == 1) {
      // Compute shape function and derivatives of shape function
      double gp_eta = gp[0];
      // double gp_w = gp[3];

      for (size_t n = 0; n < nnodes; n++) {
        B(cpt, 0) = 6 * gp_eta / (length * length);
        B(cpt, 1) = (3 * gp_eta - 1) / length;
        B(cpt, 2) = -6 * gp_eta / (length * length);
        B(cpt, 3) = (3 * gp_eta + 1) / length;
      }
      cpt++;
    }
  }
}

void siconos::mechanics::fem::FiniteElementModel::computeBeamElementaryBMatrix_direct(
    siconos::mechanics::fem::FiniteElement& fe, siconos::algebra::SiconosDenseMatrix& B,
    const std::map<int, const Material>& materials) {
  double length = fe.length();
  double length2 = length * length;
  int dim = mesh_->dim();
  // Direct computation without Gauss Integration

  if (dim == 2) {
    B(0, 0) = -1 / length;
    B(0, 3) = 1 / length;

    int integrationOrder = 3, gp_cnt = 1;
    for (const auto& gp : fe.GaussPoints(integrationOrder)) {
      double gp_eta = gp[0];
      double gp_w = gp[1];

      B(gp_cnt, 1) = sqrt(gp_w) * 6 * gp_eta / length2;
      B(gp_cnt, 2) = sqrt(gp_w) * (3 * gp_eta - 1) / length;
      B(gp_cnt, 4) = sqrt(gp_w) * (-6) * gp_eta / length2;
      B(gp_cnt, 5) = sqrt(gp_w) * (3 * gp_eta + 1) / length;
      gp_cnt++;
    }

  } else if (dim == 3) {
    B(0, 0) = -1 / length;  // axial
    B(0, 6) = 1 / length;   // axial

    B(3, 3) = -1 / length;  // twisting
    B(3, 9) = 1 / length;   // twisting

    int integrationOrder = 3, gp_cnt = 1;
    for (const auto& gp : fe.GaussPoints(integrationOrder)) {
      double gp_eta = gp[0];
      double gp_w = gp[1];

      // B(gp_cnt,1) = sqrt(gp_w)*6*gp_eta/length2);
      // B(gp_cnt+3,2) = sqrt(gp_w)*6*gp_eta/length2);
      B(gp_cnt, 2) = -sqrt(gp_w) * 6 * gp_eta / length2;
      B(gp_cnt + 3, 1) = sqrt(gp_w) * 6 * gp_eta / length2;

      B(gp_cnt, 4) = sqrt(gp_w) * (3 * gp_eta - 1) / length;
      B(gp_cnt + 3, 5) = sqrt(gp_w) * (3 * gp_eta - 1) / length;

      // B(gp_cnt,7) = sqrt(gp_w)*(-6)*gp_eta/length2);
      // B(gp_cnt+3,8) = sqrt(gp_w)*(-6)*gp_eta/length2);
      B(gp_cnt, 8) = sqrt(gp_w) * (6) * gp_eta / length2;
      B(gp_cnt + 3, 7) = sqrt(gp_w) * (-6) * gp_eta / length2;

      B(gp_cnt, 10) = sqrt(gp_w) * (3 * gp_eta + 1) / length;
      B(gp_cnt + 3, 11) = sqrt(gp_w) * (3 * gp_eta + 1) / length;

      gp_cnt++;
    }
  }
  // Transfer to global coordinates
  auto Te = fe.TeMatrix();
  B = B * Te;
}

void siconos::mechanics::fem::FiniteElementModel::computeElementaryBMatrix_direct(
    FiniteElement& fe, siconos::algebra::SiconosDenseMatrix& B, double thickness) {
  const auto& nodes = fe.nodes();
  // Direct computation without Gauss Integration
  double x1 = nodes[0]->x();
  double x2 = nodes[1]->x();
  double x3 = nodes[2]->x();

  double y1 = nodes[0]->y();
  double y2 = nodes[1]->y();
  double y3 = nodes[2]->y();

  double x21 = x2 - x1;
  double x31 = x3 - x1;
  double x32 = x3 - x2;

  double y21 = y2 - y1;
  double y31 = y3 - y1;
  double y32 = y3 - y2;

  double twoA = x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3 + x1 * y2 - x2 * y1;
  DEBUG_PRINTF("twoA = %e\n", twoA);
  //    std::shared_ptr<SiconosDenseMatrix> B = std::make_shared<SiconosDenseMatrix>(3,ndof);

  B(0, 0) = -y32;
  B(0, 2) = y31;
  B(0, 4) = -y21;

  B(1, 1) = x32;
  B(1, 3) = -x31;
  B(1, 5) = x21;

  B(2, 0) = x32;
  B(2, 1) = -y32;
  B(2, 2) = -x31;
  B(2, 3) = y31;
  B(2, 4) = x21;
  B(2, 5) = -y21;

  B = 1.0 / twoA * B;
  double coeff = sqrt(twoA / 2.0 * thickness);
  B = coeff * B;
}

void siconos::mechanics::fem::FiniteElementModel::computeBMatrix(
    siconos::algebra::SiconosSparseMatrix& B, const std::map<int, const Material>& mat) {
  DEBUG_BEGIN("siconos::mechanics::fem::native::FiniteElementModel::computeBMatrix(...)\n");
  size_t elem_cnt = 0;
  siconos::algebra::SiconosDenseMatrix Be;
  /* loop over the elements */
  for (std::shared_ptr<FiniteElement> fe : elements()) {
    if (fe->type() == FiniteElementType::B2 || fe->type() == FiniteElementType::B3) {
      double dim = mesh_->dim();
      if (dim == 2)
        Be.resize(3, 6);
      else
        Be.resize(6, 12);
      computeBeamElementaryBMatrix_direct(*fe, Be, mat);
      assembleElementary_B_Matrix(B, Be, *fe, elem_cnt);
      elem_cnt++;
    } else {
      auto dimStress = fe->ndofPerNode() * (fe->ndofPerNode() + 1) / 2;
      Be.resize(dimStress, fe->ndof());
      const Material& material = mat.at(fe->mElement()->tags(0));
      computeElementaryBMatrix_direct(*fe, Be, material.thickness());
      assembleElementary_B_Matrix(B, Be, *fe, elem_cnt);
      elem_cnt++;
    }
  }
  DEBUG_END("siconos::mechanics::fem::native::FiniteElementModel::computeBMatrix(...)\n");
}

void siconos::mechanics::fem::FiniteElementModel::computeElasticityMatrix(
    siconos::algebra::SiconosSparseMatrix& S, const std::map<int, const Material>& materials) {
  DEBUG_BEGIN(
      "siconos::mechanics::fem::native::FiniteElementModel::computeElasticityMatrix(...)\n");
  size_t elem_cnt = 0;

  //  Dinv = std::make_shared<siconos::algebra::SiconosDenseMatrix>(dimStress,dimStress);
  int dim = mesh_->dim();

  if (dim == 2) {
    /* loop over the elements */
    for (auto fe : elements()) {
      const Material& mat = materials.at(fe->mElement()->tags(0));
      auto E = mat.elasticYoungModulus();

      if (fe->type() == FiniteElementType::B2) {
        auto I = mat.momentOfInertia();
        auto A = mat.crossSectionArea();
        auto length = fe->length();
        // D is diagonal
        double c1 = 1. / (E * A * length);
        double c2 = 1. / (0.5 * E * I * length);
        siconos::algebra::SiconosMatrix33Diagonal Dinv{c1, c2, c2};
        assembleElementary_S_Matrix(S, Dinv, *fe, elem_cnt);

      } else {  // Other types of 2D elements
        siconos::algebra::SiconosMatrix33 Dmat;
        Dmat.setZero();
        auto nu = mat.Poisson_s_ratio();
        auto coef = E / ((1 + nu) * (1 - 2. * nu));
        Dmat(0, 0) = coef * (1. - nu);
        Dmat(0, 1) = coef * nu;
        Dmat(1, 0) = Dmat(0, 1);
        Dmat(1, 1) = Dmat(0, 0);
        Dmat(2, 2) = 0.5 * coef * (1.0 - 2 * nu);

        Dmat = Dmat.inverse().eval();  // We need the inverse
        assembleElementary_S_Matrix(S, Dmat, *fe, elem_cnt);
      }
      elem_cnt++;
    }
  } else if (dim == 3) {
    /* loop over the elements */
    for (auto fe : elements()) {
      const Material& mat = materials.at(fe->mElement()->tags(0));
      auto E = mat.elasticYoungModulus();

      if (fe->type() == FiniteElementType::B3) {
        auto G = mat.shearModulus();
        auto J = mat.secondMomentOfArea();
        auto I = mat.momentOfInertia();
        auto A = mat.crossSectionArea();
        auto length = fe->length();
        // D is diagonal
        double c1 = 1. / (E * A * length);
        double c2 = 2. / (E * I * length);
        double c3 = 1. / (G * J * length);
        double c4 = 2. / (E * I * length);
        siconos::algebra::SiconosMatrix66Diagonal Dinv{c1, c2, c2, c3, c4, c4};
        assembleElementary_S_Matrix(S, Dinv, *fe, elem_cnt);

      } else {  // Other types of 3D elements

        THROW_EXCEPTION(
            "computeElasticityMatrix not yet implememented in 3D for non-beam elements");
      }
      elem_cnt++;
    }
  }

  DEBUG_END(
      "siconos::mechanics::fem::native::FiniteElementModel::computeElasticityMatrix(...)\n");
}

void siconos::mechanics::fem::FiniteElementModel::applyDirichletBoundaryConditions(
    int physical_entity_tag,
    const std::vector<siconos::algebra::SiconosVector::Index>& input_dof_index,
    std::shared_ptr<siconos::modeling::BoundaryCondition> boundaryConditions,
    double imposedVelocity) {
  assert(boundaryConditions);
  for (auto elem : mesh_->elements()) {
    if (elem->tags(0) == physical_entity_tag) {
      for (auto& vertex : elem->vertices()) {
        auto node = vertexToNode_[vertex];
        auto n_dof_index = node->global_dof_index();
        for (siconos::algebra::SiconosVector::Index i : input_dof_index) {
          boundaryConditions->appendIndex(n_dof_index[i], imposedVelocity);
        }
      }
    }
  }
}

void siconos::mechanics::fem::FiniteElementModel::applyNodalForces(
    int physical_entity_tag, const siconos::algebra::SiconosVector& nodal_forces,
    Eigen::Ref<siconos::algebra::SiconosVector> forces) {
  std::vector<int> f_index(0);
  for (auto elem : mesh_->elements()) {
    if (elem->tags(0) == physical_entity_tag) {
      for (auto vertex : elem->vertices()) {
        auto node = vertexToNode_[vertex];
        // check if the node has already been handled
        auto found = find(f_index.begin(), f_index.end(), node->num());
        // if not, apply force
        if (found == f_index.end()) {
          std::cout << "Apply nodal force on node number " << node->num() << " vertex number "
                    << vertex->num() << std::endl;
          f_index.push_back(node->num());
          auto n_dof_index = node->global_dof_index();
          for (siconos::algebra::Index i = 0; i < nodal_forces.size(); i++) {
            forces(n_dof_index[i]) = nodal_forces(i);
          }
        }
      }
    }
  }
}

std::vector<std::shared_ptr<siconos::mechanics::fem::FENode>>
siconos::mechanics::fem::FiniteElementModel::contactingNodes(int contact_entity_tag) {
  DEBUG_BEGIN("siconos::mechanics::fem::FiniteElementModel::contactingNodes(...)\n");

  std::vector<std::shared_ptr<FENode>> contacting_nodes = {};

  for (auto elem : mesh_->elements()) {
    if (elem->tags(0) == contact_entity_tag) {
      for (auto v : elem->vertices()) {
        auto node = vertexToNode_[v];

        if (find(contacting_nodes.begin(), contacting_nodes.end(), node) ==
            contacting_nodes.end())  // check if the node is already existing
        {
          std::cout << "node number in contact zone " << node->num() << " vertex number "
                    << v->num() << std::endl;
          contacting_nodes.push_back(node);
        }
      }
    }
  }
  DEBUG_END("siconos::mechanics::fem::FiniteElementModel::contactingNodes(...)\n");
  return contacting_nodes;
}

void siconos::mechanics::fem::FiniteElementModel::display(bool brief) const {
  std::cout << "===== FiniteElementModel display ===== " << std::endl;
  std::cout << "- numberOfNodes : " << nodes_.size() << std::endl;
  std::cout << "- numberOfElements : " << elements_.size() << std::endl;
  for (const auto& f : elements_) {
    f->display();
  }
}
