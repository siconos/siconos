/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2023 INRIA.
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
// #include <vector>

#include "BoundaryCondition.hpp"
#include "FENode.hpp"
#include "FiniteElement.hpp"
#include "Material.hpp"
#include "Mesh.hpp"  // MVertex, MElement ...
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "op3x3.h"  // det3x3

#include "SiconosMatrixOp.hpp"  // prod ...

#define DEBUG_STDOUT
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#include "siconos_debug.h"

// #include <stdio.h>

// #include <map>

std::shared_ptr<siconos::mechanics::fem::FENode>
siconos::mechanics::fem::FiniteElementModel::vertexToNode(
    std::shared_ptr<MVertex> v) {
  if (_vertexToNode.find(v) == _vertexToNode.end()) {
    return nullptr;  // the vertex is not in the map
  }
  return _vertexToNode.at(v);
}

unsigned int siconos::mechanics::fem::FiniteElementModel::init() {
  DEBUG_BEGIN("FiniteElement::init()\n");

  //_nodes.resize(_mesh->vertices().size());

  /* Construction of the elements w.r.t type */
  unsigned int dofIdx = 0;
  unsigned int dim = _mesh->dim();
  unsigned int ndof = 0;
  std::size_t num_node = 0;

  std::vector<std::shared_ptr<MElement>> ignored_elements;

  for (auto e : _mesh->elements()) {
    /* contruction of a FE element */
    DEBUG_PRINTF("MElement->num() : %zu\n", e->num());
    std::shared_ptr<FElement> fe;
    if (dim == 2) {
      switch (e->type())  // we follow gmsh convention
      {
        case FiniteElementType::T3:  // 3-node triangle.
        {
          /* We should normally ask the user for the type of element we
           * associate with the type of MElement of gmsh */
          fe = std::make_shared<FElement>(FiniteElementType::T3, 6, e);
          /* the FE element number is equal to the MElement number */
          break;
        }
        default:
          ignored_elements.push_back(e);
          continue;
      }
    } else if (dim == 3) {
      switch (e->type())  // we follow gmsh convention
      {
        case FiniteElementType::TH4:  // 4-node tetra.
        {
          /* We should normally ask the user for the type of element we
           * associate with the type of MElement of gmsh */
          fe = std::make_shared<FElement>(FiniteElementType::TH4, 12, e);
          /* the FE element number is equal to the MElement number */
          break;
        }
        default:
          ignored_elements.push_back(e);
          continue;
      }
    }

    /* ------------- add the FE element in the elements vector */
    _elements.push_back(fe);
    _mElementTOFElement[e] = _elements.back();

    int ndofPerNode = fe->ndofPerNode();

    /* ------------- contruction of  FE nodes */
    for (auto v : e->vertices()) {
      if (_vertexToNode.find(v) ==
          _vertexToNode.end())  // check if the node is already existing
      {
        ndof += ndofPerNode;
        auto dofIndex = std::make_shared<std::vector<std::size_t>>(0);
        for (int d = 0; d < ndofPerNode; d++) dofIndex->push_back(dofIdx++);
        DEBUG_PRINTF(
            "  create node num_mode: %lu with dofIndex.size():%lu for vertex "
            "num = %lu\n",
            num_node, dofIndex->size(), v->num());
        _nodes.push_back(std::make_shared<FENode>(num_node++, v, dofIndex));
        _vertexToNode[v] = _nodes.back();
      } else {
        DEBUG_PRINTF("  node already exists for vertex : %zu \n", v->num());
      }
      // assert(_nodes[num_node]);
      // DEBUG_EXPR(vertexNode.at(v)->display(););
      fe->nodes().push_back(_vertexToNode.at(v));
    }
  }
  std::cout << "Element type not recognised or ignored : [";
  for (auto e : ignored_elements) {
    std::cout << " " << e->num();
  }
  std::cout << "]" << std::endl;

  DEBUG_PRINTF("number of nodes : %lu\n", _nodes.size());
  DEBUG_PRINTF("number of elements : %lu\n", _elements.size());

  DEBUG_END("FiniteElement::init()\n");
  return ndof;
}

void siconos::mechanics::fem::FiniteElementModel::AssembleElementaryMatrix(
    std::shared_ptr<siconos::algebra::SiconosMatrix> M,
    siconos::algebra::SimpleMatrix &Me, FElement &fe) {
  int node1_cnt = 0;
  for (auto &node1 : fe.nodes()) {
    auto &dofIndex1 = *node1->dofIndex();
    int node2_cnt = 0;
    for (auto &node2 : fe.nodes()) {
      auto &dofIndex2 = *node2->dofIndex();
      for (std::size_t i = 0; i < dofIndex1.size(); i++) {
        for (std::size_t j = 0; j < dofIndex2.size(); j++) {
          // DEBUG_PRINTF("i = %i\t j=%i,  Me.getValue(i,j) = %e\n",
          // i+node1_cnt*dofIndex1.size(), j+node2_cnt*dofIndex2.size(),
          // Me.getValue(i,j));

          M->setValue(dofIndex1[i], dofIndex2[j],
                      Me.getValue(i + node1_cnt * dofIndex1.size(),
                                  j + node2_cnt * dofIndex2.size()) +
                          M->getValue(dofIndex1[i], dofIndex2[j]));
        }
      }
      node2_cnt++;
    }
    node1_cnt++;
  }
  // DEBUG_EXPR(M->display());
}

void siconos::mechanics::fem::FiniteElementModel::AssembleElementary_B_Matrix(std::shared_ptr<siconos::algebra::SiconosMatrix> bigB,
                                                                                      siconos::algebra::SimpleMatrix& Be, FElement& fe, int elem_cnt)
{
    int node1_cnt =0;
    for(std::shared_ptr<FENode> node1 : fe.nodes())
    {
        auto& dofIndex1 = *node1->dofIndex();
        for(int i=0;i<3;i++)
            for(int j=0;j<dofIndex1.size();j++){
                bigB->setValue(3*elem_cnt+i, dofIndex1[j], Be.getValue(i,j+node1_cnt*dofIndex1.size())+ bigB->getValue(3*elem_cnt+i, dofIndex1[j]));
            }
        node1_cnt++;
    }
}

void siconos::mechanics::fem::FiniteElementModel::AssembleElementary_S_Matrix(std::shared_ptr<siconos::algebra::SiconosMatrix> S,
                                                                                      siconos::algebra::SimpleMatrix& Se, FElement& fe, int elem_cnt)
{
    int dimStress = _mesh->dim()*(_mesh->dim()+1)/2;
    int node1_cnt =0;
    for(std::shared_ptr<FENode> node1 : fe.nodes())
    {
        auto& dofIndex1 = *node1->dofIndex();
        for(int i=0;i<dimStress;i++)
            for(int j=0;j<dimStress;j++){
                S->setValue(dimStress*elem_cnt+i, dimStress*elem_cnt+j, Se.getValue(i,j)+ S->getValue(dimStress*elem_cnt+i, dimStress*elem_cnt+j));
            }
        node1_cnt++;
    }
}


void siconos::mechanics::fem::FiniteElementModel::computeElementaryMassMatrix(
    siconos::algebra::SimpleMatrix &Me, FElement &fe, double massDensity) {
  DEBUG_BEGIN(
      "siconos::mechanics::fem::FiniteElementModel::"
      "computeElementaryMassMatrix(siconos::"
      "algebra::SimpleMatrix& Me, FElement& fe, double massDensity )\n");

  Me.zero();

  auto &nodes = fe.nodes();
  int nnodes = nodes.size();

  int dim = _mesh->dim();

  std::vector<double> N(nnodes);
  std::vector<double> Nksi(nnodes);  // only in 2D
  std::vector<double> Neta(nnodes);
  std::vector<double> Nzeta(nnodes);

  std::vector<double> J(dim * dim);

  /** We perform integration by summing over the gauss points
   * this could be simplified by explicit formulae
   */

  int integrationOrder = 2;
  for (const auto &gp : fe.GaussPoints(integrationOrder)) {
    if (_mesh->dim() == 2 and
        fe.family() == FiniteElementFamily::isoparametric) {
      // Compute shape function and derivatives of shape function
      double gp_eta = gp[0];
      double gp_ksi = gp[1];
      double gp_w = gp[2];
      fe.shapeFunctionIso2D(gp_eta, gp_ksi, N, Nksi, Neta);
      // Compute element determinant
      for (int i = 0; i < 4; i++) J[i] = 0.0;
      for (int n = 0; n < nnodes; n++) {
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
      for (int i = 0; i < nnodes; i++) {
        for (int j = 0; j < nnodes; j++) {
          // DEBUG_PRINTF(" N[%i] = %e\t N[%i] = %e\t entry = %e \n", i, N[i],
          // j, N[j], coeff* N[i]*N[j]);
          Me.setValue(i * 2, j * 2,
                      coeff * N[i] * N[j] + Me.getValue(i * 2, j * 2));
          Me.setValue(i * 2 + 1, j * 2 + 1,
                      coeff * N[i] * N[j] + Me.getValue(i * 2 + 1, j * 2 + 1));
        }
      }
    } else if (_mesh->dim() == 3 and
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
      for (int n = 0; n < nnodes; n++) {
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

      double coeff = gp_w * massDensity * detJ /
                     6.0;  // we divide again by 6.0 since the reference element
                           // has volume equal to 1/6.0

      DEBUG_EXPR(std::cout << "coeff: " << coeff << std::endl;);
      /* M += (coeff * Nt N)*/
      for (int i = 0; i < nnodes; i++) {
        for (int j = 0; j < nnodes; j++) {
          // DEBUG_PRINTF(" N[%i] = %e\t N[%i] = %e\t entry = %e \n", i, N[i],
          // j, N[j], coeff* N[i]*N[j]);
          Me.setValue(i * 3, j * 3,
                      coeff * N[i] * N[j] + Me.getValue(i * 3, j * 3));
          Me.setValue(i * 3 + 1, j * 3 + 1,
                      coeff * N[i] * N[j] + Me.getValue(i * 3 + 1, j * 3 + 1));
          Me.setValue(i * 3 + 2, j * 3 + 2,
                      coeff * N[i] * N[j] +
                          Me.getValue(i * 3 + 2,
                                      j * 3 + 2));  // to be checked carefully
        }
      }
    }
  }

  DEBUG_END(
      "siconos::mechanics::fem::FiniteElementModel::"
      "computeElementaryMassMatrix(siconos::"
      "algebra::SimpleMatrix& Me, FElement& fe, double massDensity )\n");
}

void siconos::mechanics::fem::FiniteElementModel::computeMassMatrix(
    std::shared_ptr<siconos::algebra::SiconosMatrix> M,
    std::map<unsigned int, std::shared_ptr<Material>> &mat) {
  DEBUG_BEGIN(
      "siconos::mechanics::fem::FiniteElementModel::computeMassMatrix(std::"
      "shared_ptr<siconos:"
      ":algebra::SiconosMatrix> M, double massDensity )\n");
  M->zero();

  /* loop over the elements */
  for (auto &fe : elements()) {
    unsigned int ndofElement = fe->ndof();
    auto Me = std::make_shared<siconos::algebra::SimpleMatrix>(ndofElement,
                                                               ndofElement);
    // to be optimized if all the element are similar
    double massDensity = mat[fe->mElement()->tags(0)]->massDensity();
    computeElementaryMassMatrix(*Me, *fe, massDensity);
    AssembleElementaryMatrix(M, *Me, *fe);
  }
  DEBUG_END(
      "siconos::mechanics::fem::FiniteElementModel::computeMassMatrix(std::"
      "shared_ptr<siconos:"
      ":algebra::SiconosMatrix M>, double massDensity )\n");
}

void siconos::mechanics::fem::FiniteElementModel::
    computeElementaryStiffnessMatrix_direct(
        siconos::algebra::SimpleMatrix &Ke, FElement &fe,
        std::shared_ptr<siconos::algebra::SimpleMatrix> D, double thickness) {
  DEBUG_BEGIN(
      "siconos::mechanics::fem::FiniteElementModel::"
      "computeElementaryStiffnessMatrix_direct("
      "siconos::algebra::SimpleMatrix& Ke, FElement& fe, Material& mat  )\n");

  Ke.zero();

  // Compute element determinant
  int ndof = fe.ndof();

  auto &nodes = fe.nodes();

  /** We perform integration by summing over the gauss points
   * this could be simplified by explicit formulae
   */

  if (_mesh->dim() == 2 and
      fe.family() == FiniteElementFamily::isoparametric)  // Ugly
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
    auto B = std::make_shared<siconos::algebra::SimpleMatrix>(3, ndof);

    B->setValue(0, 0, -y32);
    B->setValue(0, 2, y31);
    B->setValue(0, 4, -y21);

    B->setValue(1, 1, x32);
    B->setValue(1, 3, -x31);
    B->setValue(1, 5, x21);

    B->setValue(2, 0, x32);
    B->setValue(2, 1, -y32);
    B->setValue(2, 2, -x31);
    B->setValue(2, 3, y31);
    B->setValue(2, 4, x21);
    B->setValue(2, 5, -y21);

    *B = 1.0 / twoA * *B;

    DEBUG_EXPR(B->display());

    // Compte BT D B
    auto DB = std::make_shared<siconos::algebra::SimpleMatrix>(3, ndof);
    siconos::algebra::prod(*D, *B, *DB, true);
    auto BT = std::make_shared<siconos::algebra::SimpleMatrix>(ndof, 3);
    BT->trans(*B);
    auto BTDB = std::make_shared<siconos::algebra::SimpleMatrix>(ndof, ndof);
    siconos::algebra::prod(*BT, *DB, *BTDB, true);

    Ke = (twoA / 2.0 * thickness * *BTDB);

  } else if (_mesh->dim() == 3 and
             fe.family() == FiniteElementFamily::isoparametric)  // Ugly
  {
    /* Construct the B matrix (its form is consistent with the choice
     * of the representation of strain) */
    auto B = std::make_shared<siconos::algebra::SimpleMatrix>(6, ndof);

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

    double sixV = x21 * (y31 * z41 - y41 * z31) +
                  y21 * (x41 * z31 - x31 * z41) + z21 * (x31 * y41 - x41 * y31);

    DEBUG_PRINTF("V = %e\n", sixV / 6.0);

    B->setValue(0, 0, a1);
    B->setValue(0, 3, a2);
    B->setValue(0, 6, a3);
    B->setValue(0, 9, a4);

    B->setValue(1, 1, b1);
    B->setValue(1, 4, b2);
    B->setValue(1, 7, b3);
    B->setValue(1, 10, b4);

    B->setValue(2, 2, c1);
    B->setValue(2, 5, c2);
    B->setValue(2, 8, c3);
    B->setValue(2, 11, c4);

    B->setValue(3, 0, b1);
    B->setValue(3, 1, a1);
    B->setValue(3, 3, b2);
    B->setValue(3, 4, a2);
    B->setValue(3, 6, b3);
    B->setValue(3, 7, a3);
    B->setValue(3, 9, b4);
    B->setValue(3, 10, a4);

    B->setValue(4, 1, c1);
    B->setValue(4, 2, b1);
    B->setValue(4, 4, c2);
    B->setValue(4, 5, b2);
    B->setValue(4, 7, c3);
    B->setValue(4, 8, b3);
    B->setValue(4, 10, c4);
    B->setValue(4, 11, b4);

    B->setValue(5, 0, c1);
    B->setValue(5, 2, a1);
    B->setValue(5, 3, c2);
    B->setValue(5, 5, a2);
    B->setValue(5, 6, c3);
    B->setValue(5, 8, a3);
    B->setValue(5, 9, c4);
    B->setValue(5, 11, a4);

    DEBUG_EXPR(B->display(););
    *B = 1. / sixV * *B;

    // Compte BT D B
    auto DB = std::make_shared<siconos::algebra::SimpleMatrix>(6, ndof);
    prod(*D, *B, *DB, true);
    auto BT = std::make_shared<siconos::algebra::SimpleMatrix>(ndof, 6);
    BT->trans(*B);
    auto BTDB = std::make_shared<siconos::algebra::SimpleMatrix>(ndof, ndof);
    prod(*BT, *DB, *BTDB, true);
    DEBUG_EXPR(BTDB->display(););

    Ke = sixV / 6.0 * *BTDB;
  }

  DEBUG_EXPR(Ke.display(););

  DEBUG_END(
      "siconos::mechanics::fem::FiniteElementModel::"
      "computeElementaryStiffnessMatrix_direct("
      "siconos::algebra::SimpleMatrix& Ke, FElement& fe, Material& mat  )\n");
}

void siconos::mechanics::fem::FiniteElementModel::
    computeElementaryStiffnessMatrix(
        siconos::algebra::SimpleMatrix &Ke, FElement &fe,
        std::shared_ptr<siconos::algebra::SimpleMatrix> D, double thickness) {
  DEBUG_BEGIN(
      "siconos::mechanics::fem::FiniteElementModel::"
      "computeElementaryStiffnessMatrix(siconos::"
      "algebra::SimpleMatrix& Ke, FElement& fe, Material& mat  )\n");

  Ke.zero();
  // Compute element determinant
  int ndof = fe.ndof();
  auto &nodes = fe.nodes();
  int nnodes = nodes.size();

  int dim = _mesh->dim();
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
  for (const auto &gp : fe.GaussPoints(integrationOrder)) {
    if (_mesh->dim() == 2 and
        fe.family() == FiniteElementFamily::isoparametric)  // Ugly
    {
      // Compute shape function and derivatives of shape function
      double gp_eta = gp[0];
      double gp_ksi = gp[1];
      double gp_w = gp[2];
      fe.shapeFunctionIso2D(gp_eta, gp_ksi, N, Nksi, Neta);
      // Compute element determinant
      for (int i = 0; i < 4; i++) J[i] = 0.0;
      for (int n = 0; n < nnodes; n++) {
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
      for (int n = 0; n < nnodes; n++) {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n,
        // Neta[n]);
        Nx[n] = Jinv[0] * Nksi[n] + Jinv[1] * Neta[n];
        Ny[n] = Jinv[2] * Nksi[n] + Jinv[3] * Neta[n];
        // DEBUG_PRINTF(" Nx[%i] = %e\t Ny[%i] = %e\n", n, Nx[n], n, Ny[n]);
      }

      // Construct the B matrix (its form is consistent with the choice of the
      // representation of strain)
      auto B = std::make_shared<siconos::algebra::SimpleMatrix>(3, ndof);
      B->zero();
      for (int n = 0; n < nnodes; n++) {
        B->setValue(0, 2 * n, Nx[n]);
        B->setValue(1, 2 * n, 0.0);
        B->setValue(2, 2 * n, Ny[n]);
        B->setValue(0, 2 * n + 1, 0.0);
        B->setValue(1, 2 * n + 1, Ny[n]);
        B->setValue(2, 2 * n + 1, Nx[n]);
      }

      // Compte BT D B
      auto DB = std::make_shared<siconos::algebra::SimpleMatrix>(3, ndof);
      prod(*D, *B, *DB, true);
      auto BT = std::make_shared<siconos::algebra::SimpleMatrix>(ndof, 3);
      BT->trans(*B);
      auto BTDB = std::make_shared<siconos::algebra::SimpleMatrix>(ndof, ndof);
      prod(*BT, *DB, *BTDB, true);

      double coeff = gp_w * detJ * thickness;

      Ke += (coeff * *BTDB);

      // // check with direct computation (see IFEM Chap 15 Felippa)
      // std::shared_ptr<siconos::algebra::SimpleMatrix> Ke_direct =
      // std::make_shared<siconos::algebra::SimpleMatrix>(ndof,ndof);
      // computeElementaryStiffnessMatrix_direct(*Ke_direct, fe, D, thickness );
      // std::cout << "diff " <<   (*Ke_direct- Ke).normInf() << std::endl;
    } else if (_mesh->dim() == 3 and
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
      for (int n = 0; n < nnodes; n++) {
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
      for (int n = 0; n < nnodes; n++) {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n,
        // Neta[n]);
        Nx[n] = Jinv[0] * Nksi[n] + Jinv[1] * Neta[n] + Jinv[2] * Nzeta[n];
        Ny[n] = Jinv[3] * Nksi[n] + Jinv[4] * Neta[n] + Jinv[5] * Nzeta[n];
        Nz[n] = Jinv[6] * Nksi[n] + Jinv[7] * Neta[n] + Jinv[8] * Nzeta[n];
        // DEBUG_PRINTF(" Nx[%i] = %e\t Ny[%i] = %e\n", n, Nx[n], n, Ny[n]);
      }

      /* Construct the B matrix (its form is consistent with the choice
       * of the representation of strain) */
      auto B = std::make_shared<siconos::algebra::SimpleMatrix>(6, ndof);
      B->zero();
      for (int n = 0; n < nnodes; n++) {
        B->setValue(0, 3 * n, Nx[n]);
        B->setValue(1, 3 * n, 0.0);
        B->setValue(2, 3 * n, 0.0);
        B->setValue(3, 3 * n, Ny[n]);
        B->setValue(4, 3 * n, 0.0);
        B->setValue(5, 3 * n, Nz[n]);

        B->setValue(0, 3 * n + 1, 0.0);
        B->setValue(1, 3 * n + 1, Ny[n]);
        B->setValue(2, 3 * n + 1, 0.0);
        B->setValue(3, 3 * n + 1, Nx[n]);
        B->setValue(4, 3 * n + 1, Nz[n]);
        B->setValue(5, 3 * n + 1, 0.0);

        B->setValue(0, 3 * n + 2, 0.0);
        B->setValue(1, 3 * n + 2, 0.0);
        B->setValue(2, 3 * n + 2, Nz[n]);
        B->setValue(3, 3 * n + 2, 0.0);
        B->setValue(4, 3 * n + 2, Ny[n]);
        B->setValue(5, 3 * n + 2, Nx[n]);
      }
      DEBUG_EXPR(B->display(););
      // Compte BT D B
      auto DB = std::make_shared<siconos::algebra::SimpleMatrix>(6, ndof);
      prod(*D, *B, *DB, true);
      auto BT = std::make_shared<siconos::algebra::SimpleMatrix>(ndof, 6);
      BT->trans(*B);
      auto BTDB = std::make_shared<siconos::algebra::SimpleMatrix>(ndof, ndof);
      prod(*BT, *DB, *BTDB, true);
      DEBUG_EXPR(BTDB->display(););

      double coeff = 0.0;
      coeff = gp_w * detJ / 6.0;  // we divide again by 6.0 since the reference
                                  // element has volume equal to 1/6.0
      Ke += (coeff * *BTDB);

      // // check with direct computation (see AFEM Chap 16 Felippa)
      // auto Ke_direct =
      // std::make_shared<siconos::algebra::SimpleMatrix>(ndof,ndof);
      // computeElementaryStiffnessMatrix_direct(*Ke_direct, fe, D, thickness );
      // std::cout << "diff " <<   (*Ke_direct- Ke).normInf() << std::endl;
    }
  }
  DEBUG_EXPR(Ke.display(););

  DEBUG_END(
      "siconos::mechanics::fem::FiniteElementModel::"
      "computeElementaryStiffnessMatrix(siconos::"
      "algebra::SimpleMatrix& Ke, FElement& fe, Material& mat  )\n");
}

void siconos::mechanics::fem::FiniteElementModel::computeStiffnessMatrix(
    std::shared_ptr<siconos::algebra::SiconosMatrix> K,
    std::map<unsigned int, std::shared_ptr<Material>> &materials) {
  DEBUG_BEGIN(
      "siconos::mechanics::fem::FiniteElementModel::computeStiffnessMatrix(std:"
      ":shared_ptr<"
      "siconos::algebra::SiconosMatrix K, Material>& mat )\n");
  K->zero();

  // We compute first the D matrix. Warning:  to be adpated if several
  // materials.
  std::shared_ptr<siconos::algebra::SimpleMatrix> D;

  /* loop over the elements */
  for (std::shared_ptr<FElement> fe : elements()) {
    Material &mat = *(materials[fe->mElement()->tags(0)]);
    if (_mesh->dim() == 2) {
      D = std::make_shared<siconos::algebra::SimpleMatrix>(3, 3);
      double E = mat.elasticYoungModulus();
      double nu = mat.poissonCoefficient();

      if (mat.analysisType2D() == AnalysisType2D::plane_strain) {
        double coef = E / ((1 + nu) * (1 - 2. * nu));
        (*D)(0, 0) = coef * (1. - nu);
        (*D)(0, 1) = coef * nu;
        (*D)(0, 2) = 0.0;

        (*D)(1, 0) = (*D)(0, 1);
        (*D)(1, 1) = (*D)(0, 0);
        (*D)(1, 2) = 0.0;

        (*D)(2, 0) = 0.0;
        (*D)(2, 1) = 0.0;
        (*D)(2, 2) = 0.5 * coef * (1.0 - 2 * nu);
      } else if (mat.analysisType2D() == AnalysisType2D::plane_stress) {
        double coef = E / (1 - nu * nu);
        (*D)(0, 0) = coef;
        (*D)(0, 1) = coef * nu;
        (*D)(0, 2) = 0.0;

        (*D)(1, 0) = (*D)(0, 1);
        (*D)(1, 1) = (*D)(0, 0);
        (*D)(1, 2) = 0.0;

        (*D)(2, 0) = 0.0;
        (*D)(2, 1) = 0.0;
        (*D)(2, 2) = 0.5 * coef * (1.0 - nu);
      } else
        THROW_EXCEPTION(
            "siconos::mechanics::fem::FiniteElementModel::"
            "computeStiffnessMatrix. Other type "
            "of analysis not yet implemented");

      DEBUG_EXPR(D->display(););
    } else if (_mesh->dim() == 3) {
      // Compute 3D elastic tensor.
      D = std::make_shared<siconos::algebra::SimpleMatrix>(6, 6);
      D->zero();
      double E = mat.elasticYoungModulus();
      double nu = mat.poissonCoefficient();
      double coef = E / ((1 + nu) * (1 - 2. * nu));

      (*D)(0, 0) = coef * (1. - nu);
      (*D)(0, 1) = coef * nu;
      (*D)(0, 2) = coef * nu;

      (*D)(1, 1) = (*D)(0, 0);

      (*D)(1, 0) = (*D)(0, 1);
      (*D)(1, 2) = coef * nu;

      (*D)(2, 2) = (*D)(0, 0);

      (*D)(2, 0) = (*D)(0, 2);
      (*D)(2, 1) = (*D)(1, 2);

      (*D)(3, 3) = coef * (1. - 2. * nu) / 2.;
      (*D)(4, 4) = coef * (1. - 2. * nu) / 2.;
      (*D)(5, 5) = coef * (1. - 2. * nu) / 2.;
    }

    unsigned int ndofElement = fe->ndof();
    auto Ke = std::make_shared<siconos::algebra::SimpleMatrix>(
        ndofElement,
        ndofElement);  // to be optimized if all the element are similar
    computeElementaryStiffnessMatrix(*Ke, *fe, D, mat.thickness());
    AssembleElementaryMatrix(K, *Ke, *fe);
  }
  // DEBUG_EXPR(K->display(););
  DEBUG_END(
      "siconos::mechanics::fem::FiniteElementModel::computeStiffnessMatrix(std:"
      ":shared_ptr<"
      "siconos::algebra::SiconosMatrix K, Material& mat )\n");
}

void siconos::mechanics::fem::FiniteElementModel::computeElementaryBMatrix_direct(FElement& fe, siconos::algebra::SimpleMatrix& B)
{
    int ndof = fe.ndof();
    int order = fe.order();

    std::vector<std::shared_ptr<FENode>> & nodes= fe.nodes();
    int nnodes= nodes.size();
    int dim = _mesh->dim();
    // Direct computation without Gauss Integration
    double x1 = nodes[0]->x();
    double x2 = nodes[1]->x();
    double x3 = nodes[2]->x();

    double y1 = nodes[0]->y();
    double y2 = nodes[1]->y();
    double y3 = nodes[2]->y();

    double x21 = x2-x1;
    double x31 = x3-x1;
    double x32 = x3-x2;

    double y21 = y2-y1;
    double y31 = y3-y1;
    double y32 = y3-y2;

    double twoA =
      x2*y3-x3*y2 +
      x3*y1-x1*y3 +
      x1*y2-x2*y1;
    DEBUG_PRINTF("twoA = %e\n", twoA);
//    std::shared_ptr<SimpleMatrix> B = std::make_shared<SimpleMatrix>(3,ndof);

    B.setValue(0,0, - y32);
    B.setValue(0,2,   y31);
    B.setValue(0,4, - y21);

    B.setValue(1,1,   x32);
    B.setValue(1,3, - x31);
    B.setValue(1,5,   x21);

    B.setValue(2,0,   x32);
    B.setValue(2,1, - y32);
    B.setValue(2,2, - x31);
    B.setValue(2,3,   y31);
    B.setValue(2,4,   x21);
    B.setValue(2,5, - y21);

    B = 1.0/twoA * B;

}


void siconos::mechanics::fem::FiniteElementModel::computeBMatrix(
  std::shared_ptr<siconos::algebra::SiconosMatrix> B,
  std::map<unsigned int, 	std::shared_ptr<Material> > & mat)
{
  DEBUG_BEGIN("siconos::mechanics::fem::native::FiniteElementModel::computeMassMatrix(std::shared_ptr<SiconosMatrix> M, double massDensity )\n");
  int elem_cnt=0;
  /* loop over the elements */
  for(std::shared_ptr<FElement> fe : elements())
  {
      int dimStress = fe->ndofPerNode()*(fe->ndofPerNode()+1)/2;
      std::shared_ptr<siconos::algebra::SimpleMatrix> Be = std::make_shared<siconos::algebra::SimpleMatrix>(dimStress, fe->ndof());
      computeElementaryBMatrix_direct(*fe, *Be);
      AssembleElementary_B_Matrix(B, *Be, *fe, elem_cnt);
      elem_cnt++;
  }
  DEBUG_END("siconos::mechanics::fem::native::FiniteElementModel::computeMassMatrix(std::shared_ptr<SiconosMatrix M>, double massDensity )\n");
}

void siconos::mechanics::fem::FiniteElementModel::computeSMatrix(
  std::shared_ptr<siconos::algebra::SiconosMatrix> S,
  std::map<unsigned int, 	std::shared_ptr<Material> > & mat)
{
  DEBUG_BEGIN("siconos::mechanics::fem::native::FiniteElementModel::computeMassMatrix(std::shared_ptr<SiconosMatrix> M, double massDensity )\n");
  std::cout << "Start Compute S matrix ! " << std::endl;
  int elem_cnt=0;
  int dimStress = _mesh->dim()*(_mesh->dim()+1)/2;
  Material & mate= *(mat[0]);
  std::shared_ptr<siconos::algebra::SimpleMatrix> D, Dinv;
  D = std::make_shared<siconos::algebra::SimpleMatrix>(dimStress,dimStress);
  Dinv = std::make_shared<siconos::algebra::SimpleMatrix>(dimStress,dimStress);
  std::cout << "D initialized ! " << std::endl;
  double E;
  double nu;
  double coef;
//  (*D)(0,0) = coef*(1.-nu);
//  (*D)(0,1) = coef*nu;
//  (*D)(0,2) = 0.0;

//  (*D)(1,0) = (*D)(0,1);
//  (*D)(1,1) = (*D)(0,0);
//  (*D)(1,2) = 0.0;

//  (*D)(2,0) = 0.0;
//  (*D)(2,1) = 0.0;
//  (*D)(2,2) = 0.5*coef*(1.0 - 2* nu);
//  std::cout << "before D Display ! " << std::endl;
//  D->display();
//  D->PLUFactorizationInPlace();
//  std::cout << "After PLU inverse: " << std::endl;
//  D->display();
  /* loop over the elements */
  for(std::shared_ptr<FElement> fe : elements())
  {
      Material & mate= *(mat[fe->mElement()->tags(0)]);
      E = mate.elasticYoungModulus();
      nu =  mate.poissonCoefficient();
      coef = E/((1+nu)*(1-2.*nu));
      (*D)(0,0) = coef*(1.-nu);
      (*D)(0,1) = coef*nu;
      (*D)(0,2) = 0.0;

      (*D)(1,0) = (*D)(0,1);
      (*D)(1,1) = (*D)(0,0);
      (*D)(1,2) = 0.0;

      (*D)(2,0) = 0.0;
      (*D)(2,1) = 0.0;
      (*D)(2,2) = 0.5*coef*(1.0 - 2* nu);
      std::cout << "before D Display ! " << (*D)(0,0) << std::endl;
      D->display();
      (*Dinv)(0,0) = (*D)(0,0)/((*D)(0,0)*(*D)(0,0) - (*D)(0,1)*(*D)(0,1));
      (*Dinv)(0,1) = (*D)(0,1)/((*D)(0,1)*(*D)(0,1) - (*D)(0,0)*(*D)(0,0));
      (*Dinv)(0,2) = 0.0;

      (*Dinv)(1,0) = (*Dinv)(0,1);
      (*Dinv)(1,1) = (*Dinv)(0,0);
      (*Dinv)(1,2) = 0.0;

      (*Dinv)(2,0) = 0.0;
      (*Dinv)(2,1) = 0.0;
      (*Dinv)(2,2) = 1.0/(*D)(2,2);
      std::cout << "After PLU inverse: " << std::endl;
      Dinv->display();

//      std::shared_ptr<siconos::algebra::SimpleMatrix> Se = std::make_shared<siconos::algebra::SimpleMatrix>(dimStress, dimStress);
      AssembleElementary_S_Matrix(S, *Dinv, *fe, elem_cnt);
      elem_cnt++;
  }
  DEBUG_END("siconos::mechanics::fem::native::FiniteElementModel::computeMassMatrix(std::shared_ptr<SiconosMatrix M>, double massDensity )\n");
}



void siconos::mechanics::fem::FiniteElementModel::
    applyDirichletBoundaryConditions(
        int physical_entity_tag,
        std::shared_ptr<std::vector<int>> node_dof_index,
        std::shared_ptr<siconos::modeling::BoundaryCondition>
            boundaryConditions) {
  assert(boundaryConditions);
  for (auto &e : _mesh->elements()) {
    if (e->tags(0) == physical_entity_tag) {
      for (auto &v : e->vertices()) {
        auto n = _vertexToNode[v];
        auto n_dof_index = n->dofIndex();
        for (const auto &i : *node_dof_index) {
          boundaryConditions->appendIndex((*n_dof_index)[i]);
        }
      }
    }
  }
}

void siconos::mechanics::fem::FiniteElementModel::applyNodalForces(
    int physical_entity_tag,
    std::shared_ptr<siconos::algebra::SiconosVector> nodal_forces,
    std::shared_ptr<siconos::algebra::SiconosVector> forces) {
  auto f_index = std::make_shared<std::vector<int>>(0);
  for (auto &e : _mesh->elements()) {
    if (e->tags(0) == physical_entity_tag) {
      for (auto &v : e->vertices()) {
        auto n = _vertexToNode[v];
        // check if the node is already existing
        if (find(f_index->begin(), f_index->end(), n->num()) ==
            f_index->end()) {
          std::cout << "Apply nodal force on node number " << n->num()
                    << " vertex number " << v->num() << std::endl;
          f_index->push_back(n->num());
          auto n_dof_index = n->dofIndex();
          for (unsigned int i = 0; i < nodal_forces->size(); i++) {
            forces->setValue((*n_dof_index)[i], (*nodal_forces)(i));
          }
        }
      }
    }
  }
}

std::shared_ptr<std::list<std::shared_ptr<siconos::mechanics::fem::FENode>>>
siconos::mechanics::fem::FiniteElementModel::contactingNodes(
    int contact_entity_tag) {
  DEBUG_BEGIN(
      "siconos::mechanics::fem::FiniteElementModel::applyNodalForces(int "
      "physical_entity_tag, "
      "std::shared_ptr<siconos::algebra::SiconosVector> nodal_forces, "
      "std::shared_ptr<siconos::algebra::SiconosVector> forces)\n");

  // std::shared_ptr<unisgned int> f_index  = std::make_shared<IndexInt>(0);
  auto contacting_nodes =
      std::make_shared<std::list<std::shared_ptr<FENode>>>();
  ;

  for (auto &e : _mesh->elements()) {
    if (e->tags(0) == contact_entity_tag) {
      for (auto &v : e->vertices()) {
        auto n = _vertexToNode[v];

        if (find(contacting_nodes->begin(), contacting_nodes->end(), n) ==
            contacting_nodes->end())  // check if the node is already existing
        {
          std::cout << "node number in contact zone " << n->num()
                    << " vertex number " << v->num() << std::endl;
          contacting_nodes->push_back(n);
        }
      }
    }
  }

  DEBUG_END(
      "siconos::mechanics::fem::FiniteElementModel::applyNodalForces(int "
      "physical_entity_tag, "
      "std::shared_ptr<siconos::algebra::SiconosVector> nodal_forces, "
      "std::shared_ptr<siconos::algebra::SiconosVector> forces)\n");
  return contacting_nodes;
}

void siconos::mechanics::fem::FiniteElementModel::display(bool brief) const {
  std::cout << "===== FiniteElementModel display ===== " << std::endl;
  std::cout << "- numberOfNodes : " << _nodes.size() << std::endl;
  std::cout << "- numberOfElements : " << _elements.size() << std::endl;
  for (const auto &f : _elements) {
    f->display();
  }
}
