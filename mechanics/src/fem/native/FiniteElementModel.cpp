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

#include "FiniteElementModel.hpp"
#include "SimpleMatrix.hpp"
#include "Material.hpp"
#include "SiconosAlgebraProd.hpp"
#include "SimpleMatrixFriends.hpp"
#include "op3X3.h"

// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

#include <stdio.h>
#include <iostream>
#include <map>

unsigned int siconos::mechanics::fem::native::FiniteElementModel::init()
{
  DEBUG_BEGIN("FiniteElement::init()\n");


  //_nodes.resize(_mesh->vertices().size());

  /* Construction of the elements w.r.t type */
  unsigned int numElement=0;
  unsigned int ndof=0;
  unsigned int dofIdx=0;
  unsigned int dim = _mesh->dim();

  unsigned int num_node=0;

  std::vector<MElement*> ignored_elements;

  for(MElement * e : _mesh->elements())
  {
    /* ------------- contruction of an FE element */
    DEBUG_PRINTF("MElement->num() : %zu\n", e->num());
    std::shared_ptr<FElement> fe;
    if(dim ==2)
    {
      switch(e->type()) // we follow gmsh convention
      {
      case 2: // 3-node triangle.
      {
        /* We should normally ask the user for the type of element we associate with
         * the type of MElement of gmsh */
        fe.reset(new FElement(T3, 6, e)); /* the FE element number is equal to the MElement number */
        break;
      }
      default:
        ignored_elements.push_back(e);
        continue;
      }
    }
    else if(dim ==3)
    {
      switch(e->type()) // we follow gmsh convention
      {
      case 4: // 4-node tetra.
      {
        /* We should normally ask the user for the type of element we associate with
         * the type of MElement of gmsh */
        fe.reset(new FElement(TH4, 12, e)); /* the FE element number is equal to the MElement number */
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
    for(MVertex * v : e->vertices())
    {
      if(_vertexToNode.find(v) == _vertexToNode.end())  // check if the node is already existing
      {
        ndof+=ndofPerNode;
        std::shared_ptr<Index> dofIndex = std::make_shared<Index>(0);
        for(int d =0 ; d < ndofPerNode; d++)
          dofIndex->push_back(dofIdx++);
        DEBUG_PRINTF("  create node num_mode: %lu with dofIndex.size():%lu for vertex num = %lu\n", num_node, dofIndex->size(), v->num());
        std::shared_ptr<FENode> n = std::make_shared<FENode>(num_node++, v, dofIndex);
        _nodes.push_back(n);
        _vertexToNode[v] = _nodes.back();
      }
      else
      {
        DEBUG_PRINTF("  node already exists for vertex : %zu \n", v->num());
      }
      //assert(_nodes[num_node]);
      //DEBUG_EXPR(vertexNode.at(v)->display(););
      fe->nodes().push_back(_vertexToNode.at(v));
    }
  }
  std::cout << "Element type not recognised or ignored : [" ;
  for(MElement * e: ignored_elements)
  {
    std::cout << " " << e->num() ;
  }
  std::cout << "]" << std::endl ;


  DEBUG_PRINTF("number of nodes : %i\n", _nodes.size());
  DEBUG_PRINTF("number of elements : %i\n", _elements.size());

  return ndof;
  DEBUG_END("FiniteElement::init()\n");

}
void siconos::mechanics::fem::native::FiniteElementModel::AssembleElementaryMatrix(std::shared_ptr<SiconosMatrix> M,
    SimpleMatrix& Me, FElement& fe)
{
  int node1_cnt =0;
  for(std::shared_ptr<FENode> node1 : fe.nodes())
  {
    Index& dofIndex1 = *node1->dofIndex();
    int node2_cnt =0;
    for(std::shared_ptr<FENode>  node2 : fe.nodes())
    {
      Index& dofIndex2 = *node2->dofIndex();
      for(int i = 0; i < dofIndex1.size() ; i++)
      {
        for(int j = 0; j < dofIndex2.size() ; j++)
        {
          //DEBUG_PRINTF("i = %i\t j=%i,  Me.getValue(i,j) = %e\n", i+node1_cnt*dofIndex1.size(), j+node2_cnt*dofIndex2.size(), Me.getValue(i,j));

          M->setValue(dofIndex1[i], dofIndex2[j], Me.getValue(i+node1_cnt*dofIndex1.size(),j+node2_cnt*dofIndex2.size())+ M->getValue(dofIndex1[i], dofIndex2[j]));
        }
      }
      node2_cnt++;
    }
    node1_cnt++;
  }
  //DEBUG_EXPR(M->display());
}

void siconos::mechanics::fem::native::FiniteElementModel::computeElementaryMassMatrix(SimpleMatrix& Me, FElement& fe, double massDensity)
{
  DEBUG_BEGIN("siconos::mechanics::fem::native::FiniteElementModel::computeElementaryMassMatrix(SimpleMatrix& Me, FElement& fe, double massDensity )\n");


  Me.zero();

  int ndof = fe.ndof();
  int order = fe.order();
  std::vector<std::shared_ptr<FENode>> & nodes= fe.nodes();
  int nnodes= nodes.size();

  int dim = _mesh->dim();
  double * N =(double *)malloc(nnodes *sizeof(double));

  double * Nksi =(double *)malloc(nnodes *sizeof(double)); // only in 2D
  double * Neta =(double *)malloc(nnodes *sizeof(double));
  double * Nzeta =(double *)malloc(nnodes *sizeof(double));

  double * J =(double *)malloc(dim*dim*sizeof(double));

  /** We perform integration by summing over the gauss points
   * this could be simplified by explicit formulae
   */

  int integrationOrder=2;
  for(const double* gp : fe.GaussPoints(integrationOrder))
  {
    if(_mesh->dim() ==2 and fe.family() == ISOPARAMETRIC)
    {
      // Compute shape function and derivatives of shape function
      double  gp_eta= gp[0];
      double  gp_ksi= gp[1];
      double  gp_w= gp[2];
      fe.shapeFunctionIso2D(gp_eta, gp_ksi, N, Nksi, Neta);
      // Compute element determinant
      for(int i =0; i<4; i++) J[i]=0.0;
      for(int n =0; n < nnodes; n++)
      {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n, Neta[n]);
        // DEBUG_PRINTF(" x = %e\t y = %e\n", nodes[n]->_mVertex->x(), nodes[n]->_mVertex->y());
        J[0] = J[0] + Nksi[n]*nodes[n]->_mVertex->x();
        J[1] = J[1] + Nksi[n]*nodes[n]->_mVertex->y();
        J[2] = J[2] + Neta[n]*nodes[n]->_mVertex->x();
        J[3] = J[3] + Neta[n]*nodes[n]->_mVertex->y();
      }
      double detJ = J[0]*J[3] - J[1]*J[2];
      DEBUG_PRINTF("detJ = %e\n", detJ);
      // DEBUG_EXPR(std::cout << "Gauss points : "<< gp_eta << " "  << gp_ksi << " "  << gp_w << " "   << std::endl;);

      double coeff = gp_w * massDensity * detJ;

      /* M += (coeff * Nt N)*/
      for(int i = 0; i < nnodes; i++)
      {
        for(int j= 0; j < nnodes; j++)
        {
          // DEBUG_PRINTF(" N[%i] = %e\t N[%i] = %e\t entry = %e \n", i, N[i], j, N[j], coeff* N[i]*N[j]);
          Me.setValue(i*2,j*2, coeff* N[i]*N[j] + Me.getValue(i*2,j*2));
          Me.setValue(i*2+1,j*2+1, coeff* N[i]*N[j] + Me.getValue(i*2+1,j*2+1));
        }
      }
    }
    else if(_mesh->dim() ==3 and fe.family() == ISOPARAMETRIC)  //Ugly
    {
      // Compute shape function and derivatives of shape function
      double  gp_eta  = gp[0];
      double  gp_ksi  = gp[1];
      double  gp_zeta = gp[2];
      double  gp_w    = gp[3];
      fe.shapeFunctionIso3D(gp_eta, gp_ksi, gp_zeta, N, Nksi, Neta, Nzeta);
      // Compute element determinant
      for(int i =0; i<9; i++) J[i]=0.0;
      for(int n =0; n < nnodes; n++)
      {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\t Nzeta[%i] = %e\n", n, Nksi[n], n, Neta[n], n, Nzeta[n]);
        // DEBUG_PRINTF(" x = %e\t y = %e, z = %e\n", nodes[n]->_mVertex->x(), nodes[n]->_mVertex->y(), nodes[n]->_mVertex->z());
        J[0] = J[0] + Nksi[n]*nodes[n]->_mVertex->x();
        J[1] = J[1] + Nksi[n]*nodes[n]->_mVertex->y();
        J[2] = J[2] + Nksi[n]*nodes[n]->_mVertex->z();
        J[3] = J[3] + Neta[n]*nodes[n]->_mVertex->x();
        J[4] = J[4] + Neta[n]*nodes[n]->_mVertex->y();
        J[5] = J[5] + Neta[n]*nodes[n]->_mVertex->z();
        J[6] = J[6] + Nzeta[n]*nodes[n]->_mVertex->x();
        J[7] = J[7] + Nzeta[n]*nodes[n]->_mVertex->y();
        J[8] = J[8] + Nzeta[n]*nodes[n]->_mVertex->z();
      }
      double detJ = det3x3(J);
      DEBUG_PRINTF("detJ = %e\n", detJ);
      // DEBUG_EXPR(std::cout << "Gauss points : "<< gp_eta << " "  << gp_ksi << " "  << gp_w << " "   << std::endl;);

      double coeff = gp_w * massDensity * detJ / 6.0; // we divide again by 6.0 since the reference element has volume equal to 1/6.0

      DEBUG_EXPR(std::cout << "coeff: "<< coeff << std::endl;);
      /* M += (coeff * Nt N)*/
      for(int i = 0; i < nnodes; i++)
      {
        for(int j= 0; j < nnodes; j++)
        {
          // DEBUG_PRINTF(" N[%i] = %e\t N[%i] = %e\t entry = %e \n", i, N[i], j, N[j], coeff* N[i]*N[j]);
          Me.setValue(i*3,  j*3,   coeff* N[i]*N[j] + Me.getValue(i*3,  j*3));
          Me.setValue(i*3+1,j*3+1, coeff* N[i]*N[j] + Me.getValue(i*3+1,j*3+1));
          Me.setValue(i*3+2,j*3+2, coeff* N[i]*N[j] + Me.getValue(i*3+2,j*3+2)); // to be checked carefully
        }
      }
    }
  }
  free(N);
  free(Nksi);
  free(Neta);
  free(J);

  DEBUG_END("siconos::mechanics::fem::native::FiniteElementModel::computeElementaryMassMatrix(SimpleMatrix& Me, FElement& fe, double massDensity )\n");
}



void siconos::mechanics::fem::native::FiniteElementModel::computeMassMatrix(
  std::shared_ptr<SiconosMatrix> M,
  std::map<unsigned int, 	std::shared_ptr<Material> > & mat)
{
  DEBUG_BEGIN("siconos::mechanics::fem::native::FiniteElementModel::computeMassMatrix(std::shared_ptr<SiconosMatrix> M, double massDensity )\n");
  M->zero();

  /* loop over the elements */
  for(std::shared_ptr<FElement> fe : elements())
  {
    unsigned int ndofElement  = fe->_ndof;
    std::shared_ptr<SimpleMatrix> Me = std::make_shared<SimpleMatrix>(ndofElement,ndofElement); // to be optimized if all the element are similar
    double massDensity = mat[fe->_mElement->tags(0)]->massDensity();
    computeElementaryMassMatrix(*Me, *fe, massDensity);
    AssembleElementaryMatrix(M, *Me, *fe);
  }
  DEBUG_END("siconos::mechanics::fem::native::FiniteElementModel::computeMassMatrix(std::shared_ptr<SiconosMatrix M>, double massDensity )\n");
}





void siconos::mechanics::fem::native::FiniteElementModel::computeElementaryStiffnessMatrix_direct(
  SimpleMatrix& Ke,
  FElement& fe,
  std::shared_ptr<SimpleMatrix> D,
  double thickness)
{
  DEBUG_BEGIN("siconos::mechanics::fem::native::FiniteElementModel::computeElementaryStiffnessMatrix_direct(SimpleMatrix& Ke, FElement& fe, Material& mat  )\n");

  Ke.zero();

  // Compute element determinant
  int ndof = fe.ndof();
  int order = fe.order();

  std::vector<std::shared_ptr<FENode>> & nodes= fe.nodes();
  int nnodes= nodes.size();
  int dim = _mesh->dim();


  /** We perform integration by summing over the gauss points
   * this could be simplified by explicit formulae
   */

  if(_mesh->dim() ==2 and fe.family() == ISOPARAMETRIC)  //Ugly
  {


    // Direct computation without Gauss Integration
    double x1 = nodes[0]->_mVertex->x();
    double x2 = nodes[1]->_mVertex->x();
    double x3 = nodes[2]->_mVertex->x();

    double y1 = nodes[0]->_mVertex->y();
    double y2 = nodes[1]->_mVertex->y();
    double y3 = nodes[2]->_mVertex->y();

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
    std::shared_ptr<SimpleMatrix> B = std::make_shared<SimpleMatrix>(3,ndof);

    B->setValue(0,0, - y32);
    B->setValue(0,2,   y31);
    B->setValue(0,4, - y21);

    B->setValue(1,1,   x32);
    B->setValue(1,3, - x31);
    B->setValue(1,5,   x21);

    B->setValue(2,0,   x32);
    B->setValue(2,1, - y32);
    B->setValue(2,2, - x31);
    B->setValue(2,3,   y31);
    B->setValue(2,4,   x21);
    B->setValue(2,5, - y21);

    *B = 1.0/twoA * *B;

    DEBUG_EXPR(B->display());

    // Compte BT D B
    std::shared_ptr<SimpleMatrix> DB = std::make_shared<SimpleMatrix>(3,ndof);
    prod(*D, *B, *DB, true);
    std::shared_ptr<SimpleMatrix> BT = std::make_shared<SimpleMatrix>(ndof,3);
    BT->trans(*B);
    std::shared_ptr<SimpleMatrix> BTDB = std::make_shared<SimpleMatrix>(ndof,ndof);
    prod(*BT, *DB, *BTDB, true);


    Ke = (twoA/2.0 * thickness*  *BTDB) ;


  }
  else if(_mesh->dim() ==3 and fe.family() == ISOPARAMETRIC)  //Ugly
  {

    /* Construct the B matrix (its form is consistent with the choice
     * of the representation of strain) */
    std::shared_ptr<SimpleMatrix> B = std::make_shared<SimpleMatrix>(6,ndof);

    // Direct computation without Gauss Integration
    double x1 = nodes[0]->_mVertex->x();
    double x2 = nodes[1]->_mVertex->x();
    double x3 = nodes[2]->_mVertex->x();
    double x4 = nodes[3]->_mVertex->x();

    double y1 = nodes[0]->_mVertex->y();
    double y2 = nodes[1]->_mVertex->y();
    double y3 = nodes[2]->_mVertex->y();
    double y4 = nodes[3]->_mVertex->y();

    double z1 = nodes[0]->_mVertex->z();
    double z2 = nodes[1]->_mVertex->z();
    double z3 = nodes[2]->_mVertex->z();
    double z4 = nodes[3]->_mVertex->z();

    double x21 = x2-x1;

    double x31 = x3-x1;
    double x32 = x3-x2;

    double x41 = x4-x1;
    double x42 = x4-x2;
    double x43 = x4-x3;

    double y21 = y2-y1;

    double y31 = y3-y1;
    double y32 = y3-y2;

    double y41 = y4-y1;
    double y42 = y4-y2;
    double y43 = y4-y3;

    double z21 = z2-z1;

    double z31 = z3-z1;
    double z32 = z3-z2;

    double z41 = z4-z1;
    double z42 = z4-z2;
    double z43 = z4-z3;

    double a1 =   y2 * z43 - y3 * z42 + y4 * z32;
    double a2 = - y1 * z43 + y3 * z41 - y4 * z31;
    double a3 =   y1 * z42 - y2 * z41 + y4 * z21;
    double a4 = - y1 * z32 + y2 * z31 - y3 * z21;

    double b1 = - x2 * z43 + x3 * z42 - x4 * z32;
    double b2 =   x1 * z43 - x3 * z41 + x4 * z31;
    double b3 = - x1 * z42 + x2 * z41 - x4 * z21;
    double b4 =   x1 * z32 - x2 * z31 + x3 * z21;

    double c1 =   x2 * y43 - x3 * y42 + x4 * y32;
    double c2 = - x1 * y43 + x3 * y41 - x4 * y31;
    double c3 =   x1 * y42 - x2 * y41 + x4 * y21;
    double c4 = - x1 * y32 + x2 * y31 - x3 * y21;


    double sixV =
      x21*(y31 * z41 - y41* z31) +
      y21*(x41 * z31 - x31* z41) +
      z21*(x31 * y41 - x41* y31);

    DEBUG_PRINTF("V = %e\n", sixV/6.0);

    B->setValue(0, 0,   a1);
    B->setValue(0, 3,   a2);
    B->setValue(0, 6,   a3);
    B->setValue(0, 9,   a4);

    B->setValue(1, 1,   b1);
    B->setValue(1, 4,   b2);
    B->setValue(1, 7,   b3);
    B->setValue(1, 10,  b4);

    B->setValue(2, 2,   c1);
    B->setValue(2, 5,   c2);
    B->setValue(2, 8,   c3);
    B->setValue(2, 11,  c4);

    B->setValue(3, 0,   b1);
    B->setValue(3, 1,   a1);
    B->setValue(3, 3,   b2);
    B->setValue(3, 4,   a2);
    B->setValue(3, 6,   b3);
    B->setValue(3, 7,   a3);
    B->setValue(3, 9,   b4);
    B->setValue(3, 10,  a4);

    B->setValue(4, 1,   c1);
    B->setValue(4, 2,   b1);
    B->setValue(4, 4,   c2);
    B->setValue(4, 5,   b2);
    B->setValue(4, 7,   c3);
    B->setValue(4, 8,   b3);
    B->setValue(4, 10,  c4);
    B->setValue(4, 11,  b4);

    B->setValue(5, 0,   c1);
    B->setValue(5, 2,   a1);
    B->setValue(5, 3,   c2);
    B->setValue(5, 5,   a2);
    B->setValue(5, 6,   c3);
    B->setValue(5, 8,   a3);
    B->setValue(5, 9,   c4);
    B->setValue(5, 11,  a4);

    DEBUG_EXPR(B->display(););
    *B = 1./sixV * *B ;

    // Compte BT D B
    std::shared_ptr<SimpleMatrix> DB = std::make_shared<SimpleMatrix>(6,ndof);
    prod(*D, *B, *DB, true);
    std::shared_ptr<SimpleMatrix> BT = std::make_shared<SimpleMatrix>(ndof,6);
    BT->trans(*B);
    std::shared_ptr<SimpleMatrix> BTDB = std::make_shared<SimpleMatrix>(ndof,ndof);
    prod(*BT, *DB, *BTDB, true);
    DEBUG_EXPR(BTDB->display(););

    Ke =  sixV /6.0 * *BTDB ;

  }

  DEBUG_EXPR(Ke.display(););

  DEBUG_END("siconos::mechanics::fem::native::FiniteElementModel::computeElementaryStiffnessMatrix_direct(SimpleMatrix& Ke, FElement& fe, Material& mat  )\n");
}


void siconos::mechanics::fem::native::FiniteElementModel::computeElementaryStiffnessMatrix(SimpleMatrix& Ke,
    FElement& fe,
    std::shared_ptr<SimpleMatrix> D,
    double thickness)
{
  DEBUG_BEGIN("siconos::mechanics::fem::native::FiniteElementModel::computeElementaryStiffnessMatrix(SimpleMatrix& Ke, FElement& fe, Material& mat  )\n");

  Ke.zero();
  // Compute element determinant
  int ndof = fe.ndof();
  int order = fe.order();
  std::vector<std::shared_ptr<FENode>> & nodes= fe.nodes();
  int nnodes= nodes.size();

  int dim = _mesh->dim();
  double * N =(double *)malloc(nnodes *sizeof(double));

  double * Nksi =(double *)malloc(nnodes *sizeof(double)); // only in 2D
  double * Neta =(double *)malloc(nnodes *sizeof(double));
  double * Nzeta =(double *)malloc(nnodes *sizeof(double));

  double * Nx =(double *)malloc(nnodes*sizeof(double));
  double * Ny =(double *)malloc(nnodes*sizeof(double));
  double * Nz =(double *)malloc(nnodes*sizeof(double));

  double * J =(double *)malloc(dim*dim*sizeof(double));
  double * Jinv =(double *)malloc(dim*dim*sizeof(double));
  double b[3];

  /** We perform integration by summing over the gauss points
   * this could be simplified by explicit formulae
   */

  int integrationOrder=1;
  for(const double* gp : fe.GaussPoints(integrationOrder))
  {
    if(_mesh->dim() ==2 and fe.family() == ISOPARAMETRIC)  //Ugly
    {
      // Compute shape function and derivatives of shape function
      double  gp_eta= gp[0];
      double  gp_ksi= gp[1];
      double  gp_w= gp[2];
      fe.shapeFunctionIso2D(gp_eta, gp_ksi, N, Nksi, Neta);
      // Compute element determinant
      for(int i =0; i<4; i++) J[i]=0.0;
      for(int n =0; n < nnodes; n++)
      {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n, Neta[n]);
        // DEBUG_PRINTF(" x = %e\t y = %e\n", nodes[n]->_mVertex->x(), nodes[n]->_mVertex->y());
        J[0] = J[0] + Nksi[n]*nodes[n]->_mVertex->x();
        J[1] = J[1] + Nksi[n]*nodes[n]->_mVertex->y();
        J[2] = J[2] + Neta[n]*nodes[n]->_mVertex->x();
        J[3] = J[3] + Neta[n]*nodes[n]->_mVertex->y();
      }
      double detJ = J[0]*J[3] - J[1]*J[2];
      DEBUG_PRINTF("detJ = %e\n", detJ);

      // compute inverse of the Jacobian
      Jinv[0] = J[3]/detJ;
      Jinv[1] =-J[1]/detJ;
      Jinv[2] =-J[2]/detJ;
      Jinv[3] = J[0]/detJ;

      // Compute the derivative w.r.t x and y of the shape function
      for(int n =0; n < nnodes; n++)
      {
        //DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n, Neta[n]);
        Nx[n] = Jinv[0] * Nksi[n] + Jinv[1] * Neta[n];
        Ny[n] = Jinv[2] * Nksi[n] + Jinv[3] * Neta[n];
        //DEBUG_PRINTF(" Nx[%i] = %e\t Ny[%i] = %e\n", n, Nx[n], n, Ny[n]);
      }

      // Construct the B matrix (its form is consistent with the choice of the representation of strain)
      std::shared_ptr<SimpleMatrix> B = std::make_shared<SimpleMatrix>(3,ndof);
      B->zero();
      for(int n =0; n < nnodes; n++)
      {
        B->setValue(0, 2*n,   Nx[n]);
        B->setValue(1, 2*n,   0.0);
        B->setValue(2, 2*n,   Ny[n]);
        B->setValue(0, 2*n+1, 0.0);
        B->setValue(1, 2*n+1, Ny[n]);
        B->setValue(2, 2*n+1, Nx[n]);
      }

      // Compte BT D B
      std::shared_ptr<SimpleMatrix> DB = std::make_shared<SimpleMatrix>(3,ndof);
      prod(*D, *B, *DB, true);
      std::shared_ptr<SimpleMatrix> BT = std::make_shared<SimpleMatrix>(ndof,3);
      BT->trans(*B);
      std::shared_ptr<SimpleMatrix> BTDB = std::make_shared<SimpleMatrix>(ndof,ndof);
      prod(*BT, *DB, *BTDB, true);

      double coeff =  gp_w  * detJ *  thickness;

      Ke += (coeff * *BTDB) ;

      // // check with direct computation (see IFEM Chap 15 Felippa)
      // std::shared_ptr<SimpleMatrix> Ke_direct = std::make_shared<SimpleMatrix>(ndof,ndof);
      // computeElementaryStiffnessMatrix_direct(*Ke_direct, fe, D, thickness );
      // std::cout << "diff " <<   (*Ke_direct- Ke).normInf() << std::endl;
    }
    else if(_mesh->dim() ==3 and fe.family() == ISOPARAMETRIC)  //Ugly
    {
      // Compute shape function and derivatives of shape function
      double  gp_eta  = gp[0];
      double  gp_ksi  = gp[1];
      double  gp_zeta = gp[2];
      double  gp_w    = gp[3];
      fe.shapeFunctionIso3D(gp_eta, gp_ksi, gp_zeta, N, Nksi, Neta, Nzeta);
      // Compute element determinant
      for(int i =0; i<9; i++) J[i]=0.0;
      for(int n =0; n < nnodes; n++)
      {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n, Neta[n]);
        // DEBUG_PRINTF(" x = %e\t y = %e\n", nodes[n]->_mVertex->x(), nodes[n]->_mVertex->y());
        J[0] = J[0] + Nksi[n]*nodes[n]->_mVertex->x();
        J[1] = J[1] + Nksi[n]*nodes[n]->_mVertex->y();
        J[2] = J[2] + Nksi[n]*nodes[n]->_mVertex->z();
        J[3] = J[3] + Neta[n]*nodes[n]->_mVertex->x();
        J[4] = J[4] + Neta[n]*nodes[n]->_mVertex->y();
        J[5] = J[5] + Neta[n]*nodes[n]->_mVertex->z();
        J[6] = J[6] + Nzeta[n]*nodes[n]->_mVertex->x();
        J[7] = J[7] + Nzeta[n]*nodes[n]->_mVertex->y();
        J[8] = J[8] + Nzeta[n]*nodes[n]->_mVertex->z();
      }
      double detJ = det3x3(J);
      DEBUG_PRINTF("detJ = %e\n", detJ);
      for(int j =0; j <3 ; j++)
      {
        for(int i =0; i <3 ; i++)  b[i]=0.0;
        b[j]=1.0;
        int info = solv3x3(J, &Jinv[j*3], b);
      }

      // Compute the derivative w.r.t x and y of the shape function
      for(int n =0; n < nnodes; n++)
      {
        //DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n, Neta[n]);
        Nx[n] = Jinv[0] * Nksi[n] + Jinv[1] * Neta[n] + Jinv[2] * Nzeta[n];
        Ny[n] = Jinv[3] * Nksi[n] + Jinv[4] * Neta[n] + Jinv[5] * Nzeta[n];
        Nz[n] = Jinv[6] * Nksi[n] + Jinv[7] * Neta[n] + Jinv[8] * Nzeta[n];
        //DEBUG_PRINTF(" Nx[%i] = %e\t Ny[%i] = %e\n", n, Nx[n], n, Ny[n]);
      }

      /* Construct the B matrix (its form is consistent with the choice
       * of the representation of strain) */
      std::shared_ptr<SimpleMatrix> B = std::make_shared<SimpleMatrix>(6,ndof);
      B->zero();
      for(int n =0; n < nnodes; n++)
      {
        B->setValue(0, 3*n,   Nx[n]);
        B->setValue(1, 3*n,   0.0);
        B->setValue(2, 3*n,   0.0);
        B->setValue(3, 3*n,   Ny[n]);
        B->setValue(4, 3*n,   0.0);
        B->setValue(5, 3*n,   Nz[n]);

        B->setValue(0, 3*n+1, 0.0);
        B->setValue(1, 3*n+1, Ny[n]);
        B->setValue(2, 3*n+1, 0.0);
        B->setValue(3, 3*n+1, Nx[n]);
        B->setValue(4, 3*n+1, Nz[n]);
        B->setValue(5, 3*n+1, 0.0);

        B->setValue(0, 3*n+2, 0.0);
        B->setValue(1, 3*n+2, 0.0);
        B->setValue(2, 3*n+2, Nz[n]);
        B->setValue(3, 3*n+2, 0.0);
        B->setValue(4, 3*n+2, Ny[n]);
        B->setValue(5, 3*n+2, Nx[n]);
      }
      DEBUG_EXPR(B->display(););
      // Compte BT D B
      std::shared_ptr<SimpleMatrix> DB = std::make_shared<SimpleMatrix>(6,ndof);
      prod(*D, *B, *DB, true);
      std::shared_ptr<SimpleMatrix> BT = std::make_shared<SimpleMatrix>(ndof,6);
      BT->trans(*B);
      std::shared_ptr<SimpleMatrix> BTDB = std::make_shared<SimpleMatrix>(ndof,ndof);
      prod(*BT, *DB, *BTDB, true);
      DEBUG_EXPR(BTDB->display(););

      double coeff =0.0;
      coeff = gp_w  * detJ / 6.0; // we divide again by 6.0 since the reference element has volume equal to 1/6.0
      Ke += (coeff * *BTDB) ;

      // // check with direct computation (see AFEM Chap 16 Felippa)
      // std::shared_ptr<SimpleMatrix> Ke_direct = std::make_shared<SimpleMatrix>(ndof,ndof);
      // computeElementaryStiffnessMatrix_direct(*Ke_direct, fe, D, thickness );
      // std::cout << "diff " <<   (*Ke_direct- Ke).normInf() << std::endl;
    }
  }
  DEBUG_EXPR(Ke.display(););
  free(N);
  free(Nksi);
  free(Neta);
  free(J);
  free(Jinv);


  DEBUG_END("siconos::mechanics::fem::native::FiniteElementModel::computeElementaryStiffnessMatrix(SimpleMatrix& Ke, FElement& fe, Material& mat  )\n");
}



void siconos::mechanics::fem::native::FiniteElementModel::computeStiffnessMatrix(
  std::shared_ptr<SiconosMatrix> K, std::map<unsigned int, std::shared_ptr<Material>> & materials)
{
  DEBUG_BEGIN("siconos::mechanics::fem::native::FiniteElementModel::computeStiffnessMatrix(std::shared_ptr<SiconosMatrix K, Material>& mat )\n");
  K->zero();

  // We compute first the D matrix. Warning:  to be adpated if several materials.
  std::shared_ptr<SimpleMatrix> D;

  /* loop over the elements */
  for(std::shared_ptr<FElement> fe : elements())
  {
    Material & mat= *(materials[fe->_mElement->tags(0)]);
    if(_mesh->dim() ==2)
    {
      D = std::make_shared<SimpleMatrix>(3,3);
      double E = mat.elasticYoungModulus();
      double nu =  mat.poissonCoefficient();

      if(mat.analysisType2D() == PLANE_STRAIN)
      {
        double coef = E/((1+nu)*(1-2.*nu));
        (*D)(0,0) = coef*(1.-nu);
        (*D)(0,1) = coef*nu;
        (*D)(0,2) = 0.0;

        (*D)(1,0) = (*D)(0,1);
        (*D)(1,1) = (*D)(0,0);
        (*D)(1,2) = 0.0;

        (*D)(2,0) = 0.0;
        (*D)(2,1) = 0.0;
        (*D)(2,2) = 0.5*coef*(1.0 - 2* nu);
      }
      else if(mat.analysisType2D() == PLANE_STRESS)
      {
        double coef = E/(1-nu*nu);
        (*D)(0,0) = coef;
        (*D)(0,1) = coef*nu;
        (*D)(0,2) = 0.0;

        (*D)(1,0) = (*D)(0,1);
        (*D)(1,1) = (*D)(0,0);
        (*D)(1,2) = 0.0;

        (*D)(2,0) = 0.0;
        (*D)(2,1) = 0.0;
        (*D)(2,2) = 0.5*coef*(1.0 - nu);
      }
      else
        THROW_EXCEPTION("siconos::mechanics::fem::native::FiniteElementModel::computeStiffnessMatrix. Other type of analysis not yet implemented");

      DEBUG_EXPR(D->display(););
    }
    else if(_mesh->dim() ==3)
    {
      //Compute 3D elastic tensor.
      D= std::make_shared<SimpleMatrix>(6,6);
      D->zero();
      double E = mat.elasticYoungModulus();
      double nu =  mat.poissonCoefficient();
      double coef = E/((1+nu)*(1-2.*nu));

      (*D)(0,0) = coef*(1.-nu);
      (*D)(0,1) = coef*nu;
      (*D)(0,2) = coef*nu;

      (*D)(1,1) = (*D)(0,0);

      (*D)(1,0) = (*D)(0,1);
      (*D)(1,2) = coef*nu;

      (*D)(2,2) = (*D)(0,0);

      (*D)(2,0) = (*D)(0,2);
      (*D)(2,1) = (*D)(1,2);

      (*D)(3,3) = coef*(1.-2.*nu)/2.;
      (*D)(4,4) = coef*(1.-2.*nu)/2.;
      (*D)(5,5) = coef*(1.-2.*nu)/2.;
    }


    unsigned int ndofElement  = fe->_ndof;
    std::shared_ptr<SimpleMatrix> Ke = std::make_shared<SimpleMatrix>(ndofElement,ndofElement); // to be optimized if all the element are similar
    computeElementaryStiffnessMatrix(*Ke, *fe, D, mat.thickness());
    AssembleElementaryMatrix(K, *Ke, *fe);
  }
  //DEBUG_EXPR(K->display(););
  DEBUG_END("siconos::mechanics::fem::native::FiniteElementModel::computeStiffnessMatrix(std::shared_ptr<SiconosMatrix K, Material& mat )\n");
}


void siconos::mechanics::fem::native::FiniteElementModel::applyDirichletBoundaryConditions(int physical_entity_tag,
    std::shared_ptr<IndexInt> node_dof_index,
    std::shared_ptr<BoundaryCondition> _boundaryCondition)
{
  DEBUG_BEGIN("siconos::mechanics::fem::native::FiniteElementModel::applyDirichletBoundaryConditions(int physical_entity_tag, std::shared_ptr<IndexInt node_dof_index, std::shared_ptr<BoundaryCondition>> _boundaryCondition)\n");


  std::shared_ptr<IndexInt> bdIndex = _boundaryCondition->velocityIndices();
  std::shared_ptr<SiconosVector> bdPrescribedVelocity = _boundaryCondition->prescribedVelocity();

  int  bdIndex_old_size = bdIndex->size();
  for(MElement * e : _mesh->elements())
  {
    if(e->tags(0) == physical_entity_tag)
    {
      for(MVertex * v : e->vertices())
      {
        std::shared_ptr<FENode> n = _vertexToNode[v];
        std::shared_ptr<Index> n_dof_index = n->dofIndex();
        for(unsigned int i : *node_dof_index)
        {
          // std::cout << "e->tags(0) " << e->tags(0)<< std::endl;;
          // std::cout << " i dof " << i <<std ::endl;
          // std::cout << " n_dof_index(i) " << (*n_dof_index)[i] << std ::endl;
          // std::cout << " x y z" <<  v->x() << v-> y() << v->z() << std::endl;
          if(find(bdIndex->begin(), bdIndex->end(), (*n_dof_index)[i]) == bdIndex->end())  // check if the dof is already existing
          {
            //std::cout << " n_dof_index(i) added" << (*n_dof_index)[i] << std ::endl;
            bdIndex->push_back((*n_dof_index)[i]);
          }
        }
      }
    }
  }
  int  bdIndex_added_size = bdIndex->size() - bdIndex_old_size;

  int  bdPrescribedVelocity_old_size = bdPrescribedVelocity->size();
  bdPrescribedVelocity->resize(bdPrescribedVelocity_old_size +  bdIndex_added_size);

  for(unsigned int k =0; k <bdIndex_added_size; k++)
  {
    bdPrescribedVelocity->setValue(k+bdPrescribedVelocity_old_size, 0.0);
  }

  DEBUG_END("siconos::mechanics::fem::native::FiniteElementModel::applyDirichletBoundaryConditions(int physical_entity_tag, std::shared_ptr<IndexInt node_dof_index, std::shared_ptr<BoundaryCondition >> _boundaryCondition)\n");
}

void siconos::mechanics::fem::native::FiniteElementModel::applyNodalForces(
  int physical_entity_tag, std::shared_ptr<SiconosVector> nodal_forces,
  std::shared_ptr<SiconosVector> forces)
{
  DEBUG_BEGIN("siconos::mechanics::fem::native::FiniteElementModel::applyNodalForces(int physical_entity_tag, std::shared_ptr<SiconosVector> nodal_forces, std::shared_ptr<SiconosVector> forces)\n");

  std::shared_ptr<IndexInt> f_index  = std::make_shared<IndexInt>(0);
  for(MElement * e : _mesh->elements())
  {
    if(e->tags(0) == physical_entity_tag)
    {
      for(MVertex * v : e->vertices())
      {
        std::shared_ptr<FENode> n = _vertexToNode[v];

        if(find(f_index->begin(), f_index->end(), n->num()) == f_index->end())  // check if the node is already existing
        {
          std::cout <<"Apply nodal force on node number " <<  n->num()  << " vertex number " << v->num() << std::endl;
          f_index->push_back(n->num());
          std::shared_ptr<Index> n_dof_index = n->dofIndex();
          for(unsigned int i =0; i < nodal_forces->size(); i++)
          {
            forces->setValue((*n_dof_index)[i], (*nodal_forces)(i));
          }
        }
      }
    }
  }
  DEBUG_END("siconos::mechanics::fem::native::FiniteElementModel::applyNodalForces(int physical_entity_tag, std::shared_ptr<SiconosVector> nodal_forces, std::shared_ptr<SiconosVector> forces)\n");
}


void siconos::mechanics::fem::native::FiniteElementModel::display(bool brief) const
{
  std::cout << "===== FiniteElementModel display ===== " <<std::endl;
  std::cout << "- numberOfNodes : " << _nodes.size() << std::endl;
  std::cout << "- numberOfElements : " << _elements.size() << std::endl;
  for(std::shared_ptr<FElement> f : _elements)
  {
    f->display();
  }
}
