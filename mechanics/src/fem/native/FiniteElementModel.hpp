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

/*! \file FiniteElementLinearTIDS.hpp

 */
#ifndef FINITEELEMENTMODEL_H
#define FINITEELEMENTMODEL_H
#include <vector>
#include <iostream>
#include <map>


#include "SiconosAlgebraTypeDef.hpp"          //for Index
#include "SimulationTypeDef.hpp"              //for IndexInt
#include "BoundaryCondition.hpp"

#include "Mesh.hpp"
#include "Material.hpp"
namespace siconos::mechanics::fem::native
{
// a Finite Element node
struct FENode
{
  /* node number */
  size_t _num;

  /* associated Mvertex */
  MVertex * _mVertex;

  /* associated dof number  in the global dof vector*/
  std::shared_ptr<Index> _dofIndex;

  /* */

  FENode(size_t num, MVertex *v, std::shared_ptr<Index> dofIndex): _num(num), _mVertex(v), _dofIndex(dofIndex) {};

  size_t num()
  {
    return _num;
  }

  std::shared_ptr<Index> dofIndex()
  {
    return _dofIndex;
  };

  double x()
  {
    return _mVertex->x();
  }

  double y()
  {
    return _mVertex->y();
  }

  double z()
  {
    return _mVertex->z();
  }

  void display()
  {
    std::cout << "     - Fe Node - number: " << _num
              << "               - ndof:" << _dofIndex->size()
              << "               - dofIndex: "<< _dofIndex->front() << ":" << _dofIndex->back()  ;
    std::cout << std::endl;
  };
};
DEFINE_SPTR(FENode)



enum FINITE_ELEMENT_TYPE // we follow the gmsh numbering convention.
{
  L2=1,  // 2-node line.
  T3=2,  // 3-node triangle.
  Q4=3,  // 4-node quadrangle.
  TH4=4, // 4-node tetrahedron.
  H8=5,  // 8-node hexahedron.
  P6=6,  // 6-node prism.
  PY5=7, // 5-node pyramid.
  L3=8,  // 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
  T6=9,  // 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
};
enum FINITE_ELEMENT_FAMILY
{
  ISOPARAMETRIC
};


static const double I_TET2_X [] = {0.13819660112501052, 0.13819660112501052, 0.13819660112501052, 0.58541019662496840};
static const double I_TET2_Y [] = {0.13819660112501052, 0.13819660112501052, 0.58541019662496840, 0.13819660112501052};
static const double I_TET2_Z [] = {0.13819660112501052, 0.58541019662496840, 0.13819660112501052, 0.13819660112501052};
static const double I_TET2_W [] = {0.04166666666666666, 0.04166666666666666, 0.04166666666666666, 0.04166666666666666};


static const double GP_T3_1_p1[]= {0.333333333333333333, 0.333333333333333333, 0.5 };
static  std::vector<const double* > GaussPointsT3_1 = {GP_T3_1_p1};

static const double GP_T3_2_p1[]= {0.66666666666666667, 0.16666666666666667, 0.1666666666666666 };
static const double GP_T3_2_p2[]= {0.16666666666666667, 0.66666666666666667, 0.1666666666666666 };
static const double GP_T3_2_p3[]= {0.16666666666666667, 0.16666666666666667, 0.1666666666666666 };
static  std::vector<const double* > GaussPointsT3_2 = {GP_T3_2_p1, GP_T3_2_p2, GP_T3_2_p3};

static const double GP_TH4_1_p1[]= {0.25, 0.25, 0.25,  1.0 };
static  std::vector<const double* > GaussPointsTH4_1 = {GP_TH4_1_p1};

static const double GP_TH4_2_p1[]= {0.13819660112501052, 0.13819660112501052, 0.13819660112501052,  0.04166666666666666 };
static const double GP_TH4_2_p2[]= {0.13819660112501052, 0.13819660112501052, 0.58541019662496840,  0.04166666666666666 };
static const double GP_TH4_2_p3[]= {0.13819660112501052, 0.585410196624968405,0.13819660112501052,  0.04166666666666666 };
static const double GP_TH4_2_p4[]= {0.58541019662496840, 0.13819660112501052, 0.13819660112501052,  0.04166666666666666 };
static  std::vector<const double* > GaussPointsTH4_2 = {GP_TH4_2_p1, GP_TH4_2_p2, GP_TH4_2_p3, GP_TH4_2_p4};


static  std::vector<const double* > GaussPointsEmpty = {};


// a Finite Element
struct FElement
{
  /* element number */
  size_t _num;

  /* Element type */
  FINITE_ELEMENT_TYPE _type;

  /* Element Family */
  FINITE_ELEMENT_FAMILY _family;

  /* number of dof by Element  */
  unsigned _ndof;

  /** nodes */
  std::vector<std::shared_ptr<FENode> > _nodes;

  /* associated Mesh element */
  MElement * _mElement;

  FElement(FINITE_ELEMENT_TYPE type, unsigned int ndof, MElement *e):
    _num(e->num()), _type(type), _ndof(ndof), _mElement(e), _family(ISOPARAMETRIC) {};

  unsigned int ndof()
  {
    return _ndof;
  }
  unsigned int num()
  {
    return _num;
  }

  MElement * mElement()
  {
    return _mElement;
  }

  FINITE_ELEMENT_FAMILY family()
  {
    return _family;
  }

  int order()
  {
    switch(_type)
    {
    case T3:
    case TH4:
      return 1;

      break;
    default:
      throw("FElement::order(). element type not recognized");
    }
    return 0;
  }
  int ndofPerNode()
  {
    switch(_type)
    {
    case T3:
    case Q4:
      return 2;
    case TH4:
      return 3;
    default:
      throw("FElement::ndorPernode(). element type not recognized");
    }
    return 0;
  }

  std::vector<const double*> & GaussPoints(int order)
  {
    switch(_type)
    {
    case T3:
      if(order==1)
        return  GaussPointsT3_1;
      else if(order==2)
        return  GaussPointsT3_2;
      break;
    case TH4:
      if(order==1)
        return  GaussPointsTH4_1;
      else if(order==2)
        return  GaussPointsTH4_2;
      break;

    default:
      throw("FElement::GaussPoints(). element type not recognized");
    }
    return  GaussPointsEmpty;
  }

  std::vector<std::shared_ptr<FENode> > & nodes()
  {
    return _nodes;
  }

  void shapeFunctionIso2D(double ksi, double eta, double* N, double * Nksi, double * Neta)
  {
    switch(_type)
    {
    case T3:
    {
      N[0] = 1.0 - ksi - eta;
      N[1] = ksi;
      N[2] = eta;
      Nksi[0] = -1.0;
      Nksi[1] = 1.0;
      Nksi[2] = 0.0;
      Neta[0] = -1.0;
      Neta[1] = 0.0;
      Neta[2] = 1.0;
      break;
    }
    default:
      THROW_EXCEPTION("FElement::shapeFunctionIso2D(). element type not recognized");
    }
  }
  void shapeFunctionIso3D(double ksi, double eta, double zeta,
                          double* N, double * Nksi, double * Neta, double * Nzeta)
  {
    switch(_type)
    {
    case TH4:
    {
      N[0] = 1.0 - ksi - eta - zeta;
      N[1] = ksi;
      N[2] = eta;
      N[3] = zeta;

      Nksi[0] = -1.0;
      Nksi[1] = 1.0;
      Nksi[2] = 0.0;
      Nksi[3] = 0.0;

      Neta[0] = -1.0;
      Neta[1] = 0.0;
      Neta[2] = 1.0;
      Neta[3] = 0.0;

      Nzeta[0] = -1.0;
      Nzeta[1] = 0.0;
      Nzeta[2] = 0.0;
      Nzeta[3] = 1.0;

      break;
    }
    default:
      throw("FElement::shapeFunctionIso3D(). element type not recognized");
    }
  }

  void display()
  {
    std::cout << " - FElement - number: " << _num
              << "            - type: " << _type
              << "            - ndof: " << _ndof
              << "            - number of nodes: " << _nodes.size() ;
    std::cout << std::endl;
    for(std::shared_ptr<FENode> n : _nodes)
    {
      n->display();
    }
  };
};
DEFINE_SPTR(FElement)

// a finite element model
class FiniteElementModel
{
protected :

  /** a mesh */
  std::shared_ptr<Mesh> _mesh;

  /** nodes */
  std::vector<std::shared_ptr<FENode>> _nodes;

  /** elements */
  std::vector<std::shared_ptr<FElement>> _elements;

  /** vertex to node map **/
  std::map<MVertex *, std::shared_ptr<FENode> > _vertexToNode;

  /** MElement to FElement map **/
  std::map<MElement *, std::shared_ptr<FElement> >  _mElementTOFElement;

  /** default constructor */
  FiniteElementModel() {};

public:
  FiniteElementModel(std::shared_ptr<Mesh> mesh):
    _mesh(mesh) {};

  std::vector<std::shared_ptr<FElement> > &  elements()
  {
    return _elements;
  }

  std::vector<std::shared_ptr<FENode> > &  nodes()
  {
    return _nodes;
  }
  std::shared_ptr<FENode> vertexToNode(MVertex * v)
  {
    if(_vertexToNode.find(v) == _vertexToNode.end())
    {
      std::shared_ptr<FENode> f;
      return f;
    }
    return _vertexToNode.at(v);
  }
  /* create the FEM model from the mesh and the element type
   * \return the number of dof */
  unsigned int init();

  /* Assembly method for elemetary matrix */
  void AssembleElementaryMatrix(std::shared_ptr<SiconosMatrix> M,
                                SimpleMatrix& Me, FElement& fe);

  /** compute Mass Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeMassMatrix(std::shared_ptr<SiconosMatrix>, std::map<unsigned int, std::shared_ptr<Material> > & mat);

  /** compute elementary Mass Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeElementaryMassMatrix(SimpleMatrix& Me, FElement& fe,  double massDensity);

  /** compute Stiffness Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeStiffnessMatrix(std::shared_ptr<SiconosMatrix>,  std::map<unsigned int, std::shared_ptr<Material> > & mat);

  /** compute elementary Stiffness Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeElementaryStiffnessMatrix(SimpleMatrix& Me, FElement& fe,
                                        std::shared_ptr<SimpleMatrix> D, double thickness);

  /** compute elementary Stiffness Matrix with a direct method
   * for linear element
   **/
  void computeElementaryStiffnessMatrix_direct(SimpleMatrix& Me, FElement& fe,
      std::shared_ptr<SimpleMatrix> D, double thickness);


  void applyDirichletBoundaryConditions(int physical_entity_tag, std::shared_ptr<IndexInt> node_dof_index,
                                        std::shared_ptr<BoundaryCondition> _boundaryCondition);

  void applyNodalForces(int physical_entity_tag, std::shared_ptr<SiconosVector> nodal_forces, std::shared_ptr<SiconosVector> forces);

  void display(bool brief) const;
};
DEFINE_SPTR(FiniteElementModel)

} // namespace siconos::mechanics::fem::native


#endif // FINITEELEMENTMODEL_H
