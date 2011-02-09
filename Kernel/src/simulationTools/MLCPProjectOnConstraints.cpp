/* Siconos-Kernel, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "MLCPProjectOnConstraints.hpp"
#include "MixedComplementarityConditionNSL.hpp"
#include "EqualityConditionNSL.hpp"
#include "Simulation.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonEulerRImpact.hpp"
#include "OSNSMatrixProjectOnConstraints.hpp"
using namespace std;
using namespace RELATION;
#define MLCPPROJ_DEBUG
void MLCPProjectOnConstraints::updateM()
{
  assert(0);
  // Get index set from Simulation
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());


  if (!_M)
  {
    // Creates and fills M using UR of indexSet
    _M.reset(new OSNSMatrixProjectOnConstraints(_n + _m, _n + _m, _MStorageType));
    _numerics_problem.M = &*_M->getNumericsMatrix();
    _numerics_problem.A = 0;
    _numerics_problem.B = 0;
    _numerics_problem.C = 0;
    _numerics_problem.D = 0;
    _numerics_problem.a = 0;
    _numerics_problem.b = 0;
    _numerics_problem.problemType = 0;
    _numerics_problem.n = _n;
    _numerics_problem.m = _m;
  }
  else
  {
    _M->setStorageType(_MStorageType);
    //_M->fill(indexSet);
  }
  _M->fill(indexSet);
  _sizeOutput = _M->size();
}

void MLCPProjectOnConstraints::initOSNSMatrix()
{
  _M.reset(new OSNSMatrixProjectOnConstraints(0, 0, _MStorageType));
  _n = 0;
  _m = 0;
  _curBlock = 0;
}
// Constructor from a set of data
MLCPProjectOnConstraints::MLCPProjectOnConstraints(const int newNumericsSolverId):
  MLCP(newNumericsSolverId)
{
  _levelMin = 1;
  _levelMax = 1;
}


void MLCPProjectOnConstraints::display() const
{
  cout << "======= MLCPProjectOnConstraints of size " << _sizeOutput << " with: " << endl;
  cout << "======= m " << _m << " _n " << _n << endl;
  LinearOSNS::display();
}

void MLCPProjectOnConstraints::updateUnitaryBlocks()
{
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(_levelMin);


  bool hasTopologyChanged = simulation()->model()->
                            nonSmoothDynamicalSystem()->topology()->hasChanged();
  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();

  if (hasTopologyChanged || !isLinear)
  {
    _n = 0;
    _m = 0;
    _curBlock = 0;
    UnitaryRelationsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {
      SP::UnitaryRelation UR = indexSet->bundle(*vi);
      unsigned int nslawSize = UR->getNonSmoothLawSizeProjectOnConstraints();
      if (! indexSet->properties(*vi).blockProj)
      {
        indexSet->properties(*vi).blockProj.reset(new SimpleMatrix(nslawSize, nslawSize));
      }

      computeDiagonalUnitaryBlock(*vi);
    }
  }
}

void MLCPProjectOnConstraints::computeDiagonalUnitaryBlock(const UnitaryRelationsGraph::VDescriptor& vd)
{
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());

  SP::DynamicalSystem DS1 = indexSet->properties(vd).source;
  SP::DynamicalSystem DS2 = indexSet->properties(vd).target;
  SP::UnitaryRelation UR = indexSet->bundle(vd);


  unsigned int nslawSize = UR->getNonSmoothLawSizeProjectOnConstraints();


  assert(indexSet->properties(vd).blockProj->size(0) == nslawSize);
  assert(indexSet->properties(vd).blockProj->size(1) == nslawSize);

  SP::SiconosMatrix currentUnitaryBlock = indexSet->properties(vd).blockProj;

  computeOptions(UR, UR);
  // Computes matrix _unitaryBlocks[UR1][UR2] (and allocates memory if
  // necessary) if UR1 and UR2 have commond DynamicalSystem.  How
  // _unitaryBlocks are computed depends explicitely on the type of
  // Relation of each UR.

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

  // Get the W and Theta maps of one of the Unitary Relation -
  // Warning: in the current version, if OSI!=Moreau, this fails.
  // If OSI = MOREAU, centralUnitaryBlocks = W if OSI = LSODAR,
  // centralUnitaryBlocks = M (mass matrices)
  MapOfDSMatrices centralUnitaryBlocks;
  getOSIMaps(UR, centralUnitaryBlocks);

  SP::SiconosMatrix leftUnitaryBlock, rightUnitaryBlock;

  double h = simulation()->timeDiscretisation()->currentTimeStep();

  // General form of the unitaryBlock is : unitaryBlock =
  // a*extraUnitaryBlock + b * leftUnitaryBlock * centralUnitaryBlocks
  // * rightUnitaryBlock a and b are scalars, centralUnitaryBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  RELATION::TYPES relationType = UR->getRelationType();

  currentUnitaryBlock->zero();


  // loop over the common DS
  bool endl = false;
  for (SP::DynamicalSystem ds = DS1; !endl; ds = DS2)
  {
    assert(ds == DS1 || ds == DS2);

    endl = (ds == DS2);

    if (Type::value(*ds) != Type::NewtonEulerDS)
      RuntimeException::selfThrow("MLCPProjectOnConstraints::computeUnitaryBlock - ds is not from NewtonEulerDS.");

    unsigned int sizeDS = (boost::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();

    // get _unitaryBlocks corresponding to the current DS
    // These _unitaryBlocks depends on the relation type.
    leftUnitaryBlock.reset(new SimpleMatrix(nslawSize, sizeDS));
    UR->getLeftUnitaryBlockForDSProjectOnConstraints(ds, leftUnitaryBlock);

    if (relationType == Lagrangian ||
        relationType == NewtonEuler)
    {
      // UR1 == UR2
      SP::SiconosMatrix work(new SimpleMatrix(*leftUnitaryBlock));
      //
      work->trans();
      //        cout<<"LinearOSNS : leftUBlock'\n";
      //        work->display();
      //        cout<<"LinearOSNS::computeUnitaryBlock leftUnitaryBlock"<<endl;
      //        leftUnitaryBlock->display();

      //*currentUnitaryBlock +=  *leftUnitaryBlock ** work;
      prod(*leftUnitaryBlock, *work, *currentUnitaryBlock, false);
      //      gemm(CblasNoTrans,CblasNoTrans,1.0,*leftUnitaryBlock,*work,1.0,*currentUnitaryBlock);
      //*currentUnitaryBlock *=h;
      //      cout<<"MLCPProjectOnConstraints::computeUnitaryBlock unitaryBlock "<<endl;
      //      currentUnitaryBlock->display();
    }
    else RuntimeException::selfThrow("LinearOSNS::computeUnitaryBlock not yet implemented for relation of type " + relationType);
  }
}

void MLCPProjectOnConstraints::computeUnitaryBlock(const UnitaryRelationsGraph::EDescriptor& ed)
{
  assert(0);

}

void MLCPProjectOnConstraints::computeqBlock(SP::UnitaryRelation UR, unsigned int pos)
{
  unsigned int sizeY = UR->getNonSmoothLawSizeProjectOnConstraints();
  SP::Relation R = UR->interaction()->relation();
  SP::NewtonEulerR ner = (boost::static_pointer_cast<NewtonEulerR>(R));
  for (int i = 0; i < sizeY; i++)
    _q->setValue(pos + i, ner->yProj()->getValue(UR->getRelativePosition() + i));
#ifdef MLCPPROJ_DEBUG
  printf("MLCPProjectOnConstraints::computeqBlock, _q from yProj\n");
  _q->display();
#endif

}
void MLCPProjectOnConstraints::postCompute()
{

  // This function is used to set y/lambda values using output from
  // lcp_driver (w,z).  Only UnitaryRelations (ie Interactions) of
  // indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());

  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;

  // === Loop through "active" Unitary Relations (ie present in
  // indexSets[1]) ===
  (*_z) *= 0.1;
  unsigned int pos = 0;
#ifdef MLCPPROJ_DEBUG
  printf("MLCPProjectOnConstraints::postCompute\n");
#endif
  UnitaryRelationsGraph::VIterator ui, uiend;

  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet->bundle(*ui);
    // Get the relative position of UR-unitaryBlock in the vector w
    // or z
    pos = _M->getPositionOfUnitaryBlock(ur);

    // Get Y and Lambda for the current Unitary Relation
    //y = ur->y(levelMin());
    //lambda = ur->lambda(levelMin());
    // Copy _w/_z values, starting from index pos into y/lambda.

    //      setBlock(*_w, y, y->size(), pos, 0);// Warning: yEquivalent is
    // saved in y !!
    //setBlock(*_z, lambda, lambda->size(), pos, 0);

    SP::Relation R = ur->interaction()->relation();
    Type::Siconos RType = Type::value(*R);
    if (RType != Type::NewtonEulerR)
      RuntimeException::selfThrow("MLCPProjectOnConstraints::postCompute - R is not from NewtonEulerR.");
    SP::NewtonEulerR ner = (boost::static_pointer_cast<NewtonEulerR>(R));
#ifdef MLCPPROJ_DEBUG
    printf("MLCPProjectOnConstraints::postCompute q before update\n");
    (ner->getq())->display();
#endif
    SP::SimpleVector yProj = ner->yProj();
    SP::SimpleVector aBuff(new SimpleVector(yProj->size()));
    setBlock(*_z, aBuff, yProj->size(), pos, 0);
    /*Now, update the ds's dof throw the relation*/
    SP::SiconosMatrix J = ner->jachqProj();
    SP::SimpleMatrix aux(new SimpleMatrix(*J));
    aux->trans();

    Index coord(8);
    coord[0] = 0; /*first line of aux*/
    coord[1] = aux->size(0); /*last line of aux*/
    coord[2] = 0; /*first col of aux*/
    coord[3] = yProj->size();
    coord[4] = 0;
    coord[5] = yProj->size();
    coord[6] = 0;
    coord[7] = (ner->getq())->size();
    subprod(*aux, *aBuff, *(ner->getq()), coord, false);
#ifdef MLCPPROJ_DEBUG
    printf("MLCPProjectOnConstraints::postCompute _z\n");
    _z->display();
    printf("MLCPProjectOnConstraints::postCompute q updated\n");
    (ner->getq())->display();
#endif
  }

}
void MLCPProjectOnConstraints::computeOptions(SP::UnitaryRelation UR1, SP::UnitaryRelation UR2)
{
  // Get dimension of the NonSmoothLaw (ie dim of the unitaryBlock)
  unsigned int nslawSize1 = UR1->getNonSmoothLawSizeProjectOnConstraints();
  unsigned int nslawSize2 = UR2->getNonSmoothLawSizeProjectOnConstraints();

  unsigned int equalitySize1 =  0;
  unsigned int equalitySize2 =  0;
  if (Type::value(*(UR1->interaction()->nonSmoothLaw()))
      == Type::MixedComplementarityConditionNSL)
    equalitySize1 =  MixedComplementarityConditionNSL::convert(UR1->interaction()->nonSmoothLaw())->getEqualitySize();
  else if (Type::value(*(UR1->interaction()->nonSmoothLaw()))
           == Type::EqualityConditionNSL)
    equalitySize1 = nslawSize1;

  if (Type::value(*(UR2->interaction()->nonSmoothLaw()))
      == Type::MixedComplementarityConditionNSL)
    equalitySize2 = MixedComplementarityConditionNSL::
                    convert(UR2->interaction()->nonSmoothLaw())->getEqualitySize();
  else if (Type::value(*(UR2->interaction()->nonSmoothLaw()))
           == Type::EqualityConditionNSL)
    equalitySize2 = nslawSize2;

  if (Type::value(*(UR1->interaction()->nonSmoothLaw())) == Type::NewtonImpactFrictionNSL ||
      Type::value(*(UR1->interaction()->nonSmoothLaw())) == Type::NewtonImpactNSL)
  {
    SP::NewtonEulerRImpact ri = boost::static_pointer_cast<NewtonEulerRImpact> (UR1->interaction()->relation());
    if (ri->_isOnContact)
      equalitySize1 = nslawSize1;
  }


  if (UR1 == UR2)
  {
    //UR1->getExtraUnitaryBlock(currentUnitaryBlock);
    _m += nslawSize1 - equalitySize1;
    _n += equalitySize1;
    //    _m=0;
    //_n=6;
    if (_curBlock > MLCP_NB_BLOCKS - 2)
      printf("MLCP.cpp : number of block to small, memory crach below!!!\n");
    /*add an equality block.*/
    if (equalitySize1 > 0)
    {
      _numerics_problem.blocksLine[_curBlock + 1] = _numerics_problem.blocksLine[_curBlock] + equalitySize1;
      _numerics_problem.blocksIsComp[_curBlock] = 0;
      _curBlock++;
    }
    /*add a complementarity block.*/
    if (nslawSize1 - equalitySize1 > 0)
    {
      _numerics_problem.blocksLine[_curBlock + 1] = _numerics_problem.blocksLine[_curBlock] + nslawSize1 - equalitySize1;
      _numerics_problem.blocksIsComp[_curBlock] = 1;
      _curBlock++;
    }
  }
}
