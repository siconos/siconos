/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include "NewtonEulerFrom1DLocalFrameR.hpp"
#include "OSNSMatrixProjectOnConstraints.hpp"
using namespace std;
using namespace RELATION;
//#define MLCPPROJ_DEBUG


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


  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();
  //  cout<<"isLinear: "<<isLinear<<" hasTopologyChanged: "<<hasTopologyChanged<<"hasBeUpdated: "<<_hasBeUpdated<<endl;

  if (!_hasBeUpdated || !isLinear)
  {
    if (!_hasBeUpdated)
    {
      //      printf("MLCPProjectOnConstraints::updateUnitaryBlocks must be updated.\n");
      _n = 0;
      _m = 0;
      _curBlock = 0;
    }
    UnitaryRelationsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {
      SP::UnitaryRelation UR = indexSet->bundle(*vi);
      SP::NewtonEulerR  nR = boost::static_pointer_cast<NewtonEulerR>(UR->interaction()->relation());
      unsigned int nslawSize = nR->yProj()->size();

      if (! indexSet->properties(*vi).blockProj)
      {
        indexSet->properties(*vi).blockProj.reset(new SimpleMatrix(nslawSize, nslawSize));
      }

      computeDiagonalUnitaryBlock(*vi);
    }


    UnitaryRelationsGraph::EIterator ei, eiend;
    for (boost::tie(ei, eiend) = indexSet->edges();
         ei != eiend; ++ei)
    {
      SP::UnitaryRelation UR1 = indexSet->bundle(indexSet->source(*ei));
      SP::NewtonEulerR  nR1 = boost::static_pointer_cast<NewtonEulerR>(UR1->interaction()->relation());
      unsigned int nslawSize1 = nR1->yProj()->size();

      SP::UnitaryRelation UR2 = indexSet->bundle(indexSet->target(*ei));
      SP::NewtonEulerR  nR2 = boost::static_pointer_cast<NewtonEulerR>(UR2->interaction()->relation());
      unsigned int nslawSize2 = nR2->yProj()->size();

      // Memory allocation if needed
      // unsigned int nslawSize1 = UR1->getNonSmoothLawSizeProjectOnConstraints();
      // unsigned int nslawSize2 = UR2->getNonSmoothLawSizeProjectOnConstraints();
      unsigned int isrc = indexSet->index(indexSet->source(*ei));
      unsigned int itar = indexSet->index(indexSet->target(*ei));

      if (itar > isrc) // upper block
      {
        if (! indexSet->properties(*ei).upper_blockProj)
        {
          indexSet->properties(*ei).upper_blockProj.reset(new SimpleMatrix(nslawSize1, nslawSize2));
        }
      }
      else  // lower block
      {
        if (! indexSet->properties(*ei).lower_blockProj)
        {
          indexSet->properties(*ei).lower_blockProj.reset(new SimpleMatrix(nslawSize1, nslawSize2));
        }
      }


      computeUnitaryBlock(*ei);


      // allocation for transposed block
      // should be avoided

      if (itar > isrc) // upper block has been computed
      {
        // if (!indexSet->properties(*ei).lower_blockProj)
        //   {
        //     indexSet->properties(*ei).lower_blockProj.
        //  reset(new SimpleMatrix(indexSet->properties(*ei).upper_blockProj->size(1),
        //             indexSet->properties(*ei).upper_blockProj->size(0)));
        //   }
        indexSet->properties(*ei).lower_blockProj.reset(new SimpleMatrix(*(indexSet->properties(*ei).upper_blockProj)));
        indexSet->properties(*ei).lower_blockProj->trans();
        //          indexSet->properties(*ei).lower_blockProj->trans(*indexSet->properties(*ei).upper_blockProj);
      }
      else
      {
        assert(itar < isrc);    // lower block has been computed
        // if (!indexSet->properties(*ei).upper_blockProj)
        //   {
        //     indexSet->properties(*ei).upper_blockProj.
        //  reset(new SimpleMatrix(indexSet->properties(*ei).lower_blockProj->size(1),
        //             indexSet->properties(*ei).lower_blockProj->size(0)));
        //   }
        indexSet->properties(*ei).upper_blockProj.
        reset(new SimpleMatrix(*(indexSet->properties(*ei).lower_blockProj)));
        indexSet->properties(*ei).upper_blockProj->trans();
      }
      //      printf("MLCPP upper: %i %i\n",indexSet->properties(*ei).upper_blockProj->size(0),indexSet->properties(*ei).upper_blockProj->size(1));
      //      printf("MLCPP lower: %i %i\n",indexSet->properties(*ei).lower_blockProj->size(0),indexSet->properties(*ei).lower_blockProj->size(1));

    }

  }
}

void MLCPProjectOnConstraints::computeDiagonalUnitaryBlock(const UnitaryRelationsGraph::VDescriptor& vd)
{
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());

  SP::DynamicalSystem DS1 = indexSet->properties(vd).source;
  SP::DynamicalSystem DS2 = indexSet->properties(vd).target;
  SP::UnitaryRelation UR = indexSet->bundle(vd);
  SP::NewtonEulerR  nR = boost::static_pointer_cast<NewtonEulerR>(UR->interaction()->relation());
  unsigned int nslawSize = nR->yProj()->size();


  // unsigned int nslawSize = UR->getNonSmoothLawSizeProjectOnConstraints();


  assert(indexSet->properties(vd).blockProj->size(0) == nslawSize);
  assert(indexSet->properties(vd).blockProj->size(1) == nslawSize);

  SP::SiconosMatrix currentUnitaryBlock = indexSet->properties(vd).blockProj;
  if (!_hasBeUpdated)
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

  SP::SiconosMatrix leftUnitaryBlock, rightUnitaryBlock, leftUnitaryBlock1;


  // General form of the unitaryBlock is : unitaryBlock =
  // a*extraUnitaryBlock + b * leftUnitaryBlock * centralUnitaryBlocks
  // * rightUnitaryBlock a and b are scalars, centralUnitaryBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  RELATION::TYPES relationType = UR->getRelationType();
  if (relationType != NewtonEuler)
  {
    RuntimeException::selfThrow("MLCPProjectOnConstraints::computeUnitaryBlock - relation is not from NewtonEulerR.");
  }

  currentUnitaryBlock->zero();

  cout << "compute DiagUnitary block" << endl;
  // loop over the common DS
  bool endl = false;
  for (SP::DynamicalSystem ds = DS1; !endl; ds = DS2)
  {
    assert(ds == DS1 || ds == DS2);

    endl = (ds == DS2);

    if (Type::value(*ds) != Type::NewtonEulerDS)
      RuntimeException::selfThrow("MLCPProjectOnConstraints::computeUnitaryBlock - ds is not from NewtonEulerDS.");

    SP::NewtonEulerDS neds = (boost::static_pointer_cast<NewtonEulerDS>(ds));
#ifdef MLCPPROJ_WITH_CT

    unsigned int sizeDS = neds->getDim();
    SP::SimpleMatrix T = neds->T();
    SP::SimpleMatrix workT(new SimpleMatrix(*T));
    workT->trans();
    SP::SimpleMatrix workT2(new SimpleMatrix(6, 6));
    prod(*workT, *T, *workT2, true);
    leftUnitaryBlock.reset(new SimpleMatrix(nslawSize, sizeDS));
    UR->getLeftUnitaryBlockForDS(ds, leftUnitaryBlock);
    SP::SiconosMatrix work(new SimpleMatrix(*leftUnitaryBlock));
    cout << "LinearOSNS : leftUBlock\n";
    work->display();
    work->trans();
    cout << "LinearOSNS::computeUnitaryBlock workT2" << endl;
    workT2->display();
    workT2->PLUForwardBackwardInPlace(*work);
    prod(*leftUnitaryBlock, *work, *currentUnitaryBlock, false);
#else
    if (0) //(boost::static_pointer_cast<NewtonEulerR> UR->interaction()->relation())->_isConstact){
    {
      unsigned int sizeDS = neds->getDim();
      SP::SimpleMatrix T = neds->T();
      SP::SimpleMatrix workT(new SimpleMatrix(*T));
      workT->trans();
      SP::SimpleMatrix workT2(new SimpleMatrix(6, 6));
      prod(*workT, *T, *workT2, true);
      leftUnitaryBlock1.reset(new SimpleMatrix(nslawSize, sizeDS));
      UR->getLeftUnitaryBlockForDS(ds, leftUnitaryBlock);
      leftUnitaryBlock.reset(new SimpleMatrix(1, sizeDS));
      for (unsigned int ii = 0; ii < sizeDS; ii++)
        leftUnitaryBlock->setValue(1, ii, leftUnitaryBlock1->getValue(1, ii));

      SP::SiconosMatrix work(new SimpleMatrix(*leftUnitaryBlock));
      //cout<<"LinearOSNS : leftUBlock\n";
      //work->display();
      work->trans();
      //cout<<"LinearOSNS::computeUnitaryBlock workT2"<<endl;
      //workT2->display();
      workT2->PLUForwardBackwardInPlace(*work);
      prod(*leftUnitaryBlock, *work, *currentUnitaryBlock, false);
    }
    else
    {
      unsigned int sizeDS = (boost::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();
      leftUnitaryBlock.reset(new SimpleMatrix(nslawSize, sizeDS));
      UR->getLeftUnitaryBlockForDSProjectOnConstraints(ds, leftUnitaryBlock);
      SP::SiconosMatrix work(new SimpleMatrix(*leftUnitaryBlock));
      //cout<<"LinearOSNS nslawSize="<<nslawSize<<": leftUBlock\n";
      //work->display();
      work->trans();
      prod(*leftUnitaryBlock, *work, *currentUnitaryBlock, false);
    }

#endif

    cout << "MLCPProjectOnConstraints::computeUnitaryBlock unitaryBlock " << endl;
    currentUnitaryBlock->display();
  }
}

void MLCPProjectOnConstraints::computeUnitaryBlock(const UnitaryRelationsGraph::EDescriptor& ed)
{

  // Computes matrix _unitaryBlocks[UR1][UR2] (and allocates memory if
  // necessary) if UR1 and UR2 have commond DynamicalSystem.  How
  // _unitaryBlocks are computed depends explicitely on the type of
  // Relation of each UR.

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

  // Get dimension of the NonSmoothLaw (ie dim of the unitaryBlock)
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(_levelMin);

  SP::DynamicalSystem ds = indexSet->bundle(ed);
  SP::UnitaryRelation UR1 = indexSet->bundle(indexSet->source(ed));
  SP::UnitaryRelation UR2 = indexSet->bundle(indexSet->target(ed));

  unsigned int index1 = indexSet->index(indexSet->source(ed));
  unsigned int index2 = indexSet->index(indexSet->target(ed));

  SP::NewtonEulerR  nR1 = boost::static_pointer_cast<NewtonEulerR>(UR1->interaction()->relation());
  unsigned int nslawSize1 = nR1->yProj()->size();

  SP::NewtonEulerR  nR2 = boost::static_pointer_cast<NewtonEulerR>(UR2->interaction()->relation());
  unsigned int nslawSize2 = nR2->yProj()->size();


  // unsigned int nslawSize1 = UR1->getNonSmoothLawSizeProjectOnConstraints();
  // unsigned int nslawSize2 = UR2->getNonSmoothLawSizeProjectOnConstraints();

  /*
    DynamicalSystemsSet commonDS;
    intersection(*UR1->dynamicalSystems(),*UR2->dynamicalSystems(), commonDS);
    assert (!commonDS.isEmpty()) ;
    for (DSIterator itDS = commonDS.begin(); itDS!=commonDS.end(); itDS++)
    {
    assert (*itDS == ds);
    }
  */

  SP::SiconosMatrix currentUnitaryBlock;

  assert(index1 != index2);

  if (index2 > index1) // upper block
  {
    //     if (! indexSet->properties(ed).upper_block)
    //     {
    //       indexSet->properties(ed).upper_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
    //     }

    assert(indexSet->properties(ed).upper_blockProj->size(0) == nslawSize1);
    assert(indexSet->properties(ed).upper_blockProj->size(1) == nslawSize2);

    currentUnitaryBlock = indexSet->properties(ed).upper_blockProj;
  }
  else  // lower block
  {
    //     if (! indexSet->properties(ed).lower_block)
    //     {
    //       indexSet->properties(ed).lower_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
    //     }

    assert(indexSet->properties(ed).lower_blockProj->size(0) == nslawSize1);
    assert(indexSet->properties(ed).lower_blockProj->size(1) == nslawSize2);

    currentUnitaryBlock = indexSet->properties(ed).lower_blockProj;
  }


  SP::SiconosMatrix leftUnitaryBlock, rightUnitaryBlock;

  RELATION::TYPES relationType1, relationType2;

  // General form of the unitaryBlock is : unitaryBlock =
  // a*extraUnitaryBlock + b * leftUnitaryBlock * centralUnitaryBlocks
  // * rightUnitaryBlock a and b are scalars, centralUnitaryBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  relationType1 = UR1->getRelationType();
  relationType2 = UR2->getRelationType();
  if (relationType1 != NewtonEuler ||
      relationType2 != NewtonEuler)
    RuntimeException::selfThrow("LinearOSNS::computeUnitaryBlock not yet implemented for relation of type " + relationType1);

  // ==== First Order Relations - Specific treatment for diagonal
  // _unitaryBlocks ===
  assert(UR1 != UR2);
  currentUnitaryBlock->zero();
#ifdef MLCPPROJ_WITH_CT
  unsigned int sizeDS = (boost::static_pointer_cast<NewtonEulerDS>(ds))->getDim();
  leftUnitaryBlock.reset(new SimpleMatrix(nslawSize1, sizeDS));
  UR1->getLeftUnitaryBlockForDS(ds, leftUnitaryBlock);
  SP::NewtonEulerDS neds = (boost::static_pointer_cast<NewtonEulerDS>(ds));
  SP::SimpleMatrix T = neds->T();
  SP::SimpleMatrix workT(new SimpleMatrix(*T));
  workT->trans();
  SP::SimpleMatrix workT2(new SimpleMatrix(6, 6));
  prod(*workT, *T, *workT2, true);
  rightUnitaryBlock.reset(new SimpleMatrix(nslawSize2, sizeDS));
  UR2->getLeftUnitaryBlockForDS(ds, rightUnitaryBlock);
  rightUnitaryBlock->trans();
  workT2->PLUForwardBackwardInPlace(*rightUnitaryBlock);
  prod(*leftUnitaryBlock, *rightUnitaryBlock, *currentUnitaryBlock, false);

#else

  unsigned int sizeDS = (boost::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();
  leftUnitaryBlock.reset(new SimpleMatrix(nslawSize1, sizeDS));
  UR1->getLeftUnitaryBlockForDSProjectOnConstraints(ds, leftUnitaryBlock);
  SP::NewtonEulerDS neds = (boost::static_pointer_cast<NewtonEulerDS>(ds));
  rightUnitaryBlock.reset(new SimpleMatrix(nslawSize2, sizeDS));
  UR2->getLeftUnitaryBlockForDSProjectOnConstraints(ds, rightUnitaryBlock);
  rightUnitaryBlock->trans();
  prod(*leftUnitaryBlock, *rightUnitaryBlock, *currentUnitaryBlock, false);

#endif

}

void MLCPProjectOnConstraints::computeqBlock(SP::UnitaryRelation UR, unsigned int pos)
{

  SP::NewtonEulerR  nR = boost::static_pointer_cast<NewtonEulerR>(UR->interaction()->relation());
  unsigned int sizeY = nR->yProj()->size();

  // unsigned int sizeY = UR->getNonSmoothLawSizeProjectOnConstraints();
  SP::Relation R = UR->interaction()->relation();
  SP::NewtonEulerR ner = (boost::static_pointer_cast<NewtonEulerR>(R));
  for (unsigned int i = 0; i < sizeY; i++)
    _q->setValue(pos + i, ner->yProj()->getValue(UR->getRelativePosition() + i));
#ifdef MLCPPROJ_DEBUG
  printf("MLCPProjectOnConstraints::computeqBlock, _q from yProj\n");
  _q->display();
#endif

}
void MLCPProjectOnConstraints::postCompute()
{
  _hasBeUpdated = true;
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
  /** We chose to do a small step in view of stabilized the algorithm.*/
  (*_z) *= 0.2;
  unsigned int pos = 0;
#ifdef MLCPPROJ_DEBUG
  printf("MLCPProjectOnConstraints::postCompute _z\n");
  _z->display();
  display();
#endif



  UnitaryRelationsGraph::VIterator ui, uiend;
#ifdef MLCPPROJ_DEBUG
  printf("MLCPProjectOnConstraints::postCompute BEFORE UPDATE:\n");
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet->bundle(*ui);
    SP::Relation R = ur->interaction()->relation();
    SP::NewtonEulerR ner = (boost::static_pointer_cast<NewtonEulerR>(R));
    ner->computeh(0);
    ner->getq()->display();
  }
#endif


  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {

    SP::UnitaryRelation ur = indexSet->bundle(*ui);
    // Get the relative position of UR-unitaryBlock in the vector w
    // or z
    pos = _M->getPositionOfUnitaryBlock(ur);
    SP::NewtonEulerR  nR = boost::static_pointer_cast<NewtonEulerR>(ur->interaction()->relation());
    unsigned int nslawSize = nR->yProj()->size();

    // unsigned int nslawSize = ur->getNonSmoothLawSizeProjectOnConstraints();


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
    printf("MLCPProjectOnConstraints::postCompute ner->jachq \n");
    ner->jachq()->display();
    ner->jachqT()->display();
#endif

#ifdef MLCPPROJ_DEBUG
    printf("MLCPProjectOnConstraints::postCompute q before update\n");
    (ner->getq())->display();
#endif
    SP::SimpleVector yProj = ner->yProj();
    SP::SimpleVector aBuff(new SimpleVector(yProj->size()));
    setBlock(*_z, aBuff, yProj->size(), pos, 0);
    /*Now, update the ds's dof throw the relation*/
    //proj_with_q SP::SiconosMatrix J=ner->jachqProj();
    printf("MLCPP lambda of ur is pos =%i :", pos);
    aBuff->display();

    unsigned int k = 0;
    unsigned int NumDS = 0;
    DSIterator itDS;
    Index coord(8);
#ifdef MLCPPROJ_WITH_CT
    SP::SiconosMatrix J = ner->jachqT();
    //printf("J\n");
    //J->display();

    SP::SimpleMatrix aux(new SimpleMatrix(*J));
    aux->trans();

    itDS = ur->interaction()->dynamicalSystemsBegin();
    while (itDS != ur->interaction()->dynamicalSystemsEnd())
    {
      SP::NewtonEulerDS neds = (boost::static_pointer_cast<NewtonEulerDS>(*itDS));
      k += neds->getqDim();
      itDS++;
    }
    SP::SimpleVector deltaqG(new SimpleVector(k));
    // prod(*aux,*aBuff,*vel,true);
    k = 0;
    deltaqG->zero();

    itDS = ur->interaction()->dynamicalSystemsBegin();
    while (itDS != ur->interaction()->dynamicalSystemsEnd())
    {
      Type::Siconos dsType = Type::value(**itDS);
      if (dsType != Type::NewtonEulerDS)
        RuntimeException::selfThrow("MLCPProjectOnConstaraint::postCompute- ds is not from NewtonEulerDS.");
      SP::NewtonEulerDS neds = (boost::static_pointer_cast<NewtonEulerDS>(*itDS));

      SP::SimpleVector vel(new SimpleVector(neds->getDim()));
      SP::SimpleVector deltaQ(new SimpleVector(neds->getqDim()));
      coord[0] = k;
      coord[1] = k + neds->getDim();
      coord[2] = 0;
      coord[3] = aux->size(1);
      coord[4] = 0;
      coord[5] = aux->size(1);
      coord[6] = 0;
      coord[7] = neds->getDim();
      subprod(*aux   , *aBuff, *vel, coord, true);
      SP::SimpleMatrix T = neds->T();
      SP::SimpleMatrix workT(new SimpleMatrix(*T));
      workT->trans();
      SP::SimpleMatrix workT2(new SimpleMatrix(6, 6));
      prod(*workT, *T, *workT2, true);

      workT2->PLUForwardBackwardInPlace(*vel);
      prod(*T, *vel, *deltaQ);

      for (unsigned int ii = 0; ii < neds->getqDim(); ii++)
      {
        ner->getq()->setValue(k + ii,
                              deltaQ->getValue(ii) +
                              ner->getq()->getValue(k + ii));
        deltaqG->setValue(k + ii, deltaQ->getValue(ii) + deltaqG->getValue(k + ii));
      }

      k += neds->getDim();
      itDS++;

      printf("MLCPProjectOnConstraints::postCompute ds %d, q updated\n", NumDS);
      (ner->getq())->display();

      NumDS++;

    }
    printf("MLCPProjectOnConstraints deltaqG:\n");
    deltaqG->display();
#else

    SP::SiconosMatrix J = ner->jachq();
    SP::SimpleMatrix aux(new SimpleMatrix(*J));
    aux->trans();
    prod(*aux, *aBuff, *(ner->getq()), false);
#endif


#ifdef MLCPPROJ_DEBUG
    printf("MLCPProjectOnConstraints::postCompute _z\n");
    _z->display();
    printf("MLCPProjectOnConstraints::postCompute q updated\n");
    (ner->getq())->display();
#endif
  }
#ifdef MLCPPROJ_DEBUG
  printf("MLCPProjectOnConstraints::postCompute AFTER UPDATE:\n");
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet->bundle(*ui);
    SP::Relation R = ur->interaction()->relation();
    SP::NewtonEulerR ner = (boost::static_pointer_cast<NewtonEulerR>(R));
    ner->computeh(0);
    ner->getq()->display();
  }
#endif
}
void MLCPProjectOnConstraints::computeOptions(SP::UnitaryRelation UR1, SP::UnitaryRelation UR2)
{
  //  printf("MLCPProjectOnConstraints::computeOptions\n");
  // Get dimension of the NonSmoothLaw (ie dim of the unitaryBlock)


  SP::NewtonEulerR  nR1 = boost::static_pointer_cast<NewtonEulerR>(UR1->interaction()->relation());
  unsigned int nslawSize1 = nR1->yProj()->size();

  SP::NewtonEulerR  nR2 = boost::static_pointer_cast<NewtonEulerR>(UR2->interaction()->relation());
  unsigned int nslawSize2 = nR2->yProj()->size();


  // unsigned int nslawSize1 = UR1->getNonSmoothLawSizeProjectOnConstraints();
  // unsigned int nslawSize2 = UR2->getNonSmoothLawSizeProjectOnConstraints();



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
    SP::NewtonEulerFrom1DLocalFrameR ri = boost::static_pointer_cast<NewtonEulerFrom1DLocalFrameR> (UR1->interaction()->relation());
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
