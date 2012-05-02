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
#include "LinearOSNS.hpp"
#include "OneStepNSProblemXML.hpp"
#include "Simulation.hpp"
#include "Topology.hpp"
#include "Model.hpp"
#include "Moreau.hpp"
#include "Lsodar.hpp"
#include "NewtonEulerR.hpp"
#include "FirstOrderLinearR.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "LagrangianLinearTIR.hpp"
#include "NewtonImpactNSL.hpp"
#include "MultipleImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "LagrangianRheonomousR.hpp"
#include "LagrangianScleronomousR.hpp"
#include "LagrangianLinearTIDS.hpp"


using namespace std;
using namespace RELATION;
//#define LINEAROSNS_DEBUG
LinearOSNS::LinearOSNS(): _MStorageType(0), _keepLambdaAndYState(false)
{
}
// xml constructor
LinearOSNS::LinearOSNS(SP::OneStepNSProblemXML onestepnspbxml,
                       const string & name) :
  OneStepNSProblem(onestepnspbxml), _MStorageType(0), _keepLambdaAndYState(false)
{
  // Read storage type if given (optional , default = dense)
  if (onestepnspbxml->hasStorageType())
    _MStorageType = onestepnspbxml->getStorageType();
}

// Constructor from a set of data
LinearOSNS::LinearOSNS(const int numericsSolverId,  const string& name,
                       const string& newId):
  OneStepNSProblem(newId, numericsSolverId), _MStorageType(0), _keepLambdaAndYState(false)
{}

// Setters

void LinearOSNS::setW(const SiconosVector& newValue)
{
  assert(_sizeOutput == newValue.size() &&
         "LinearOSNS: setW, inconsistent size between given velocity size and problem size. You should set sizeOutput before");
  setObject<SimpleVector, SP::SiconosVector, SiconosVector>(_w, newValue);
}

void LinearOSNS::setz(const SiconosVector& newValue)
{
  assert(_sizeOutput == newValue.size() &&
         "LinearOSNS: setz, inconsistent size between given velocity size and problem size. You should set sizeOutput before");
  setObject<SimpleVector, SP::SiconosVector, SiconosVector>(_z, newValue);
}

void LinearOSNS::setM(const OSNSMatrix& newValue)
{
  // Useless ?
  RuntimeException::selfThrow("LinearOSNS: setM, forbidden operation. Try setMPtr().");
}

void LinearOSNS::initVectorsMemory()
{
  // Memory allocation for _w, M, z and q.
  // If one of them has already been allocated, nothing is done.
  // We suppose that user has chosen a correct size.
  if (! _w)
    _w.reset(new SimpleVector(maxSize()));
  else
  {
    if (_w->size() != maxSize())
      _w->resize(maxSize());
  }

  if (! _z)
    _z.reset(new SimpleVector(maxSize()));
  else
  {
    if (_z->size() != maxSize())
      _z->resize(maxSize());
  }

  if (! _q)
    _q.reset(new SimpleVector(maxSize()));
  else
  {
    if (_q->size() != maxSize())
      _q->resize(maxSize());
  }
}
void LinearOSNS::initOSNSMatrix()
{
  // Default size for M = maxSize()
  if (! _M)
  {
    if (_MStorageType == 0)
      _M.reset(new OSNSMatrix(maxSize(), _MStorageType));
    else // if(_MStorageType == 1) size = number of _interactionBlocks
      // = number of Interactionin the largest considered indexSet
      if (levelMin() != LEVELMAX)
      {
        _M.reset(new OSNSMatrix(simulation()->indexSet(levelMin())->size(), _MStorageType));
      }
      else
      {
        _M.reset(new OSNSMatrix(1, _MStorageType));
      }
  }
}
void LinearOSNS::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (_M,q,_w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize(sim);

  initVectorsMemory();

  // get topology
  SP::Topology topology = simulation()->model()
                          ->nonSmoothDynamicalSystem()->topology();

  // Note that _interactionBlocks is up to date since updateInteractionBlocks
  // has been called during OneStepNSProblem::initialize()

  initOSNSMatrix();

}


void LinearOSNS::computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd)
{

  // Computes matrix _interactionBlocks[inter1][inter1] (and allocates memory if
  // necessary) one or two DS are concerned by inter1 .  How
  // _interactionBlocks are computed depends explicitely on the type of
  // Relation of each Interaction.

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

  // Get dimension of the NonSmoothLaw (ie dim of the interactionBlock)
  SP::InteractionsGraph indexSet = simulation()->indexSet(_levelMin);
  SP::Interaction inter = indexSet->bundle(vd);

  // At most 2 DS are linked by an Interaction
  SP::DynamicalSystem DS1;
  SP::DynamicalSystem DS2;

  if (indexSet->properties(vd).source != indexSet->properties(vd).target)
  {
    // a two DS Interaction
    DS1 = indexSet->properties(vd).source;
    DS2 = indexSet->properties(vd).target;
  }
  else
  {
    DS1 = indexSet->properties(vd).source;
    DS2 = DS1;


    // #ifdef LINEAROSNS_DEBUG
    //   cout<<"\nLinearOSNS::computeDiagonalInteractionBlock"<<endl;
    //   std::cout << "levelMin()" << levelMin()<<std::endl;
    //   std::cout << "indexSet :"<< indexSet << std::endl;
    //   std::cout << "vd :"<< vd << std::endl;
    //   indexSet->display();
    //   std::cout << "DS1 :" << std::endl;
    //   DS1->display();
    //   std::cout << "DS2 :" << std::endl;
    //   DS2->display();
    // #endif


    InteractionsGraph::OEIterator oei, oeiend;
    for (boost::tie(oei, oeiend) = indexSet->out_edges(vd);
         oei != oeiend; ++oei)
    {
      // note : at most 4 edges
      DS2 = indexSet->bundle(*oei);
      if (DS2 != DS1) break;
    }
  }
  assert(DS1);
  assert(DS2);

  /*
    SP::DynamicalSystemsSet commonDS = inter->dynamicalSystems();
    assert (!commonDS->isEmpty()) ;
    for (DSIterator itDS = commonDS->begin(); itDS!=commonDS->end(); itDS++)
    {
    assert (*itDS == DS1 || *itDS == DS2);
    }
  */

  unsigned int nslawSize = inter->getNonSmoothLawSize();

  //   if (! indexSet->properties(vd).block)
  //   {
  //     indexSet->properties(vd).block.reset(new SimpleMatrix(nslawSize, nslawSize));
  //   }

  assert(indexSet->properties(vd).block->size(0) == inter->getNonSmoothLawSize());
  assert(indexSet->properties(vd).block->size(1) == inter->getNonSmoothLawSize());

  SP::SiconosMatrix currentInteractionBlock = indexSet->properties(vd).block;

  // Get the W and Theta maps of one of the Interaction -
  // Warning: in the current version, if OSI!=Moreau, this fails.
  // If OSI = MOREAU, centralInteractionBlocks = W if OSI = LSODAR,
  // centralInteractionBlocks = M (mass matrices)
  MapOfDSMatrices centralInteractionBlocks;
  getOSIMaps(inter , centralInteractionBlocks);

  SP::SiconosMatrix leftInteractionBlock, rightInteractionBlock;

  RELATION::TYPES relationType;
  double h = simulation()->timeDiscretisation()->currentTimeStep();

  // General form of the interactionBlock is : interactionBlock =
  // a*extraInteractionBlock + b * leftInteractionBlock * centralInteractionBlocks
  // * rightInteractionBlock a and b are scalars, centralInteractionBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  relationType = inter->getRelationType();

  inter->getExtraInteractionBlock(currentInteractionBlock);


  // loop over the DS
  bool endl = false;
  for (SP::DynamicalSystem ds = DS1; !endl; ds = DS2)
  {
    assert(ds == DS1 || ds == DS2);

    endl = (ds == DS2);
    unsigned int sizeDS = ds->getDim();

    // get _interactionBlocks corresponding to the current DS
    // These _interactionBlocks depends on the relation type.
    leftInteractionBlock.reset(new SimpleMatrix(nslawSize, sizeDS));
    inter->getLeftInteractionBlockForDS(ds, leftInteractionBlock);

    // Computing depends on relation type -> move this in Interaction method?
    if (relationType == FirstOrder)
    {

      rightInteractionBlock.reset(new SimpleMatrix(sizeDS, nslawSize));

      inter->getRightInteractionBlockForDS(ds, rightInteractionBlock);
      SP::OneStepIntegrator Osi = simulation()->integratorOfDS(ds);
      OSI::TYPES  osiType = Osi->getType();

      if (osiType == OSI::MOREAU)
      {
        if ((boost::static_pointer_cast<Moreau> (Osi))->useGamma() || (boost::static_pointer_cast<Moreau> (Osi))->useGammaForRelation())
        {
          *rightInteractionBlock *= (boost::static_pointer_cast<Moreau> (Osi))->gamma();
        }
      }


      // centralInteractionBlock contains a lu-factorized matrix and we solve
      // centralInteractionBlock * X = rightInteractionBlock with PLU
      centralInteractionBlocks[ds->number()]->PLUForwardBackwardInPlace(*rightInteractionBlock);

      //      integration of r with theta method removed
      //      *currentInteractionBlock += h *Theta[*itDS]* *leftInteractionBlock * (*rightInteractionBlock); //left = C, right = W.B
      //gemm(h,*leftInteractionBlock,*rightInteractionBlock,1.0,*currentInteractionBlock);
      *leftInteractionBlock *= h;

      prod(*leftInteractionBlock, *rightInteractionBlock, *currentInteractionBlock, false);
      //left = C, right = inv(W).B

    }
    else if (relationType == Lagrangian ||
             relationType == NewtonEuler)
    {

      Type::Siconos dsType = Type::value(*ds);

      if (dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
      {
        SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

        if (d->boundaryConditions())
        {
          for (vector<unsigned int>::iterator itindex =
                 d->boundaryConditions()->velocityIndices()->begin() ;
               itindex != d->boundaryConditions()->velocityIndices()->end();
               ++itindex)
          {
            // (nslawSize,sizeDS));
            SP::SiconosVector coltmp(new SimpleVector(nslawSize));
            coltmp->zero();
            leftInteractionBlock->setCol(*itindex, *coltmp);
          }
        }
      }

      // (inter1 == inter2)
      SP::SiconosMatrix work(new SimpleMatrix(*leftInteractionBlock));
      //
      //        cout<<"LinearOSNS : leftUBlock\n";
      //        work->display();
      work->trans();
      //        cout<<"LinearOSNS::computeInteractionBlock leftInteractionBlock"<<endl;
      //        leftInteractionBlock->display();
      centralInteractionBlocks[ds->number()]->PLUForwardBackwardInPlace(*work);
      //*currentInteractionBlock +=  *leftInteractionBlock ** work;


      prod(*leftInteractionBlock, *work, *currentInteractionBlock, false);
      //      gemm(CblasNoTrans,CblasNoTrans,1.0,*leftInteractionBlock,*work,1.0,*currentInteractionBlock);
      //*currentInteractionBlock *=h;
      //        cout<<"LinearOSNS::computeInteractionBlock interactionBlock"<<endl;
      //        currentInteractionBlock->display();
#ifdef LINEAROSNS_DEBUG
      std::cout << "LinearOSNS::computeDiagonalInteractionBlock : DiagInteractionBlock" << std::endl;
      currentInteractionBlock->display();
#endif
    }


    else RuntimeException::selfThrow("LinearOSNS::computeInteractionBlock not yet implemented for relation of type " + relationType);
  }

}

void LinearOSNS::computeInteractionBlock(const InteractionsGraph::EDescriptor& ed)
{

  // Computes matrix _interactionBlocks[inter1][inter2] (and allocates memory if
  // necessary) if inter1 and inter2 have commond DynamicalSystem.  How
  // _interactionBlocks are computed depends explicitely on the type of
  // Relation of each Interaction.

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

  // Get dimension of the NonSmoothLaw (ie dim of the interactionBlock)
  SP::InteractionsGraph indexSet = simulation()->indexSet(_levelMin);

  SP::DynamicalSystem ds = indexSet->bundle(ed);
  SP::Interaction inter1 = indexSet->bundle(indexSet->source(ed));
  SP::Interaction inter2 = indexSet->bundle(indexSet->target(ed));

  unsigned int index1 = indexSet->index(indexSet->source(ed));
  unsigned int index2 = indexSet->index(indexSet->target(ed));

  unsigned int nslawSize1 = inter1->getNonSmoothLawSize();
  unsigned int nslawSize2 = inter2->getNonSmoothLawSize();

  /*
    DynamicalSystemsSet commonDS;
    intersection(*inter1->dynamicalSystems(),*inter2->dynamicalSystems(), commonDS);
    assert (!commonDS.isEmpty()) ;
    for (DSIterator itDS = commonDS.begin(); itDS!=commonDS.end(); itDS++)
    {
    assert (*itDS == ds);
    }
  */

  SP::SiconosMatrix currentInteractionBlock;

  assert(index1 != index2);

  if (index2 > index1) // upper block
  {
    //     if (! indexSet->properties(ed).upper_block)
    //     {
    //       indexSet->properties(ed).upper_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
    //     }

    assert(indexSet->properties(ed).upper_block->size(0) == nslawSize1);
    assert(indexSet->properties(ed).upper_block->size(1) == nslawSize2);

    currentInteractionBlock = indexSet->properties(ed).upper_block;
  }
  else  // lower block
  {
    //     if (! indexSet->properties(ed).lower_block)
    //     {
    //       indexSet->properties(ed).lower_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
    //     }

    assert(indexSet->properties(ed).lower_block->size(0) == nslawSize1);
    assert(indexSet->properties(ed).lower_block->size(1) == nslawSize2);

    currentInteractionBlock = indexSet->properties(ed).lower_block;
  }

  // Get the W and Theta maps of one of the Interaction -
  // Warning: in the current version, if OSI!=Moreau, this fails.
  // If OSI = MOREAU, centralInteractionBlocks = W if OSI = LSODAR,
  // centralInteractionBlocks = M (mass matrices)
  MapOfDSMatrices centralInteractionBlocks;
  getOSIMaps(inter1, centralInteractionBlocks);

  SP::SiconosMatrix leftInteractionBlock, rightInteractionBlock;

  RELATION::TYPES relationType1, relationType2;
  double h = simulation()->timeDiscretisation()->currentTimeStep();

  // General form of the interactionBlock is : interactionBlock =
  // a*extraInteractionBlock + b * leftInteractionBlock * centralInteractionBlocks
  // * rightInteractionBlock a and b are scalars, centralInteractionBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  relationType1 = inter1->getRelationType();
  relationType2 = inter2->getRelationType();


  // ==== First Order Relations - Specific treatment for diagonal
  // _interactionBlocks ===
  assert(inter1 != inter2);
  //  currentInteractionBlock->zero();

  // loop over the common DS
  unsigned int sizeDS = ds->getDim();

  // get _interactionBlocks corresponding to the current DS
  // These _interactionBlocks depends on the relation type.
  leftInteractionBlock.reset(new SimpleMatrix(nslawSize1, sizeDS));
  inter1->getLeftInteractionBlockForDS(ds, leftInteractionBlock);

  // Computing depends on relation type -> move this in Interaction method?
  if (relationType1 == FirstOrder && relationType2 == FirstOrder)
  {

    rightInteractionBlock.reset(new SimpleMatrix(sizeDS, nslawSize2));

    inter2->getRightInteractionBlockForDS(ds, rightInteractionBlock);
    // centralInteractionBlock contains a lu-factorized matrix and we solve
    // centralInteractionBlock * X = rightInteractionBlock with PLU
    centralInteractionBlocks[ds->number()]->PLUForwardBackwardInPlace(*rightInteractionBlock);

    //      integration of r with theta method removed
    //      *currentInteractionBlock += h *Theta[*itDS]* *leftInteractionBlock * (*rightInteractionBlock); //left = C, right = W.B
    //gemm(h,*leftInteractionBlock,*rightInteractionBlock,1.0,*currentInteractionBlock);
    *leftInteractionBlock *= h;

    prod(*leftInteractionBlock, *rightInteractionBlock, *currentInteractionBlock, false);
    //left = C, right = inv(W).B


  }
  else if (relationType1 == Lagrangian ||
           relationType2 == Lagrangian ||
           relationType1 == NewtonEuler ||
           relationType2 == NewtonEuler)
  {

    Type::Siconos dsType = Type::value(*ds);

    if (dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
    {
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      if (d->boundaryConditions())
      {
        for (vector<unsigned int>::iterator itindex =
               d->boundaryConditions()->velocityIndices()->begin() ;
             itindex != d->boundaryConditions()->velocityIndices()->end();
             ++itindex)
        {
          // (nslawSize1,sizeDS));
          SP::SiconosVector coltmp(new SimpleVector(nslawSize1));
          coltmp->zero();
          leftInteractionBlock->setCol(*itindex, *coltmp);
        }
      }
    }

    // inter1 != inter2
    rightInteractionBlock.reset(new SimpleMatrix(nslawSize2, sizeDS));
    inter2->getLeftInteractionBlockForDS(ds, rightInteractionBlock);

    // Warning: we use getLeft for Right interactionBlock
    // because right = transpose(left) and because of
    // size checking inside the getBlock function, a
    // getRight call will fail.
    rightInteractionBlock->trans();
    centralInteractionBlocks[ds->number()]->PLUForwardBackwardInPlace(*rightInteractionBlock);
    //*currentInteractionBlock +=  *leftInteractionBlock ** work;
    prod(*leftInteractionBlock, *rightInteractionBlock, *currentInteractionBlock, false);
#ifdef LINEAROSNS_DEBUG
    std::cout << "LinearOSNS::computeInteractionBlock : currentInteractionBlock" << std::endl;
    currentInteractionBlock->display();
#endif
  }
  else RuntimeException::selfThrow("LinearOSNS::computeInteractionBlock not yet implemented for relation of type " + relationType1);

}



void LinearOSNS::computeqBlock(SP::Interaction inter, unsigned int pos)
{

  SP::DynamicalSystem ds = *(inter->dynamicalSystemsBegin());
  SP::OneStepIntegrator Osi = simulation()->integratorOfDS(ds);
  OSI::TYPES  osiType = Osi->getType();


  unsigned int sizeY = inter->getNonSmoothLawSize();
  Index coord(8);
  unsigned int relativePosition = 0;

  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = pos;
  coord[7] = pos + sizeY;

  SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();



  if (osiType == OSI::MOREAU ||
      osiType == OSI::LSODAR ||
      osiType == OSI::D1MINUSLINEAR ||
      osiType == OSI::SCHATZMANPAOLI)
  {
    Osi->computeFreeOutput(inter , this);
    setBlock(*inter->yp(), _q, sizeY , 0, pos);
  }
  else if (osiType == OSI::MOREAU2)
  {

  }
  else
    RuntimeException::selfThrow("LinearOSNS::computeqBlock not yet implemented for OSI of type " + osiType);

  // Add "non-smooth law effect" on q only for the case LCP at velocity level and with the NewtonImpactNSL
  if (osiType != OSI::LSODAR || ((*allOSNS)[SICONOS_OSNSP_ED_IMPACT]).get() == this) // added by Son Nguyen
  {
    if (inter->getRelationType() == Lagrangian || inter->getRelationType() == NewtonEuler)
    {
      //    SP::SiconosVisitor nslEffectOnSim(new _NSLEffectOnSim(this, inter, pos));
      //       simulation()->accept(*nslEffectOnSim);
    }
  }
}

void LinearOSNS::computeq(double time)
{
  if (_q->size() != _sizeOutput)
    _q->resize(_sizeOutput);
  _q->zero();

  // === Get index set from Simulation ===
  SP::InteractionsGraph indexSet =
    simulation()->indexSet(levelMin());
  // === Loop through "active" Interactions (ie present in
  // indexSets[level]) ===

  unsigned int pos = 0;
  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet->bundle(*ui);

    // Compute q, this depends on the type of non smooth problem, on
    // the relation type and on the non smooth law
    pos = _M->getPositionOfInteractionBlock(inter);
    computeqBlock(inter, pos); // free output is saved in y
  }
}



void LinearOSNS::preCompute(double time)
{
  // This function is used to prepare data for the
  // LinearComplementarityProblem

  // - computation of M and q
  // - set _sizeOutput
  // - check dim. for _z,_w

  // If the topology is time-invariant, only q needs to be computed at
  // each time step.  M, _sizeOutput have been computed in initialize
  // and are uptodate.

  // Get topology
  SP::Topology topology = simulation()->model()
                          ->nonSmoothDynamicalSystem()->topology();
  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();

  //  std::cout << "!b || !isLinear :"  << boolalpha <<  (!b || !isLinear) <<  std::endl;

  // nothing to do
  if (_levelMin == LEVELMAX)
    return;

  if (!_hasBeUpdated || !isLinear)
  {
    // Computes new _interactionBlocks if required
    updateInteractionBlocks();

    // Updates matrix M
    SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin());
    //    _M->fill(indexSet);
    _M->fill(indexSet, !_hasBeUpdated);
    //      updateOSNSMatrix();
    _sizeOutput = _M->size();

    // Checks z and _w sizes and reset if necessary
    if (_z->size() != _sizeOutput)
    {
      _z->resize(_sizeOutput, false);
      _z->zero();
    }

    if (_w->size() != _sizeOutput)
    {
      _w->resize(_sizeOutput);
      _w->zero();
    }

    // _w and _z <- old values. Note : sizeOuput can be unchanged,
    // but positions may have changed
    if (_keepLambdaAndYState)
    {
      InteractionsGraph::VIterator ui, uiend;
      for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
      {
        SP::Interaction inter = indexSet->bundle(*ui);
        // Get the relative position of inter-interactionBlock in the vector w
        // or z
        unsigned int pos = _M->getPositionOfInteractionBlock(inter);

        SPC::SiconosVector yOld = inter->yOld(levelMin());
        SPC::SiconosVector lambdaOld = inter->lambdaOld(levelMin());

        setBlock(*yOld, _w, yOld->size(), 0, pos);
        setBlock(*lambdaOld, _z, lambdaOld->size(), 0, pos);
      }
    }
  }
  else
  {
    // nothing to do (IsLinear and not changed)
  }


  // Computes q of LinearOSNS
  computeq(time);

}

void LinearOSNS::postCompute()
{
  // This function is used to set y/lambda values using output from
  // lcp_driver (w,z).  Only Interactions (ie Interactions) of
  // indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin());

  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;

  // === Loop through "active" Interactions (ie present in
  // indexSets[1]) ===

  unsigned int pos = 0;

  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet->bundle(*ui);
    // Get the relative position of inter-interactionBlock in the vector w
    // or z
    pos = _M->getPositionOfInteractionBlock(inter);

    // Get Y and Lambda for the current Interaction
    y = inter->y(levelMin());
    lambda = inter->lambda(levelMin());
    // Copy _w/_z values, starting from index pos into y/lambda.

    setBlock(*_w, y, y->size(), pos, 0);// Warning: yEquivalent is
    // saved in y !!
    setBlock(*_z, lambda, lambda->size(), pos, 0);
  }

}

void LinearOSNS::display() const
{
  cout << "_M  ";
  if (_M) _M->display();
  else cout << "-> NULL" << endl;
  cout << endl << " q : " ;
  if (_q) _q->display();
  else cout << "-> NULL" << endl;
  cout << "w  ";
  if (_w) _w->display();
  else cout << "-> NULL" << endl;
  cout << endl << "z : " ;
  if (_z) _z->display();
  else cout << "-> NULL" << endl;
  cout << "==========================" << endl;
}

void LinearOSNS::saveNSProblemToXML()
{
  //   if(onestepnspbxml != NULL)
  //     {
  // //       (boost::static_pointer_cast<LinearOSNSXML>(onestepnspbxml))->setM(*_M);
  //       (boost::static_pointer_cast<LinearOSNSXML>(onestepnspbxml))->setQ(*q);
  //     }
  //   else RuntimeException::selfThrow("LinearOSNS::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  RuntimeException::selfThrow("LinearOSNS::saveNSProblemToXML - Not yet implemented.");
}

