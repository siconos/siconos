/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "Simulation.hpp"
#include "Topology.hpp"
#include "Model.hpp"
#include "MoreauJeanOSI.hpp"
#include "EulerMoreauOSI.hpp"
#include "LsodarOSI.hpp"
#include "NewMarkAlphaOSI.hpp"
#include "ZeroOrderHoldOSI.hpp"
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
#include "NewtonEulerDS.hpp"
#include "OSNSMatrix.hpp"

#include "Tools.hpp"

using namespace RELATION;
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

LinearOSNS::LinearOSNS(): OneStepNSProblem(), _MStorageType(0), _keepLambdaAndYState(true)
{}

// Constructor from a set of data
LinearOSNS::LinearOSNS(const int numericsSolverId):
  OneStepNSProblem(numericsSolverId), _MStorageType(0), _keepLambdaAndYState(true)
{}

// Setters

void LinearOSNS::setW(const SiconosVector& newValue)
{
  assert(_sizeOutput == newValue.size() &&
         "LinearOSNS: setW, inconsistent size between given velocity size and problem size. You should set sizeOutput before");
  setObject<SiconosVector, SP::SiconosVector, SiconosVector>(_w, newValue);
}

void LinearOSNS::setz(const SiconosVector& newValue)
{
  assert(_sizeOutput == newValue.size() &&
         "LinearOSNS: setz, inconsistent size between given velocity size and problem size. You should set sizeOutput before");
  setObject<SiconosVector, SP::SiconosVector, SiconosVector>(_z, newValue);
}

void LinearOSNS::initVectorsMemory()
{
  // Memory allocation for _w, M, z and q.
  // If one of them has already been allocated, nothing is done.
  // We suppose that user has chosen a correct size.
  if (! _w)
    _w.reset(new SiconosVector(maxSize()));
  else
  {
    if (_w->size() != maxSize())
      _w->resize(maxSize());
  }

  if (! _z)
    _z.reset(new SiconosVector(maxSize()));
  else
  {
    if (_z->size() != maxSize())
      _z->resize(maxSize());
  }

  if (! _q)
    _q.reset(new SiconosVector(maxSize()));
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
      if (indexSetLevel() != LEVELMAX && simulation()->model()->nonSmoothDynamicalSystem()->topology()->indexSetsSize() > indexSetLevel())
      {
        _M.reset(new OSNSMatrix(simulation()->indexSet(indexSetLevel())->size(), _MStorageType));
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

  // Note that _interactionBlocks is up to date since updateInteractionBlocks
  // has been called during OneStepNSProblem::initialize()

  initOSNSMatrix();

}

void LinearOSNS::computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd)
{
  DEBUG_PRINT("LinearOSNS::computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd)\n");

  // Computes matrix _interactionBlocks[inter1][inter1] (and allocates memory if
  // necessary) one or two DS are concerned by inter1 .  How
  // _interactionBlocks are computed depends explicitely on the type of
  // Relation of each Interaction.

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

  // Get dimension of the NonSmoothLaw (ie dim of the interactionBlock)
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vd);
  // Get osi property from interaction
  // We assume that all ds in vertex_inter have the same osi.
  SP::OneStepIntegrator Osi = indexSet->properties(vd).osi;
  //SP::OneStepIntegrator Osi = simulation()->integratorOfDS(ds);
  OSI::TYPES  osiType = Osi->getType();


  // At most 2 DS are linked by an Interaction
  SP::DynamicalSystem DS1;
  SP::DynamicalSystem DS2;
  unsigned int pos1, pos2;
  // --- Get the dynamical system(s) (edge(s)) connected to the current interaction (vertex) ---
  if (indexSet->properties(vd).source != indexSet->properties(vd).target)
  {
    DEBUG_PRINT("a two DS Interaction\n");
    DS1 = indexSet->properties(vd).source;
    DS2 = indexSet->properties(vd).target;
  }
  else
  {
    DEBUG_PRINT("a single DS Interaction\n");
    DS1 = indexSet->properties(vd).source;
    DS2 = DS1;
    InteractionsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = indexSet->out_edges(vd);
         oei != oeiend; ++oei)
    {
      // note : at most 4 edges
      DS2 = indexSet->bundle(*oei);
      if (DS2 != DS1)
      {
        assert(false);
        break;
      }
    }
  }
  assert(DS1);
  assert(DS2);
  pos1 = indexSet->properties(vd).source_pos;
  pos2 = indexSet->properties(vd).target_pos;

  // --- Check block size ---
  assert(indexSet->properties(vd).block->size(0) == inter->nonSmoothLaw()->size());
  assert(indexSet->properties(vd).block->size(1) == inter->nonSmoothLaw()->size());

  // --- Compute diagonal block ---
  // Block to be set in OSNS Matrix, corresponding to
  // the current interaction
  SP::SiconosMatrix currentInteractionBlock = indexSet->properties(vd).block;
  SP::SiconosMatrix leftInteractionBlock, rightInteractionBlock;

  RELATION::TYPES relationType;
  double h = simulation()->currentTimeStep();

  // General form of the interactionBlock is : interactionBlock =
  // a*extraInteractionBlock + b * leftInteractionBlock * centralInteractionBlocks
  // * rightInteractionBlock a and b are scalars, centralInteractionBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  relationType = inter->relation()->getType();
  VectorOfSMatrices& workMInter = *indexSet->properties(vd).workMatrices;

  inter->getExtraInteractionBlock(currentInteractionBlock, workMInter);

  unsigned int nslawSize = inter->nonSmoothLaw()->size();
  // loop over the DS connected to the interaction.
  bool endl = false;
  unsigned int pos = pos1;
  for (SP::DynamicalSystem ds = DS1; !endl; ds = DS2)
  {
    assert(ds == DS1 || ds == DS2);
    endl = (ds == DS2);
    unsigned int sizeDS = ds->getDim();
    // get _interactionBlocks corresponding to the current DS
    // These _interactionBlocks depends on the relation type.
    leftInteractionBlock.reset(new SimpleMatrix(nslawSize, sizeDS));
    inter->getLeftInteractionBlockForDS(pos, leftInteractionBlock, workMInter);
    DEBUG_EXPR(leftInteractionBlock->display(););
    // Computing depends on relation type -> move this in Interaction method?
    if (relationType == FirstOrder)
    {

      rightInteractionBlock.reset(new SimpleMatrix(sizeDS, nslawSize));

      inter->getRightInteractionBlockForDS(pos, rightInteractionBlock, workMInter);

      if (osiType == OSI::EULERMOREAUOSI)
      {
        if ((std11::static_pointer_cast<EulerMoreauOSI> (Osi))->useGamma() || (std11::static_pointer_cast<EulerMoreauOSI> (Osi))->useGammaForRelation())
        {
          *rightInteractionBlock *= (std11::static_pointer_cast<EulerMoreauOSI> (Osi))->gamma();
        }
      }

      // for ZOH, we have a different formula ...
      if (osiType == OSI::ZOHOSI && indexSet->properties(vd).forControl)
      {
        *rightInteractionBlock = std11::static_pointer_cast<ZeroOrderHoldOSI>(Osi)->Bd(ds);
        prod(*leftInteractionBlock, *rightInteractionBlock, *currentInteractionBlock, false);
      }
      else
      {
        // centralInteractionBlock contains a lu-factorized matrix and we solve
        // centralInteractionBlock * X = rightInteractionBlock with PLU
        SP::SiconosMatrix centralInteractionBlock = getOSIMatrix(Osi, ds);
        centralInteractionBlock->PLUForwardBackwardInPlace(*rightInteractionBlock);
        inter->computeKhat(*rightInteractionBlock, workMInter, h); // if K is non 0

        //      integration of r with theta method removed
        //      *currentInteractionBlock += h *Theta[*itDS]* *leftInteractionBlock * (*rightInteractionBlock); //left = C, right = W.B
        //gemm(h,*leftInteractionBlock,*rightInteractionBlock,1.0,*currentInteractionBlock);
        *leftInteractionBlock *= h;
        prod(*leftInteractionBlock, *rightInteractionBlock, *currentInteractionBlock, false);
        //left = C, right = inv(W).B
      }

    }
    else if (relationType == Lagrangian ||
             relationType == NewtonEuler)
    {

      SP::BoundaryCondition bc;
      Type::Siconos dsType = Type::value(*ds);
      if (dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
      {
        SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
        if (d->boundaryConditions()) bc = d->boundaryConditions();
      }
      else if (dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
        if (d->boundaryConditions()) bc = d->boundaryConditions();
      }
      if (bc)
      {
        for (std::vector<unsigned int>::iterator itindex = bc->velocityIndices()->begin() ;
             itindex != bc->velocityIndices()->end();
             ++itindex)
        {
          // (nslawSize,sizeDS));
          SP::SiconosVector coltmp(new SiconosVector(nslawSize));
          coltmp->zero();
          leftInteractionBlock->setCol(*itindex, *coltmp);
        }
      }
      DEBUG_PRINT("leftInteractionBlock after application of boundary conditions\n");
      DEBUG_EXPR(leftInteractionBlock->display(););


      // (inter1 == inter2)
      DEBUG_EXPR(leftInteractionBlock->display(););
      SP::SiconosMatrix work(new SimpleMatrix(*leftInteractionBlock));
      work->trans();
      SP::SiconosMatrix centralInteractionBlock = getOSIMatrix(Osi, ds);
      DEBUG_EXPR(centralInteractionBlock->display(););
      DEBUG_EXPR_WE(std::cout <<  std::boolalpha << " centralInteractionBlock->isPLUFactorized() = "<< centralInteractionBlock->isPLUFactorized() << std::endl;);
      centralInteractionBlock->PLUForwardBackwardInPlace(*work);
      //*currentInteractionBlock +=  *leftInteractionBlock ** work;
      prod(*leftInteractionBlock, *work, *currentInteractionBlock, false);
      //      gemm(CblasNoTrans,CblasNoTrans,1.0,*leftInteractionBlock,*work,1.0,*currentInteractionBlock);
      //*currentInteractionBlock *=h;
    }
    else RuntimeException::selfThrow("LinearOSNS::computeInteractionBlock not yet implemented for relation of type " + relationType);
    // Set pos for next loop.
    pos = pos2;
  }
  DEBUG_PRINT("LinearOSNS::computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd) ends \n");
}

void LinearOSNS::computeInteractionBlock(const InteractionsGraph::EDescriptor& ed)
{
  
  DEBUG_PRINT("LinearOSNS::computeInteractionBlock(const InteractionsGraph::EDescriptor& ed)\n");

  // Computes matrix _interactionBlocks[inter1][inter2] (and allocates memory if
  // necessary) if inter1 and inter2 have commond DynamicalSystem.  How
  // _interactionBlocks are computed depends explicitely on the type of
  // Relation of each Interaction.

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

  // Get dimension of the NonSmoothLaw (ie dim of the interactionBlock)
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  SP::DynamicalSystem ds = indexSet->bundle(ed);
  SP::Interaction inter1 = indexSet->bundle(indexSet->source(ed));
  SP::Interaction inter2 = indexSet->bundle(indexSet->target(ed));
  // Once again we assume that inter1 and inter2 have the same osi ...
  SP::OneStepIntegrator Osi = indexSet->properties(indexSet->source(ed)).osi;

  // For the edge 'ds', we need to find relative position of this ds
  // in inter1 and inter2 relation matrices (--> pos1 and pos2 below)
  // - find if ds is source or target in inter_i
  InteractionsGraph::VDescriptor vertex_inter;
  // - get the corresponding position
  unsigned int pos1, pos2;
  // source of inter1 :
  vertex_inter = indexSet->source(ed);
  VectorOfSMatrices& workMInter1 = *indexSet->properties(vertex_inter).workMatrices;
  SP::DynamicalSystem tmpds = indexSet->properties(vertex_inter).source;
  if (tmpds == ds)
    pos1 =  indexSet->properties(vertex_inter).source_pos;
  else
  {
    tmpds  = indexSet->properties(vertex_inter).target;
    pos1 =  indexSet->properties(vertex_inter).target_pos;
  }
  // now, inter2
  vertex_inter = indexSet->target(ed);
  VectorOfSMatrices& workMInter2 = *indexSet->properties(vertex_inter).workMatrices;
  tmpds = indexSet->properties(vertex_inter).source;
  if (tmpds == ds)
    pos2 =  indexSet->properties(vertex_inter).source_pos;
  else
  {
    tmpds  = indexSet->properties(vertex_inter).target;
    pos2 =  indexSet->properties(vertex_inter).target_pos;
  }

  unsigned int index1 = indexSet->index(indexSet->source(ed));
  unsigned int index2 = indexSet->index(indexSet->target(ed));
  unsigned int nslawSize1 = inter1->nonSmoothLaw()->size();
  unsigned int nslawSize2 = inter2->nonSmoothLaw()->size();

  SP::SiconosMatrix currentInteractionBlock;

  assert(index1 != index2);

  if (index2 > index1) // upper block
  {
    assert(indexSet->properties(ed).upper_block->size(0) == nslawSize1);
    assert(indexSet->properties(ed).upper_block->size(1) == nslawSize2);

    currentInteractionBlock = indexSet->properties(ed).upper_block;
  }
  else  // lower block
  {
    assert(indexSet->properties(ed).lower_block->size(0) == nslawSize1);
    assert(indexSet->properties(ed).lower_block->size(1) == nslawSize2);

    currentInteractionBlock = indexSet->properties(ed).lower_block;
  }

  SP::SiconosMatrix leftInteractionBlock, rightInteractionBlock;

  RELATION::TYPES relationType1, relationType2;
  double h = simulation()->currentTimeStep();

  // General form of the interactionBlock is : interactionBlock =
  // a*extraInteractionBlock + b * leftInteractionBlock * centralInteractionBlocks
  // * rightInteractionBlock a and b are scalars, centralInteractionBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  relationType1 = inter1->relation()->getType();
  relationType2 = inter2->relation()->getType();

  // ==== First Order Relations - Specific treatment for diagonal
  // _interactionBlocks ===
  assert(inter1 != inter2);
  //  currentInteractionBlock->zero();


  // loop over the common DS
  unsigned int sizeDS = ds->getDim();

  // get _interactionBlocks corresponding to the current DS
  // These _interactionBlocks depends on the relation type.
  leftInteractionBlock.reset(new SimpleMatrix(nslawSize1, sizeDS));
  inter1->getLeftInteractionBlockForDS(pos1, leftInteractionBlock, workMInter1);

  // Computing depends on relation type -> move this in Interaction method?
  if (relationType1 == FirstOrder && relationType2 == FirstOrder)
  {

    rightInteractionBlock.reset(new SimpleMatrix(sizeDS, nslawSize2));

    inter2->getRightInteractionBlockForDS(pos2, rightInteractionBlock, workMInter2);
    // centralInteractionBlock contains a lu-factorized matrix and we solve
    // centralInteractionBlock * X = rightInteractionBlock with PLU
    SP::SiconosMatrix centralInteractionBlock = getOSIMatrix(Osi, ds);
    centralInteractionBlock->PLUForwardBackwardInPlace(*rightInteractionBlock);

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

    //Type::Siconos dsType = Type::value(*ds);


    // if (d->boundaryConditions())
    // {
    //   for (std::vector<unsigned int>::iterator itindex =
    //          d->boundaryConditions()->velocityIndices()->begin() ;
    //        itindex != d->boundaryConditions()->velocityIndices()->end();
    //        ++itindex)
    //   {
    //     // (nslawSize1,sizeDS));
    //     SP::SiconosVector coltmp(new SiconosVector(nslawSize1));
    //     coltmp->zero();
    //     leftInteractionBlock->setCol(*itindex, *coltmp);
    //   }
    // }

    SP::BoundaryCondition bc;
    Type::Siconos dsType = Type::value(*ds);
    if (dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
      if (d->boundaryConditions()) bc = d->boundaryConditions();
    }
    else if (dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
      if (d->boundaryConditions()) bc = d->boundaryConditions();
    }
    if (bc)
    {
      for (std::vector<unsigned int>::iterator itindex = bc->velocityIndices()->begin() ;
           itindex != bc->velocityIndices()->end();
           ++itindex)
      {
        // (nslawSize,sizeDS));
        SP::SiconosVector coltmp(new SiconosVector(nslawSize1));
        coltmp->zero();
        leftInteractionBlock->setCol(*itindex, *coltmp);
      }
    }

    // inter1 != inter2
    rightInteractionBlock.reset(new SimpleMatrix(nslawSize2, sizeDS));
    inter2->getLeftInteractionBlockForDS(pos2, rightInteractionBlock, workMInter2);

    // Warning: we use getLeft for Right interactionBlock
    // because right = transpose(left) and because of
    // size checking inside the getBlock function, a
    // getRight call will fail.
    rightInteractionBlock->trans();
    SP::SimpleMatrix centralInteractionBlock = getOSIMatrix(Osi, ds);
    centralInteractionBlock->PLUForwardBackwardInPlace(*rightInteractionBlock);
    //*currentInteractionBlock +=  *leftInteractionBlock ** work;
    prod(*leftInteractionBlock, *rightInteractionBlock, *currentInteractionBlock, false);
  }
  else RuntimeException::selfThrow("LinearOSNS::computeInteractionBlock not yet implemented for relation of type " + relationType1);

}

void LinearOSNS::computeqBlock(InteractionsGraph::VDescriptor& vertex_inter, unsigned int pos)
{
  DEBUG_PRINT("LinearOSNS::computeqBlock(SP::Interaction inter, unsigned int pos)\n");
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());

  // Get the osi for the concerned interaction. Note FP: this must be a property
  // of ds rather than interaction?
  SP::OneStepIntegrator Osi = indexSet->properties(vertex_inter).osi;
  //SP::OneStepIntegrator Osi = simulation()->integratorOfDS(ds);
  SP::Interaction inter = indexSet->bundle(vertex_inter);
  OSI::TYPES  osiType = Osi->getType();
  unsigned int sizeY = inter->nonSmoothLaw()->size();

  if (osiType == OSI::EULERMOREAUOSI ||
      osiType == OSI::MOREAUJEANOSI ||
      osiType == OSI::MOREAUDIRECTPROJECTIONOSI ||
      osiType == OSI::LSODAROSI ||
      osiType == OSI::NEWMARKALPHAOSI ||
      osiType == OSI::D1MINUSLINEAROSI ||
      osiType == OSI::SCHATZMANPAOLIOSI ||
      osiType == OSI::ZOHOSI)
  {
    Osi->computeFreeOutput(vertex_inter, this);
    setBlock(*inter->yForNSsolver(), _q, sizeY , 0, pos);
    DEBUG_EXPR(_q->display());
  }
  else if (osiType == OSI::MOREAUJEANOSI2)
  {

  }
  else
    RuntimeException::selfThrow("LinearOSNS::computeqBlock not yet implemented for OSI of type " + osiType);

}

void LinearOSNS::computeq(double time)
{
  if (_q->size() != _sizeOutput)
    _q->resize(_sizeOutput);
  _q->zero();

  // === Get index set from Simulation ===
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  // === Loop through "active" Interactions (ie present in
  // indexSets[level]) ===

  unsigned int pos = 0;
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    Interaction& inter = *indexSet->bundle(*ui);

    // Compute q, this depends on the type of non smooth problem, on
    // the relation type and on the non smooth law
    pos = _M->getPositionOfInteractionBlock(inter);
    computeqBlock(*ui, pos); // free output is saved in y
  }
}



bool LinearOSNS::preCompute(double time)
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
  SP::Topology topology = simulation()->model()->nonSmoothDynamicalSystem()->topology();
  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();

  //   std::cout << "!b || !isLinear :"  << boolalpha <<  (!b || !isLinear) <<  std::endl;

  // nothing to do
  if (indexSetLevel() == LEVELMAX)
    return false;

  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  assert(indexSet);
  if (indexSet->size() == 0)
    return false;

  if (!_hasBeenUpdated || !isLinear)
  {
    // Computes new _interactionBlocks if required
    updateInteractionBlocks();

    //    _M->fill(indexSet);
    _M->fill(indexSet, !_hasBeenUpdated);
    DEBUG_EXPR(_M->display(););

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
      for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
      {
        Interaction& inter = *indexSet->bundle(*ui);
        // Get the relative position of inter-interactionBlock in the vector w
        // or z
        unsigned int pos = _M->getPositionOfInteractionBlock(inter);

        SiconosVector& yOutputOld = *inter.yOld(inputOutputLevel());
        SiconosVector& lambdaOld = *inter.lambdaOld(inputOutputLevel());


        if (_sizeOutput >= yOutputOld.size() + pos)
        {
          setBlock(yOutputOld, _w, yOutputOld.size(), 0, pos);
          setBlock(lambdaOld, _z, lambdaOld.size(), 0, pos);
        }
        else
        {
          //std::cout << std::endl;
          //std::cout << "LinearOSNS::preCompute FIXME: Old variables are bigger than the LinearOSNS variables w and z !" << std::endl;
        }

      }
    }
  }
  else
  {
    // nothing to do (IsLinear and not changed)
  }


  // Computes q of LinearOSNS
  computeq(time);

  return true;

}

void LinearOSNS::postCompute()
{
  // This function is used to set y/lambda values using output from
  // lcp_driver (w,z).  Only Interactions (ie Interactions) of
  // indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());

  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;

  // === Loop through "active" Interactions (ie present in
  // indexSets[1]) ===

  unsigned int pos = 0;

  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    Interaction& inter = *indexSet->bundle(*ui);
    // Get the relative position of inter-interactionBlock in the vector w
    // or z
    pos = _M->getPositionOfInteractionBlock(inter);

    // Get Y and Lambda for the current Interaction
    y = inter.y(inputOutputLevel());
    lambda = inter.lambda(inputOutputLevel());
    // Copy _w/_z values, starting from index pos into y/lambda.

    //setBlock(*_w, y, y->size(), pos, 0);// Warning: yEquivalent is
    // saved in y !!
    setBlock(*_z, lambda, lambda->size(), pos, 0);
    DEBUG_EXPR(lambda->display(););
  }

}

void LinearOSNS::display() const
{
  std::cout << "==========================" <<std::endl;
  std::cout << "_M  ";
  if (_M) _M->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout <<std::endl << "q : " ;
  if (_q) _q->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << std::endl;
  std::cout << "w : ";
  if (_w) _w->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout <<std::endl << "z : " ;
  if (_z) _z->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << std::endl;
  std::cout << "The linearOSNSP works on the index set of level  " << _indexSetLevel<< std::endl;
  std::cout << "==========================" <<std::endl;
}
