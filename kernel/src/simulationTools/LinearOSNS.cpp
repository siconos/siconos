/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, softwareﬁ
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "LinearOSNS.hpp"
#include "Simulation.hpp"
#include "Topology.hpp"
#include "Model.hpp"
#include "MoreauJeanOSI.hpp"
#include "MoreauJeanBilbaoOSI.hpp"
#include "D1MinusLinearOSI.hpp"
#include "EulerMoreauOSI.hpp"
#include "SchatzmanPaoliOSI.hpp"
#include "LsodarOSI.hpp"
#include "NewMarkAlphaOSI.hpp"
#include "ZeroOrderHoldOSI.hpp"
#include "NewtonEulerR.hpp"
#include "FirstOrderLinearR.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "LagrangianLinearTIR.hpp"
#include "LagrangianCompliantLinearTIR.hpp"
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
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

LinearOSNS::LinearOSNS(): OneStepNSProblem(), _numericsMatrixStorageType(NM_DENSE), _keepLambdaAndYState(true)
{
}

// Constructor from a set of data
LinearOSNS::LinearOSNS(const int numericsSolverId):
  OneStepNSProblem(numericsSolverId), _numericsMatrixStorageType(0), _keepLambdaAndYState(true)
{}

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
      switch (_numericsMatrixStorageType)
	{
	case NM_DENSE:
	case NM_SPARSE:
	  {
	    _M.reset(new OSNSMatrix(maxSize(), _numericsMatrixStorageType));
	    break;
	  }
	case NM_SPARSE_BLOCK:
	  {
	    // = number of Interactionin the largest considered indexSet
	    if (indexSetLevel() != LEVELMAX && simulation()->nonSmoothDynamicalSystem()->topology()->indexSetsSize() > indexSetLevel())
	      {
		_M.reset(new OSNSMatrix(simulation()->indexSet(indexSetLevel())->size(), _numericsMatrixStorageType));
	      }
	    else
	      {
		_M.reset(new OSNSMatrix(1, _numericsMatrixStorageType));
	      }
	    break;
	  }
	  {
	    default:
	      RuntimeException::selfThrow("LinearOSNS::initOSNSMatrix unknown _storageType");
	  }
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
  DEBUG_BEGIN("LinearOSNS::computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd)\n");

  // Compute diagonal block (graph property) for a given interaction.
  // - How  blocks are computed depends explicitely on the nonsmooth law, the relation and the dynamical system types.
  // - No memory allocation there : blocks should be previously (updateInteractionBlocks) and properly allocated.

  // -- Get index sets from simulation --
  // Active interactions
  InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());
  // All declared interactions, where properties are saved.
  InteractionsGraph& parentSet = *simulation()->indexSet(0);
  // --- Get the dynamical system(s) (edge(s)) connected to the current interaction (vertex) ---
  // Remark : ds1 == ds2 is allowed.
  SP::DynamicalSystem ds1 = indexSet.properties(vd).source;
  SP::DynamicalSystem ds2 = indexSet.properties(vd).target;
  assert(ds1);
  assert(ds2);
  // Get dynamical systems positions in the interaction
  unsigned int pos1 = indexSet.properties(vd).source_pos;
  unsigned int pos2 = indexSet.properties(vd).target_pos;
  
  // --- Diagonal block will be computed into 'block' property from current vertex. ---
  // Interaction corresponding to current vertex 
  SP::Interaction inter = indexSet.bundle(vd);
  
  // and the vertex for this interaction in parent set
  InteractionsGraph::VDescriptor vd0 = parentSet.descriptor(inter);

  // Get block from index set and check if shape is consistent with the interaction.
  SiconosMatrix& currentInteractionBlock = *parentSet.properties(vd0).block;
  unsigned int nslawSize = inter->nonSmoothLaw()->size();
  assert(currentInteractionBlock.size(0) == nslawSize);
  assert(currentInteractionBlock.size(1) == nslawSize);
  
  SP::SiconosMatrix leftInteractionBlock, rightInteractionBlock;

  // Get interaction properties (relation type)
  RELATION::TYPES relationType = inter->relation()->getType();
  RELATION::SUBTYPES relationSubType= inter->relation()->getSubType();
  double h = simulation()->currentTimeStep();

  // General form of the diagonal block is :
  // block = a * extraInteractionBlock
  //       + b * leftInteractionBlock * centralInteractionBlocks * rightInteractionBlock
  // - a, b : scalars
  // - centralInteractionBlocks : matrix which depends on integrator and dynamical system
  // - leftInteractionBlock, rightInteractionBlock : matrices depending on the interaction and on the dynamical system.
  

  VectorOfSMatrices& workMInter = *indexSet.properties(vd).workMatrices;
  // Get extraInteractionBlock. May be null.
  // Saved in currentInteractionBlock.
  // workMinter is some internal buffer for interaction properties. Read only
  inter->getExtraInteractionBlock(currentInteractionBlock, workMInter);


  // Dynamical systems graph. Required to get iteration matrix from the osi.
  DynamicalSystemsGraph& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

  // loop over the DS connected to the interaction.
  bool endl = false;
  unsigned int pos = pos1;
  for (SP::DynamicalSystem ds = ds1; !endl; ds = ds2)
    {
      assert(ds == ds1 || ds == ds2);
      endl = (ds == ds2);
      OneStepIntegrator& osi = *DSG0.properties(DSG0.descriptor(ds)).osi;
      OSI::TYPES osiType = osi.getType();
      unsigned int sizeDS = ds->dimension();

      // TEMP BUFFER for leftInteractionBlock
      leftInteractionBlock.reset(new SimpleMatrix(nslawSize, sizeDS));

      // Block in relation matrix corresponding current DS is saved in leftInteractionBlock.
      inter->getLeftInteractionBlockForDS(pos, leftInteractionBlock, workMInter);

      // Computing depends on relation type -> move this in Interaction method?
      if (relationType == FirstOrder)
	{

	  // TEMP BUFFER for rightInteractionBlock
	  rightInteractionBlock.reset(new SimpleMatrix(sizeDS, nslawSize));
	  inter->getRightInteractionBlockForDS(pos, rightInteractionBlock, workMInter);

	  if (osiType == OSI::EULERMOREAUOSI)
	    {
	      if ((static_cast<EulerMoreauOSI&> (osi)).useGamma() || (static_cast<EulerMoreauOSI&> (osi)).useGammaForRelation())
		{
		  *rightInteractionBlock *= (static_cast<EulerMoreauOSI&> (osi)).gamma();
		}
	    }

	  // for ZOH, we have a different formula ...
	  if ((osiType == OSI::ZOHOSI) && indexSet.properties(vd).forControl)
	    {
	      *rightInteractionBlock = static_cast<ZeroOrderHoldOSI&>(osi).Bd(ds);
	      prod(*leftInteractionBlock, *rightInteractionBlock, currentInteractionBlock, false);
	    }
	  else
	    {
	      // centralInteractionBlock contains a lu-factorized matrix and we solve
	      // centralInteractionBlock * X = rightInteractionBlock with PLU
	      SiconosMatrix& centralInteractionBlock = *getOSIMatrix(osi, ds);
	      centralInteractionBlock.PLUForwardBackwardInPlace(*rightInteractionBlock);
	      inter->computeKhat(*rightInteractionBlock, workMInter, h); // if K is non 0

	      //      integration of r with theta method removed
	      //      *currentInteractionBlock += h *Theta[*itDS]* *leftInteractionBlock * (*rightInteractionBlock); //left = C, right = W.B
	      //gemm(h,*leftInteractionBlock,*rightInteractionBlock,1.0,*currentInteractionBlock);
	      *leftInteractionBlock *= h;
	      prod(*leftInteractionBlock, *rightInteractionBlock, currentInteractionBlock, false);
	      //left = C, right = inv(W).B
	    }

	}
      else if (relationType == Lagrangian ||
	       relationType == NewtonEuler)
	{

	  SP::BoundaryCondition bc;
	  Type::Siconos dsType = Type::value(*ds);
	  if (dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS
	      || dsType == Type::LagrangianLinearDiagonalDS)
	    {
	      LagrangianDS& d = static_cast<LagrangianDS&> (*ds);
	      if (d.boundaryConditions()) bc = d.boundaryConditions();
	    }
	  else if (dsType == Type::NewtonEulerDS)
	    {
	      NewtonEulerDS& d = static_cast<NewtonEulerDS&> (*ds);
	      if (d.boundaryConditions()) bc = d.boundaryConditions();
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
	  if(osiType == OSI::MOREAUJEANBILBAOOSI || dsType == Type::LagrangianLinearDiagonalDS)
	    {
	      // TEMP buffer, once again ...
	      SimpleMatrix work(nslawSize, sizeDS);

	      // Get inverse of the iteration matrix
	      // Remind that for diagonal ds, the inverse of the iteration matrix has been saved.
	      SiconosMatrix& inv_iteration_matrix = *getOSIMatrix(osi, ds);
	      // Compute work = leftInteractionBlock * inv_iteration_matrix = Ht.W-1
	      //axpy_prod(*leftInteractionBlock, inv_iteration_matrix, work, true);
	      prod(*leftInteractionBlock, inv_iteration_matrix, work, true);
	      // transpose in place --> leftInteractionBlock = leftInteractionBlock.transpose = rightInteractionBlock
	      leftInteractionBlock->trans();
	      // currentInteractionBlock += work * leftInteractionBlock i.e. current = (leftInteractionBlock * inv_W * rightInteractionBlock)
	      prod(work,* leftInteractionBlock, currentInteractionBlock, true);
	      if (relationSubType == CompliantLinearTIR)
		{
		  if (osiType == OSI::MOREAUJEANOSI)
		    {
		      // note FP : move this to extraInteractionBlock??
		      currentInteractionBlock *= (static_cast<MoreauJeanOSI&> (osi)).theta() ;
		      currentInteractionBlock +=  *std11::static_pointer_cast<LagrangianCompliantLinearTIR>(inter->relation())->D()/simulation()->timeStep() ;
		    }
		}
	    }
	  else
	    {
	      // TEMP buffer, once again ...
	      SP::SiconosMatrix work(new SimpleMatrix(*leftInteractionBlock));
	      work->trans();
	      SiconosMatrix& centralInteractionBlock = *getOSIMatrix(osi, ds);
	      // solve iteration_matrix * X = work for X, result in work. 
	      centralInteractionBlock.PLUForwardBackwardInPlace(*work);
	      // currentInteractionBlock +=  leftInteractionBlock * work i.e current += H.W-1.Ht 
	      prod(*leftInteractionBlock, *work, currentInteractionBlock, false);
	      //  gemm(CblasNoTrans,CblasNoTrans,1.0,*leftInteractionBlock,*work,1.0,*currentInteractionBlock);
	      //*currentInteractionBlock *=h;
	      if (relationSubType == CompliantLinearTIR)
		{
		  if (osiType == OSI::MOREAUJEANOSI)
		    {
		      // note FP : move this to extraInteractionBlock??
		      currentInteractionBlock *= (static_cast<MoreauJeanOSI&> (osi)).theta() ;
		      currentInteractionBlock +=  *std11::static_pointer_cast<LagrangianCompliantLinearTIR>(inter->relation())->D()/simulation()->timeStep() ;
		    }
		}
	    }
	}
	  else RuntimeException::selfThrow("LinearOSNS::computeDiagonalInteractionBlock not yet implemented for relation of type " + relationType);
	  // Set pos for next dynamical system.
	  pos = pos2;
    }
  DEBUG_END("LinearOSNS::computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd) ends \n");
}

  void LinearOSNS::computeInteractionBlock(const InteractionsGraph::EDescriptor& ed)
  {

    DEBUG_PRINT("LinearOSNS::computeInteractionBlock(const InteractionsGraph::EDescriptor& ed)\n");
    // Assume:
    // - edge properties (lower/upper blocks) are properly allocated

  
    // Index set of active interactions
    InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());
    // All declared interactions, where properties are saved.
    InteractionsGraph& parentSet = *simulation()->indexSet(0);

    // Current dynamical system (from edge)
    // Required to find/get :
    //    - ds type (block computation)
    //    - osi type and integrator matrix (osi.W) for ds
    SP::DynamicalSystem ds = indexSet.bundle(ed);
    DynamicalSystemsGraph& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
    OneStepIntegrator& osi = *DSG0.properties(DSG0.descriptor(ds)).osi;
    OSI::TYPES osiType = osi.getType();

    // Source and target of current edge, i.e. connected interactions.
    // Required to find/get :
    //  - relations types and nslaw sizes (i.e. block sizes)
    //  - left/right blocks for current ds
    InteractionsGraph::VDescriptor vertex_source, vertex_target;
    vertex_source = indexSet.source(ed);
    vertex_target = indexSet.target(ed);
    Interaction& inter_source = *indexSet.bundle(vertex_source);
    Interaction& inter_target = *indexSet.bundle(vertex_target);

    // -- Check ds 'position' in connected interactions --
    // For the edge 'ds', we need to find relative position of this ds
    // in inter_source and inter_target relation matrices.
    unsigned int pos_source, pos_target;
    // Get 'source' ds of inter_source
    SP::DynamicalSystem tmpds = indexSet.properties(vertex_source).source;
    if(tmpds == ds) // is ds source of inter_source?
      pos_source = indexSet.properties(vertex_source).source_pos;
    else
      pos_source = indexSet.properties(vertex_source).target_pos;

    // Repeat process for inter_target
    tmpds = indexSet.properties(vertex_target).source;
    if(tmpds == ds) // is ds source of inter_target?
      pos_target = indexSet.properties(vertex_target).source_pos;
    else
      pos_target = indexSet.properties(vertex_target).target_pos;

    // Get indices of source and target
    unsigned int index_source = indexSet.index(vertex_source);
    unsigned int index_target = indexSet.index(vertex_target);
    unsigned int nslawSize1 = inter_source.nonSmoothLaw()->size();
    unsigned int nslawSize2 = inter_target.nonSmoothLaw()->size();
    bool compute_upper = (index_target > index_source);
    SP::SiconosMatrix currentInteractionBlock;

    InteractionsGraph::EDescriptor ed_parent = indexSet.parent_edge(ed);

    if (compute_upper)
      currentInteractionBlock = parentSet.properties(ed_parent).upper_block;
    else  // lower block
      currentInteractionBlock = indexSet.properties(ed_parent).lower_block;

    // std::cout << "COMPUTER LOWER, DEAL WITH EDGE " << indexSet.index(ed) << ", ds number : " << indexSet.bundle(*ei)->number()<< std::endl;
    // std::cout << "BEtween (index/inter number)" << isource << "/" << source->number() << " and " << itarget << "/" << target->number() << std:: 
    // source of inter_source :

    VectorOfSMatrices& workMInter1 = *indexSet.properties(vertex_source).workMatrices;
    VectorOfSMatrices& workMInter2 = *indexSet.properties(vertex_target).workMatrices;
    SP::SiconosMatrix leftInteractionBlock;

    RELATION::TYPES relationType1, relationType2;
    double h = simulation()->currentTimeStep();

    // General form of the interactionBlock is : interactionBlock =
    // a*extraInteractionBlock + b * leftInteractionBlock * centralInteractionBlocks
    // * rightInteractionBlock a and b are scalars, centralInteractionBlocks a
    // matrix depending on the integrator (and on the DS), the
    // simulation type ...  left, right and extra depend on the relation
    // type and the non smooth law.
    relationType1 = inter_source.relation()->getType();
    relationType2 = inter_target.relation()->getType();

    // ==== First Order Relations - Specific treatment for diagonal
    // _interactionBlocks ===
    // assert(inter_source != inter_target);
    //  currentInteractionBlock->zero();


    // loop over the common DS
    unsigned int sizeDS = ds->dimension();

    // get _interactionBlocks corresponding to the current DS
    // These _interactionBlocks depends on the relation type.
    leftInteractionBlock.reset(new SimpleMatrix(nslawSize1, sizeDS));
    inter_source.getLeftInteractionBlockForDS(pos_source, leftInteractionBlock, workMInter1);

    // Computing depends on relation type -> move this in Interaction method?
    if (relationType1 == FirstOrder && relationType2 == FirstOrder)
      {
    
	SP::SiconosMatrix rightInteractionBlock(new SimpleMatrix(sizeDS, nslawSize2));
	inter_target.getRightInteractionBlockForDS(pos_target, rightInteractionBlock, workMInter2);
	// centralInteractionBlock contains a lu-factorized matrix and we solve
	// centralInteractionBlock * X = rightInteractionBlock with PLU
	SiconosMatrix& centralInteractionBlock = *getOSIMatrix(osi, ds);
	centralInteractionBlock.PLUForwardBackwardInPlace(*rightInteractionBlock);

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
	SP::BoundaryCondition bc;
	Type::Siconos dsType = Type::value(*ds);
	if (dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearDiagonalDS)
	  {
	    LagrangianDS& d = static_cast<LagrangianDS&> (*ds);
	    if (d.boundaryConditions()) bc = d.boundaryConditions();
	  }
	else if (dsType == Type::NewtonEulerDS)
	  {
	    NewtonEulerDS& d = static_cast<NewtonEulerDS&> (*ds);
	    if (d.boundaryConditions()) bc = d.boundaryConditions();
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

	if(osiType == OSI::MOREAUJEANBILBAOOSI || dsType == Type::LagrangianLinearDiagonalDS)
	  {
	    // Rightinteractionblock used first as buffer to save left * W-1
	    SP::SiconosMatrix rightInteractionBlock(new SimpleMatrix(nslawSize2, sizeDS));
	    // Get inverse of the iteration matrix
	    SiconosMatrix& inv_iteration_matrix = *getOSIMatrix(osi, ds);
	    // remind that W contains the inverse of the iteration matrix
	    axpy_prod(*leftInteractionBlock, inv_iteration_matrix, *rightInteractionBlock, true);
	    // Then save block corresponding to the 'right' interaction into leftInteractionBlock
	    inter_target.getLeftInteractionBlockForDS(pos_target, leftInteractionBlock, workMInter2);
	    leftInteractionBlock->trans();
	    // and compute LW-1R == rightInteractionBlock * leftInteractionBlock into currentInteractionBlock
	    prod(*rightInteractionBlock, *leftInteractionBlock, *currentInteractionBlock, false);
	  }
	else
	  {
	    // inter_source != inter_target
	    SP::SiconosMatrix rightInteractionBlock(new SimpleMatrix(nslawSize2, sizeDS));
	    inter_target.getLeftInteractionBlockForDS(pos_target, rightInteractionBlock, workMInter2);
	    rightInteractionBlock->trans();
	    // Warning: we use getLeft for Right interactionBlock
	    // because right = transpose(left) and because of
	    // size checking inside the getBlock function, a
	    // getRight call will fail.
	    SimpleMatrix& centralInteractionBlock = *getOSIMatrix(osi, ds);
	    centralInteractionBlock.PLUForwardBackwardInPlace(*rightInteractionBlock);
	    //*currentInteractionBlock +=  *leftInteractionBlock ** work;
	    prod(*leftInteractionBlock, *rightInteractionBlock, *currentInteractionBlock, false);
	  }
      }
    else RuntimeException::selfThrow("LinearOSNS::computeInteractionBlock not yet implemented for relation of type " + relationType1);

  }

void LinearOSNS::computeqBlock(InteractionsGraph::VDescriptor& vertex_inter, unsigned int pos)
{
  DEBUG_PRINT("LinearOSNS::computeqBlock(SP::Interaction inter, unsigned int pos)\n");
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());

  // At most 2 DS are linked by an Interaction
  SP::DynamicalSystem ds1 = indexSet->properties(vertex_inter).source;
  SP::DynamicalSystem ds2 = indexSet->properties(vertex_inter).target;
  assert(ds1);
  assert(ds2);

  DynamicalSystemsGraph& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

  OneStepIntegrator& osi1 = *DSG0.properties(DSG0.descriptor(ds1)).osi;
  OneStepIntegrator& osi2 = *DSG0.properties(DSG0.descriptor(ds2)).osi;

  OSI::TYPES osi1Type = osi1.getType();
  OSI::TYPES osi2Type = osi2.getType();

  SP::Interaction inter = indexSet->bundle(vertex_inter);
  unsigned int sizeY = inter->nonSmoothLaw()->size();

  // We assume that the osi of ds1 (osi1) is integrating the interaction
  DEBUG_EXPR(display());
  if ((osi1Type == OSI::EULERMOREAUOSI && osi2Type == OSI::EULERMOREAUOSI) ||
      (osi1Type == OSI::ZOHOSI && osi2Type == OSI::ZOHOSI))
    {
      osi1.computeFreeOutput(vertex_inter, this);
      SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[FirstOrderR::osnsp_rhs];
      setBlock(osnsp_rhs, _q, sizeY , 0, pos);
    }
  else if ((osi1Type == OSI::MOREAUJEANOSI  && osi2Type == OSI::MOREAUJEANOSI  )||
           (osi1Type == OSI::MOREAUDIRECTPROJECTIONOSI && osi2Type == OSI::MOREAUDIRECTPROJECTIONOSI))
    {
      osi1.computeFreeOutput(vertex_inter, this);
      SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[MoreauJeanOSI::OSNSP_RHS];
      setBlock(osnsp_rhs, _q, sizeY , 0, pos);
    }
  else if ((osi1Type == OSI::MOREAUJEANBILBAOOSI && osi2Type == OSI::MOREAUJEANBILBAOOSI ))
    {
      osi1.computeFreeOutput(vertex_inter, this);
      SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[MoreauJeanBilbaoOSI::OSNSP_RHS];
      setBlock(osnsp_rhs, _q, sizeY , 0, pos);
    }
  else if ((osi1Type == OSI::LSODAROSI && osi2Type == OSI::LSODAROSI  ) )
    {
      osi1.computeFreeOutput(vertex_inter, this);
      SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[LsodarOSI::OSNSP_RHS];
      setBlock(osnsp_rhs, _q, sizeY , 0, pos);
    }
  else if ((osi1Type == OSI::NEWMARKALPHAOSI && osi2Type == OSI::NEWMARKALPHAOSI  ))
    {
      osi1.computeFreeOutput(vertex_inter, this);
      SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[NewMarkAlphaOSI::OSNSP_RHS];
      setBlock(osnsp_rhs, _q, sizeY , 0, pos);
    }
  else if ((osi1Type == OSI::SCHATZMANPAOLIOSI && osi2Type == OSI::SCHATZMANPAOLIOSI ) )
    {
      osi1.computeFreeOutput(vertex_inter, this);
      SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[SchatzmanPaoliOSI::OSNSP_RHS];
      setBlock(osnsp_rhs, _q, sizeY , 0, pos);
    }
  else if ((osi1Type == OSI::D1MINUSLINEAROSI && osi2Type == OSI::D1MINUSLINEAROSI  ))
    {
      osi1.computeFreeOutput(vertex_inter, this);
      SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[D1MinusLinearOSI::OSNSP_RHS];
      setBlock(osnsp_rhs, _q, sizeY , 0, pos);
    }

  else if (osi1Type == OSI::MOREAUJEANGOSI && osi2Type == OSI::MOREAUJEANGOSI)
    {

    }
  else
    RuntimeException::selfThrow("LinearOSNS::computeqBlock not yet implemented for OSI1 and OSI2 of type " + osi1Type  + osi2Type);
  DEBUG_EXPR(_q->display());
}

void LinearOSNS::computeq(double time)
{
  if (_q->size() != _sizeOutput)
    _q->resize(_sizeOutput);
  _q->zero();

  // === Get index set from Simulation ===
  InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());
  // === Loop through "active" Interactions (ie present in
  // indexSets[level]) ===

  unsigned int pos = 0;
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
    {
      // Compute q, this depends on the type of non smooth problem, on
      // the relation type and on the non smooth law
      pos = indexSet.properties(*ui).absolute_position;
      computeqBlock(*ui, pos); // free output is saved in y
    }
}


void LinearOSNS::updateOperators()
{
  InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());
  InteractionsGraph& parentSet = *simulation()->indexSet(0);

  // Update OSNSpb matrix, according to indexSet current state
  _M->fillW(indexSet, parentSet, !_hasBeenUpdated);
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

  // Reset _w and _z with previous values of y and lambda
  // (i.e. val saved in yOutputOld and lambdaOld of the interaction).
  // Note : sizeOuput can be unchanged, but positions may have changed. (??)
  if (_keepLambdaAndYState)
    {
      InteractionsGraph::VIterator ui, uiend;
      for (std11::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
	{
	  Interaction& inter = *indexSet.bundle(*ui);
	  // Get the position of inter-interactionBlock in the vector w
	  // or z
	  unsigned int pos = indexSet.properties(*ui).absolute_position;
	  SiconosVector& yOutputOld = *inter.yOld(inputOutputLevel());
	  SiconosVector& lambdaOld = *inter.lambdaOld(inputOutputLevel());

	  if (_sizeOutput >= yOutputOld.size() + pos)
	    {
	      setBlock(yOutputOld, _w, yOutputOld.size(), 0, pos);
	      setBlock(lambdaOld, _z, lambdaOld.size(), 0, pos);
	    }
	}
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
  SP::Topology topology = simulation()->nonSmoothDynamicalSystem()->topology();
  bool isLinear = simulation()->nonSmoothDynamicalSystem()->isLinear();

  //   std::cout << "!b || !isLinear :"  << boolalpha <<  (!b || !isLinear) <<  std::endl;

  // nothing to do
  if (indexSetLevel() == LEVELMAX)
    return false;

  InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());

  if (indexSet.size() == 0)
    return false;

  // 
  if(!_hasBeenUpdated) 
    {
      // Should be called even for linear nsds.
      // It's up to updateInteractionBlocks to switch between linear/non linear
      // cases.
      updateInteractionBlocks();
      updateOperators();
    }
  // else
  // nothing to do (IsLinear and not changed)

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
  InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());

  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;

  // === Loop through "active" Interactions (ie present in
  // indexSets[1]) ===

  unsigned int pos = 0;

  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
    {
      Interaction& inter = *indexSet.bundle(*ui);
      // Get the  position of inter-interactionBlock in the vector w
      // or z
      pos = indexSet.properties(*ui).absolute_position;

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
