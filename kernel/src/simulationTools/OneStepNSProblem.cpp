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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "OneStepNSProblem.hpp"
#include "NonSmoothDynamicalSystem.hpp"
//#include "Interaction.hpp"
#include "Interaction.hpp"
#include "Topology.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "EulerMoreauOSI.hpp"
#include "MoreauJeanOSI.hpp"
#include "MoreauJeanBilbaoOSI.hpp"
#include "SchatzmanPaoliOSI.hpp"
#include "NewMarkAlphaOSI.hpp"
#include "LagrangianDS.hpp"
#include "NewtonEulerDS.hpp"
#include "ZeroOrderHoldOSI.hpp"
#include "NonSmoothLaw.hpp"
#include "Simulation.hpp"

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"
#include "numerics_verbose.h" // numerics to set verbose mode ...


OneStepNSProblem::OneStepNSProblem():
  _indexSetLevel(0), _inputOutputLevel(0), _maxSize(0), _hasBeenUpdated(false)
{
  _numerics_solver_options.reset(new SolverOptions);
  _numerics_solver_options->iWork = NULL;   _numerics_solver_options->callback = NULL;
  _numerics_solver_options->dWork = NULL;
}
// --- CONSTRUCTORS/DESTRUCTOR ---


// Constructor with given simulation and a pointer on Solver (Warning, solver is an optional argument)
OneStepNSProblem::OneStepNSProblem(int numericsSolverId):
  _numerics_solver_id(numericsSolverId), _sizeOutput(0),
  _indexSetLevel(0), _inputOutputLevel(0), _maxSize(0), _hasBeenUpdated(false)
{

  _numerics_solver_options.reset(new SolverOptions);
  _numerics_solver_options->iWork = NULL;   _numerics_solver_options->callback = NULL;
  _numerics_solver_options->dWork = NULL;
  _numerics_solver_options->solverId = numericsSolverId;
}

bool OneStepNSProblem::hasInteractions() const
{
  return _simulation->nonSmoothDynamicalSystem()->topology()->indexSet(_indexSetLevel)->size() > 0 ;
}


void OneStepNSProblem::updateInteractionBlocks()
{
  DEBUG_PRINT("OneStepNSProblem::updateInteractionBlocks() starts\n");
  // -- Get index sets from simulation --
  // Active interactions
  InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());
  // All declared interactions, where properties are saved.
  InteractionsGraph& parentSet = *simulation()->indexSet(0);

  bool isLinear = simulation()->nonSmoothDynamicalSystem()->isLinear();

  // -- Compute diagonal blocks (vertices property) --
  //  Linear case:
  //   - Loop over all vertices
  //       if block exists (isallocated) : nothing
  //       else allocate and compute
  //
  //  Nonlinear case:
  //   - Loop over all vertices
  //       if block exists (isallocated) : compute
  //       else allocate and compute

  InteractionsGraph::VIterator vi, viend;
  SP::Interaction inter;
  InteractionsGraph::VDescriptor vd_parent;
  unsigned int nslawSize;

  // Loop over all active interactions (I(level) vertices)
  for (std11::tie(vi, viend) = indexSet.vertices(); vi != viend; ++vi)
    {
      inter = indexSet.bundle(*vi);
      //std::cout << "START INTER LOOP, inter number : " << inter->number() <<std::endl;
      nslawSize = inter->nonSmoothLaw()->size();
      // Get corresponding descriptor in parent graph
      vd_parent = parentSet.descriptor(inter);
      // Allocate and compute only if inter->block is null in parent set
      if(isLinear)
	{
	  if(! parentSet.properties(vd_parent).block)
	    {
	      // Allocate ...
	      // std::cout << "Allocate diagonal block " << parentSet.index(vd_parent) << std::endl;
	      parentSet.properties(vd_parent).block.reset(new SimpleMatrix(nslawSize, nslawSize));
	      // std::cout << "compute diagonal block " << std::endl;
	      // Compute ...
	      computeDiagonalInteractionBlock(*vi);
	    }
	}
      else
	{
	  if(! parentSet.properties(vd_parent).block)
	    {
	      // Allocate (at first touch)
	      // std::cout << "Allocate diagonal block " << std::endl;
	      parentSet.properties(vd_parent).block.reset(new SimpleMatrix(nslawSize, nslawSize));
	    }
	  // Reset and compute : every time
	  // std::cout << "compute diagonal block " << std::endl;
	  parentSet.properties(vd_parent).block->zero(); // check if this is really required.
	  computeDiagonalInteractionBlock(*vi);
	}
    }
    
  // -- Compute extra-diagonal blocks (edges property) --
  //
  InteractionsGraph::VDescriptor vd_source, vd_target;
  InteractionsGraph::EDescriptor ed_parent;
  SP::Interaction source, target;
  InteractionsGraph::EIterator ei, eiend;
  unsigned int isource, itarget;
  unsigned int source_size, target_size;
  bool compute_upper;
  // Loop over all edges (dynamical systems!) in active interactions set
  for (std11::tie(ei, eiend) = indexSet.edges(); ei != eiend; ++ei)
    {
      // Get edge in parent set corresponding to current edge in active set
      ed_parent = indexSet.parent_edge(*ei);
      //std::cout << "START EDGES LOOP " << std::endl;
      // For each edge (ds) get connected interactions (2 exactly, there can't be self loops in IG)
      vd_source = indexSet.source(*ei);
      vd_target = indexSet.target(*ei);
      source = indexSet.bundle(vd_source);
      target = indexSet.bundle(vd_target);
      // their sizes, 
      source_size = source->nonSmoothLaw()->size();
      target_size = source->nonSmoothLaw()->size();
      // and their indices.
      isource = indexSet.index(vd_source);
      itarget = indexSet.index(vd_target);
      // Set case (position in OSNSMatrix, upper or lower, depending on indices).
      compute_upper = (itarget > isource);
      //std::cout << "DEAL WITH EDGE " << indexSet.index(*ei) << ", ds number : " << indexSet.bundle(*ei)->number()<< std::endl;
      //std::cout << "BEtween (index/inter number)" << isource << "/" << source->number() << " and " << itarget << "/" << target->number() << std::endl;
      if(isLinear)
	{
	  if (compute_upper) // upper block
	    {
	      if (! parentSet.properties(ed_parent).upper_block)
		{
		  // Allocate and compute at first touch
		  // std::cout << "UPPER, DEAL WITH EDGE " << indexSet.index(*ei) << ", ds number : " << indexSet.bundle(*ei)->number()<< std::endl;
		  // std::cout << "BEtween (index/inter number)" << isource << "/" << source->number() << " and " << itarget << "/" << target->number() << std::endl;
		  //std::cout << "ED lin upper block allocation " << isource << " " << itarget << " " << indexSet.bundle(*ei)->number() << std::endl;
		  parentSet.properties(ed_parent).upper_block.reset(new SimpleMatrix(source_size, target_size));
		  computeInteractionBlock(*ei);
		  // TEMP to deal with non-sym case
		  // Allocate lower block and transpose in place.
		  parentSet.properties(ed_parent).lower_block.reset(new SimpleMatrix(*parentSet.properties(ed_parent).upper_block));
		  parentSet.properties(ed_parent).lower_block->trans();	     
		}
	    }
	  else // lower block
	    {
	      if (! parentSet.properties(ed_parent).lower_block)
		{
		  // Allocate and compute at first touch
		  // std::cout << "LOWER, DEAL WITH EDGE " << indexSet.index(*ei) << ", ds number : " << indexSet.bundle(*ei)->number()<< std::endl;
		  // std::cout << "BEtween (index/inter number)" << isource << "/" << source->number() << " and " << itarget << "/" << target->number() << std::endl;
		  parentSet.properties(ed_parent).lower_block.reset(new SimpleMatrix(source_size, target_size));
		  computeInteractionBlock(*ei);
		  // TEMP to deal with non-sym case
		  parentSet.properties(ed_parent).upper_block.reset(new SimpleMatrix(*parentSet.properties(ed_parent).lower_block));
		  parentSet.properties(ed_parent).upper_block->trans();	     
		}
	    }
	}
      else
	{
	  if (compute_upper) // upper block
	    {
	      // Allocate at first touch
	      if (! parentSet.properties(ed_parent).upper_block)
		{
		  // std::cout << "ED nonlin upper block allocation " <<  isource << " " << itarget<< std::endl;	       
		  parentSet.properties(ed_parent).upper_block.reset(new SimpleMatrix(source_size, target_size));
		}
	      assert(parentSet.properties(ed_parent).upper_block->size(0) == source_size);
	      assert(parentSet.properties(ed_parent).upper_block->size(1) == target_size);
	    }
	  else // lower block
	    {
	      // Allocate at first touch
	      if (! parentSet.properties(ed_parent).lower_block)
		{
		  // std::cout << "ED nonlin lower block allocation " <<  isource << " " << itarget<< std::endl;	       
		  parentSet.properties(ed_parent).lower_block.reset(new SimpleMatrix(source_size, target_size));
		}
	      assert(parentSet.properties(ed_parent).lower_block->size(0) == source_size);
	      assert(parentSet.properties(ed_parent).lower_block->size(1) == target_size);
	    }
	  // Compute every time
	  // std::cout << "ED nonlin  block comp "<< std::endl;	       
	  computeInteractionBlock(*ei);
	}
    }
}

void OneStepNSProblem::displayBlocks(InteractionsGraph& indexSet)
{

  std::cout <<  "OneStepNSProblem::displayBlocks(InteractionsGraph& indexSet) " << std::endl;
  InteractionsGraph::VIterator vi, viend;
  for (std11::tie(vi, viend) = indexSet.vertices();
       vi != viend; ++vi)
  {
    std::cout << " vertex ..." << indexSet.index(*vi) << std::endl;
    if (indexSet.properties(*vi).block)
    {
      indexSet.properties(*vi).block->display();
    }

    InteractionsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = indexSet.out_edges(*vi);
         oei != oeiend; ++oei)
    {
      InteractionsGraph::EDescriptor ed1, ed2;
      std11::tie(ed1, ed2) = indexSet.edges(indexSet.source(*oei), indexSet.target(*oei));
      std::cout << " edge ..." << indexSet.index(ed1) << std::endl;

      if (indexSet.properties(ed1).upper_block)
      {
        indexSet.properties(ed1).upper_block->display();
      }
      if (indexSet.properties(ed1).lower_block)
      {
        indexSet.properties(ed1).lower_block->display();
      }
      if (indexSet.properties(ed2).upper_block)
      {
        indexSet.properties(ed2).upper_block->display();
      }
      if (indexSet.properties(ed2).lower_block)
      {
        indexSet.properties(ed2).lower_block->display();
      }
    }

  }
}

void OneStepNSProblem::initialize(SP::Simulation sim)
{
  // Link with the simulation that owns this osnsp

  assert(sim && "OneStepNSProblem::initialize(sim), sim is null.");

  _simulation = sim;

  // The maximum size of the problem (for example, the dim. of M in
  // LCP or Friction problems).  Set to the number of possible scalar
  // constraints declared in the topology.
  if (_maxSize == 0) // if maxSize not set explicitely by user before
    _maxSize = simulation()->nonSmoothDynamicalSystem()->topology()->numberOfConstraints();
}

SP::SimpleMatrix OneStepNSProblem::getOSIMatrix(OneStepIntegrator& Osi, SP::DynamicalSystem ds)
{
  // Returns the integration matrix from one-step integrator and dynamical system.

  // Matrix depends on OSI type.
  SP::SimpleMatrix block;
  OSI::TYPES osiType; // type of the current one step integrator
  Type::Siconos dsType; // type of the current Dynamical System

  osiType = Osi.getType();
  dsType = Type::value(*ds);

  if (osiType == OSI::MOREAUJEANOSI
      || osiType == OSI::MOREAUDIRECTPROJECTIONOSI)
  {
      block = (static_cast<MoreauJeanOSI&> (Osi)).W(ds); // get its W matrix ( pointer link!)
  }
  else if (osiType == OSI::MOREAUJEANBILBAOOSI)
  {
    block = (static_cast<MoreauJeanBilbaoOSI&> (Osi)).iteration_matrix(ds); // get its W matrix ( pointer link!)
  }
  else if (osiType == OSI::SCHATZMANPAOLIOSI)
  {
      block = (static_cast<SchatzmanPaoliOSI&> (Osi)).W(ds); // get its W matrix ( pointer link!)
  }
  else if (osiType == OSI::EULERMOREAUOSI)
  {
    block = (static_cast<EulerMoreauOSI&>(Osi)).W(ds); // get its W matrix ( pointer link!)
  }
  else if (osiType == OSI::LSODAROSI) // Warning: LagrangianDS only at the time !!!
  {
    if (dsType != Type::LagrangianDS && dsType != Type::LagrangianLinearTIDS)
      RuntimeException::selfThrow("OneStepNSProblem::getOSIMatrix not yet implemented for LsodarOSI Integrator with dynamical system of type " + dsType);

    // get lu-factorized mass
    block = (std11::static_pointer_cast<LagrangianDS>(ds))->inverseMass();
  }
  else if (osiType == OSI::NEWMARKALPHAOSI)
  {
    if (dsType != Type::LagrangianDS && dsType != Type::LagrangianLinearTIDS)
    {
      RuntimeException::selfThrow("OneStepNSProblem::getOSIMatrix not yet implemented for NewmarkAlphaOSI Integrator with dynamical system of type " + dsType);
    }
    //
    SP::OneStepNSProblems  allOSNS  = Osi.simulation()->oneStepNSProblems();
    // If LCP at acceleration level
    if (((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == this)
    {
      block = (std11::static_pointer_cast<LagrangianDS>(ds))->inverseMass();
    }
    else // It LCP at position level
    {
      block = (static_cast<NewMarkAlphaOSI&>(Osi)).W(ds);
    }
  } // End Newmark OSI
  else if (osiType == OSI::D1MINUSLINEAROSI)
  {
    DEBUG_PRINT("OneStepNSProblem::getOSIMatrix  for osiType   OSI::D1MINUSLINEAR\n");
    /** \warning V.A. 30/052013 for implicit D1Minus it will not be the mass matrix for all OSNSP*/
    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      // SP::SimpleMatrix Mold;
      // Mold.reset(new SimpleMatrix(*(std11::static_pointer_cast<LagrangianDS>(ds))->mass()));
      // DEBUG_EXPR(Mold->display(););
      // DEBUG_EXPR_WE(std::cout <<  std::boolalpha << " Mold->isPLUFactorized() = "<< Mold->isPLUFactorized() << std::endl;);
      //(std11::static_pointer_cast<LagrangianDS>(ds))->computeMass();
      SP::SiconosMatrix Mass = ((std11::static_pointer_cast<LagrangianDS>(ds))->mass()) ;
      DEBUG_EXPR(Mass->display(););
      DEBUG_EXPR_WE(std::cout <<  std::boolalpha << " Mass->isPLUFactorized() = "<< Mass->isPLUFactorized() << std::endl;);

      //DEBUG_EXPR(std::cout << (*Mass-*Mold).normInf() << std::endl;);
      /*Copy of the current mass matrix. */
      block.reset(new SimpleMatrix(*Mass));
    }
    else if (dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
      //   d->computeMass();
      //   d->mass()->resetLU();
      DEBUG_EXPR(d->mass()->display(););
      block.reset(new SimpleMatrix(*(d->mass())));
    }
    else
      RuntimeException::selfThrow("OneStepNSProblem::getOSIMatrix not yet implemented for D1MinusLinearOSI integrator with dynamical system of type " + dsType);
  }
  // for ZeroOrderHoldOSI, the central block is Ad = \int exp{As} ds over t_k, t_{k+1}
  else if (osiType == OSI::ZOHOSI)
  {
    if (!block)
      block.reset(new SimpleMatrix((static_cast<ZeroOrderHoldOSI&>(Osi)).Ad(ds)));
    else
      *block = (static_cast<ZeroOrderHoldOSI&>(Osi)).Ad(ds);
  }
  else
    RuntimeException::selfThrow("OneStepNSProblem::getOSIMatrix not yet implemented for Integrator of type " + osiType);
  return block;
}

void OneStepNSProblem::setSolverId(int solverId)
{
  RuntimeException::selfThrow("OneStepNSProblem::setSolverId - this virtual method should be implemented in all derived classes!");
}

void OneStepNSProblem::setNumericsVerboseMode(bool vMode)
{
  numerics_set_verbose(vMode);
}
