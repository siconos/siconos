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
  // The present functions checks various conditions and possibly
  // compute interactionBlocks matrices.
  //
  // Let interi and interj be two Interactions.
  //
  // Things to be checked are:
  //  1 - is the topology time invariant?
  //  2 - does interactionBlocks[interi][interj] already exists (ie has been
  //  computed in a previous time step)?
  //  3 - do we need to compute this interactionBlock? An interactionBlock has
  //  to be computed if interi and interj are in IndexSet1 AND if interi and
  //  interj have common DynamicalSystems.
  //

  //
  //  - If 1 and 2 are true then it does nothing. 3 is not checked.
  //  - If 1 == true, 2 == false, 3 == false, it does nothing.
  //  - If 1 == true, 2 == false, 3 == true, it computes the interactionBlock.
  //  - If 1==false, 2 is not checked, and the interactionBlock is computed if 3==true.
  //

  // Get index set from Simulation
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());

  bool isLinear = simulation()->nonSmoothDynamicalSystem()->isLinear();

  // we put diagonal informations on vertices
  // self loops with bgl are a *nightmare* at the moment
  // (patch 65198 on standard boost install)
  if (indexSet->properties().symmetric)
  {
    DEBUG_PRINT("OneStepNSProblem::updateInteractionBlocks(). Symmetric case");
    InteractionsGraph::VIterator vi, viend;
    for (std11::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {
      SP::Interaction inter = indexSet->bundle(*vi);
      unsigned int nslawSize = inter->nonSmoothLaw()->size();
      if (! indexSet->properties(*vi).block)
      {
        indexSet->properties(*vi).block.reset(new SimpleMatrix(nslawSize, nslawSize));
      }

      if (!isLinear || !_hasBeenUpdated)
      {
        computeDiagonalInteractionBlock(*vi);
      }
    }

    /* interactionBlock must be zeroed at init */
    std::vector<bool> initialized;
    initialized.resize(indexSet->edges_number());
    std::fill(initialized.begin(), initialized.end(), false);

    InteractionsGraph::EIterator ei, eiend;
    for (std11::tie(ei, eiend) = indexSet->edges();
         ei != eiend; ++ei)
    {
      SP::Interaction inter1 = indexSet->bundle(indexSet->source(*ei));
      SP::Interaction inter2 = indexSet->bundle(indexSet->target(*ei));

      /* on adjoint graph there is at most 2 edges between source and target */
      InteractionsGraph::EDescriptor ed1, ed2;
      std11::tie(ed1, ed2) = indexSet->edges(indexSet->source(*ei), indexSet->target(*ei));

      assert(*ei == ed1 || *ei == ed2);

      /* the first edge has the lower index */
      assert(indexSet->index(ed1) <= indexSet->index(ed2));

      // Memory allocation if needed
      unsigned int nslawSize1 = inter1->nonSmoothLaw()->size();
      unsigned int nslawSize2 = inter2->nonSmoothLaw()->size();
      unsigned int isrc = indexSet->index(indexSet->source(*ei));
      unsigned int itar = indexSet->index(indexSet->target(*ei));

      SP::SiconosMatrix currentInteractionBlock;

      if (itar > isrc) // upper block
      {
        if (! indexSet->properties(ed1).upper_block)
        {
          indexSet->properties(ed1).upper_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
          if (ed2 != ed1)
            indexSet->properties(ed2).upper_block = indexSet->properties(ed1).upper_block;
        }
        currentInteractionBlock = indexSet->properties(ed1).upper_block;
      }
      else  // lower block
      {
        if (! indexSet->properties(ed1).lower_block)
        {
          indexSet->properties(ed1).lower_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
          if (ed2 != ed1)
            indexSet->properties(ed2).lower_block = indexSet->properties(ed1).lower_block;
        }
        currentInteractionBlock = indexSet->properties(ed1).lower_block;
      }

      if (!initialized[indexSet->index(ed1)])
      {
        initialized[indexSet->index(ed1)] = true;
        currentInteractionBlock->zero();
      }
      if (!isLinear || !_hasBeenUpdated)
      {
        {
          computeInteractionBlock(*ei);
        }

        // allocation for transposed block
        // should be avoided

        if (itar > isrc) // upper block has been computed
        {
          if (!indexSet->properties(ed1).lower_block)
          {
            indexSet->properties(ed1).lower_block.
            reset(new SimpleMatrix(indexSet->properties(ed1).upper_block->size(1),
                                   indexSet->properties(ed1).upper_block->size(0)));
          }
          indexSet->properties(ed1).lower_block->trans(*indexSet->properties(ed1).upper_block);
          indexSet->properties(ed2).lower_block = indexSet->properties(ed1).lower_block;
        }
        else
        {
          assert(itar < isrc);    // lower block has been computed
          if (!indexSet->properties(ed1).upper_block)
          {
            indexSet->properties(ed1).upper_block.
            reset(new SimpleMatrix(indexSet->properties(ed1).lower_block->size(1),
                                   indexSet->properties(ed1).lower_block->size(0)));
          }
          indexSet->properties(ed1).upper_block->trans(*indexSet->properties(ed1).lower_block);
          indexSet->properties(ed2).upper_block = indexSet->properties(ed1).upper_block;
        }
      }
    }
  }
  else // not symmetric => follow out_edges for each vertices
  {
    DEBUG_PRINT("OneStepNSProblem::updateInteractionBlocks(). Non symmetric case\n");

    InteractionsGraph::VIterator vi, viend;
    for (std11::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {
      DEBUG_PRINT("OneStepNSProblem::updateInteractionBlocks(). Computation of diaganal block\n");
      SP::Interaction inter = indexSet->bundle(*vi);
      unsigned int nslawSize = inter->nonSmoothLaw()->size();
      if (! indexSet->properties(*vi).block)
      {
        indexSet->properties(*vi).block.reset(new SimpleMatrix(nslawSize, nslawSize));
      }

      if (!isLinear || !_hasBeenUpdated)
      {
        computeDiagonalInteractionBlock(*vi);
      }

      /* on a undirected graph, out_edges gives all incident edges */
      InteractionsGraph::OEIterator oei, oeiend;
      /* interactionBlock must be zeroed at init */
      std::map<SP::SiconosMatrix, bool> initialized;
      for (std11::tie(oei, oeiend) = indexSet->out_edges(*vi);
           oei != oeiend; ++oei)
      {
        /* on adjoint graph there is at most 2 edges between source and target */
        InteractionsGraph::EDescriptor ed1, ed2;
        std11::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));
        if (indexSet->properties(ed1).upper_block)
        {
          initialized[indexSet->properties(ed1).upper_block] = false;
        }
        // if(indexSet->properties(ed2).upper_block)
        // {
        //   initialized[indexSet->properties(ed2).upper_block] = false;
        // }

        if (indexSet->properties(ed1).lower_block)
        {
          initialized[indexSet->properties(ed1).lower_block] = false;
        }
        // if(indexSet->properties(ed2).lower_block)
        // {
        //   initialized[indexSet->properties(ed2).lower_block] = false;
        // }

      }

      for (std11::tie(oei, oeiend) = indexSet->out_edges(*vi);
           oei != oeiend; ++oei)
      {
        DEBUG_PRINT("OneStepNSProblem::updateInteractionBlocks(). Computation of extra-diaganal block\n");

        /* on adjoint graph there is at most 2 edges between source and target */
        InteractionsGraph::EDescriptor ed1, ed2;
        std11::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));

        assert(*oei == ed1 || *oei == ed2);

        /* the first edge as the lower index */
        assert(indexSet->index(ed1) <= indexSet->index(ed2));

        SP::Interaction inter1 = indexSet->bundle(indexSet->source(*oei));
        SP::Interaction inter2 = indexSet->bundle(indexSet->target(*oei));

        // Memory allocation if needed
        unsigned int nslawSize1 = inter1->nonSmoothLaw()->size();
        unsigned int nslawSize2 = inter2->nonSmoothLaw()->size();
        unsigned int isrc = indexSet->index(indexSet->source(*oei));
        unsigned int itar = indexSet->index(indexSet->target(*oei));

        SP::SiconosMatrix currentInteractionBlock;

        if (itar > isrc) // upper block
        {
          if (! indexSet->properties(ed1).upper_block)
          {
            indexSet->properties(ed1).upper_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
            initialized[indexSet->properties(ed1).upper_block] = false;
            if (ed2 != ed1)
              indexSet->properties(ed2).upper_block = indexSet->properties(ed1).upper_block;
          }
          currentInteractionBlock = indexSet->properties(ed1).upper_block;

        }
        else  // lower block
        {
          if (! indexSet->properties(ed1).lower_block)
          {
            indexSet->properties(ed1).lower_block.reset(new SimpleMatrix(nslawSize1, nslawSize2));
            initialized[indexSet->properties(ed1).lower_block] = false;
            if (ed2 != ed1)
              indexSet->properties(ed2).lower_block = indexSet->properties(ed1).lower_block;
          }
          currentInteractionBlock = indexSet->properties(ed1).lower_block;
        }


        if (!initialized[currentInteractionBlock])
        {
          initialized[currentInteractionBlock] = true;
          currentInteractionBlock->zero();
        }

        if (!isLinear || !_hasBeenUpdated)
        {
          if (isrc != itar)
            computeInteractionBlock(*oei);
        }

      }
    }
  }


  DEBUG_EXPR(displayBlocks(indexSet););

  DEBUG_PRINT("OneStepNSProblem::updateInteractionBlocks() ends\n");


}

void OneStepNSProblem::displayBlocks(SP::InteractionsGraph indexSet)
{

  std::cout <<  "OneStepNSProblem::displayBlocks(SP::InteractionsGraph indexSet) " << std::endl;
  InteractionsGraph::VIterator vi, viend;
  for (std11::tie(vi, viend) = indexSet->vertices();
       vi != viend; ++vi)
  {
    SP::Interaction inter = indexSet->bundle(*vi);
    if (indexSet->properties(*vi).block)
    {
      indexSet->properties(*vi).block->display();
    }

    InteractionsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = indexSet->out_edges(*vi);
         oei != oeiend; ++oei)
    {
      InteractionsGraph::EDescriptor ed1, ed2;
      std11::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));

      if (indexSet->properties(ed1).upper_block)
      {
        indexSet->properties(ed1).upper_block->display();
      }
      if (indexSet->properties(ed1).lower_block)
      {
        indexSet->properties(ed1).lower_block->display();
      }
      if (indexSet->properties(ed2).upper_block)
      {
        indexSet->properties(ed2).upper_block->display();
      }
      if (indexSet->properties(ed2).lower_block)
      {
        indexSet->properties(ed2).lower_block->display();
      }
    }

  }
}

void OneStepNSProblem::initialize(SP::Simulation sim)
{
  // Link with the simulation that owns this osnsp

  assert(sim && "OneStepNSProblem::initialize(sim), sim is null.");

  _simulation = sim;

  // === Adds this in the simulation set of OneStepNSProblem === First
  // checks the id if required.  An id is required if there is more
  // than one OneStepNSProblem in the simulation

  // The maximum size of the problem (for example, the dim. of M in
  // LCP or Friction problems).  Set to the number of possible scalar
  // constraints declared in the topology.
  if (_maxSize == 0) // if maxSize not set explicitely by user before
    // initialize
    _maxSize = simulation()->nonSmoothDynamicalSystem()->topology()->numberOfConstraints();
}

SP::SimpleMatrix OneStepNSProblem::getOSIMatrix(OneStepIntegrator& Osi, SP::DynamicalSystem ds)
{
  // Connect block to the OSI matrix of a dynamical system for the current simulation.
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
