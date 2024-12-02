/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "EulerMoreauOSI.hpp"
#include "Interaction.hpp"
#include "LsodarOSI.hpp"
#include "MoreauJeanOSI.hpp"
#include "NewMarkAlphaOSI.hpp"
#include "OneStepIntegrator.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "Simulation.hpp"
#include "Topology.hpp"
#include "MoreauJeanBilbaoOSI.hpp"
#include "SchatzmanPaoliOSI.hpp"
#include "LagrangianDS.hpp"
#include "NewtonEulerDS.hpp"
#include "NonSmoothLaw.hpp"
#include "NumericsSolversNamespace.h"  // for SolverOptions tools
#include "NumericsToolsNamespace.h"    // for verbose mode
#include "Tools.hpp"                   // enum_to_string
#include "ZeroOrderHoldOSI.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

// --- CONSTRUCTORS/DESTRUCTOR ---

bool siconos::nonsmooth_formulations::OneStepNSProblem::hasInteractions() const
{
  return _simulation->nonSmoothDynamicalSystem()
             ->topology()
             ->indexSet(_indexSetLevel)
             ->size() > 0;
}

void siconos::nonsmooth_formulations::OneStepNSProblem::updateInteractionBlocks()
{
  DEBUG_PRINT("siconos::nonsmooth_formulations::OneStepNSProblem::updateInteractionBlocks() starts\n");
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
  auto indexSet = simulation()->indexSet(indexSetLevel());

  bool isLinear = simulation()->nonSmoothDynamicalSystem()->isLinear();

  // we put diagonal information on vertices
  // self loops with bgl are a *nightmare* at the moment
  // (patch 65198 on standard boost install)
  if (indexSet->properties().symmetric) {
    DEBUG_PRINT(
        "siconos::nonsmooth_formulations::OneStepNSProblem::updateInteractionBlocks(). Symmetric case");
    siconos::graphs::InteractionsGraph::VIterator vi, viend;
    for (std::tie(vi, viend) = indexSet->vertices(); vi != viend; ++vi) {
      std::shared_ptr<siconos::modeling::Interaction> inter = indexSet->bundle(*vi);
      auto nslawSize = inter->nonSmoothLaw()->size();
      if (!indexSet->properties(*vi).block) {
        indexSet->properties(*vi).block =
            std::make_shared<siconos::algebra::SimpleMatrix>(nslawSize, nslawSize);
      }

      if (!isLinear || !_hasBeenUpdated) {
        computeDiagonalInteractionBlock(*vi);
      }
    }

    /* interactionBlock must be zeroed at init */
    std::vector<bool> initialized;
    initialized.resize(indexSet->edges_number());
    std::fill(initialized.begin(), initialized.end(), false);

    siconos::graphs::InteractionsGraph::EIterator ei, eiend;
    for (std::tie(ei, eiend) = indexSet->edges(); ei != eiend; ++ei) {
      std::shared_ptr<siconos::modeling::Interaction> inter1 =
          indexSet->bundle(indexSet->source(*ei));
      std::shared_ptr<siconos::modeling::Interaction> inter2 =
          indexSet->bundle(indexSet->target(*ei));

      /* on adjoint graph there is at most 2 edges between source and target */
      siconos::graphs::InteractionsGraph::EDescriptor ed1, ed2;
      std::tie(ed1, ed2) = indexSet->edges(indexSet->source(*ei), indexSet->target(*ei));

      assert(*ei == ed1 || *ei == ed2);

      /* the first edge has the lower index */
      assert(indexSet->index(ed1) <= indexSet->index(ed2));

      // Memory allocation if needed
      auto nslawSize1 = inter1->nonSmoothLaw()->size();
      auto nslawSize2 = inter2->nonSmoothLaw()->size();
      auto isrc = indexSet->index(indexSet->source(*ei));
      auto itar = indexSet->index(indexSet->target(*ei));

      std::shared_ptr<siconos::algebra::SiconosMatrix> currentInteractionBlock;

      if (itar > isrc)  // upper block
      {
        if (!indexSet->properties(ed1).upper_block) {
          indexSet->properties(ed1).upper_block =
              std::make_shared<siconos::algebra::SimpleMatrix>(nslawSize1, nslawSize2);
          if (ed2 != ed1)
            indexSet->properties(ed2).upper_block = indexSet->properties(ed1).upper_block;
        }
        currentInteractionBlock = indexSet->properties(ed1).upper_block;
      }
      else  // lower block
      {
        if (!indexSet->properties(ed1).lower_block) {
          indexSet->properties(ed1).lower_block =
              std::make_shared<siconos::algebra::SimpleMatrix>(nslawSize1, nslawSize2);
          if (ed2 != ed1)
            indexSet->properties(ed2).lower_block = indexSet->properties(ed1).lower_block;
        }
        currentInteractionBlock = indexSet->properties(ed1).lower_block;
      }

      if (!initialized[indexSet->index(ed1)]) {
        initialized[indexSet->index(ed1)] = true;
        currentInteractionBlock->zero();
      }
      if (!isLinear || !_hasBeenUpdated) {
        {
          computeInteractionBlock(*ei);
        }

        // allocation for transposed block
        // should be avoided

        if (itar > isrc)  // upper block has been computed
        {
          if (!indexSet->properties(ed1).lower_block) {
            indexSet->properties(ed1).lower_block =
                std::make_shared<siconos::algebra::SimpleMatrix>(
                    indexSet->properties(ed1).upper_block->size(1),
                    indexSet->properties(ed1).upper_block->size(0));
          }
          indexSet->properties(ed1).lower_block->trans(*indexSet->properties(ed1).upper_block);
          indexSet->properties(ed2).lower_block = indexSet->properties(ed1).lower_block;
        }
        else {
          assert(itar < isrc);  // lower block has been computed
          if (!indexSet->properties(ed1).upper_block) {
            indexSet->properties(ed1).upper_block =
                std::make_shared<siconos::algebra::SimpleMatrix>(
                    indexSet->properties(ed1).lower_block->size(1),
                    indexSet->properties(ed1).lower_block->size(0));
          }
          indexSet->properties(ed1).upper_block->trans(*indexSet->properties(ed1).lower_block);
          indexSet->properties(ed2).upper_block = indexSet->properties(ed1).upper_block;
        }
      }
    }
  }
  else  // not symmetric => follow out_edges for each vertices
  {
    DEBUG_PRINT(
        "siconos::nonsmooth_formulations::OneStepNSProblem::updateInteractionBlocks(). Non symmetric "
        "case\n");

    siconos::graphs::InteractionsGraph::VIterator vi, viend;
    for (std::tie(vi, viend) = indexSet->vertices(); vi != viend; ++vi) {
      DEBUG_PRINT(
          "siconos::nonsmooth_formulations::OneStepNSProblem::updateInteractionBlocks(). Computation of "
          "diaganal block\n");
      std::shared_ptr<siconos::modeling::Interaction> inter = indexSet->bundle(*vi);
      auto nslawSize = inter->nonSmoothLaw()->size();
      if (!indexSet->properties(*vi).block) {
        indexSet->properties(*vi).block =
            std::make_shared<siconos::algebra::SimpleMatrix>(nslawSize, nslawSize);
      }

      if (!isLinear || !_hasBeenUpdated) {
        computeDiagonalInteractionBlock(*vi);
      }

      /* on a undirected graph, out_edges gives all incident edges */
      siconos::graphs::InteractionsGraph::OEIterator oei, oeiend;
      /* interactionBlock must be zeroed at init */
      std::map<std::shared_ptr<siconos::algebra::SiconosMatrix>, bool> initialized;
      for (std::tie(oei, oeiend) = indexSet->out_edges(*vi); oei != oeiend; ++oei) {
        /* on adjoint graph there is at most 2 edges between source and target */
        siconos::graphs::InteractionsGraph::EDescriptor ed1, ed2;
        std::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));
        if (indexSet->properties(ed1).upper_block) {
          initialized[indexSet->properties(ed1).upper_block] = false;
        }
        // if(indexSet->properties(ed2).upper_block)
        // {
        //   initialized[indexSet->properties(ed2).upper_block] = false;
        // }

        if (indexSet->properties(ed1).lower_block) {
          initialized[indexSet->properties(ed1).lower_block] = false;
        }
        // if(indexSet->properties(ed2).lower_block)
        // {
        //   initialized[indexSet->properties(ed2).lower_block] = false;
        // }
      }

      for (std::tie(oei, oeiend) = indexSet->out_edges(*vi); oei != oeiend; ++oei) {
        DEBUG_PRINT(
            "siconos::nonsmooth_formulations::OneStepNSProblem::updateInteractionBlocks(). Computation of "
            "extra-diaganal block\n");

        /* on adjoint graph there is at most 2 edges between source and target */
        siconos::graphs::InteractionsGraph::EDescriptor ed1, ed2;
        std::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));

        assert(*oei == ed1 || *oei == ed2);

        /* the first edge as the lower index */
        assert(indexSet->index(ed1) <= indexSet->index(ed2));

        auto inter1 = indexSet->bundle(indexSet->source(*oei));
        auto inter2 = indexSet->bundle(indexSet->target(*oei));

        // Memory allocation if needed
        auto nslawSize1 = inter1->nonSmoothLaw()->size();
        auto nslawSize2 = inter2->nonSmoothLaw()->size();
        auto isrc = indexSet->index(indexSet->source(*oei));
        auto itar = indexSet->index(indexSet->target(*oei));

        std::shared_ptr<siconos::algebra::SiconosMatrix> currentInteractionBlock;

        if (itar > isrc)  // upper block
        {
          if (!indexSet->properties(ed1).upper_block) {
            indexSet->properties(ed1).upper_block =
                std::make_shared<siconos::algebra::SimpleMatrix>(nslawSize1, nslawSize2);
            initialized[indexSet->properties(ed1).upper_block] = false;
            if (ed2 != ed1)
              indexSet->properties(ed2).upper_block = indexSet->properties(ed1).upper_block;
          }
          currentInteractionBlock = indexSet->properties(ed1).upper_block;
        }
        else  // lower block
        {
          if (!indexSet->properties(ed1).lower_block) {
            indexSet->properties(ed1).lower_block =
                std::make_shared<siconos::algebra::SimpleMatrix>(nslawSize1, nslawSize2);
            initialized[indexSet->properties(ed1).lower_block] = false;
            if (ed2 != ed1)
              indexSet->properties(ed2).lower_block = indexSet->properties(ed1).lower_block;
          }
          currentInteractionBlock = indexSet->properties(ed1).lower_block;
        }

        if (!initialized[currentInteractionBlock]) {
          initialized[currentInteractionBlock] = true;
          currentInteractionBlock->zero();
        }

        if (!isLinear || !_hasBeenUpdated) {
          if (isrc != itar) computeInteractionBlock(*oei);
        }
      }
    }
  }

  DEBUG_EXPR(displayBlocks(indexSet););

  DEBUG_PRINT("siconos::nonsmooth_formulations::OneStepNSProblem::updateInteractionBlocks() ends\n");
}

void siconos::nonsmooth_formulations::OneStepNSProblem::displayBlocks(
    std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet)
{
  std::cout << "siconos::nonsmooth_formulations::OneStepNSProblem::displayBlocks(std::shared_ptr<siconos::"
               "graphs::InteractionsGraph> indexSet) "
            << std::endl;
  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  for (std::tie(vi, viend) = indexSet->vertices(); vi != viend; ++vi) {
    std::shared_ptr<siconos::modeling::Interaction> inter = indexSet->bundle(*vi);
    if (indexSet->properties(*vi).block) {
      indexSet->properties(*vi).block->display();
    }

    siconos::graphs::InteractionsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = indexSet->out_edges(*vi); oei != oeiend; ++oei) {
      siconos::graphs::InteractionsGraph::EDescriptor ed1, ed2;
      std::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));

      if (indexSet->properties(ed1).upper_block) {
        indexSet->properties(ed1).upper_block->display();
      }
      if (indexSet->properties(ed1).lower_block) {
        indexSet->properties(ed1).lower_block->display();
      }
      if (indexSet->properties(ed2).upper_block) {
        indexSet->properties(ed2).upper_block->display();
      }
      if (indexSet->properties(ed2).lower_block) {
        indexSet->properties(ed2).lower_block->display();
      }
    }
  }
}

void siconos::nonsmooth_formulations::OneStepNSProblem::initialize(
    std::shared_ptr<siconos::simulation::Simulation> sim)
{
  // Link with the simulation that owns this osnsp

  assert(sim && "siconos::nonsmooth_formulations::OneStepNSProblem::initialize(sim), sim is null.");

  _simulation = sim;

  // The maximum size of the problem (for example, the dim. of M in
  // LCP or Friction problems).  Set to the number of possible scalar
  // constraints declared in the topology.
  if (_maxSize == 0)  // if maxSize not set explicitely by user before
    _maxSize = simulation()->nonSmoothDynamicalSystem()->topology()->numberOfConstraints();
}

std::shared_ptr<siconos::algebra::SimpleMatrix>
siconos::nonsmooth_formulations::OneStepNSProblem::getOSIMatrix(
    siconos::integrators::OneStepIntegrator& Osi,
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  // Returns the integration matrix from one-step integrator and dynamical system.

  // Matrix depends on OSI type.
  std::shared_ptr<siconos::algebra::SimpleMatrix> block;

  auto osiType = Osi.getType();
  // auto dsType = Type::value(*ds);

  if (osiType == siconos::integrators::IntegratorType::MOREAUJEANOSI ||
      osiType == siconos::integrators::IntegratorType::MOREAUDIRECTPROJECTIONOSI) {
    block = (static_cast<siconos::integrators::MoreauJeanOSI&>(Osi))
                .W(ds);  // get its W matrix ( pointer link!)
  }
  else if (osiType == siconos::integrators::IntegratorType::MOREAUJEANBILBAOOSI) {
    block =
        (static_cast<siconos::integrators::MoreauJeanBilbaoOSI&>(Osi)).iteration_matrix(ds);
  }
  else if (osiType == siconos::integrators::IntegratorType::SCHATZMANPAOLIOSI) {
    block = (static_cast<siconos::integrators::SchatzmanPaoliOSI&>(Osi)).W(ds);
  }
  else if (osiType == siconos::integrators::IntegratorType::EULERMOREAUOSI) {
    block = (static_cast<siconos::integrators::EulerMoreauOSI&>(Osi)).W(ds);
  }
  else if (osiType == siconos::integrators::IntegratorType::LSODAROSI)
  // Warning: LagrangianDS only at
  // the time !!!
  {
    // get lu-factorized mass
    if (auto lds = dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      block = lds->inverseMass();
    }
    else
      THROW_EXCEPTION(
          "siconos::nonsmooth_formulations::OneStepNSProblem::getOSIMatrix is implemented for LsodarOSI "
          "only with LagrangianDS systems.");
  }
  else if (osiType == siconos::integrators::IntegratorType::NEWMARKALPHAOSI) {
    if (auto lds = dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto allOSNS = Osi.simulation()->oneStepNSProblems();
      // If LCP at acceleration level
      if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == this) {
        block = lds->inverseMass();
      }
      else  // It LCP at position level
      {
        block = (static_cast<siconos::integrators::NewMarkAlphaOSI&>(Osi)).W(ds);
      }
    }
    else {
      THROW_EXCEPTION(
          "siconos::nonsmooth_formulations::OneStepNSProblem::getOSIMatrix is implemented for "
          "NewmarkAlphaOSI only with LagrangianDS systems.");
    }
  }  // End Newmark OSI
  else if (osiType == siconos::integrators::IntegratorType::D1MINUSLINEAROSI) {
    DEBUG_PRINT(
        "siconos::nonsmooth_formulations::OneStepNSProblem::getOSIMatrix  for osiType "
        "siconos::integrators::IntegratorType::D1MINUSLINEAR\n");
    /** \warning V.A. 30/052013 for implicit D1Minus it will not be the mass matrix for all
    OSNSP*/
    if (auto lds = dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto Mass = ((std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds))->mass());
      DEBUG_EXPR(Mass->display(););
      DEBUG_EXPR_WE(std::cout << std::boolalpha
                              << " Mass->isFactorized() = " << Mass->isFactorized() << "\n";);

      // DEBUG_EXPR(std::cout << (*Mass-*Mold).normInf() << std::endl;);
      /*Copy of the current mass matrix. */
      block = std::make_shared<siconos::algebra::SimpleMatrix>(*Mass);
    }
    else if (auto d = dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      //   d->computeMass();
      //   d->mass()->resetFactorizationFlags();
      DEBUG_EXPR(d->mass()->display(););
      block = std::make_shared<siconos::algebra::SimpleMatrix>(*(d->mass()));
    }
    else
      THROW_EXCEPTION(
          "siconos::nonsmooth_formulations::OneStepNS::getOSIMatrix for D1Minus, only implemented for "
          "Lagrangian or NewtonEuler");
  }
  // for ZeroOrderHoldOSI, the central block is Ad = \int exp{As} ds over t_k, t_{k+1}
  else if (osiType == siconos::integrators::IntegratorType::ZOHOSI) {
    if (!block)
      block = std::make_shared<siconos::algebra::SimpleMatrix>(
          (static_cast<siconos::integrators::ZeroOrderHoldOSI&>(Osi)).Ad(ds));
    else
      *block = (static_cast<siconos::integrators::ZeroOrderHoldOSI&>(Osi)).Ad(ds);
  }
  else
    THROW_EXCEPTION(
        "siconos::nonsmooth_formulations::OneStepNSProblem::getOSIMatrix not yet implemented for "
        "Integrator of type " +
        siconos::tools::enum_to_string(osiType));
  return block;
}

void siconos::nonsmooth_formulations::OneStepNSProblem::setSolverId(int solverId)
{
  // And create a new one, with default parameters values.
  _numerics_solver_options.reset(solver_options_create(solverId),
                                 solver_options_delete);
}

void siconos::nonsmooth_formulations::OneStepNSProblem::setNumericsVerboseMode(bool vMode)
{
  numerics_set_verbose(vMode);
}

void siconos::nonsmooth_formulations::OneStepNSProblem::setNumericsVerboseLevel(int level)
{
  numerics_set_verbose(level);
}
