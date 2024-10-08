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
#include "TimeSteppingCombinedProjection.hpp"

#include "Interaction.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "MLCPProjectOnConstraints.hpp"
#include "MoreauJeanOSI.hpp"
#include "NewtonEulerDS.hpp"
#include "NonSmoothLaw.hpp"
#include "OneStepNSProblem.hpp"
#include "Relation.hpp"
#include "SiconosVector.hpp"
#include "Topology.hpp"
// #define TSPROJ_DEBUG_LEVEL1
// #define TSPROJ_WITHOUT_PROJECTION
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

siconos::simulation::TimeSteppingCombinedProjection::TimeSteppingCombinedProjection(
    std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
    std::shared_ptr<TimeDiscretisation> td,
    std::shared_ptr<siconos::integrators::OneStepIntegrator> osi,
    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb_velo,
    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb_pos,
    unsigned int level)
    : TimeStepping{nsds, td, osi, osnspb_velo}, _indexSetLevelForProjection{level} {
  (*_allNSProblems).resize(SICONOS_NB_OSNSP_TSP);
  insertNonSmoothProblem(osnspb_pos, siconos::simulation::SICONOS_OSNSP_TS_POS);

  if (_indexSetLevelForProjection != 2) {
    THROW_EXCEPTION(
        "siconos::simulation::TimeSteppingCombinedProjection::"
        "TimeSteppingCombinedProjection "
        "level not equal to 2 is not yet implemented.  ");
  }
}

// struct
// siconos::simulation::TimeSteppingCombinedProjection::_SimulationEffectOnOSNSP
//   : public ssiconos::internal::SiconosVisitor {

//   TimeSteppingCombinedProjection* _parent{nullptr};

//   _SimulationEffectOnOSNSP(TimeSteppingCombinedProjection* p) : _parent(p)
//   {
//     std::cout << "hello\n";
//   };

//   void visit(MLCPProjectOnConstraints& onsnsp) const override
//   {
//     bool toto = (bool)_parent->doCombinedProjOnEquality();
//     onsnsp.setDoProjOnEquality(toto);
//   }
//   void visit(MLCPProjectOnConstraints* onsnsp) { std::cout << "hello\n"; }
//   void visit(MLCPProjectOnConstraints onsnsp) { std::cout << "hello\n"; }
// };

void siconos::simulation::TimeSteppingCombinedProjection::initializeOneStepNSProblem() {
  updateIndexSets();
  TimeStepping::initializeOneStepNSProblem();

  auto osnspb_pos = (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_POS];

  osnspb_pos->setIndexSetLevel(_indexSetLevelForProjection);
  osnspb_pos->setInputOutputLevel(0);

  (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->setIndexSetLevel(1);
  (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->setInputOutputLevel(1);

  if (auto mlcp =
          std::dynamic_pointer_cast<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>(
              osnspb_pos)) {
    mlcp->setDoProjOnEquality(_doCombinedProjOnEquality);
  }
}

void siconos::simulation::TimeSteppingCombinedProjection::advanceToEvent() {
  DEBUG_PRINT("================================================");
  DEBUG_PRINT("siconos::simulation::TimeSteppingCombinedProjection::advanceToEvent()");
  DEBUG_PRINT("================================================\n");

  initialize();

  _isIndexSetsStable = false;
  _maxViolationUnilateral = 0.0;
  _maxViolationEquality = 0.0;

  // Update interactions if a manager was provided
  updateInteractions();

  auto topo = _nsds->topology();
  if (topo->numberOfIndexSet() > _indexSetLevelForProjection) {
    auto indexSet1 = topo->indexSet(1);
    auto indexSet2 = topo->indexSet(2);
    assert(indexSet1);
    assert(indexSet2);

    siconos::graphs::InteractionsGraph::VIterator ui, uiend, vnext;

    // zeroing the lambda of indexSet1
    std::tie(ui, uiend) = indexSet1->vertices();
    for (vnext = ui; ui != uiend; ui = vnext) {
      ++vnext;
      auto inter1 = indexSet1->bundle(*ui);
      inter1->lambda(1)->zero();
      indexSet1->eraseProperties(*ui);
      siconos::graphs::InteractionsGraph::OEIterator oei, oeiend;
      for (std::tie(oei, oeiend) = indexSet1->out_edges(*ui); oei != oeiend; ++oei) {
        siconos::graphs::InteractionsGraph::EDescriptor ed1, ed2;
        std::tie(ed1, ed2) =
            indexSet1->edges(indexSet1->source(*oei), indexSet1->target(*oei));
        if (ed2 != ed1) {
          indexSet1->eraseProperties(ed1);
          indexSet1->eraseProperties(ed2);
        } else {
          indexSet1->eraseProperties(ed1);
        }
      }
    }

    indexSet1->clear();

    std::tie(ui, uiend) = indexSet2->vertices();
    for (vnext = ui; ui != uiend; ui = vnext) {
      ++vnext;
      indexSet2->eraseProperties(*ui);
      siconos::graphs::InteractionsGraph::OEIterator oei, oeiend;
      for (std::tie(oei, oeiend) = indexSet2->out_edges(*ui); oei != oeiend; ++oei) {
        siconos::graphs::InteractionsGraph::EDescriptor ed1, ed2;
        std::tie(ed1, ed2) =
            indexSet2->edges(indexSet2->source(*oei), indexSet2->target(*oei));
        if (ed2 != ed1) {
          indexSet2->eraseProperties(ed1);
          indexSet2->eraseProperties(ed2);
        } else {
          indexSet2->eraseProperties(ed1);
        }
      }
    }

    indexSet2->clear();
  }

  _nbIndexSetsIteration = 0;
  _cumulatedNewtonNbIterations = 0;
  _nbCumulatedProjectionIteration = 0;

  while (!_isIndexSetsStable) {
    _nbIndexSetsIteration++;
    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    DEBUG_PRINTF("====================== _nbIndexSetsIteration = %i \n ",
                 _nbIndexSetsIteration);

#ifdef TSPROJ_DEBUG_LEVEL1
    auto indexSet0 = topo->indexSet(0);
    std::cout << "indexSet0->size() " << indexSet0->size() << std::endl;
    unsigned int level;
    std::shared_ptr<siconos::modeling::Interaction> inter;
#endif

#ifdef TSPROJ_DEBUG_LEVEL2

    if (topo->numberOfIndexSet() > _indexSetLevelForProjection) {
      auto indexSet1 = topo->indexSet(1);
      auto indexSet2 = topo->indexSet(2);
      std::cout << "indexSet1->size() " << indexSet1->size() << std::endl;
      std::cout << "indexSet2->size() " << indexSet2->size() << std::endl;
    }

    for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
      inter = indexSet0->bundle(*ui);

      std::cout << "inter->number()" << inter->number() << std::endl;

      inter->computeOutput(getTkp1(), indexSet0->properties(*ui), 0);
      inter->computeOutput(getTkp1(), indexSet0->properties(*ui), 1);

      //  inter->swapInMemory();

      level = 0;

      assert(inter->lowerLevelForOutput() <= level);
      assert(inter->upperLevelForOutput() >= level);

      // std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS()     <<
      // std::endl;
      std::cout << "inter->y(" << level << ")\n";
      inter->y(level)->display();

      std::cout << "inter->y_k(" << level << ")\n";
      inter->y_k(level)->display();

      level = 1;
      assert(inter->lowerLevelForOutput() <= level);
      assert(inter->upperLevelForOutput() >= level);
      // std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS()     <<
      // std::endl;
      std::cout << "inter->y(" << level << ")\n";
      inter->y(level)->display();
      std::cout << "inter->y_k(" << level << ")\n";
      inter->y_k(level)->display();
    }

#endif

    if (_nbIndexSetsIteration > _kIndexSetMax) {
      THROW_EXCEPTION(
          "siconos::simulation::TimeSteppingCombinedProjection::"
          "TimeSteppingCombinedProjection _nbIndexSetsIteration >  "
          "_kIndexSetMax ");
    }

    /** First step, Solve the standard velocity formulation.*/
    TimeStepping::newtonSolve(_newtonTolerance, _newtonMaxIteration);
    _cumulatedNewtonNbIterations += getNewtonNbIterations();

#ifdef TSPROJ_DEBUG_LEVEL1

    for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
      inter = indexSet0->bundle(*ui);
      std::cout << "inter->number()" << inter->number() << std::endl;

      level = 0;

      assert(inter->lowerLevelForOutput() <= level);
      assert(inter->upperLevelForOutput() >= level);

      //   std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS()     <<
      //   std::endl;
      std::cout << "inter->y(" << level << ")\n";
      inter->y(level)->display();
      std::cout << "inter->y_k(" << level << ")\n";
      inter->y_k(level)->display();

      level = 1;
      assert(inter->lowerLevelForOutput() <= level);
      assert(inter->upperLevelForOutput() >= level);
      // std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS()     <<
      // std::endl;
      std::cout << "inter->y(" << level << ")\n";
      inter->y(level)->display();
      std::cout << "inter->y_k(" << level << ")\n";
      inter->y_k(level)->display();
    }
#endif

    int info = 0;

    // Zeroing Lambda Muliplier of indexSet()

    auto indexSet = _nsds->topology()->indexSet(0);
    for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
      auto inter = indexSet->bundle(*ui);
      inter->lambda(0)->zero();
    }
    _nsds->updateInput(nextTime(), 0);

#ifdef TSPROJ_WITHOUT_PROJECTION

#else
    /** Second step, Perform the projection on constraints.*/
    DEBUG_PRINT(
        "siconos::simulation::TimeSteppingCombinedProjection::newtonSolve "
        "begin "
        "projection:\n");

    auto dsGraph = _nsds->dynamicalSystems();

    bool runningProjection = false;
    _nbProjectionIteration = 0;

    if (_nsds->topology()->numberOfIndexSet() > _indexSetLevelForProjection) {
      updateIndexSet(2);
      computeCriteria(&runningProjection);

#ifdef TSPROJ_DEBUG_LEVEL1

      auto indexSet2 = topo->indexSet(2);
      auto indexSet1 = topo->indexSet(1);
      if (indexSet2->size() > 1) {
        printf("indexSet2->size() = %i >1 \n", (int)indexSet2->size());
      }
      if (indexSet1->size() > 0) {
        printf("indexSet1->size() = %i >0 \n", (int)indexSet1->size());
      }

#endif
    }

    if (!runningProjection) {
      // Zeroing Lambda Muliplier of indexSet()

      auto indexSet = _nsds->topology()->indexSet(0);
      siconos::graphs::InteractionsGraph::VIterator ui, uiend;
      for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
        auto inter = indexSet->bundle(*ui);
        inter->lambda(0)->zero();
      }
      _nsds->updateInput(nextTime(), 0);
    }

    // Store the q vector of each DS.

    for (auto aVi2 : *dsGraph) {  // = dsGraph->begin(); aVi2 != dsGraph->end(); ++aVi2) {
      auto ds = dsGraph->bundle(aVi2);
      auto& workVectors = *dsGraph->properties(aVi2).workVectors;

      if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
        *workVectors[siconos::integrators::MoreauJeanOSI::QTMP] = *neds->q();
      } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
        *workVectors[siconos::integrators::MoreauJeanOSI::QTMP] = *d->q();
      } else
        THROW_EXCEPTION(
            "siconos::simulation::TimeSteppingCombinedProjection::"
            "advanceToEvent() :: - Ds is "
            "not from NewtonEulerDS neither from LagrangianDS.");
    }

    _nbProjectionIteration = 0;

    while ((runningProjection && _nbProjectionIteration < _projectionMaxIteration) &&
           _doCombinedProj) {
      _nbProjectionIteration++;

      DEBUG_PRINTF("Projection iteration number   %d\t", _nbProjectionIteration);
      DEBUG_PRINT("================================================\n");

      // Zeroing Lambda Muliplier of indexSet()

      auto indexSet = _nsds->topology()->indexSet(0);
      siconos::graphs::InteractionsGraph::VIterator ui, uiend;
      for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
        auto inter = indexSet->bundle(*ui);
        inter->lambda(0)->zero();
      }

      info = 0;
#ifdef TSPROJ_DEBUG_LEVEL1
      std::cout << "TimeSteppingCombinedProjection compute OSNSP POS.\n";
#endif
      info = computeOneStepNSProblem(siconos::simulation::SICONOS_OSNSP_TS_POS);

      if (info) {
        std::cout << " siconos::simulation::TimeSteppingCombinedProjection::"
                     "advanceToEvent() "
                     "project on constraints. solver failed."
                  << std::endl;
        return;
      }

      _nsds->updateInput(nextTime(), 0);

      for (auto aVi2 : *dsGraph) {
        auto ds = dsGraph->bundle(aVi2);
        auto& workVectors = *dsGraph->properties(aVi2).workVectors;

        if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
          auto q = neds->q();
          auto qtmp = workVectors[siconos::integrators::MoreauJeanOSI::QTMP];
          if (neds->p(0)) {
            //*q = * qtmp +  *neds->p(0);
            *q += *neds->p(0);
          }
          neds->normalizeq();
          // siconos::modeling::newton_euler::computeT);

#ifdef TSPROJ_DEBUG_LEVEL1
          neds->display();
#endif
        } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
          auto q = d->q();
          auto qtmp = workVectors[siconos::integrators::MoreauJeanOSI::QTMP];
          if (d->p(0)) {
            //*q = * qtmp +  *d->p(0);
            *q += *d->p(0);
          }
#ifdef TSPROJ_DEBUG_LEVEL1
          std::cout << " q=\n";
          q->display();
          std::cout << " p(0)=\n";
          d->p(0)->display();
          std::cout << " p(1)=\n";
          d->p(1)->display();
#endif
        } else
          THROW_EXCEPTION(
              "siconos::simulation::TimeSteppingCombinedProjection::"
              "advanceToEvent() - Ds is "
              "not from NewtonEulerDS neither from LagrangianDS.");
      }

      updateWorldFromDS();

      computeCriteria(&runningProjection);

      // cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      // Z:"<<endl;
      //(*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_POS]->display();
      //(std::static_pointer_cast<LinearOSNS>((*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_POS]))->z()->display();

#ifdef TSPROJ_DEBUG_LEVEL1

      // auto indexSet1 = _nsds->topology()->indexSet(1);
      // std ::cout << "lambda(1) in IndexSet1\n";
      // for (std::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui)
      // {
      //   auto inter = indexSet1->bundle(*ui);
      //   inter->lambda(1)->display();
      // }
      auto indexSet2 = _nsds->topology()->indexSet(2);
      std ::cout << "lambda(0) in indexSet2\n";
      for (std::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui) {
        auto inter = indexSet2->bundle(*ui);
        inter->lambda(0)->display();
      }

#endif

      // cout<<"during projection before normalizing of q:\n";
      // for (InteractionsIterator it = allInteractions->begin(); it !=
      // allInteractions->end(); it++)
      //{
      //   (*it)->relation()->computeh(getTkp1());
      // }
    }  // end while ((runningProjection && _nbProjectionIteration <
       // _projectionMaxIteration) && _doCombinedProj)

    DEBUG_PRINTF(
        "siconos::simulation::TimeSteppingCombinedProjection::Projection end : "
        "Number of "
        "iterations= %i\n",
        _nbProjectionIteration);

    _nbCumulatedProjectionIteration += _nbProjectionIteration;
    if (_nbProjectionIteration == _projectionMaxIteration) {
      std::cout << "siconos::simulation::TimeSteppingCombinedProjection::"
                   "advanceToEvent() Max "
                   "number of projection iterations reached ("
                << _nbProjectionIteration << ")\n";
      printf("              max criteria equality =  %e.\n", _maxViolationEquality);
      printf("              max criteria unilateral =  %e.\n", _maxViolationUnilateral);
    }

#endif  // TSPROJ_WITHOUT_PROJECTION

    DEBUG_PRINT(
        "siconos::simulation::TimeSteppingCombinedProjection::newtonSolve end "
        "projection:\n");

    // We update forces to start the Newton Loop the next tiem step with a
    // correct value in swap
    for (auto aVi2 : *dsGraph) {
      auto ds = dsGraph->bundle(aVi2);
      if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
        double time = nextTime();
        neds->computeWrench(neds->twist_read(), neds->q_read(), time);
      } else if (std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(
                     ds)) {  // nothing ...
      } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
        auto time = nextTime();
        d->computeTotalForces(d->velocity_read(), d->q_read(), time);
      } else
        THROW_EXCEPTION(
            "siconos::simulation::TimeSteppingCombinedProjection::"
            "advanceToEvent() - Ds is "
            "not from NewtonEulerDS neither from LagrangianDS.");
    }

    if (_nsds->topology()->numberOfIndexSet() > _indexSetLevelForProjection) {
      updateIndexSet(1);
    } else {
      _isIndexSetsStable = true;
    }
#ifdef TSPROJ_DEBUG_LEVEL1

    if (topo->numberOfIndexSet() > _indexSetLevelForProjection) {
      auto indexSet1 = topo->indexSet(1);
      auto indexSet2 = topo->indexSet(2);
      std::cout << "indexSet1->size() " << indexSet1->size() << std::endl;
      std::cout << "indexSet2->size() " << indexSet2->size() << std::endl;
    }

    level = 0;
    for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
      inter = indexSet0->bundle(*ui);
      assert(inter->lowerLevelForOutput() <= level);
      assert(inter->upperLevelForOutput() >= level);
      inter->computeOutput(getTkp1(), indexSet0->properties(*ui), 0);
      inter->computeOutput(getTkp1(), indexSet0->properties(*ui), 1);

      std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS() << std::endl;
      std::cout << "inter->y(" << level << ")\n";
      inter->y(level)->display();
      std::cout << "inter->y_k(" << level << ")\n";
      inter->y_k(level)->display();
    }
    level = 1;
    for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
      inter = indexSet0->bundle(*ui);
      assert(inter->lowerLevelForOutput() <= level);
      assert(inter->upperLevelForOutput() >= level);
      std::cout << "inter->getSizeOfDS()" << inter->getSizeOfDS() << std::endl;
      std::cout << "inter->y(" << level << ")\n";
      inter->y(level)->display();
      std::cout << "inter->y_k(" << level << ")\n";
      inter->y_k(level)->display();
    }

#endif

  }  // end  while (!_isIndexSetsStable)

  DEBUG_PRINTF(
      "siconos::simulation::TimeSteppingCombinedProjection::indexset stable "
      "end : Number of "
      "iterations= %i \n ",
      _nbIndexSetsIteration);

  return;
}

void siconos::simulation::TimeSteppingCombinedProjection::computeCriteria(
    bool* runningProjection) {
  DEBUG_PRINT(
      "siconos::simulation::TimeSteppingCombinedProjection::computeCriteria("
      "bool * "
      "runningProjection)\n");
  // auto indexSet = _nsds->topology()->indexSet(_indexSetLevelForProjection);
  auto indexSet = _nsds->topology()->indexSet(_indexSetLevelForProjection);

  siconos::graphs::InteractionsGraph::VIterator aVi, viend;

  double maxViolationEquality = -1e24;
  double maxViolationUnilateral = -1e24;

  *runningProjection = false;

  for (std::tie(aVi, viend) = indexSet->vertices(); aVi != viend; ++aVi) {
    auto interac = indexSet->bundle(*aVi);

    interac->computeOutput(getTkp1(), 0);
    interac->relation()->computeJach(getTkp1(), *interac);

    if (siconos::types::type_value(*(interac->nonSmoothLaw())) ==
            siconos::modeling::Type::NewtonImpactFrictionNSL ||
        siconos::types::type_value(*(interac->nonSmoothLaw())) ==
            siconos::modeling::Type::NewtonImpactNSL) {
#ifdef TSPROJ_DEBUG_LEVEL1
      printf(
          "  "
          "siconos::simulation::TimeSteppingCombinedProjection::"
          "computeCriteria  "
          "Unilateral "
          "interac->y(0)->getValue(0) %e.\n",
          interac->y(0)->getValue(0));
#endif
      if (!_doCombinedProjOnEquality) {
        if (maxViolationUnilateral > _constraintTolUnilateral) {
          double criteria = std::max(0.0, -interac->y(0)->getValue(0));
          if (criteria > maxViolationUnilateral) maxViolationUnilateral = criteria;

          *runningProjection = true;
#ifdef TSPROJ_DEBUG_LEVEL1
          printf("TSProj newton criteria unilateral true %e.\n", criteria);
#endif
        }
      } else {
        auto criteria = interac->y(0)->getValue(0);
        if (criteria > maxViolationUnilateral) maxViolationUnilateral = criteria;

        if (std::abs(criteria) >= _constraintTolUnilateral) {
          *runningProjection = true;
#ifdef TSPROJ_DEBUG_LEVEL1
          printf("TSProj newton criteria unilateral true %e.\n", criteria);
#endif
        }
      }
    } else {
      DEBUG_EXPR(interac->y(0)->display(););
      if (interac->y(0)->normInf() > maxViolationEquality) {
        DEBUG_EXPR(interac->y(0)->display(););
        maxViolationEquality = interac->y(0)->normInf();
      }
      if (interac->y(0)->normInf() > _constraintTol) {
        *runningProjection = true;
#ifdef TSPROJ_DEBUG_LEVEL1
        printf("TSProj  newton criteria equality true %e.\n", interac->y(0)->normInf());
#endif
      }
    }
  }

  _maxViolationUnilateral = maxViolationUnilateral;
  _maxViolationEquality = maxViolationEquality;

  DEBUG_PRINTF("              max criteria equality =  %e.\n", _maxViolationEquality);
  DEBUG_PRINTF("              max criteria unilateral =  %e.\n", _maxViolationUnilateral);

#ifdef TSPROJ_DEBUG_LEVEL1
  printf("TSProj newton min/max criteria projection\n");
  std::cout << "                 runningProjection  " << *runningProjection << std::endl;
  printf("              max criteria equality =  %e.\n", maxViolationEquality);
  printf("              max criteria unilateral =  %e.\n", maxViolationUnilateral);
  //  printf("              min criteria unilateral =
  //  %e.\n",minViolationUnilateral);
#endif
}

void siconos::simulation::TimeSteppingCombinedProjection::updateIndexSet(unsigned int i) {
  // To update IndexSet i: add or remove Interactions from
  // this set, depending on y values.
  // boost::default_color_type is used to organize update in InteractionsGraph:
  // - white_color : undiscovered vertex (Interaction)
  // - gray_color : discovered vertex (Interactions) but searching descendants
  // - black_color : discovered vertex (Interaction) together with the
  // descendants

  assert(_nsds);
  assert(_nsds->topology());

  auto topo = _nsds->topology();

  assert(i < topo->indexSetsSize() &&
         "TimeStepping::updateIndexSet(i), indexSets[i] does not exist.");
  // IndexSets[0] must not be updated in simulation, since it belongs to
  // Topology.
  assert(i > 0 && "TimeStepping::updateIndexSet(i=0), indexSets[0] cannot be updated.");

  // For all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i].

  auto indexSet0 = topo->indexSet(0);
  auto indexSet1 = topo->indexSet(1);
  auto indexSet2 = topo->indexSet(2);
  auto& DSG0 = *nonSmoothDynamicalSystem()->dynamicalSystems();
  assert(indexSet0);
  assert(indexSet1);
  assert(indexSet2);

  topo->setHasChanged(false);

  DEBUG_PRINTF("update indexSets start : indexSet0 size : %i\n", (int)(indexSet0->size()));

  // Check indexSet1

  if (i == 1) {
    siconos::graphs::InteractionsGraph::VIterator ui1, ui1end, v1next;

    std::tie(ui1, ui1end) = indexSet1->vertices();
    _isIndexSetsStable = true;

    DEBUG_PRINTF("update IndexSets start : indexSet1 size : %i\n", (int)(indexSet1->size()));
    // indexSet1->display();
    // Remove interactions from the indexSet1
    for (v1next = ui1; ui1 != ui1end; ui1 = v1next) {
      ++v1next;
      auto inter1 = indexSet1->bundle(*ui1);
      if (indexSet0->is_vertex(inter1)) {
        auto ur1_descr0 = indexSet0->descriptor(inter1);
        assert((indexSet0->color(ur1_descr0) == boost::white_color));
        indexSet0->color(ur1_descr0) = boost::gray_color;
      }
      // else
      // {
      //   // Interactions is not in indexSet0 anymore.
      //   // ui1 becomes invalid
      //   indexSet1->remove_vertex(inter1);
      //   topo->setHasChanged(true);
      //   _isIndexSetsStable=false;
      // }
    }

    // indexSet0\indexSet1 scan
    siconos::graphs::InteractionsGraph::VIterator ui0, ui0end;
    // Add interaction in indexSet1
    for (std::tie(ui0, ui0end) = indexSet0->vertices(); ui0 != ui0end; ++ui0) {
      if (indexSet0->color(*ui0) == boost::black_color) {
        // reset
        indexSet0->color(*ui0) = boost::white_color;
      } else {
        if (indexSet0->color(*ui0) == boost::gray_color) {
          // reset
          indexSet0->color(*ui0) = boost::white_color;

          assert(indexSet1->is_vertex(indexSet0->bundle(*ui0)));
          /*assert( { !predictorDeactivate(indexSet0->bundle(*ui0),i) ||
            siconos::types::type_value(*(indexSet0->bundle(*ui0)->interaction()->nonSmoothLaw()))
            == siconos::modeling::Type::EqualityConditionNSL ;
            });*/
        } else {
          assert(indexSet0->color(*ui0) == boost::white_color);

          auto inter0 = indexSet0->bundle(*ui0);
          assert(!indexSet1->is_vertex(inter0));

          bool activate = true;
          if (siconos::types::type_value(*(inter0->nonSmoothLaw())) !=
              siconos::modeling::Type::EqualityConditionNSL) {
            // auto Osi = indexSet0->properties(*ui0).osi;
            //  We assume that the integrator of the ds1 drive the update of the
            //  index set
            auto ds1 = indexSet1->properties(*ui0).source;
            auto& osi = *DSG0.properties(DSG0.descriptor(ds1)).osi;

            activate = osi.addInteractionInIndexSet(inter0, i);
          }
          if (activate) {
            assert(!indexSet1->is_vertex(inter0));

            // vertex and edges insertion in indexSet1
            indexSet1->copy_vertex(inter0, *indexSet0);
            topo->setHasChanged(true);
            _isIndexSetsStable = false;
            assert(indexSet1->is_vertex(inter0));
          }
        }
      }
    }
    indexSet1->update_vertices_indices();
    indexSet1->update_edges_indices();
    assert(indexSet1->size() <= indexSet0->size());
    DEBUG_PRINTF("update indexSets end : indexSet0 size : %i\n", (int)(indexSet0->size()));
    DEBUG_PRINTF("update IndexSets end : indexSet1 size : %i\n", (int)(indexSet1->size()));
  }  // i==1

  if (i == 2) {
    siconos::graphs::InteractionsGraph::VIterator ui1, ui1end, v1next;
    std::tie(ui1, ui1end) = indexSet2->vertices();

    for (v1next = ui1; ui1 != ui1end; ui1 = v1next) {
      ++v1next;
      indexSet2->eraseProperties(*ui1);
      siconos::graphs::InteractionsGraph::OEIterator oei, oeiend;
      for (std::tie(oei, oeiend) = indexSet2->out_edges(*ui1); oei != oeiend; ++oei) {
        siconos::graphs::InteractionsGraph::EDescriptor ed1, ed2;
        std::tie(ed1, ed2) =
            indexSet2->edges(indexSet2->source(*oei), indexSet2->target(*oei));
        if (ed2 != ed1) {
          indexSet2->eraseProperties(ed1);
          indexSet2->eraseProperties(ed2);
        } else {
          indexSet2->eraseProperties(ed1);
        }
      }
    }

    indexSet2->clear();
    DEBUG_PRINTF("update IndexSets start : indexSet2 size : %i\n", (int)(indexSet2->size()));

    // Scan indexSet1
    std::tie(ui1, ui1end) = indexSet1->vertices();
    for (v1next = ui1; ui1 != ui1end; ui1 = v1next) {
      ++v1next;
      auto inter1 = indexSet1->bundle(*ui1);
      bool activate = true;
      if (siconos::types::type_value(*(inter1->nonSmoothLaw())) !=
          siconos::modeling::Type::EqualityConditionNSL) {
        // auto Osi = indexSet1->properties(*ui1).osi;
        //  We assume that the integrator of the ds1 drive the update of the
        //  index set
        auto ds1 = indexSet1->properties(*ui1).source;
        auto& osi = *DSG0.properties(DSG0.descriptor(ds1)).osi;

        activate = osi.addInteractionInIndexSet(inter1, i);
      }
      if (activate) {
        assert(!indexSet2->is_vertex(inter1));

        // vertex and edges insertion in indexSet2
        indexSet2->copy_vertex(inter1, *indexSet1);
        topo->setHasChanged(true);
        assert(indexSet2->is_vertex(inter1));
      }
    }
    DEBUG_PRINTF("update IndexSets end : indexSet0 size : %i\n", (int)(indexSet0->size()));
    DEBUG_PRINTF("update IndexSets end : indexSet2 size : %i\n", (int)(indexSet2->size()));
    indexSet2->update_vertices_indices();
    indexSet2->update_edges_indices();
  }
}
