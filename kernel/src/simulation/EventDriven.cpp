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
#include "EventDriven.hpp"

#include <memory>

#include "EventsManager.hpp"
#include "FirstOrderNonLinearDS.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianSparseDS.hpp"
#include "LsodarOSI.hpp"
#include "NewMarkAlphaOSI.hpp"
#include "NonSmoothLaw.hpp"
#include "OneStepNSProblem.hpp"
#include "Relation.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Topology.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

/** defaut constructor
 *  \param a pointer to a timeDiscretisation (linked to the model that owns this
 * simulation)
 */
siconos::simulation::EventDriven::EventDriven(
    std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
    std::shared_ptr<TimeDiscretisation> td)
    : Simulation{nsds, td} {
  (*_allNSProblems).resize(_numberOfOneStepNSproblems);
}

siconos::simulation::EventDriven::EventDriven(
    std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
    std::shared_ptr<TimeDiscretisation> td, std::size_t nb)
    : Simulation{nsds, td}, _numberOfOneStepNSproblems{0} {
  (*_allNSProblems).resize(nb);
}

void siconos::simulation::EventDriven::insertIntegrator(
    std::shared_ptr<siconos::integrators::OneStepIntegrator> osi) {
  Simulation::insertIntegrator(osi);
  // Determine the number of OneStepNSproblems depending on the
  // OneStepIntegrator type
  auto osiType = osi->getType();
  if (osiType ==
      siconos::integrators::IntegratorType::NEWMARKALPHAOSI)  // EventDrivent asscociated with
                                                              // NewMarkAlpha
  {
    _numberOfOneStepNSproblems = 3;
    if (_allNSProblems->size() != 3) {
      (*_allNSProblems).resize(_numberOfOneStepNSproblems);
    }
  }
}

void siconos::simulation::EventDriven::updateIndexSet(
    siconos::simulation::Topology::size_type i) {
  DEBUG_BEGIN("siconos::simulation::EventDriven::updateIndexSet(i)\n");
  DEBUG_PRINTF("with i = %i\n", i);
  assert(_nsds);
  assert(_nsds->topology());
  auto topo = _nsds->topology();

  assert(i < topo->indexSetsSize() &&
         "siconos::simulation::EventDriven::updateIndexSet(i), indexSets[i] "
         "does not exist.");
  // IndexSets[0] must not be updated in simulation, since it belongs to
  // Topology.
  assert(i > 0 &&
         "siconos::simulation::EventDriven::updateIndexSet(i=0), indexSets[0] "
         "cannot be "
         "updated.");

  // For all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i].
  auto indexSet1 = topo->indexSet(1);
  auto indexSet2 = topo->indexSet(2);
  assert(_indexSet0);
  assert(indexSet1);
  assert(indexSet2);

  DEBUG_PRINTF("update indexSets start : _indexSet0 size : %ld\n", _indexSet0->size());
  DEBUG_PRINTF("update IndexSets start : indexSet1 size : %ld\n", indexSet1->size());
  DEBUG_PRINTF("update IndexSets start : indexSet2 size : %ld\n", indexSet2->size());

  siconos::graphs::InteractionsGraph::VIterator uibegin, uipend, uip;
  std::tie(uibegin, uipend) = _indexSet0->vertices();
  // loop over all vertices of the indexSet[i-1]
  for (uip = uibegin; uip != uipend; ++uip) {
    auto inter = _indexSet0->bundle(*uip);
    if (i == 1)  // IndexSet[1]
    {
      // if indexSet[1]=>getYRef(0): output y
      // if indexSet[2]=>getYRef(1): output ydot
      auto y = (*inter->y(0))(0);  // output to define the IndexSets at this Interaction
      if (y < -tolerance_)         // y[0] < 0
      {
        inter->display();
        std::cout << "y = " << y << " < -tolerance_ =  " << -tolerance_ << "\n";
        THROW_EXCEPTION(
            "siconos::simulation::EventDriven::updateIndexSet, output of level "
            "0 must be "
            "positive!!! ");
      }
      // 1 - If the Interaction is not yet in the set
      if (!indexSet1->is_vertex(inter))  // Interaction is not yet in the indexSet[i]
      {
        if (fabs(y) <= tolerance_) {
          // vertex and edges insertions
          indexSet1->copy_vertex(inter, *_indexSet0);
        }
      } else  // if the Interaction was already in the set
      {
        if (fabs(y) > tolerance_) {
          indexSet1->remove_vertex(inter);  // remove the Interaction from IndexSet[1]
          inter->lambda(1)->setZero();      // reset the lambda[1] to zero
        }
      }
    } else if (i == 2)  // IndexSet[2]
    {
      if (indexSet1->is_vertex(inter))  // Interaction is in the indexSet[1]
      {
        auto y = (*inter->y(1))(0);        // output of level 1 at this Interaction
        if (!indexSet2->is_vertex(inter))  // Interaction is not yet in the indexSet[2]
        {
          if (fabs(y) <= tolerance_) {
            // vertex and edges insertions
            indexSet2->copy_vertex(inter, *_indexSet0);
          }
        } else  // if the Interaction was already in the set
        {
          if (fabs(y) > tolerance_) {
            indexSet2->remove_vertex(inter);  // remove the Interaction from IndexSet[1]
            inter->lambda(2)->setZero();      // reset the lambda[i] to zero
          }
        }
      } else  // Interaction is not in the indexSet[1]
      {
        if (indexSet2->is_vertex(inter))  // Interaction is in the indexSet[2]
        {
          indexSet2->remove_vertex(inter);  // remove the Interaction from IndexSet[2]
          inter->lambda(2)->setZero();      // reset the lambda[i] to zero
        }
      }
    } else {
      THROW_EXCEPTION(
          "siconos::simulation::EventDriven::updateIndexSet, IndexSet[i > 2] "
          "doesn't exist");
    }
  }

  DEBUG_PRINTF("update indexSets end : _indexSet0 size : %ld\n", _indexSet0->size());
  DEBUG_PRINTF("update IndexSets end : indexSet1 size : %ld\n", indexSet1->size());
  DEBUG_PRINTF("update IndexSets end : indexSet2 size : %ld\n", indexSet2->size());
  DEBUG_END("siconos::simulation::EventDriven::updateIndexSet(i)\n");
}

void siconos::simulation::EventDriven::updateIndexSetsWithDoubleCondition() {
  assert(_nsds);
  assert(_nsds->topology());

  // for all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i]

  auto topo = _nsds->topology();

  auto indexSet2 = topo->indexSet(2);

  siconos::graphs::InteractionsGraph::VIterator ui, uiend, vnext;
  std::tie(ui, uiend) = indexSet2->vertices();

  for (vnext = ui; ui != uiend; ui = vnext) {
    ++vnext;

    auto inter = indexSet2->bundle(*ui);
    auto gamma = (*inter->y(2))(0);
    auto F = (*inter->lambda(2))(0);
    DEBUG_PRINTF("ED 1 update with double condition%f\n", F);
    DEBUG_PRINTF("ED 2 update with double condition %f\n", gamma);
    DEBUG_PRINTF("ED 3 update with double condition%f\n", tolerance_);

    if (fabs(F) < tolerance_)
      indexSet2->remove_vertex(inter);
    else if ((gamma < -tolerance_) || (F < -tolerance_)) {
      THROW_EXCEPTION(
          "siconos::simulation::EventDriven::"
          "updateIndexSetsWithDoubleCondition(), output[2] "
          "and lambda[2] for Interaction of indexSet[2] must be nonnegative.");
    } else if (((fabs(gamma) > tolerance_) && (fabs(F) > tolerance_))) {
      THROW_EXCEPTION(
          "siconos::simulation::EventDriven::"
          "updateIndexSetsWithDoubleCondition(), something "
          "is wrong for the LCP resolution.");
    } else {
    };
    DEBUG_PRINTF("End update with double condition %f\n", tolerance_);
  }
}
void siconos::simulation::EventDriven::initializeOneStepNSProblem() {
  DEBUG_BEGIN("siconos::simulation::EventDriven::initializeOneStepNSProblem()\n");
  assert(_nsds);
  assert(_nsds->topology());
  // for all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i]
  // Note that interactions set may be empty.
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto topo = _nsds->topology();

  // === update all index sets ===
  // updateIndexSets();
  initOSIRhs();

  if (!_allNSProblems->empty())  // ie if at least a non smooth problem has been built.
  {
    auto osiType = (*_allOSI->begin())->getType();
    if (osiType == siconos::integrators::IntegratorType::LSODAROSI ||
        osiType == siconos::integrators::IntegratorType::HEM5OSI)  // EventDriven
                                                                   // associated with
                                                                   // LsodarOSI OSI
    {
    } else if (osiType ==
               siconos::integrators::IntegratorType::NEWMARKALPHAOSI)  // EventDriven
                                                                       // associated with
                                                                       // NewMarkAlpha
    {
      if (_allNSProblems->size() != 3)
        THROW_EXCEPTION(
            " siconos::simulation::EventDriven::initialize, \n an EventDriven "
            "simulation "
            "associated with NewMarkAlphaOSI must have three non smooth "
            "problems.\n Here, "
            "there are " +
            std::to_string(_allNSProblems->size()));
      // Initialize OSNSP at position level
      (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_POS]->setInputOutputLevel(
          2);
      (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_POS]->setIndexSetLevel(2);
      (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_POS]->initialize(
          shared_from_this());
    } else {
      THROW_EXCEPTION(
          " siconos::simulation::EventDriven::initialize, OSI not yet "
          "implemented.");
    }

    if (!((*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_IMPACT])) /* ie if the impact
                                                                             * problem does not
                                                                             *  exist */
      THROW_EXCEPTION(
          "siconos::simulation::EventDriven::initialize, an EventDriven "
          "simulation must have "
          "an 'impact' non smooth problem.");

    if (!((*_allNSProblems)
              [siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC])) /* ie if the
                                                                    * acceleration-level
                                                                    * problem does
                                                                    * not exist */
      THROW_EXCEPTION(
          "siconos::simulation::EventDriven::initialize, an EventDriven "
          "simulation must have "
          "an 'acceleration' non smooth problem.");

    // Initialize OSNSP for impact problem and at the acceleration level
    // WARNING: only for Lagrangian systems - To be reviewed for other ones.
    (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_IMPACT]->setInputOutputLevel(1);
    (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_IMPACT]->setIndexSetLevel(1);
    (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_IMPACT]->initialize(
        shared_from_this());

    (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]->setInputOutputLevel(
        2);
    (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]->setIndexSetLevel(2);
    (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]->initialize(
        shared_from_this());
    //
    // // Detect NonSmoothEvent at the beginning of the simulation
    // if(topo->indexSetsSize() > 1)
    // {
    //   auto indexSet1 = _nsds->topology()->indexSet(1);
    //   if(indexSet1->size() != 0)  // There is one non-smooth event to be
    //   added
    //   {
    //     DEBUG_PRINT("Schedule an event at starting time\n");
    //     _eventsManager->scheduleNonSmoothEvent(*this,
    //     _eventsManager->startingTime(), false);
    //   };
    // }
  }
  DEBUG_END("EventDriven::initializeOneStepNSProblem()\n");
}

void siconos::simulation::EventDriven::initOSIs() {}

void siconos::simulation::EventDriven::initOSIRhs() {
  DEBUG_BEGIN("void siconos::simulation::EventDriven::initOSIRhs()\n")
  // === initialization for OneStepIntegrators ===
  auto osiType = (*_allOSI->begin())->getType();
  for (auto itosi : *_allOSI) {
    // Check whether OSIs used are of the same type
    if ((itosi)->getType() != osiType) THROW_EXCEPTION("OSIs used must be of the same type");

    // perform the initialization
    siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
    auto osiDSGraph = (itosi)->dynamicalSystemsGraph();
    for (std::tie(dsi, dsend) = osiDSGraph->vertices(); dsi != dsend; ++dsi) {
      if (!itosi->checkOSI(dsi)) continue;

      auto ds = osiDSGraph->bundle(*dsi);
      // Initialize right-hand side
      ds->initRhs(startingTime());
      // DEBUG_EXPR(ds->display());
    }
  }
  DEBUG_END("void siconos::simulation::EventDriven::initOSIRhs()\n")
}

void siconos::simulation::EventDriven::firstInitialize() {
  if (!_isInitialized) {
    DEBUG_PRINT(" - 6 - First initialization of the simulation\n");

    _T = _nsds->finalT();

    // === Events manager initialization ===
    _eventsManager->initialize(_T);
    _tinit = _eventsManager->startingTime();

    // Process events at time _tinit. Useful to save values in memories
    // for example.  Warning: can not be called during
    // eventsManager->initialize, because it needs the initialization of
    // OSI, OSNS ...
    // _eventsManager->preUpdate(*this);

    _tend = _eventsManager->nextTime();

    // End of initialize:

    //  - all OSI and OSNS (ie DS and Interactions) states are computed
    //  - for time _tinit and saved into memories.
    //  - Sensors or related objects are updated for t=_tinit.
    //  - current time of the model is equal to t1, time of the first
    //  - event after _tinit.
    //  - currentEvent of the simu. corresponds to _tinit and nextEvent
    //  - to _tend.

    auto topo = _nsds->topology();
    // Detect NonSmoothEvent at the beginning of the simulation
    if (topo->indexSetsSize() > 1) {
      auto indexSet1 = _nsds->topology()->indexSet(1);
      if (indexSet1->size() != 0)  // There is one non-smooth event to be added
      {
        DEBUG_PRINT("Schedule an event at starting time\n");
        _eventsManager->scheduleNonSmoothEvent(*this, _eventsManager->startingTime(), false);
      };
    }

    _isInitialized = true;
  }
}

void siconos::simulation::EventDriven::initialize() {
  DEBUG_BEGIN("void siconos::simulation::EventDriven::initialize()\n");

  // Initialization for Simulation
  _indexSet0 = _nsds->topology()->indexSet(0);
  _DSG0 = _nsds->topology()->dSG(0);

  Simulation::initialize();

  // Initialization for all OneStepIntegrators
  // initOSIs();
  initOSIRhs();
  DEBUG_END("void siconos::simulation::EventDriven::initialize()\n")
}

void siconos::simulation::EventDriven::computef(siconos::integrators::OneStepIntegrator& osi,
                                                int* sizeOfX, double* time, double* x,
                                                double* xdot) {
  DEBUG_BEGIN(
      "siconos::simulation::EventDriven::computef(OneStepIntegrator& osi, "
      "int * sizeOfX, "
      "double * time, double * x, "
      "double * xdot)\n");
  // computeF is supposed to fill xdot in, using the definition of the
  // dynamical systems belonging to the osi

  // Check osi type: only lsodar is allowed.
  assert(osi.getType() == siconos::integrators::IntegratorType::LSODAROSI);

  auto& lsodar = static_cast<siconos::integrators::LsodarOSI&>(osi);
  // fill in xWork vector (ie all the x of the ds of this osi) with x
  assert(*sizeOfX >= 0);
  std::size_t sizex = static_cast<std::size_t>(*sizeOfX);
  lsodar.fillXWork(sizex, x);

  double t = *time;
  // Update Jacobian matrices at all interactions
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = _indexSet0->vertices(); ui != uiend; ++ui) {
    auto& inter = *_indexSet0->bundle(*ui);
    inter.relation()->computeJach(t, inter);
  }

  // solve a LCP at "acceleration" level if required
  if (!_allNSProblems->empty()) {
    if (((*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->hasInteractions())) {
      // Update the state of the DS
      (*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->compute(t);
      _nsds->updateInput(t, 2);  // Necessary to compute DS state below
    }
    // Compute the right-hand side ( xdot = f + r in DS) for all the
    // ds, with the new value of input.  lsodar->computeRhs(t);
  }

  // update the DS of the OSI.
  lsodar.computeRhs(t);
  //  for the DS state, ie the ones computed by lsodar (x above)
  // Update Index sets? No !!

  // Get the required value, ie xdot for output.
  siconos::algebra::Index pos = 0;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  auto osiDSGraph = lsodar.dynamicalSystemsGraph();
  for (std::tie(dsi, dsend) = osiDSGraph->vertices(); dsi != dsend; ++dsi) {
    if (!(lsodar.checkOSI(dsi))) continue;

    auto ds = osiDSGraph->bundle(*dsi);
    if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto qdot = lds->velocity_read();
      auto acc = lds->acceleration_read();
      assert((pos + acc.size() + qdot.size()) <= siconos::algebra::to_index(*sizeOfX) &&
             "Destination buffer too small!");
      std::copy(qdot.data(), qdot.data() + qdot.size(), &xdot[pos]);
      pos += qdot.size();
      std::copy(acc.data(), acc.data() + acc.size(), &xdot[pos]);
      pos += acc.size();
    } else if (auto lds =
                   std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseDS>(ds)) {
      auto qdot = lds->velocity_read();
      auto acc = lds->acceleration_read();
      assert((pos + acc.size() + qdot.size()) <= siconos::algebra::to_index(*sizeOfX) &&
             "Destination buffer too small!");
      std::copy(qdot.data(), qdot.data() + qdot.size(), &xdot[pos]);
      pos += qdot.size();
      std::copy(acc.data(), acc.data() + acc.size(), &xdot[pos]);
      pos += acc.size();
    } else {
      auto rhs = ds->rhs_read();
      assert(pos + rhs.size() <= *sizeOfX && "Destination buffer too small!");
      std::copy(rhs.data(), rhs.data() + rhs.size(), &xdot[pos]);
      pos += rhs.size();
    }
  }
  DEBUG_END(
      "siconos::simulation::EventDriven::computef(OneStepIntegrator& osi, "
      "int * sizeOfX, "
      "double * time, double * x, "
      "double * xdot)\n");
}

void siconos::simulation::EventDriven::computeJacobianfx(
    siconos::integrators::OneStepIntegrator& osi, int* sizeOfX, double* time, double* x,
    double* jacob) {
  assert(osi.getType() == siconos::integrators::IntegratorType::LSODAROSI);

  auto& lsodar = static_cast<siconos::integrators::LsodarOSI&>(osi);

  // Remark A: according to DLSODAR doc, each call to jacobian is
  // preceded by a call to f with the same arguments NEQ, T, and Y.
  // Thus to gain some efficiency, intermediate quantities shared by
  // both calculations may be saved in class members?  fill in xWork
  // vector (ie all the x of the ds of this osi) with x fillXWork(x);
  // -> copy
  // Maybe this step is not necessary?  because of
  // remark A above

  // Compute the jacobian of the vector field according to x for the
  // current ds
  double t = *time;
  lsodar.computeJacobianRhs(t, *_DSG0);

  // Save jacobianX values from dynamical system into current jacob
  // (in-out parameter)

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  auto osiDSGraph = lsodar.dynamicalSystemsGraph();
  for (std::tie(dsi, dsend) = osiDSGraph->vertices(); dsi != dsend; ++dsi) {
    if (!(lsodar.checkOSI(dsi))) continue;

    auto ds = osiDSGraph->bundle(*dsi);
    const auto& jacotmp = ds->jacobianRhsOver_x();
    std::copy(jacotmp.begin(), jacotmp.end(), jacob);
  }
}

void siconos::simulation::EventDriven::computeg(
    std::shared_ptr<siconos::integrators::OneStepIntegrator> osi, int* sizeOfX, double* time,
    double* x, int* ng, double* gOut) {
  assert(_nsds);
  assert(_nsds->topology());
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto topo = _nsds->topology();
  auto indexSet2 = topo->indexSet(2);
  siconos::algebra::Index nsLawSize, k = 0;
  std::shared_ptr<siconos::algebra::SiconosVector> y, ydot, yddot, lambda;
  auto lsodar = std::dynamic_pointer_cast<siconos::integrators::LsodarOSI>(osi);
  assert(lsodar);
  // Solve LCP at acceleration level to calculate the lambda[2] at Interaction
  // of indexSet[2]
  assert(*sizeOfX >= 0);
  std::size_t sizex = static_cast<std::size_t>(*sizeOfX);
  lsodar->fillXWork(sizex, x);
  //
  double t = *time;
  if (!_allNSProblems->empty()) {
    if (((*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]
             ->hasInteractions())) {
      (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]->compute(t);
    }
  };
  /*
     double * xdottmp = (double *)malloc(*sizeOfX*sizeof(double));
     computef(osi, sizeOfX,time,x,xdottmp);
     free(xdottmp);
     */
  // Update the output from level 0 to level 1
  _nsds->updateOutput(t, 0);
  _nsds->updateOutput(t, 1);
  _nsds->updateOutput(t, 2);
  //
  for (std::tie(ui, uiend) = _indexSet0->vertices(); ui != uiend; ++ui) {
    auto inter = _indexSet0->bundle(*ui);
    nsLawSize = inter->nonSmoothLaw()->size();
    y = inter->y(0);     // output y at this Interaction
    ydot = inter->y(1);  // output of level 1 at this Interaction
    yddot = inter->y(2);
    lambda = inter->lambda(2);           // input of level 2 at this Interaction
    if (!(indexSet2->is_vertex(inter)))  // if Interaction is not in the indexSet[2]
    {
      for (decltype(nsLawSize) i = 0; i < nsLawSize; ++i) {
        if ((*y)(i) > tolerance_) {
          gOut[k] = (*y)(i);
        } else {
          if ((*ydot)(i) > -tolerance_) {
            gOut[k] = 100 * tolerance_;
          } else {
            gOut[k] = (*y)(i);
          }
        }
        k++;
      }
    } else  // If Interaction is in the indexSet[2]
    {
      for (decltype(nsLawSize) i = 0; i < nsLawSize; ++i) {
        if ((*lambda)(i) > tolerance_) {
          gOut[k] = (*lambda)(i);  // g = lambda[2]
        } else {
          if ((*yddot)(i) > tolerance_) {
            gOut[k] = (*lambda)(i);
          } else {
            gOut[k] = 100 * tolerance_;
          }
        }
        k++;
      }
    }
  }
}
void siconos::simulation::EventDriven::updateImpactState() {
  // Compute input = R(lambda[1])
  _nsds->updateInput(nextTime(), 1);

  // Compute post-impact velocity
  for (auto itOSI : *_allOSI) itOSI->updateState(1);
}

void siconos::simulation::EventDriven::updateSmoothState() {
  // Update input of level 2
  _nsds->updateInput(nextTime(),
                     2);  // Note FP : Probably already up to date? (previous
                          // call of updateInput in simu)
  // Compute acceleration
  for (auto itOSI : *_allOSI) itOSI->updateState(2);
}

void siconos::simulation::EventDriven::updateState(unsigned int levelInput) {
  assert(levelInput <= 2);
  if (levelInput == 1) {
    updateImpactState();
  } else {
    updateSmoothState();
  }
}

void siconos::simulation::EventDriven::updateOutput(unsigned int levelInput) {
  // Update output (y)
  _nsds->updateOutput(nextTime(), levelInput);
  // Warning: index sets are not updated in this function !!
}

void siconos::simulation::EventDriven::advanceToEvent() {
  DEBUG_BEGIN("siconos::simulation::EventDriven::advanceToEvent()\n");

  initialize();

  // Update interactions if a manager was provided
  updateInteractions();

  _tinit = _eventsManager->startingTime();
  _tend = _eventsManager->nextTime();
  _tout = _tend;

  DEBUG_PRINTF("_tinit = %g, _tend = %g \n", _tinit, _tend);

  bool isNewEventOccur = false;  // set to true if a new event occur during integration
  auto osiType = (*_allOSI->begin())->getType();  // Type of OSIs
  double _minConstraint = 0.0;

  // Initialize lambdas of all interactions.
  siconos::graphs::InteractionsGraph::VIterator ui, uiend, vnext;
  std::tie(ui, uiend) = _indexSet0->vertices();
  for (vnext = ui; ui != uiend; ui = vnext) {
    ++vnext;
    _indexSet0->bundle(*ui)->resetAllLambda();
  }

  if (osiType == siconos::integrators::IntegratorType::NEWMARKALPHAOSI) {
    // if the time to next event if too small, we skip the integration
    // it may happen if a first event is of time nonsmooth at the initial time
    // with _tout= _tend
    if (fabs(_tend - _tinit) >= 10 * siconos::internal::MACHINE_PREC) {
      newtonSolve(_newtonTolerance, _newtonMaxIteration);
    } else {
      auto indexSet2 = _nsds->topology()->indexSet(2);
      if (indexSet2->size() != 0)  // if indexSet2 is not empty, solve LCP to
                                   // determine contact forces
      {
        int info = computeOneStepNSProblem(siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_POS);
        if (info != 0) {
          std::cout << "Warning!!!In siconos::simulation::EventDriven::newtonSolve: "
                       "LCP solver may "
                       "fail"
                    << "\n";
        }
      }
      updateInput(2);
    }
    // Update after Newton iteration
    // Update input of level 2 >>> has already been done in newtonSolve
    // Update state of all Dynamicall Systems >>>  has already been done in
    // newtonSolve Update outputs of levels 0, 1, 2
    updateOutput(0);
    updateOutput(1);
    updateOutput(2);
    // Detect whether or not some events occur during the integration step
    _minConstraint = detectEvents();
    //
#ifdef DEBUG_MESSAGES
    std::cout << "========== siconos::simulation::EventDriven::advanceToEvent "
                 "============="
              << "\n";
    std::cout.precision(15);
    std::cout << "Istate: " << _istate << "\n";
    std::cout << "Maximum value of constraint functions: " << _minConstraint << "\n";
#endif
    //
    if (_istate != 2)  // some events occur
    {
      std::cout << "In siconos::simulation::EventDriven::advanceToEvent, some "
                   "events are detected!!!"
                << "\n";
      if (std::abs(_minConstraint) <
          tolerance_)  // events occur at the end of the integration step
      {
        isNewEventOccur = true;
      } else  // events need to be localized
      {
        isNewEventOccur = true;
        LocalizeFirstEvent();
      }
      // add new event to the list to be handled
      std::cout << "A new event occurs at time: " << _tout << "\n";
      _eventsManager->scheduleNonSmoothEvent(*this, _tout);
    }
  } else if (osiType == siconos::integrators::IntegratorType::LSODAROSI ||
             osiType == siconos::integrators::IntegratorType::HEM5OSI) {
    // WARNING: this is supposed to work for only one OSI, including all
    // the DS.  To be reviewed for multiple OSI case (if it has sense?).

    // ---> Step 1: integrrate the smooth dynamics from current event to
    // next event; Starting event = last accessed event.  Next event =
    // next time step or first root of the 'g' function found by
    // integrator (LsodarOSI)

    // if _istate == 1 => first call. It this case we suppose that _tinit
    // and _tend have been initialized before
    // if(_istate == 2 || _istate == 3)
    //  {
    //    _tinit = _eventsManager->startingTime();
    //    _tend =  _eventsManager->nextTime();
    //  }
    // _tout = _tend;
    // call integrate method for each OSI, between _tinit and _tend.
    for (auto it : *_allOSI) {
      it->resetAllNonSmoothParts();

      //====================================================================================
      //     std::cout << " Start of LsodarOSI integration" << "\n";
      it->integrate(_tinit, _tend, _tout, _istate);  // integrate must

      //  std::cout << " End of LsodarOSI integration" << "\n";
      // std::shared_ptr<siconos::integrators::LsodarOSI> lsodar =
      // std::static_pointer_cast<LsodarOSI>(it);
      // auto iwork = lsodar->getIwork();
      // auto rwork = lsodar->getRwork();
      //  std::cout << "Number of steps used: " << iwork[10] <<"\n";
      //  std::cout << "Method order last used: " << iwork[13] <<"\n";
      //  std::cout << "Step size last used: " << rwork[10] <<"\n";
      // return a flag (_istate) telling if _tend has been  reached or not.
      //====================================================================================

      if (_printStat) {
        statOut << " =================> Results after advanceToEvent "
                   "<================= "
                << "\n";
        statOut << " Starting time: " << _tinit << "\n";
        statOut << " _istate " << _istate << "\n";
      }
      if (_istate == 3)  // ie if _tout is not equal to _tend: one or more roots
                         // have been found.
      {
        DEBUG_PRINTF("An event has been found at time = %g", _tout);
        isNewEventOccur = true;
        // Add an event into the events manager list
        _eventsManager->scheduleNonSmoothEvent(*this, _tout);
        if (_printStat)
          statOut << " -----------> New non-smooth event at time " << _tout << "\n";
      }
    }
    // Set model time to _tout
    // update output[0], output[1], output[2]
    updateOutput(0);
    updateOutput(1);
    updateOutput(2);
    // update lambda[2], input[2] and indexSet[2] with double consitions for the
    // case there is no new event added during time integration, otherwise, this
    //  update is done when the new event is processed
    if (!isNewEventOccur) {
      if (!_allNSProblems->empty()) {
        // Solve LCP at acceleration level
        if (((*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->hasInteractions())) {
          auto indexSet2 = _nsds->topology()->indexSet(2);
          if (indexSet2->size() != 0) {
            (*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->compute(_tout);
            _nsds->updateInput(_tout, 2);
            // update indexSet[2] with double condition
            // updateIndexSetsWithDoubleCondition();
          }
        }
      }
    }
  } else {
    THROW_EXCEPTION(
        "In siconos::simulation::EventDriven::advanceToEvent, this type of "
        "OneStepIntegrator "
        "does not exist for Event-Driven scheme!!!");
  }
  DEBUG_END("siconos::simulation::EventDriven::advanceToEvent()\n");
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double siconos::simulation::EventDriven::computeResiduConstraints() {
  // Make sure that the state of all Dynamical Systems was updated
  double t = nextTime();  // time at the end of the step
  auto indexSet2 = _nsds->topology()->indexSet(2);
  double _y;
  // Loop over all interactions of indexSet2
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  double _maxResiduGap = 0.0;
  for (auto itosi : *_allOSI) {
    auto osi_NewMark = std::dynamic_pointer_cast<siconos::integrators::NewMarkAlphaOSI>(itosi);
    assert(osi_NewMark);

    bool _flag = osi_NewMark->handleVelocityConstraints();
    for (std::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui) {
      auto& inter = *indexSet2->bundle(*ui);
      if (!_flag)  // constraints at the position level
      {
        inter.computeOutput(t, 0);  // compute y[0] for the interaction at the end time
        _y = (*inter.y(0))(0);
      } else  // constraints at the velocity level
      {
        inter.computeOutput(t, 1);  // compute y[1] for the interaction at the end time
        _y = (*inter.y(1))(0);
      }

      if (_maxResiduGap < abs(_y)) {
        _maxResiduGap = abs(_y);
      }
      DEBUG_PRINTF("Constraint residu: =  %e \n", _y);
    }
  }

  DEBUG_PRINTF("Maximum constraint residu = %e \n", _maxResiduGap);
  return _maxResiduGap;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::simulation::EventDriven::prepareNewtonIteration() {
  DEBUG_BEGIN("siconos::simulation::EventDriven::prepareNewtonIteration()\n");
  // At this stage, we do
  // (1) compute iteration matrix W for all DSs belonging to all OSIs
  // (2) compute free residu for all DSs belonging to all OSIs and get maximum
  // residu (3) compute free state for all DSs belonging to all OSIs (4) compute
  // maximum gap residu over all interactions of indexSet 2
  _newtonResiduDSMax = 0.0;
  _newtonResiduYMax = 0.0;
  double _maxResidu;
  // Update input of level 2
  _nsds->updateInput(nextTime(), 2);
  // Loop over all OSIs

  for (auto itosi : *_allOSI) {
    auto osiType = itosi->getType();
    if (osiType != siconos::integrators::IntegratorType::NEWMARKALPHAOSI) {
      THROW_EXCEPTION(
          "In siconos::simulation::EventDriven::prepareNewtonIteration, the "
          "current OSI is "
          "not NewMarkAlpha scheme!!!");
    }
    // Compute iteration matrix W
    itosi->prepareNewtonIteration(nextTime());  // preparation for each OSI
    // Compute free residus, maximum residu
    _maxResidu = itosi->computeResidu();
    if (_newtonResiduDSMax < _maxResidu) {
      _newtonResiduDSMax = _maxResidu;
    }
    // Compute free state
    itosi->computeFreeState();
  }
  // Compute maximum gap residu
  _newtonResiduYMax = computeResiduConstraints();
  DEBUG_END("siconos::simulation::EventDriven::prepareNewtonIteration()\n");
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool siconos::simulation::EventDriven::newtonCheckConvergence(double _tol) {
  bool checkConvergence = true;
  if ((_newtonResiduDSMax > _tol) || (_newtonResiduYMax > _tol)) {
    checkConvergence = false;
  }
  return checkConvergence;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::simulation::EventDriven::predictionNewtonIteration() {
  // Prediction of the state for all Dynamical Systems before Newton iteration
  for (auto itosi : *_allOSI) {
    auto osi_NewMark = std::dynamic_pointer_cast<siconos::integrators::NewMarkAlphaOSI>(itosi);
    assert(osi_NewMark);
    osi_NewMark->prediction();
  }
  // Prediction of the output and lambda for all Interactions before Newton
  // iteration
  double t = nextTime();  // time at the end of the step
  // Loop over all interactions
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = _indexSet0->vertices(); ui != uiend; ++ui) {
    auto& inter = *_indexSet0->bundle(*ui);
    inter.computeOutput(t,
                        0);      // compute y[0] for the interaction at the end time
                                 // with the state predicted for Dynamical Systems
    inter.lambda(2)->setZero();  // reset lambda[2] to zero
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::simulation::EventDriven::correctionNewtonIteration() {
  // Update the input of level 2 for all Dynamical Systems after each iteration
  _nsds->updateInput(nextTime(), 2);
  // Correction
  for (auto itosi : *_allOSI) {
    auto osi_NewMark = std::dynamic_pointer_cast<siconos::integrators::NewMarkAlphaOSI>(itosi);
    assert(osi_NewMark);
    osi_NewMark->correction();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::simulation::EventDriven::newtonSolve(double criterion, unsigned int maxStep) {
  DEBUG_BEGIN(
      "siconos::simulation::EventDriven::newtonSolve(double criterion, "
      "unsigned int "
      "maxStep)\n");
  _isNewtonConverge = false;
  _newtonNbIterations = 0;  // number of Newton iterations
  int info = 0;
  _istate = 1;  // beginning of time integration
  // Prediction
  predictionNewtonIteration();
  while (1 != 0) {
    _newtonNbIterations++;
    // Prepare for iteration
    prepareNewtonIteration();
    // Check convergence
    _isNewtonConverge = newtonCheckConvergence(_newtonTolerance);
    //

    DEBUG_PRINTF("Iteration:  %i \n", _newtonNbIterations);
    DEBUG_PRINTF("Convergence:  %s \n", (_isNewtonConverge) ? "true" : "false");

    //
    if (_isNewtonConverge) {
      break;
    }
    if (_newtonNbIterations > maxStep) {
      std::cout << "Warning!!!In "
                   "siconos::simulation::EventDriven::newtonSolve: Number of "
                   "iterations is greater than the maximum value "
                << maxStep << "\n";
    }
    // If no convergence, proceed iteration
    auto indexSet2 = _nsds->topology()->indexSet(2);
    if (indexSet2->size() !=
        0)  // if indexSet2 is not empty, solve LCP to determine contact forces
    {
      info = computeOneStepNSProblem(SICONOS_OSNSP_ED_SMOOTH_POS);
      if (info != 0) {
        std::cout << "Warning!!!In siconos::simulation::EventDriven::newtonSolve: "
                     "LCP solver may "
                     "fail"
                  << "\n";
      }
    }
    // Correction of the state of all Dynamical Systems
    correctionNewtonIteration();
  }
  DEBUG_END(
      "siconos::simulation::EventDriven::newtonSolve(double criterion, "
      "unsigned int "
      "maxStep)\n");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double siconos::simulation::EventDriven::detectEvents(bool updateIstate) {
  DEBUG_BEGIN(
      "double siconos::simulation::EventDriven::detectEvents(bool "
      "updateIstate)\n")
  double _minResiduOutput = 0.0;  // maximum of g_i with i running over all
                                  // activated or deactivated contacts
  // Loop over all interactions to detect whether some constraints are activated
  // or deactivated
  bool _IsContactClosed = false;
  bool _IsContactOpened = false;
  bool _IsFirstTime = true;
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  std::shared_ptr<siconos::algebra::SiconosVector> y, ydot, lambda;
  auto topo = _nsds->topology();
  auto indexSet2 = topo->indexSet(2);
  for (std::tie(ui, uiend) = _indexSet0->vertices(); ui != uiend; ++ui) {
    auto inter = _indexSet0->bundle(*ui);
    auto nsLawSize = inter->nonSmoothLaw()->size();
    if (nsLawSize != 1) {
      THROW_EXCEPTION(
          "In siconos::simulation::EventDriven::detectEvents, the interaction "
          "size > 1 has "
          "not been implemented yet!!!");
    }
    y = inter->y(0);                     // output y at this Interaction
    ydot = inter->y(1);                  // output of level 1 at this Interaction
    lambda = inter->lambda(2);           // input of level 2 at this Interaction
    if (!(indexSet2->is_vertex(inter)))  // if Interaction is not in the indexSet[2]
    {
      if ((*y)(0) < tolerance_)  // gap at the current interaction <= 0
      {
        _IsContactClosed = true;
      }

      if (_IsFirstTime) {
        _minResiduOutput = (*y)(0);
        _IsFirstTime = false;
      } else {
        if (_minResiduOutput > (*y)(0)) {
          _minResiduOutput = (*y)(0);
        }
      }
    } else  // If interaction is in the indexSet[2]
    {
      if ((*lambda)(0) < tolerance_)  // normal force at the current interaction <= 0
      {
        _IsContactOpened = true;
      }

      if (_IsFirstTime) {
        _minResiduOutput = (*lambda)(0);
        _IsFirstTime = false;
      } else {
        if (_minResiduOutput > (*lambda)(0)) {
          _minResiduOutput = (*lambda)(0);
        }
      }
    }
    //
    DEBUG_EXPR_WE(std::cout.precision(15);
                  std::cout << "Contact number: " << inter->number() << "\n";
                  std::cout << "Contact gap: " << (*y)(0) << "\n";
                  std::cout << "Contact force: " << (*lambda)(0) << "\n";
                  std::cout << "Is contact is closed: " << _IsContactClosed << "\n";
                  std::cout << "Is contact is opened: " << _IsContactOpened << "\n";);
    //
  }
  //
  if (updateIstate) {
    if ((!_IsContactClosed) && (!_IsContactOpened)) {
      _istate = 2;  // no event is detected
    } else if ((_IsContactClosed) && (!_IsContactOpened)) {
      _istate = 3;  // Only some contacts are closed
    } else if ((!_IsContactClosed) && (_IsContactOpened)) {
      _istate = 4;  // Only some contacts are opened
    } else {
      _istate = 5;  // Some contacts are closed AND some contacts are opened
    }
  }
  //
  DEBUG_END(
      "double siconos::simulation::EventDriven::detectEvents(bool "
      "updateIstate)\n")
  return _minResiduOutput;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void siconos::simulation::EventDriven::LocalizeFirstEvent() {
  // We localize the first event occuring during the integration step when the
  // flag _istate = 3 or 4 Compute the coefficients of the dense output
  // polynomial for all DSs
  for (auto itosi : *_allOSI) {
    auto osi_NewMark = std::dynamic_pointer_cast<siconos::integrators::NewMarkAlphaOSI>(itosi);
    assert(osi_NewMark);
    osi_NewMark->prepareEventLocalization();
  }
  //
  double t_a = startingTime();
  double t_b = nextTime();
  double _minConstraint = 0.0;
  bool found = false;
  bool _IsupdateIstate = false;
  unsigned int _numIter = 0;
  while (!found) {
    _numIter++;
    double t_i = (t_b + t_a) / 2.0;  // mid-time of the current interval
    // set t_i as the current time
    // Generate dense output for all DSs at the time t_i
    for (auto itosi : *_allOSI) {
      auto osi_NewMark =
          std::dynamic_pointer_cast<siconos::integrators::NewMarkAlphaOSI>(itosi);
      assert(osi_NewMark);
      osi_NewMark->DenseOutputallDSs(t_i);
    }
    // If _istate = 3 or 5, i.e. some contacts are closed, we need to compute
    // y[0] for all interactions
    if ((_istate == 3) || (_istate == 5))  // some contacts are closed
    {
      _nsds->updateOutput(t_i, 0);
    }
    // If _istate = 4 or 5, i.e. some contacts are detached, we need to solve
    // LCP at the acceleration level to compute contact forces
    if ((_istate == 4) || (_istate == 5))  // some contacts are opened
    {
      if (!_allNSProblems->empty()) {
        (*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->compute(t_i);
      }
    }
    // Check whether or not some events occur in the interval [t_a, t_i]
    _minConstraint = detectEvents(_IsupdateIstate);
    if (std::abs(_minConstraint) < tolerance_)  // first event is found
    {
      _tout = t_i;
      found = true;
    }
    // if some events are detected in the interval [t_a, t_i] (if _istate != 2),
    // set t_b = t_i
    if (_minConstraint < -tolerance_) {
      t_b = t_i;
    } else  // if no event is detected in [t_a, t_i], then we have to detect
            // events in the interval [t_i, t_b]
    {
      t_a = t_i;
    }
    //
    if (_numIter > _localizeEventMaxIter) {
      THROW_EXCEPTION(
          "In siconos::simulation::EventDriven::LocalizeFirstEvent, the "
          "numbner of iterations "
          "performed is too large!!!");
    }
  }
}
