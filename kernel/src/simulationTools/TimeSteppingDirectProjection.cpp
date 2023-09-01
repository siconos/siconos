/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "TimeSteppingDirectProjection.hpp"

#include "Interaction.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "MoreauJeanDirectProjectionOSI.hpp"
#include "MoreauJeanOSI.hpp"
#include "NewtonEuler1DR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NonSmoothLaw.hpp"
#include "OneStepIntegrator.hpp"
#include "OneStepNSProblem.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "TimeStepping.hpp"
#include "Tools.hpp"
#include "Topology.hpp"

static siconos::simulation::CheckSolverFPtr
    checkSolverOutputProjectOnConstraints = nullptr;
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"
// #define CORRECTIONSVELOCITIES

siconos::simulation::TimeSteppingDirectProjection::TimeSteppingDirectProjection(
    std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
    std::shared_ptr<TimeDiscretisation> td,
    std::shared_ptr<siconos::integrators::OneStepIntegrator> osi,
    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>
        osnspb_velo,
    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>
        osnspb_pos,
    unsigned int level)
    : TimeStepping{nsds, td, osi, osnspb_velo},
      _indexSetLevelForProjection{level} {
  if (not std::dynamic_pointer_cast<
          siconos::integrators::MoreauJeanDirectProjectionOSI>(osi))
    THROW_EXCEPTION(
        "siconos::simulation::TimeSteppingDirectProjection::"
        "TimeSteppingDirectProjection.  "
        "wrong type of OneStepIntegrator");

  (*_allNSProblems).resize(SICONOS_NB_OSNSP_TSP);
  insertNonSmoothProblem(osnspb_pos, siconos::simulation::SICONOS_OSNSP_TS_POS);
}

void siconos::simulation::TimeSteppingDirectProjection::initOSNS() {
  TimeStepping::initOSNS();

  (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_POS]
      ->setIndexSetLevel(_indexSetLevelForProjection);
  (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_POS]
      ->setInputOutputLevel(0);

  (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]
      ->setIndexSetLevel(1);
  (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]
      ->setInputOutputLevel(1);
}

void siconos::simulation::TimeSteppingDirectProjection::nextStep() {
  TimeStepping::nextStep();

  // Zeroing Lambda Muliplier of indexSet()

  auto indexSet = _nsds->topology()->indexSet(0);
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    auto inter = indexSet->bundle(*ui);
    inter->lambda(0)->zero();
  }
}

void siconos::simulation::TimeSteppingDirectProjection::advanceToEvent() {
  initialize();

  /** First step, Solve the standard velocity formulation.*/

  DEBUG_BEGIN("TimeStepping::newtonSolve\n");

  if (!_doOnlyProj)
    TimeStepping::newtonSolve(_newtonTolerance, _newtonMaxIteration);
  else
    updateInteractions();

  DEBUG_EXPR_WE(
      std::cout << "TimeStepping::newtonSolve end : Number of iterations="
                << getNewtonNbIterations() << "\n";
      std::cout << "                              : newtonResiduDSMax="
                << newtonResiduDSMax() << "\n";
      std::cout << "                              : newtonResiduYMax="
                << newtonResiduYMax() << "\n";
      std::cout << "                              : newtonResiduRMax="
                << newtonResiduRMax() << "\n";);

  if (!_doProj) return;
  int info = 0;

  /** Second step, Perform the projection on constraints.*/

  DEBUG_PRINT(
      "siconos::simulation::TimeSteppingDirectProjection::newtonSolve begin "
      "projection:\n");

  auto dsGraph = _nsds->dynamicalSystems();

#ifdef TSPROJ_CORRECTIONVELOCITIES
  for (auto& vi : *dsGraph) {
    auto ds = dsGraph->bundle(vi);
    if (auto neds =
            std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds))
      *(neds->deltaq()) = *(neds->q());
    else
      THROW_EXCEPTION("TS:: - ds is not from NewtonEulerDS.");
  }
#endif

  bool runningProjection = false;
  _nbProjectionIteration = 0;
  // for (InteractionsIterator it = allInteractions->begin(); it !=
  // allInteractions->end(); it++){
  //   double criteria = (*it)->relation()->y(0)->getValue(0);
  //   if (Type::value(*((*it)->nonSmoothLaw())) ==
  //   siconos::modeling::Type::NewtonImpactFrictionNSL ||
  //  Type::value(*((*it)->nonSmoothLaw())) ==
  //  siconos::modeling::Type::NewtonImpactNSL){
  //     auto ri = std::static_pointer_cast<NewtonEuler1DR> ((*it)->relation());
  //     if (criteria < -1e-7){
  //  ri->_isOnContact=true;
  //     }else{
  //  ri->_isOnContact=false;
  //     }
  //   }
  //   if (criteria < -_constraintTol)
  //     runningNewton=true;
  // }
  if (_nsds->topology()->numberOfIndexSet() > _indexSetLevelForProjection)
    computeCriteria(&runningProjection);
  // Zeroing Lambda Muliplier of indexSet()

  auto indexSet = _nsds->topology()->indexSet(0);
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    auto inter = indexSet->bundle(*ui);
    inter->lambda(0)->zero();
  }
  _nsds->updateInput(nextTime(), 0);

  // Store the q vector of each DS.

  for (auto aVi2 = dsGraph->begin(); aVi2 != dsGraph->end(); ++aVi2) {
    auto ds = dsGraph->bundle(*aVi2);
    auto& workVectors = *dsGraph->properties(*aVi2).workVectors;
    if (auto neds =
            std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      *workVectors[siconos::integrators::MoreauJeanOSI::QTMP] = *neds->q();
    } else if (auto d =
                   std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(
                       ds)) {
      *workVectors[siconos::integrators::MoreauJeanOSI::QTMP] = *d->q();
    } else
      THROW_EXCEPTION(
          "siconos::simulation::TimeSteppingDirectProjection::advanceToEvent() "
          ":: - Ds is not "
          "from NewtonEulerDS neither from LagrangianDS.");
  }

  while (runningProjection &&
         _nbProjectionIteration < _projectionMaxIteration) {
    _nbProjectionIteration++;
    DEBUG_PRINTF("TimeSteppingDirectProjection projection step = %d\n",
                 _nbProjectionIteration);

    auto indexSet = _nsds->topology()->indexSet(0);
    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
      auto inter = indexSet->bundle(*ui);
      inter->lambda(0)->zero();
    }
    _nsds->updateInput(nextTime(), 0);
    info = 0;

    DEBUG_PRINT("TimeSteppingProjectOnConstraint compute OSNSP.\n");

    info = computeOneStepNSProblem(siconos::simulation::SICONOS_OSNSP_TS_POS);

    DEBUG_PRINTF("IndexSet0->size() = %i\n",
                 (int)_nsds->topology()->indexSet(0)->size());
    DEBUG_PRINTF("IndexSet1->size() = %i\n",
                 (int)_nsds->topology()->indexSet(1)->size());
    DEBUG_EXPR(
        oneStepNSProblem(siconos::simulation::SICONOS_OSNSP_TS_POS)->display());

    if (info && _newtonWarningOnNonConvergence) {
      std::cout << "[kernel] "
                   "siconos::simulation::TimeSteppingDirectProjection::"
                   "advanceToEvent() project on constraints. solver failed.\n";
    }
    _nsds->updateInput(nextTime(), 0);

    DEBUG_EXPR_WE(
        std ::cout << "After update input\n";
        auto indexSet1 = _nsds->topology()->indexSet(1);
        std ::cout << "lamda(1) in IndexSet1\n";
        for (std::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet1->bundle(*ui);
          inter->lambda(1)->display();
        } auto indexSet0 = _nsds->topology()->indexSet(0);
        std ::cout << "lamda(0) in indexSet0\n";
        for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet0->bundle(*ui);
          inter->lambda(0)->display();
        });

    // This part should be in
    // MoreauJeanOSIProjectOnConstraintsOS::updateState(level =0)
    for (auto aVi2 : *dsGraph) {
      auto ds = dsGraph->bundle(aVi2);
      auto& workVectors = *dsGraph->properties(aVi2).workVectors;

      if (auto neds =
              std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
        auto q = neds->q();
        auto qtmp = workVectors[siconos::integrators::MoreauJeanOSI::QTMP];

        DEBUG_EXPR_WE(std ::cout << "qtmp before  update \n"; qtmp->display();
                      std ::cout << "p(0) before  update \n";
                      neds->p(0)->display(););

        if (neds->p(0)) {
          //*q = * qtmp +  *neds->p(0);
          *q += *neds->p(
              0);  // Why it works like that and not with the previous line ?
        }

        DEBUG_EXPR_WE(std ::cout << "q after  update \n"; q->display(););

        neds->normalizeq();
        neds->computeT();
      } else if (auto d =
                     std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(
                         ds)) {
        auto q = d->q();
        auto qtmp = workVectors[siconos::integrators::MoreauJeanOSI::QTMP];

        if (d->p(0)) {
          //*q = * qtmp +  *d->p(0);
          *q += *d->p(0);
        }
      } else
        THROW_EXCEPTION(
            "siconos::simulation::TimeSteppingDirectProjection::advanceToEvent("
            ") :: - Ds is "
            "not from NewtonEulerDS neither from LagrangianDS.");
    }

    computeCriteria(&runningProjection);

    // cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    // Z:"<<endl;
    //(*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_POS]->display();
    //(std::static_pointer_cast<LinearOSNS>((*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_POS]))->z()->display();

    DEBUG_EXPR_WE(
        std::cout << "siconos::simulation::TimeSteppingDirectProjection::"
                     "Projection end : "
                     "Number of iterations="
                  << _nbProjectionIteration << "\n";
        std ::cout << "After update state in position\n";
        std ::cout << "lamda(1) in IndexSet1\n";
        auto indexSet1 = _nsds->topology()->indexSet(1);
        auto indexSet0 = _nsds->topology()->indexSet(0);

        for (std::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet1->bundle(*ui);
          inter->lambda(1)->display();
        }

        std ::cout
        << "lamda(0) in indexSet0\n";
        for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet0->bundle(*ui);
          inter->lambda(0)->display();
        } std ::cout
        << "y(1) in IndexSet1\n";
        for (std::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet1->bundle(*ui);
          inter->y(1)->display();
        } std ::cout
        << "y(0) in indexSet0\n";
        for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet0->bundle(*ui);
          inter->y(0)->display();
        });

    // cout<<"during projection before normalizing of q:\n";
    // for (InteractionsIterator it = allInteractions->begin(); it !=
    // allInteractions->end(); it++)
    //{
    //   (*it)->relation()->computeh(getTkp1());
    // }
  }  // end while(runningProjection && _nbProjectionIteration <
     // _projectionMaxIteration)

  // We update forces to start the Newton Loop the next tiem step with a correct
  // value in swap
  for (auto aVi2 : *dsGraph) {
    auto ds = dsGraph->bundle(aVi2);
    if (auto neds =
            std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      auto time = nextTime();
      neds->computeForces(time, neds->q(), neds->twist());
    } else if (auto d =
                   std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(
                       ds)) {
      auto time = nextTime();
      d->computeForces(time, d->q(), d->velocity());
    } else if (std::dynamic_pointer_cast<
                   siconos::modeling::LagrangianLinearTIDS>(ds)) {
    } else
      THROW_EXCEPTION(
          "TimeSteppingCombinedProjection::advanceToEvent() - Ds is not from "
          "NewtonEulerDS "
          "neither from LagrangianDS.");
  }

  if (_nbProjectionIteration == _projectionMaxIteration &&
      _newtonWarningOnNonConvergence) {
    std::cout << "[kernel] "
                 "siconos::simulation::TimeSteppingDirectProjection::"
                 "advanceToEvent() Max "
                 "number of projection iterations reached ("
              << _nbProjectionIteration << ")" << std::endl;
    printf("[kernel]                max criteria equality =  %e.\n",
           _maxViolationEquality);
    printf("[kernel]                max criteria unilateral =  %e.\n",
           _maxViolationUnilateral);
  }

  DEBUG_END(
      "siconos::simulation::TimeSteppingDirectProjection::newtonSolve()\n");

  return;
  // #ifdef TSPROJ_CORRECTIONVELOCITIES
  //    /*The following reduces the velocity because the position step increase
  //    the energy of the system. This formulation works only with simple
  //    systems.To activate it, comment the next line.*/

  //   for(DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi !=
  //   dsGraph->end(); ++vi)
  //   {
  //     auto ds = dsGraph->bundle(*vi);
  //     Type::Siconos dsType = Type::value(*ds);
  //     if(dsType !=Type::NewtonEulerDS)
  //       THROW_EXCEPTION("TS:: - ds is not from NewtonEulerDS.");
  //     // auto  dotq = neds->dotq();
  //     // auto  q = neds->q();
  //     // auto  qold = neds->qMemory()->getSiconosVector(0);
  //     // double h = timeStep();
  //     // dotq->setValue(3,(q->getValue(3)-qold->getValue(3))/h);
  //     // dotq->setValue(4,(q->getValue(4)-qold->getValue(4))/h);
  //     // dotq->setValue(5,(q->getValue(5)-qold->getValue(5))/h);
  //     // dotq->setValue(6,(q->getValue(6)-qold->getValue(6))/h);

  //     /*compute the new velocity seeing the work of fext*/
  //     auto neds = std::static_pointer_cast<NewtonEulerDS>(ds);
  //     *(neds->deltaq())-=*(neds->q());
  //     DEBUG_EXPR(printf("TSProj NewtonSolve :deltaq:");
  //     (neds->deltaq())->display(););
  //     //continue;
  //     double  n2q=neds->deltaq()->norm2();
  //     double n2=0.0;
  //     if(neds->fExt())
  //       n2=neds->fExt()->norm2();
  //     if(n2 > 1e-7 && n2q > 1e-14)
  //     {
  //       //if (n2q < 1e-14)
  //       // continue;

  //       auto  FextNorm =
  //       std::make_shared<siconos::algebra::SiconosVector>(3));
  //       FextNorm->setValue(0,neds->fExt()->getValue(0));
  //       FextNorm->setValue(1,neds->fExt()->getValue(1));
  //       FextNorm->setValue(2,neds->fExt()->getValue(2));
  // DEBUG_EXPR_WE(
  //       std::cout<<"siconos::simulation::TimeSteppingDirectProjection::newtonSolve
  //       deltaQ
  //       :\n"; neds->deltaq()->display();
  //       std::cout<<"siconos::simulation::TimeSteppingDirectProjection::newtonSolve
  //       Fext
  //       :\n"; FextNorm->display();
  //       );

  //       (*FextNorm)*=(1./n2);
  //       /*work of external forces.*/
  //       double workFext= neds->fExt()->getValue(0) *
  //       neds->deltaq()->getValue(0)+
  //                        neds->fExt()->getValue(1) *
  //                        neds->deltaq()->getValue(1)+
  //                        neds->fExt()->getValue(2) *
  //                        neds->deltaq()->getValue(2);
  //       //workFext*=2.0;
  //       double VkFNorm=FextNorm->getValue(0)*neds->velocity()->getValue(0)+
  //                      FextNorm->getValue(1)*neds->velocity()->getValue(1)+
  //                      FextNorm->getValue(2)*neds->velocity()->getValue(2);
  //       double VkcFNorm =VkFNorm;
  //       VkcFNorm= VkFNorm*VkFNorm - 2*fabs(workFext)/(neds->massValue());
  //       if(VkcFNorm >0)
  //       {
  //         if(VkFNorm>0)
  //           VkcFNorm=sqrt(VkcFNorm);
  //         else
  //           VkcFNorm=-sqrt(VkcFNorm);
  //       }
  //       else
  //         VkcFNorm=0;
  //       // if (VkFNorm >= 0 && workFext >0){
  //       //   ;//VkcFNorm=sqrt
  //       (2*workFext/(neds->massValue())+VkFNorm*VkFNorm);
  //       // }else if (VkFNorm <= 0 && workFext < 0){
  //       //   ;//VkcFNorm=-sqrt
  //       (fabs(2*workFext/(neds->massValue())+VkFNorm*VkFNorm));
  //       // }else if (VkFNorm > 0 && workFext <0){
  //       //   VkcFNorm= VkFNorm*VkFNorm + 2*workFext/(neds->massValue());
  //       //   if (VkcFNorm >0)
  //       //     VkcFNorm=sqrt(VkcFNorm);
  //       //   else
  //       //     VkcFNorm=0;
  //       // }else if (VkFNorm < 0 && workFext > 0){
  //       //   VkcFNorm= VkFNorm*VkFNorm - 2*workFext/(neds->massValue());
  //       //   if (VkcFNorm >0)
  //       //     VkcFNorm=-sqrt(VkcFNorm);
  //       //   else
  //       //     VkcFNorm=0;
  //       // }
  // DEBUG_EXPR_WE(
  //       printf("siconos::simulation::TimeSteppingDirectProjection::newtonSolve
  //       velocity before update(prevComp=%e, newComp=%e)\n",VkFNorm,VkcFNorm);
  //       printf("VELOCITY1 "); neds->velocity()->display();
  // );
  //       neds->velocity()->setValue(0,neds->velocity()->getValue(0)+(VkcFNorm
  //       - VkFNorm)*FextNorm->getValue(0));
  //       neds->velocity()->setValue(1,neds->velocity()->getValue(1)+(VkcFNorm
  //       - VkFNorm)*FextNorm->getValue(1));
  //       neds->velocity()->setValue(2,neds->velocity()->getValue(2)+(VkcFNorm
  //       - VkFNorm)*FextNorm->getValue(2));
  // DEBUG_EXPR_WE(
  //       std::cout<<"siconos::simulation::TimeSteppingDirectProjection::newtonSolve
  //       velocity updated\n"; printf("VELOCITY2 ");
  //       neds->velocity()->display();
  // )
  //     }
  //     auto T = neds->T();
  //     auto  dotq = neds->dotq();
  //     siconos::algebra::prod(*T,*neds->velocity(),*dotq,true);
  //     if(!_allNSProblems->empty())
  //     {
  //       for(unsigned int level = _levelMinForOutput;
  //           level < _levelMaxForOutput;
  //           level++)
  //         updateOutput(level);
  //     }
  //   }
  // #endif
}

void siconos::simulation::TimeSteppingDirectProjection::computeCriteria(
    bool* runningProjection) {
  auto indexSet = _nsds->topology()->indexSet(_indexSetLevelForProjection);
  siconos::graphs::InteractionsGraph::VIterator aVi, viend;

  double maxViolationEquality = -1e24;
  double minViolationEquality = +1e24;
  double maxViolationUnilateral = -1e24;
  // double minViolationUnilateral = +1e24;

  *runningProjection = false;

  for (std::tie(aVi, viend) = indexSet->vertices(); aVi != viend; ++aVi) {
    auto inter = indexSet->bundle(*aVi);
    inter->computeOutput(getTkp1(), 0);
    inter->relation()->computeJach(getTkp1(), *inter);

    if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
            siconos::modeling::Type::NewtonImpactFrictionNSL ||
        siconos::types::type_value(*(inter->nonSmoothLaw())) ==
            siconos::modeling::Type::NewtonImpactNSL) {
      double criteria = std::max(0.0, -inter->y(0)->getValue(0));
      DEBUG_PRINTF("Unilateral inter->y(0)->getValue(0) %e.\n",
                   inter->y(0)->getValue(0));
      if (criteria > maxViolationUnilateral) maxViolationUnilateral = criteria;
      // if (criteria < minViolationUnilateral) minViolationUnilateral=criteria;
      if (maxViolationUnilateral > _constraintTolUnilateral) {
        *runningProjection = true;

        DEBUG_PRINTF("TSProj newton criteria unilateral true %e.\n", criteria);
      }
    } else {
      DEBUG_PRINTF("Equality inter->y(0)->normInf() %e.\n",
                   inter->y(0)->normInf());
      if (inter->y(0)->normInf() > maxViolationEquality)
        maxViolationEquality = inter->y(0)->normInf();
      if (inter->y(0)->normInf() < minViolationEquality)
        minViolationEquality = inter->y(0)->normInf();
      if (inter->y(0)->normInf() > _constraintTol) {
        *runningProjection = true;
        DEBUG_PRINTF("TSProj  newton criteria equality true %e.\n",
                     inter->y(0)->normInf());
      }
    }
    _maxViolationUnilateral = maxViolationUnilateral;
    _maxViolationEquality = maxViolationEquality;
  }

  DEBUG_PRINT("TSProj newton min/max criteria projection\n");
  DEBUG_EXPR(std::cout << "             runningProjection " << std::boolalpha
                       << *runningProjection << std::endl;);
  DEBUG_PRINTF("              min criteria equality =  %e.\n",
               minViolationEquality);
  DEBUG_PRINTF("              max criteria equality =  %e.\n",
               maxViolationEquality);
  DEBUG_PRINTF("              max criteria unilateral =  %e.\n",
               maxViolationUnilateral);
  // DEBUG_PRINTF("              min criteria unilateral =
  // %e.\n",minViolationUnilateral);
}

void siconos::simulation::TimeSteppingDirectProjection::newtonSolve(
    double criterion, unsigned int maxStep) {
  bool isNewtonConverge = false;
  _newtonNbIterations = 0;  // number of Newton iterations
  int info = 0;
  // cout<<"||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||
  // BEGIN NEWTON IT
  // "<<endl;
  bool isLinear = _nsds->isLinear();
  auto indexSet = _nsds->topology()->indexSet(0);
  initializeNewtonLoop();

  if ((_newtonOptions == TimeSteppingType::LINEAR ||
       _newtonOptions == TimeSteppingType::LINEAR_IMPLICIT) ||
      isLinear) {
    _newtonNbIterations++;
    prepareNewtonIteration();
    computeFreeState();
    // updateOutput(0);
    // updateIndexSets();
    if (!_allNSProblems->empty() && indexSet->size() > 0)
      info = computeOneStepNSProblem(
          siconos::simulation::SICONOS_OSNSP_TS_VELOCITY);
    // Check output from solver (convergence or not ...)
    if (!checkSolverOutputProjectOnConstraints)
      DefaultCheckSolverOutput(info);
    else
      checkSolverOutputProjectOnConstraints(info, this);

    update();

    // isNewtonConverge = newtonCheckConvergence(criterion);
  }

  else if (_newtonOptions == TimeSteppingType::NONLINEAR) {
    while ((!isNewtonConverge) && (_newtonNbIterations < maxStep) && (!info)) {
      _newtonNbIterations++;
      prepareNewtonIteration();
      computeFreeState();
      // updateOutput(0);
      // updateIndexSets();

      // if there is not any Interaction at
      // the beginning of the simulation _allNSProblems may not be
      // empty here (check with SpaceFilter and one disk not on
      // the ground : MultiBodyTest::t2)

      // if((*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->simulation())
      // is also relevant here.
      if (!_allNSProblems->empty() && indexSet->size() > 0) {
        info = computeOneStepNSProblem(
            siconos::simulation::SICONOS_OSNSP_TS_VELOCITY);
      }
      // Check output from solver (convergence or not ...)
      if (!checkSolverOutputProjectOnConstraints)
        DefaultCheckSolverOutput(info);
      else
        checkSolverOutputProjectOnConstraints(info, this);

      updateAllInput();
      updateState();
      isNewtonConverge = newtonCheckConvergence(criterion);
      if (!isNewtonConverge && !info) {
        updateOutput();
      }
    }
    if (!isNewtonConverge) {
      if (_newtonWarningOnNonConvergence) {
        std::cout << "TimeStepping::newtonSolve -- Newton process stopped: "
                     "max. number of steps ("
                  << maxStep << ") reached." << std::endl;
      }
    } else if (info && _newtonWarningOnNonConvergence) {
      std::cout << "TimeStepping::newtonSolve -- Newton process stopped: "
                   "solver failed."
                << std::endl;
    }
    //    else
    //      std::cout << "TimeStepping::newtonSolve succed
    //      nbit="<<_newtonNbIterations<<"maxStep="<<maxStep<<endl;
  } else
    THROW_EXCEPTION(
        "TimeStepping::NewtonSolve failed. Unknown newtonOptions: " +
        siconos::tools::enum_to_string(_newtonOptions));
}
