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

#include "TimeSteppingDirectProjection.hpp"

#include "LagrangianLinearTIDS.hpp"
#include "LagrangianSparseLinearTIDS.hpp"
#include "MoreauJeanDirectProjectionOSI.hpp"
#include "MoreauJeanOSI.hpp"
#include "NewtonEulerDS.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OneStepIntegrator.hpp"
#include "OneStepNSProblem.hpp"
#include "Relation.hpp"  // IWYU pragma: keep
#include "SiconosException.hpp"
#include "TimeStepping.hpp"
#include "Tools.hpp"

static siconos::simulation::CheckSolverFPtr checkSolverOutputProjectOnConstraints = nullptr;
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"
// #define CORRECTIONSVELOCITIES

siconos::simulation::TimeSteppingDirectProjection::TimeSteppingDirectProjection(
    std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
    std::shared_ptr<TimeDiscretisation> td,
    std::shared_ptr<siconos::integrators::OneStepIntegrator> osi,
    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb_velo,
    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osnspb_pos,
    unsigned int level)
    : TimeStepping{nsds, td, osi, osnspb_velo}, _indexSetLevelForProjection{level} {
  if (not std::dynamic_pointer_cast<siconos::integrators::MoreauJeanDirectProjectionOSI>(osi))
    THROW_EXCEPTION(
        "siconos::simulation::TimeSteppingDirectProjection::"
        "TimeSteppingDirectProjection.  "
        "wrong type of OneStepIntegrator");

  (*_allNSProblems).resize(SICONOS_NB_OSNSP_TSP);
  insertNonSmoothProblem(osnspb_pos, siconos::simulation::SICONOS_OSNSP_TS_POS);
}

void siconos::simulation::TimeSteppingDirectProjection::initializeOneStepNSProblem() {
  updateIndexSets();
  TimeStepping::initializeOneStepNSProblem();

  (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_POS]->setIndexSetLevel(
      _indexSetLevelForProjection);
  (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_POS]->setInputOutputLevel(0);

  (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->setIndexSetLevel(1);
  (*_allNSProblems)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]->setInputOutputLevel(1);
}

void siconos::simulation::TimeSteppingDirectProjection::nextStep() {
  TimeStepping::nextStep();

  // Zeroing Lambda Muliplier of indexSet()

  auto indexSet = nonSmoothDynamicalSystem()->topology()->indexSet(0);
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    auto inter = indexSet->bundle(*ui);
    inter->lambda(0)->setZero();
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

  DEBUG_EXPR_WE(std::cout << "TimeStepping::newtonSolve end : Number of iterations="
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
  auto nsds = nonSmoothDynamicalSystem();
  auto dsGraph = nsds->dynamicalSystems();

#ifdef TSPROJ_CORRECTIONVELOCITIES
  for (auto& vi : *dsGraph) {
    auto ds = dsGraph->bundle(vi);
    if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds))
      *(neds->deltaq()) = *(neds->q());
    else
      THROW_EXCEPTION("TS:: - ds is not from NewtonEulerDS.");
  }
#endif

  bool runningProjection = false;
  _nbProjectionIteration = 0;
  // for (InteractionsIterator it = allInteractions->begin(); it !=
  // allInteractions->end(); it++){
  //   double criteria = (*(*it)->relation()->y(0))(0);
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
  if (nsds->topology()->numberOfIndexSet() > _indexSetLevelForProjection)
    computeCriteria(&runningProjection);
  // Zeroing Lambda Muliplier of indexSet()

  auto indexSet = nsds->topology()->indexSet(0);
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    auto inter = indexSet->bundle(*ui);
    inter->lambda(0)->setZero();
  }
  nsds->updateInput(nextTime(), 0);

  // Store the q vector of each DS.

  // for (auto aVi2 = dsGraph->begin(); aVi2 != dsGraph->end(); ++aVi2) {
  // auto ds = dsGraph->bundle(*aVi2);
  //  auto& workVectors = *dsGraph->properties(*aVi2).workVectors;
  //  if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
  //    *workVectors[siconos::integrators::MoreauJeanOSI::QTMP] = *neds->q();
  //  } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
  //    *workVectors[siconos::integrators::MoreauJeanOSI::QTMP] = *d->q();
  //  } else if (auto d =
  //  std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseDS>(ds)) {
  //    *workVectors[siconos::integrators::MoreauJeanOSI::QTMP] = *d->q();
  //  } else
  //    THROW_EXCEPTION(
  //        "siconos::simulation::TimeSteppingDirectProjection::advanceToEvent() "
  //        ":: - Ds is not "
  //        "from NewtonEulerDS neither from LagrangianDS.");
  //}

  while (runningProjection && _nbProjectionIteration < _projectionMaxIteration) {
    _nbProjectionIteration++;
    DEBUG_PRINTF("TimeSteppingDirectProjection projection step = %d\n",
                 _nbProjectionIteration);

    auto indexSet = nsds->topology()->indexSet(0);
    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
      auto inter = indexSet->bundle(*ui);
      inter->lambda(0)->setZero();
    }
    nsds->updateInput(nextTime(), 0);
    info = 0;

    DEBUG_PRINT("TimeSteppingProjectOnConstraint compute OSNSP.\n");

    info = computeOneStepNSProblem(siconos::simulation::SICONOS_OSNSP_TS_POS);

    DEBUG_PRINTF("IndexSet0->size() = %i\n", (int)nsds->topology()->indexSet(0)->size());
    DEBUG_PRINTF("IndexSet1->size() = %i\n", (int)nsds->topology()->indexSet(1)->size());
    DEBUG_EXPR(oneStepNSProblem(siconos::simulation::SICONOS_OSNSP_TS_POS)->display());

    if (info && _newtonWarningOnNonConvergence) {
      std::cout << "[kernel] "
                   "siconos::simulation::TimeSteppingDirectProjection::"
                   "advanceToEvent() project on constraints. solver failed.\n";
    }
    nsds->updateInput(nextTime(), 0);

    DEBUG_EXPR_WE(
        std ::cout << "After update input\n"; auto indexSet1 = nsds->topology()->indexSet(1);
        std ::cout << "lamda(1) in IndexSet1\n";
        for (std::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet1->bundle(*ui);
          siconos::algebra::print(*inter->lambda(1));
        } auto indexSet0 = nsds->topology()->indexSet(0);
        std ::cout << "lamda(0) in indexSet0\n";
        for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet0->bundle(*ui);
          siconos::algebra::print(*inter->lambda(0));
        });

    // This part should be in
    // MoreauJeanOSIProjectOnConstraintsOS::updateState(level =0)
    for (auto aVi2 : *dsGraph) {
      auto ds = dsGraph->bundle(aVi2);
      // auto& workVectors = *dsGraph->properties(aVi2).workVectors;

      if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
        // auto qtmp = workVectors[siconos::integrators::MoreauJeanOSI::QTMP];

        // DEBUG_EXPR_WE(std ::cout << "qtmp before  update \n";
        // siconos::algebra::print(*qtmp);
        //               std ::cout << "p(0) before  update \n";
        //               siconos::algebra::print(*neds->p(0)););

        if (neds->p(0)) {
          //*q = * qtmp +  *neds->p(0);
          *neds->q() +=
              neds->p_read(0);  // Why it works like that and not with the previous line ?
        }

        DEBUG_EXPR_WE(std ::cout << "q after  update \n"; siconos::algebra::print(*q););

        neds->normalizeq();
        neds->computeT(neds->q_read());
      } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
        // auto qtmp = workVectors[siconos::integrators::MoreauJeanOSI::QTMP];

        if (d->p(0)) {
          //*q = * qtmp +  *d->p(0);
          *d->q() += d->p_read(0);
        }
      } else if (auto d =
                     std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseDS>(ds)) {
        // auto qtmp = workVectors[siconos::integrators::MoreauJeanOSI::QTMP];

        if (d->p(0)) {
          //*q = * qtmp +  *d->p(0);
          *d->q() += d->p_read(0);
        }
      } else
        THROW_EXCEPTION(
            "siconos::simulation::TimeSteppingDirectProjection::advanceToEvent("
            ") :: - Ds is "
            "not from NewtonEulerDS neither from LagrangianDS.");
    }

    computeCriteria(&runningProjection);

    DEBUG_EXPR_WE(
        std::cout << "siconos::simulation::TimeSteppingDirectProjection::"
                     "Projection end : "
                     "Number of iterations="
                  << _nbProjectionIteration << "\n";
        std ::cout << "After update state in position\n";
        std ::cout << "lamda(1) in IndexSet1\n";
        auto indexSet1 = nsds->topology()->indexSet(1);
        auto indexSet0 = nsds->topology()->indexSet(0);

        for (std::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet1->bundle(*ui);
          siconos::algebra::print(*inter->lambda(1));
        }

        std ::cout
        << "lamda(0) in indexSet0\n";
        for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet0->bundle(*ui);
          siconos::algebra::print(*inter->lambda(0));
        } std ::cout
        << "y(1) in IndexSet1\n";
        for (std::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet1->bundle(*ui);
          siconos::algebra::print(*inter->y(1));
        } std ::cout
        << "y(0) in indexSet0\n";
        for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
          auto inter = indexSet0->bundle(*ui);
          siconos::algebra::print(*inter->y(0));
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
    if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      auto time = nextTime();
      neds->computeWrench(neds->twist_read(), neds->q_read(), time);
    } else if (std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearTIDS>(ds)) {
    } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto time = nextTime();
      d->computeTotalForces(d->velocity_read(), d->q_read(), time);
    } else if (std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseLinearTIDS>(ds)) {
    } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseDS>(ds)) {
      auto time = nextTime();
      d->computeTotalForces(d->velocity_read(), d->q_read(), time);
    } else
      THROW_EXCEPTION(
          "TimeSteppingCombinedProjection::advanceToEvent() - Ds is not from "
          "NewtonEulerDS "
          "neither from LagrangianDS.");
  }

  if (_nbProjectionIteration == _projectionMaxIteration && _newtonWarningOnNonConvergence) {
    std::cout << "[kernel] "
                 "siconos::simulation::TimeSteppingDirectProjection::"
                 "advanceToEvent() Max "
                 "number of projection iterations reached ("
              << _nbProjectionIteration << ")" << std::endl;
    printf("[kernel]                max criteria equality =  %e.\n", _maxViolationEquality);
    printf("[kernel]                max criteria unilateral =  %e.\n",
           _maxViolationUnilateral);
  }

  DEBUG_END("siconos::simulation::TimeSteppingDirectProjection::newtonSolve()\n");

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
  //     // (*dotq)(3) = ((*q)(3)-(*qold)(3))/h;
  //     // (*dotq)(4) = ((*q)(4)-(*qold)(4))/h;
  //     // (*dotq)(5) = ((*q)(5)-(*qold)(5))/h;
  //     // (*dotq)(6) = ((*q)(6)-(*qold)(6))/h;

  //     /*compute the new velocity seeing the work of fext*/
  //     auto neds = std::static_pointer_cast<NewtonEulerDS>(ds);
  //     *(neds->deltaq())-=*(neds->q());
  //     DEBUG_EXPR(printf("TSProj NewtonSolve :deltaq:");
  //     siconos::algebra::print(*(neds->deltaq())););
  //     //continue;
  //     double  n2q=neds->deltaq()->norm();
  //     double n2=0.0;
  //     if(neds->fext())
  //       n2=neds->fext()->norm();
  //     if(n2 > 1e-7 && n2q > 1e-14)
  //     {
  //       //if (n2q < 1e-14)
  //       // continue;

  //       auto  FextNorm =
  //       std::make_shared<siconos::algebra::SiconosVector>(3));
  //       (*FextNorm)(0) = neds->fext()(0);
  //       (*FextNorm)(1) = neds->fext()(1);
  //       (*FextNorm)(2) = neds->fext()(2);
  // DEBUG_EXPR_WE(
  //       std::cout<<"siconos::simulation::TimeSteppingDirectProjection::newtonSolve
  //       deltaQ
  //       :\n"; siconos::algebra::print(*neds->deltaq());
  //       std::cout<<"siconos::simulation::TimeSteppingDirectProjection::newtonSolve
  //       Fext
  //       :\n"; siconos::algebra::print(*FextNorm);
  //       );

  //       (*FextNorm)*=(1./n2);
  //       /*work of external forces.*/
  //       double workFext= (*neds->fext())(0) *
  //       (*neds->deltaq())(0)+
  //                        (*neds->fext())(1) *
  //                        (*neds->deltaq())(1)+
  //                        (*neds->fext())(2) *
  //                        (*neds->deltaq())(2);
  //       //workFext*=2.0;
  //       double VkFNorm=(*FextNorm)(0)*(*neds->velocity())(0)+
  //                      (*FextNorm)(1)*(*neds->velocity())(1)+
  //                      (*FextNorm)(2)*(*neds->velocity())(2);
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
  //       printf("VELOCITY1 "); siconos::algebra::print(*neds->velocity());
  // );
  //       (*neds->velocity()->setValue(0,neds->velocity())(0)+(VkcFNorm
  //       - VkFNorm)*(*FextNorm)(0));
  //       (*neds->velocity()->setValue(1,neds->velocity())(1)+(VkcFNorm
  //       - VkFNorm)*(*FextNorm)(1));
  //       (*neds->velocity()->setValue(2,neds->velocity())(2)+(VkcFNorm
  //       - VkFNorm)*(*FextNorm)(2));
  // DEBUG_EXPR_WE(
  //       std::cout<<"siconos::simulation::TimeSteppingDirectProjection::newtonSolve
  //       velocity updated\n"; printf("VELOCITY2 ");
  //       siconos::algebra::print(*neds->velocity());
  // )
  //     }
  //     auto T = neds->T();
  //     auto  dotq = neds->dotq();
  //     dotq->noalias() = *T **neds->velocity();
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
  auto indexSet =
      nonSmoothDynamicalSystem()->topology()->indexSet(_indexSetLevelForProjection);
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
      double criteria = std::max(0.0, -(*inter->y(0))(0));
      DEBUG_PRINTF("Unilateral inter->y(0)(0) %e.\n", (*inter->y(0))(0));
      if (criteria > maxViolationUnilateral) maxViolationUnilateral = criteria;
      // if (criteria < minViolationUnilateral) minViolationUnilateral=criteria;
      if (maxViolationUnilateral > _constraintTolUnilateral) {
        *runningProjection = true;

        DEBUG_PRINTF("TSProj newton criteria unilateral true %e.\n", criteria);
      }
    } else {
      DEBUG_PRINTF("Equality siconos::algebra::normInf(*inter->y(0)) %e.\n",
                   siconos::algebra::normInf(*inter->y(0)));
      if (siconos::algebra::normInf(*inter->y(0)) > maxViolationEquality)
        maxViolationEquality = siconos::algebra::normInf(*inter->y(0));
      if (siconos::algebra::normInf(*inter->y(0)) < minViolationEquality)
        minViolationEquality = siconos::algebra::normInf(*inter->y(0));
      if (siconos::algebra::normInf(*inter->y(0)) > _constraintTol) {
        *runningProjection = true;
        DEBUG_PRINTF("TSProj  newton criteria equality true %e.\n",
                     siconos::algebra::normInf(*inter->y(0)));
      }
    }
    _maxViolationUnilateral = maxViolationUnilateral;
    _maxViolationEquality = maxViolationEquality;
  }

  DEBUG_PRINT("TSProj newton min/max criteria projection\n");
  DEBUG_EXPR(std::cout << "             runningProjection " << std::boolalpha
                       << *runningProjection << std::endl;);
  DEBUG_PRINTF("              min criteria equality =  %e.\n", minViolationEquality);
  DEBUG_PRINTF("              max criteria equality =  %e.\n", maxViolationEquality);
  DEBUG_PRINTF("              max criteria unilateral =  %e.\n", maxViolationUnilateral);
  // DEBUG_PRINTF("              min criteria unilateral =
  // %e.\n",minViolationUnilateral);
}

void siconos::simulation::TimeSteppingDirectProjection::newtonSolve(double criterion,
                                                                    unsigned int maxStep) {
  bool isNewtonConverge = false;
  _newtonNbIterations = 0;  // number of Newton iterations
  int info = 0;
  auto nsds = nonSmoothDynamicalSystem();
  bool isLinear = nsds->isLinear();
  auto indexSet = nsds->topology()->indexSet(0);
  initializeNewtonSolve();

  if ((_newtonOptions == TimeSteppingType::LINEAR ||
       _newtonOptions == TimeSteppingType::LINEAR_IMPLICIT) ||
      isLinear) {
    _newtonNbIterations++;
    prepareNewtonIteration();
    computeFreeState();
    // updateOutput(0);
    // updateIndexSets();
    if (!_allNSProblems->empty() && indexSet->size() > 0)
      info = computeOneStepNSProblem(siconos::simulation::SICONOS_OSNSP_TS_VELOCITY);
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
        info = computeOneStepNSProblem(siconos::simulation::SICONOS_OSNSP_TS_VELOCITY);
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
    THROW_EXCEPTION("TimeStepping::NewtonSolve failed. Unknown newtonOptions: " +
                    siconos::tools::enum_to_string(_newtonOptions));
}
