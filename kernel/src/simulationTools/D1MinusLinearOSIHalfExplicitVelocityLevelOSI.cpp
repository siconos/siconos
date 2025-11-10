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

#include "D1MinusLinearOSI.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "OneStepNSProblem.hpp"
#include "Relation.hpp"
#include "SiconosAlgebraAddons.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "Tools.hpp"  // for enum_to_string

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

double siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel() {
  DEBUG_BEGIN(
      "siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel(), "
      "starts\n");

  double t = _simulation->nextTime();               // end of the time step
  double told = _simulation->startingTime();        // beginning of the time step
  double h = _simulation->timeStep();               // time step length
  auto allOSNS = _simulation->oneStepNSProblems();  // all OSNSP
  auto topo = _simulation->nonSmoothDynamicalSystem()->topology();
  auto indexSet1 = topo->indexSet(1);

  /******************************************************************************************
   *  Step 1-  solve a LCP at velocity level for lambda^+_{k} for the last set indices
   *   if index1 is empty we should skip this step
   ******************************************************************************************/

  DEBUG_PRINT("\nEVALUATE LEFT HAND SIDE\n");

  // DEBUG_EXPR(std::cout<< "allOSNS->empty()   " << std::boolalpha << allOSNS->empty() <<
  // std::endl << std::endl); DEBUG_EXPR(std::cout<< "allOSNS->size()   "  << allOSNS->size()
  // << std::endl << std::endl);

  // -- LEFT SIDE --

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    std::shared_ptr<siconos::algebra::SiconosVector> work_tdg;

    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto& vFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];
      /* POINTER CONSTRUCTOR : vFree will contain the velocity without contact force */
      vFree.setZero();

      // get left state from memory
      auto qold = d->qMemory().getSiconosVector(0);
      auto vold = d->velocityMemory().getSiconosVector(0);  // right limit

      DEBUG_EXPR(siconos::algebra::print(*qold));
      DEBUG_EXPR(siconos::algebra::print(*vold));

      /* compute the force and store in vFree */
      d->computeTotalForces(vold, qold, told);
      vFree += d->totalForces();

      /* Compute the acceleration due to the external force */
      /* vFree contains left (right limit) acceleration without contact force */
      if (d->LUMass()) {
        d->init_lu_mass();
        vFree = d->LUMass()->solve(vFree);
      }

      work_tdg = workVectors[siconos::integrators::D1MinusLinearOSI::FREE_TDG];
      *work_tdg = vFree;
      /*Compute the right limit of the (free) velocity at  t^+_k with contact force :  */
      vFree *= h;
      vFree += vold;

      /* Compute a prediction for q  that will serve for computing new values
         of the Jacobian of the constraints */

      *d->q() = qold + h * vold;
      DEBUG_PRINT(
          "vFree contains the right limit of the (free) velocity at  t^+_k with contact force "
          ":\n");
      DEBUG_EXPR(siconos::algebra::print(vFree));
      DEBUG_PRINT(
          "work_tdg contains left (right limit) acceleration without contact forcework :\n");
      DEBUG_EXPR(siconos::algebra::print(*work_tdg));
    } else if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      auto& vFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];
      vFree.setZero();

      // get left state from memory
      auto qold = neds->qMemory().getSiconosVector(0);
      auto vold = neds->wrenchMemory().getSiconosVector(0);  // right limit

      DEBUG_EXPR(siconos::algebra::print(vFree));
      DEBUG_EXPR(siconos::algebra::print(*qold));
      DEBUG_EXPR(siconos::algebra::print(*vold));

      neds->computeWrench(vold, qold, told);
      vFree += neds->wrench();

      if (neds->LUMass()) {
        vFree = neds->LUMass()->solve(vFree);
        // contains left (right limit) acceleration without contact force
      }

      work_tdg = workVectors[siconos::integrators::D1MinusLinearOSI::FREE_TDG];
      ;
      work_tdg->setZero();
      DEBUG_EXPR(siconos::algebra::print(*work_tdg));
      *work_tdg = vFree;

      /*Compute the right limit of the (free) velocity at  t^+_k with contact force :  */
      vFree *= h;
      vFree += vold;

      // * q = * qold + h* *vold ; to be written consistently for Newton Euler DS
      DEBUG_PRINT(
          "vFree contains the right limit of the (free) velocity at  t^+_k with contact force "
          ":\n");
      DEBUG_EXPR(siconos::algebra::print(vFree));
      DEBUG_PRINT(
          "work_tdg contains left (right limit) acceleration without contact forcework :\n");
      DEBUG_EXPR(siconos::algebra::print(*work_tdg));
    } else {
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel "
          "- only implemented for "
          "Lagrangian or Newton Euler dynamical systems.");
    }
  }

  if (!allOSNS->empty()) {
    if (indexSet1->size() > 0) {
      //_simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(t, *indexSet1);
      siconos::graphs::InteractionsGraph::VIterator ui, uiend;
      std::shared_ptr<siconos::modeling::Interaction> inter;
      for (std::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui) {
        inter = indexSet1->bundle(*ui);
        inter->relation()->computeJach(t, *inter);
        inter->relation()->computeJacg(told, *inter);
      }

      if (_simulation->nonSmoothDynamicalSystem()->topology()->hasChanged()) {
        for (auto osns : *allOSNS) {
          osns->setHasBeenUpdated(false);
        }
      }
      assert((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]);

      if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]
               ->hasInteractions())) /* it should be equivalent to indexSet1 */
      {
        DEBUG_PRINT("We compute lambda^+_{k} \n");
        (*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->compute(told);
        DEBUG_EXPR(siconos::algebra::print(
            *(*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]));
      }
    }
  }

  // Note Franck : at the time this results in a call to swapInMem of all Interactions of the
  // NSDS So let the simu do this.
  //(*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->saveInMemory(); // we push
  // y and lambda in Memories
  _simulation->nonSmoothDynamicalSystem()->pushInteractionsInMemory();
  _simulation->nonSmoothDynamicalSystem()->updateInput(_simulation->nextTime(), 2);

  /**************************************************************************************************************
   *  Step 2 -  compute v_{k,1}
   **************************************************************************************************************/

  DEBUG_PRINT("\n PREDICT RIGHT HAND SIDE\n");

  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    /* \warning the following conditional statement should be removed with a MechanicalDS class
     */
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      // contains residu without nonsmooth effect
      auto& residuFree = *workVectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];
      auto& vFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];
      auto work_tdg = workVectors[siconos::integrators::D1MinusLinearOSI::FREE_TDG];

      // initialize *it->residuFree and predicted right velocity (left limit)
      auto& v = *d->velocity();  // contains velocity v_{k+1}^- and not free velocity
      v.setZero();

      auto& p2 = *d->p(2);
      siconos::algebra::SiconosVector dummy = d->p_read(2);  // value = contact force
      DEBUG_EXPR(siconos::algebra::print(p2));
      /* we homogenize p(2) to a force for the user output   */
      p2 /= h;
      if (d->LUMass()) {
        dummy = d->LUMass()->solve(dummy);
        DEBUG_EXPR(siconos::algebra::print(*d->LUMass()););
      }
      DEBUG_EXPR(siconos::algebra::print(vFree));
      DEBUG_EXPR(siconos::algebra::print(dummy));

      v = vFree + dummy;
      DEBUG_EXPR(siconos::algebra::print(v));

      // get left state from memory
      const auto& qold = d->qMemory().getSiconosVector(0);
      const auto& vold = d->velocityMemory().getSiconosVector(0);
      auto& q = *d->q();  // POINTER CONSTRUCTOR : contains position q_{k+1}
      q = qold;
      q += 0.5 * h * (vold + v);
      DEBUG_EXPR(siconos::algebra::print(*q));

      residuFree.setZero();
      residuFree -= 0.5 * (h * *work_tdg) + 0.5 * dummy;
      DEBUG_EXPR(siconos::algebra::print(residuFree));
    } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      auto& vFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];
      auto& residuFree =
          *workVectors[siconos::integrators::D1MinusLinearOSI::
                           RESIDU_FREE];  // contains residu without nonsmooth effect
      auto work_tdg = workVectors[siconos::integrators::D1MinusLinearOSI::FREE_TDG];

      // get left state from memory
      const auto& qold = d->qMemory().getSiconosVector(0);
      DEBUG_EXPR(const auto& vold = d->twistMemory().getSiconosVector(0););

      // initialize *it->residuFree and predicted right velocity (left limit)
      auto& v = *d->twist();  // contains velocity v_{k+1}^- and not free velocity
      v.setZero();

      auto& p2 = *d->p(2);
      siconos::algebra::SiconosVector dummy(p2);  // value = contact force
      DEBUG_EXPR(siconos::algebra::print(p2));
      /* we homogenize p(2) to a force for the user output   */
      p2 /= h;

      if (d->LUMass()) {
        dummy = d->LUMass()->solve(dummy);
      }

      DEBUG_EXPR(siconos::algebra::print(vFree));
      DEBUG_EXPR(siconos::algebra::print(qold));
      DEBUG_EXPR(siconos::algebra::print(vold));

      v = vFree + dummy;

      DEBUG_EXPR(siconos::algebra::print(v));

      // first step consists in computing  \dot q.
      // second step consists in updating q.
      //
      auto& dotq = *d->dotq();
      dotq = d->T() * v;

      const auto& dotqold = d->dotqMemory().getSiconosVector(0);

      // *d->q()contains position q_{k+1}
      *d->q() = 0.5 * h * (dotqold + dotq);
      *d->q() += qold;
      DEBUG_PRINT("new q before normalizing\n");
      DEBUG_EXPR(siconos::algebra::print(q));
      // q[3:6] must be normalized
      d->normalizeq();
      d->computeT(d->q_read());
      DEBUG_PRINT("new q after normalizing\n");
      DEBUG_EXPR(siconos::algebra::print(q));

      residuFree.setZero();
      residuFree -= 0.5 * (h * *work_tdg) + 0.5 * dummy;
      DEBUG_EXPR(siconos::algebra::print(residuFree));
    } else
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel "
          "- only implemented for "
          "Lagrangian or Newton Euler dynamical systems.");

    /** At this step, we obtain
     * \f[
     * \begin{cases}
     * v_{k,0} = \mbox{\tt vold} \\
     * q_{k,0} = qold \\
     * F_{k,+} = F(told,qold,vold) \\
     * Work_{freefree} =  M^{-1}_k (F^+_{k})  \mbox{stored in work_tdg} \\
     * Work_{free} = vold + h* M^{-1}_k (P^+_{2,k}+F^+_{k})  \mbox{stored in vFree} \\
     * R_{free} = -h/2 * M^{-1}_k (P^+_{2,k}+F^+_{k})  \mbox{stored in ResiduFree} \\
     * v_{k,1} = v_{k,0} + h * M^{-1}_k (P^+_{2,k}+F^+_{k})  \mbox{stored in v} \\
     * q_{k,1} = q_{k,0} + \frac{h}{2} (v_{k,0} + v_{k,1}) \mbox{stored in q} \\
     * \end{cases}
     * \f]
     **/
  }

  DEBUG_PRINT("\n DECIDE STRATEGY\n");
  /** Decide of the strategy impact or smooth multiplier.
   *  Compute _isThereImpactInTheTimeStep
   */
  _isThereImpactInTheTimeStep = false;
  if (!allOSNS->empty()) {
    for (auto level = levelMinForOutput(); level < levelMaxForOutput(); level++) {
      _simulation->nonSmoothDynamicalSystem()->updateOutput(_simulation->nextTime(), level);
    }
    _simulation->updateIndexSets();

    auto topo = _simulation->nonSmoothDynamicalSystem()->topology();
    auto indexSet2 = topo->indexSet(2);

    if (indexSet2->size() > 0) {
      _isThereImpactInTheTimeStep = true;
    } else {
      _isThereImpactInTheTimeStep = false;
    }
  }

  /* If _isThereImpactInTheTimeStep = true;
   * we recompute residuFree by removing the contribution of the nonimpulsive contact forces.
   * We add the contribution of the external forces at the end
   * of the time--step
   * If _isThereImpactInTheTimeStep = false;
   * we recompute residuFree by adding   the contribution of the external forces at the end
   * and the contribution of the nonimpulsive contact forces that are computed by solving the
   * osnsp.
   */
  if (_isThereImpactInTheTimeStep) {
    DEBUG_PRINT(
        "There is an impact in the step. indexSet1->size() > 0.  _isThereImpactInTheTimeStep "
        "= true\n");
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);
      auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

      // type of the current DS
      /* \warning the following conditional statement should be removed with a MechanicalDS
       * class */
      if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
        auto& residuFree = *workVectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];
        auto work_tdg = workVectors[siconos::integrators::D1MinusLinearOSI::FREE_TDG];
        assert(work_tdg);
        DEBUG_EXPR(siconos::algebra::print(*work_tdg));
        residuFree = -0.5 * h * *work_tdg;
        d->computeTotalForces(d->velocity_read(), d->q_read(), t);
        DEBUG_EXPR(siconos::algebra::print(*d->totalForces()););
        *work_tdg = d->totalForces();

        if (d->LUMass()) {
          d->init_lu_mass();
          *work_tdg = d->LUMass()->solve(*work_tdg);
          // contains right (left limit) acceleration without contact force
        }

        DEBUG_EXPR(siconos::algebra::print(*work_tdg));
        residuFree -= 0.5 * h * *work_tdg;
        DEBUG_EXPR(siconos::algebra::print(residuFree));
      } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
        auto& residuFree = *workVectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];
        d->twist()->setZero();  // why?
        auto work_tdg = workVectors[siconos::integrators::D1MinusLinearOSI::FREE_TDG];
        assert(work_tdg);
        residuFree = 0.5 * h * *work_tdg;
        work_tdg->setZero();

        d->computeWrench(d->twist_read(), d->q_read(), t);
        *work_tdg += d->wrench();

        if (d->LUMass()) {
          *work_tdg = d->LUMass()->solve(*work_tdg);
          // contains right (left limit) acceleration without contact force
        }

        residuFree -= 0.5 * h * *work_tdg;
        DEBUG_EXPR(siconos::algebra::print(residuFree));
      } else
        THROW_EXCEPTION(
            "siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel "
            "- only implemented for "
            "Lagrangian or Newton Euler dynamical systems.");
    }
  } else {
    DEBUG_PRINT(
        "There is no  impact in the step. indexSet1->size() = 0. _isThereImpactInTheTimeStep "
        "= false;\n");
    // -- RIGHT SIDE --
    // calculate acceleration without contact force
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);
      auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

      // type of the current DS
      /* \warning the following conditional statement should be removed with a MechanicalDS
       * class */
      if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
        auto& vFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];
        vFree.setZero();
        // get right state from memory
        // d->q();         // contains position q_{k+1}
        // d->velocity();  // contains velocity v_{k+1}^- and not free velocity
        d->computeTotalForces(d->velocity_read(), d->q_read(), t);
        DEBUG_EXPR(siconos::algebra::print(*d->totalForces()));
        vFree += d->totalForces();

        if (d->LUMass()) {
          d->init_lu_mass();
          vFree = d->LUMass()->solve(vFree);
        }
        /* vFree contains right (left limit) acceleration without contact force */
        auto& residuFree = *workVectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];
        residuFree += -0.5 * h * vFree;
        DEBUG_EXPR(siconos::algebra::print(residuFree));
        /* Compute the right limit of the (free) velocity at  t^+_k with contact force : */
        const auto& vold = d->velocityMemory().getSiconosVector(0);
        DEBUG_EXPR(siconos::algebra::print(vold));

        vFree = vold - residuFree;

        DEBUG_PRINT(
            "vFree contains the right  limit of the (free) velocity at  t^-_{k+1} without "
            "contact force :\n");
        DEBUG_EXPR(siconos::algebra::print(vFree));
      } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
        auto& vFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];
        vFree.setZero();
        // get right state from memory
        auto q = d->q();      // contains position q_{k+1}
        auto v = d->twist();  // contains velocity v_{k+1}^- and not free velocity

        DEBUG_EXPR(siconos::algebra::print(*q));
        DEBUG_EXPR(siconos::algebra::print(*v));

        d->computeWrench(d->twist_read(), d->q_read(), t);
        vFree += d->wrench();

        if (d->LUMass()) {
          vFree = d->LUMass()->solve(vFree);
        }
        /* work_tdg contains right (left limit) acceleration without contact force */
        auto& residuFree = *workVectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];

        residuFree += -0.5 * h * vFree;

        /*  Compute the right limit of the (free) velocity at  t^+_k with contact force : */
        const auto& vold = d->twistMemory().getSiconosVector(0);
        vFree = vold - residuFree;
        DEBUG_PRINT(
            "vFree contains the  left limit of the (free) velocity at  t^-_{k+1} without "
            "contact force :\n");
        DEBUG_EXPR(siconos::algebra::print(vFree));
      } else
        THROW_EXCEPTION(
            "siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel "
            "- only implemented for "
            "Lagrangian or Newton Euler dynamical systems.");
    }

    // solve a LCP at acceleration level only for contacts which have been active at the
    // beginning of the time-step
    if (!allOSNS->empty()) {
      // for (unsigned int level = _simulation->levelMinForOutput(); level <
      // _simulation->levelMaxForOutput(); level++)
      // {
      //   _simulation->updateOutput(level);
      // }
      // _simulation->updateIndexSets();
      DEBUG_PRINT("We compute lambda^-_{k+1} \n");
      DEBUG_PRINTF("indexSet1->size() = %i\n", (int)indexSet1->size());

      _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(t, *indexSet1);

      if (_simulation->nonSmoothDynamicalSystem()->topology()->hasChanged()) {
        for (auto osns : *allOSNS) {
          osns->setHasBeenUpdated(false);
        }
      }
      if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]
               ->hasInteractions())) {
        (*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->compute(t);
        DEBUG_EXPR(siconos::algebra::print(
                       *(*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]););
      }
    }

    _simulation->nonSmoothDynamicalSystem()->updateInput(_simulation->nextTime(), 2);
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);
      auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
      // type of the current DS
      /* \warning the following conditional statement should be removed with a MechanicalDS
       * class */
      if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
        auto& residuFree = *workVectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];
        DEBUG_EXPR(siconos::algebra::print(residuFree););

        if (d->p(2)) {
          // get right state from memory
          auto& p2 = *d->p(2);
          siconos::algebra::SiconosVector dummy = d->p_read(2);
          DEBUG_EXPR(siconos::algebra::print(p2));
          /* we homogenize p(2) to a force for the user output   */
          p2 *= 2.0 / h;

          if (d->LUMass()) {
            d->init_lu_mass();
            dummy = d->LUMass()->solve(dummy);
          }
          residuFree -= dummy;
        }
        DEBUG_EXPR(siconos::algebra::print(residuFree));
      } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
        auto& residuFree = *workVectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];

        if (d->p(2)) {
          // get right state from memory
          auto& p2 = *d->p(2);
          siconos::algebra::SiconosVector dummy(*(d->p(2)));  // value = contact force
          DEBUG_EXPR(siconos::algebra::print(p2));
          /* we homogenize p(2) to a force for the user output   */
          p2 *= 2.0 / h;
          if (d->LUMass()) {
            dummy = d->LUMass()->solve(dummy);
          }
          residuFree -= dummy;
        }
      } else
        THROW_EXCEPTION(
            "siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel "
            "- only implemented for "
            "Lagrangian or Newton Euler dynamical systems.");

      /**
       * \f[
       * \begin{cases}
       * v_{k,0} = \mbox{vold} \\
       * q_{k,0} = \mbox{qold} \\
       * F^+_{k} = \mbox{F(told,qold,vold)} \\
       * v_{k,1} = v_{k,0} + h M^{-1}_k (P^+_{2,k}+F^+_{k}) \\[2mm]
       * q_{k,1} = q_{k,0} + \frac{h}{2} (v_{k,0} + v_{k,1})  \\[2mm]
       * F^-_{k+1} = F(t^{-}_{k+1},q_{k,1},v_{k,1}) \\[2mm]
       * R_{free} = - \frac{h}{2}  M^{-1}_k (P^+_{2,k}+F^+_{k}) -  \frac{h}{2}  M^{-1}_{k+1}
       *(P^-_{2,k+1}+F^-_{k+1}) \\[2mm] \end{cases} \f]
       **/
    }

  }  // No impact

  DEBUG_END(
      "siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel(), "
      "ends\n");

  return 0.;  // there is no Newton iteration and the residuum is assumed to vanish
}

void siconos::integrators::D1MinusLinearOSI::computeFreeOutputHalfExplicitVelocityLevel(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) {
  DEBUG_PRINT(
      "\n siconos::integrators::D1MinusLinearOSI::computeFreeOutputHalfExplicitVelocityLevel "
      "starts\n");
  auto allOSNS = _simulation->oneStepNSProblems();  // all OSNSP
  auto indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  auto inter = indexSet->bundle(vertex_inter);
  auto& DSlink = inter->linkToDSVariables();
  auto& workBlockV = *indexSet->properties(vertex_inter).workBlockVectors;
  // get relation and non smooth law information
  auto relationType = inter->relation()->getType();  // relation
  auto relationSubType = inter->relation()->getSubType();

  std::shared_ptr<siconos::algebra::SiconosMatrix>
      C;  // Jacobian of Relation with respect to degree of freedom
  std::shared_ptr<siconos::algebra::BlockVector> Xfree;  // free degree of freedom
  auto& osnsp_rhs = *(*indexSet->properties(vertex_inter)
                           .workVectors)[siconos::integrators::D1MinusLinearOSI::OSNSP_RHS];

  DEBUG_PRINT("osnsp_rhs before\n");
  DEBUG_EXPR(siconos::algebra::print(osnsp_rhs););
  // define Xfree for velocity and acceleration level
  if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp) {
    /* get the current velocity  of the aggregated ds */
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      Xfree = DSlink[tools::enum_to_index(modeling::LagrangianR::WorkDS::q1)];
      DEBUG_PRINT("Xfree = DSlink[siconos::modeling::LagrangianR::q1];\n");
    } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
      Xfree = DSlink[tools::enum_to_index(modeling::NewtonEulerR::WorkDS::velocity)];
      DEBUG_PRINT("Xfree = DSlink[siconos::modeling::NewtonEulerR::velocity];\n");
    } else
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::computeFreeOutput - unknown relation "
          "type.");
    DEBUG_PRINT("Xfree contains the current velocity of the aggregated ds];\n");
    DEBUG_EXPR(siconos::algebra::print(*Xfree));
  } else if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp) {
    /* get the free velocity of the aggregated ds */
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      Xfree = workBlockV[siconos::integrators::D1MinusLinearOSI::xfree];
    } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
      Xfree = workBlockV[siconos::integrators::D1MinusLinearOSI::xfree];
    } else
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::computeFreeOutput - unknown relation "
          "type.");
    assert(Xfree);
    DEBUG_PRINT("Xfree = DSlink[Lagrangian/NewtonEulerR::xfree];\n");
    DEBUG_PRINT("Xfree contains the free velocity of the aggregated ds];\n");
    DEBUG_EXPR(siconos::algebra::print(*Xfree));
  } else
    THROW_EXCEPTION(
        "siconos::integrators::D1MinusLinearOSI::computeFreeOutput - OSNSP neither on "
        "velocity nor on acceleration level.");

  // calculate data of interaction
  auto mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());

  // only Lagrangian Systems
  if (relationType == siconos::modeling::RelationType::Lagrangian) {
    auto lagr =
        std::static_pointer_cast<siconos::modeling::LagrangianR>(mainInteraction->relation());
    // in osnsp_rhs the linear part of velocity or acceleration relation will be saved
    auto C = lagr->jacobianhOver_q();

    DEBUG_EXPR(siconos::algebra::print(*C););
    assert(Xfree);
    siconos::algebra::matrixBlockVector_prod(C, *Xfree, osnsp_rhs, true);

    DEBUG_EXPR(siconos::algebra::print(osnsp_rhs););
    /*  explicit time dependence -> partial time derivative has to be added */
    if (relationSubType == siconos::modeling::RelationSubType::RheonomousR) {
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::computeFreeOutput is not implemented  at "
          "velocity level for LagrangianRheonomousR.");
    }

    /* add the contribution due to the coefficient of restitution*/
    if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp) {
      auto nslEffectOnFreeOutput = std::make_shared<_NSLEffectOnFreeOutput>(
          osnsp, inter, indexSet->properties(vertex_inter));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }

    if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp) {
      /*Do nothing*/
    }
    DEBUG_EXPR(siconos::algebra::print(osnsp_rhs););
  }
  /*Newton-Euler */
  else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
    auto ner =
        std::static_pointer_cast<siconos::modeling::NewtonEulerR>(mainInteraction->relation());
    auto CT = ner->H_NE_prod_T();
    DEBUG_EXPR(siconos::algebra::print(*CT));
    assert(Xfree);
    // creates a POINTER link between workX[ds] (xfree) and the
    // corresponding interactionBlock in each Interaction for each ds of the
    // current Interaction.
    // XXX Big quirks !!! -- xhub
    siconos::algebra::matrixBlockVector_prod(CT, *Xfree, osnsp_rhs, true);

    /* add the contribution due to the coefficient of restitution*/
    if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp) {
      auto nslEffectOnFreeOutput = std::make_shared<_NSLEffectOnFreeOutput>(
          osnsp, inter, indexSet->properties(vertex_inter));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }

    if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp) {
      /* Do nothing*/
    }
    DEBUG_EXPR(siconos::algebra::print(osnsp_rhs););
  }

  else
    THROW_EXCEPTION(
        "siconos::integrators::D1MinusLinearOSI::computeFreeOutput - not implemented for "
        "Relation of type " +
        siconos::tools::enum_to_string(relationType));

  DEBUG_EXPR(siconos::algebra::print(osnsp_rhs););

  DEBUG_PRINT(
      "siconos::integrators::D1MinusLinearOSI::computeFreeOutputHalfExplicitVelocityLevel "
      "ends\n");
}

bool siconos::integrators::D1MinusLinearOSI::addInteractionInIndexSetHalfExplicitVelocityLevel(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  DEBUG_PRINT(
      "siconos::integrators::D1MinusLinearOSI::"
      "addInteractionInIndexSetHalfExplicitVelocityLevel.\n");

  if (i == 1) {
    DEBUG_PRINT(" level == 1\n");
    if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
        siconos::modeling::Type::EqualityConditionNSL) {
      return true;
    }
    /* Active for impulses calculations? Contacts that are closed */
    double y = (*(inter->y(0)))(0);  // current position
    DEBUG_PRINTF("y= %18.14e\n", y);
    return (y <= DEFAULT_TOL_D1MINUS);
  } else if (i == 2) {
    if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
        siconos::modeling::Type::EqualityConditionNSL) {
      return false;
    }
    /*  Contacts which have been closing in the last time step */
    DEBUG_PRINT(" level == 2\n");
    double y = (*(inter->y(0)))(0);   // current position
    double y_k = (inter->y_k(0))(0);  // old position
    DEBUG_PRINTF("y= %18.14e\n", y);
    DEBUG_PRINTF("y_k= %18.14e\n", y_k);
    /* if Interaction has not been active in the previous calculation
       and now becomes active */
    return (y <= DEFAULT_TOL_D1MINUS && y_k > DEFAULT_TOL_D1MINUS);
  } else
    THROW_EXCEPTION(
        "siconos::integrators::D1MinusLinearOSI::"
        "addInteractionInIndexSetHalfExplicitVelocityLevel, IndexSet[i > 2] does not exist.");
  return false;
}

bool siconos::integrators::D1MinusLinearOSI::
    removeInteractionFromIndexSetHalfExplicitVelocityLevel(
        std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  DEBUG_PRINT(
      "siconos::integrators::D1MinusLinearOSI::"
      "removeInteractionFromIndexSetHalfExplicitVelocityLevel.\n");

  return !(addInteractionInIndexSetHalfExplicitVelocityLevel(inter, i));
}
