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
#include "LagrangianScleronomousR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "NonSmoothLaw.hpp"
#include "OneStepNSProblem.hpp"
#include "Relation.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // block prod
#include "SiconosPointers.hpp"        // For createSPtr
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // for scal
#include "Simulation.hpp"
#include "Tools.hpp"  // for enum_to_string
#include "Topology.hpp"
#include "TypeName.hpp"  // for siconos::types visitors
//
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

double siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitAccelerationLevel() {
  DEBUG_BEGIN(
      "\n "
      "siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitAccelerationLevel()"
      "\n");

  double t = _simulation->nextTime();               // end of the time step
  double told = _simulation->startingTime();        // beginning of the time step
  double h = _simulation->timeStep();               // time step length
  auto allOSNS = _simulation->oneStepNSProblems();  // all OSNSP
  auto topo = _simulation->nonSmoothDynamicalSystem()->topology();
  auto indexSet2 = topo->indexSet(2);

  /**************************************************************************************************************
   *  Step 1-  solve a LCP at acceleration level for lambda^+_{k} for the last set indices
   *   if index2 is empty we should skip this step
   **************************************************************************************************************/

  DEBUG_PRINT("\nEVALUATE LEFT HAND SIDE\n");

  DEBUG_EXPR(std::cout << "allOSNS->empty()   " << std::boolalpha << allOSNS->empty()
                       << std::endl
                       << std::endl);
  DEBUG_EXPR(std::cout << "allOSNS->size()   " << allOSNS->size() << std::endl << std::endl);

  // -- LEFT SIDE --

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    DEBUG_EXPR(ds->display());
    std::shared_ptr<siconos::algebra::SiconosVector> work_tdg;

    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto& accFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];
      /* POINTER CONSTRUCTOR : will contain the acceleration without contact force */
      accFree.setZero();

      // get left state from memory
      auto qold = d->qMemory().getSiconosVector(0);
      auto vold = d->velocityMemory().getSiconosVector(0);  // right limit
      DEBUG_EXPR(accFree.display());
      DEBUG_EXPR(qold.display());
      DEBUG_EXPR(vold.display());

      /* compute the force and store in accFree */
      d->computeTotalForces(vold, qold, told);

      DEBUG_EXPR(d->totalForces()->display());
      accFree += d->totalForces();

      /* Compute the acceleration due to the external force */
      /* accFree contains left (right limit) acceleration without contact force */
      if (d->LUMass()) {
        d->update_lu_mass();
        accFree = d->LUMass()->solve(accFree);
      }

      /* Store the value of accFree in workspace(::FREE_TDG) */
      work_tdg = workVectors[siconos::integrators::D1MinusLinearOSI::FREE_TDG];
      work_tdg->setZero();
      *work_tdg = accFree;  // store the value in WorkFreeFree

      DEBUG_PRINT(
          "accFree contains right limit acceleration at  t^+_k with contact force :\n");
      DEBUG_EXPR(accFree.display());
      DEBUG_PRINT(
          "work_tdg contains right limit acceleration at t^+_k without contact force :\n");
      DEBUG_EXPR(work_tdg->display());
    } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      auto& accFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];

      // get left state from memory
      auto qold = d->qMemory().getSiconosVector(0);
      auto vold = d->twistMemory().getSiconosVector(0);  // right limit

      DEBUG_EXPR(accFree.display());
      DEBUG_EXPR(qold.display());
      DEBUG_EXPR(vold.display());

      work_tdg = workVectors[siconos::integrators::D1MinusLinearOSI::FREE_TDG];
      work_tdg->setZero();
      DEBUG_EXPR(work_tdg->display());

      d->computeWrench(vold, qold, told);
      accFree += d->wrench();

      if (d->LUMass()) {
        accFree = d->LUMass()->solve(accFree);
        // contains left (right limit) acceleration without contact force
      }
      *work_tdg = accFree;  // store the value in WorkFreeFree

      DEBUG_PRINT(
          "accFree contains right limit acceleration at  t^+_k with contact force :\n");
      DEBUG_EXPR(accFree.display());
      DEBUG_PRINT(
          "work_tdg contains right limit acceleration at t^+_k without contact force :\n");
      DEBUG_EXPR(work_tdg->display());
    } else {
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitAccelerationLevel "
          "- only implemented for "
          "Lagrangian or Newton Euler dynamical systems.");
    }
  }

  if (!allOSNS->empty()) {
    if (indexSet2->size() > 0) {
      siconos::graphs::InteractionsGraph::VIterator ui, uiend;
      std::shared_ptr<siconos::modeling::Interaction> inter;
      for (std::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui) {
        inter = indexSet2->bundle(*ui);
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
               ->hasInteractions()))  // it should be equivalent to indexSet2
      {
        DEBUG_PRINT("We compute lambda^+_{k} \n");
        (*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->compute(told);
        DEBUG_EXPR((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->display());
      }

      // Note Franck : at the time this results in a call to swapInMem of all Interactions of
      // the NSDS So let the simu do this.
      //(*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->saveInMemory(); // we
      // push y and lambda in Memories
      _simulation->nonSmoothDynamicalSystem()->pushInteractionsInMemory();
      _simulation->nonSmoothDynamicalSystem()->updateInput(_simulation->nextTime(), 2);

      for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
        if (!checkOSI(dsi)) continue;
        auto ds = _dynamicalSystemsGraph->bundle(*dsi);

        auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

        if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
          auto& accFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];

          siconos::algebra::SiconosVector dummy = *(d->p(2));  // copy
          if (d->LUMass()) {
            d->update_lu_mass();
            dummy = d->LUMass()->solve(dummy);
          }
          accFree += dummy;

          DEBUG_EXPR(d->p(2)->display());
        } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
          auto& accFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];
          siconos::algebra::SiconosVector dummy = *(d->p(2));  // copy
          if (d->LUMass()) {
            dummy = d->LUMass()->solve(dummy);
          }
          accFree += dummy;

          DEBUG_EXPR(d->p(2)->display());
        } else
          THROW_EXCEPTION(
              "siconos::integrators::D1MinusLinearOSI::"
              "computeResiduHalfExplicitAccelerationLevel - only implemented for "
              "Lagrangian or Newton Euler dynamical systems.");
      }
    }
  }

  /**************************************************************************************************************
   *  Step 2 -  compute v_{k,1}
   **************************************************************************************************************/

  DEBUG_PRINT("\n PREDICT RIGHT HAND SIDE\n");

  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    // type of the current DS
    /* \warning the following conditional statement should be removed with a MechanicalDS class
     */
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      // contains residu without nonsmooth effect
      auto& residuFree = *workVectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];
      auto& accFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];

      residuFree.setZero();

      DEBUG_EXPR(accFree.display());
      DEBUG_EXPR(qold.display());
      DEBUG_EXPR(vold.display());

      residuFree -= 0.5 * h * accFree;

      const auto& vold = d->velocityMemory().getSiconosVector(0);

      *d->velocity() = h * accFree;
      *d->velocity() += vold;

      DEBUG_EXPR(residuFree.display());
      DEBUG_EXPR(v.display());

      *d->q() = 0.5 * h * vold;
      *d->q() += d->qMemory().getSiconosVector(0);

    } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      auto& residuFree =
          *workVectors[siconos::integrators::D1MinusLinearOSI::
                           RESIDU_FREE];  // contains residu without nonsmooth effect
      auto& accFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];
      // get left state from memory
      const auto& qold = d->qMemory().getSiconosVector(0);
      const auto& vold = d->twistMemory().getSiconosVector(0);

      // initialize *it->residuFree and predicted right velocity (left limit)
      residuFree.setZero();

      DEBUG_EXPR(accFree.display());
      DEBUG_EXPR(qold.display());
      DEBUG_EXPR(vold.display());

      residuFree -= 0.5 * h * accFree;

      *d->twist() = h * accFree;
      *d->twist() += vold;

      DEBUG_EXPR(residuFree.display());
      DEBUG_EXPR(v.display());

      // first step consists in computing  \dot q.
      // second step consists in updating q.
      //
      *d->dotq() = d->T() * d->twist_read();
      const auto& dotqold = d->dotqMemory().getSiconosVector(0);

      *d->q() = 0.5 * h * (dotqold + d->dotq_read());
      *d->q() += qold;

      // q[3:6] must be normalized
      DEBUG_EXPR(std::cout << d->q_read() << "\n");
      d->normalizeq();

      d->computeT(d->q_read());
      DEBUG_PRINT("new q after normalizing\n");
      DEBUG_EXPR(std::cout << d->q_read() << "\n");
    } else
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::"
          "computeResiduHalfExplicitAccelerationLevel - only implemented for "
          "Lagrangian or Newton Euler dynamical systems.");

    /** At this step, we obtain
     * \f[
     * \begin{cases}
     * v_{k,0} = \mbox{\tt vold} \\
     * q_{k,0} = qold \\
     * F_{k,+} = F(told,qold,vold) \\
     * Work_{freefree} =  M^{-1}_k (F^+_{k})  \mbox{stored in work_tdg} \\
     * Work_{free} =  M^{-1}_k (P^+_{2,k}+F^+_{k})  \mbox{stored in accFree} \\
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
    auto indexSet3 = topo->indexSet(3);

    if (indexSet3->size() > 0) {
      _isThereImpactInTheTimeStep = true;
      DEBUG_PRINT(
          "There is an impact in the step. indexSet3->size() > 0. _isThereImpactInTheTimeStep "
          "= true;\n");
    } else {
      _isThereImpactInTheTimeStep = false;
      DEBUG_PRINT(
          "There is no  impact in the step. indexSet3->size() = 0. "
          "_isThereImpactInTheTimeStep = false;\n");
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
        "There is an impact in the step. indexSet3->size() > 0.  _isThereImpactInTheTimeStep "
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
        residuFree = -0.5 * h * *work_tdg;

        d->computeTotalForces(d->velocity_read(), d->q_read(), t);
        *work_tdg = d->totalForces();
        if (d->LUMass()) {
          d->update_lu_mass();
          *work_tdg = d->LUMass()->solve(*work_tdg);
        }
        residuFree -= 0.5 * h * *work_tdg;
        DEBUG_EXPR(residuFree.display());
      } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
        auto& residuFree = *workVectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];
        const auto qold = d->qMemory().getSiconosVector(0);
        const auto vold = d->twistMemory().getSiconosVector(0);  // right limit
        d->twist()->setZero();                                   // FP: why?
        auto work_tdg = workVectors[siconos::integrators::D1MinusLinearOSI::FREE_TDG];
        assert(work_tdg);
        residuFree = 0.5 * h * *work_tdg;
        d->computeWrench(d->twist_read(), d->q_read(), t);
        *work_tdg = d->wrench();

        if (d->LUMass()) {
          *work_tdg = d->LUMass()->solve(*work_tdg);
        }
        residuFree -= 0.5 * h * *work_tdg;
        DEBUG_EXPR(residuFree.display());
      } else
        THROW_EXCEPTION(
            "siconos::integrators::D1MinusLinearOSI::"
            "computeResiduHalfExplicitAccelerationLevel - only implemented for "
            "Lagrangian or Newton Euler dynamical systems.");
    }
  } else {
    DEBUG_PRINT(
        "There is no  impact in the step. indexSet3->size() = 0. _isThereImpactInTheTimeStep "
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
        auto& accFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];
        accFree.setZero();
        // get right state from memory
        DEBUG_EXPR(accFree.display());
        DEBUG_EXPR(q->display());
        DEBUG_EXPR(v->display());
        // d->q()  contains position q_{k+1}
        // d->velocity() contains velocity v_{k+1}^- and not free velocity

        d->computeTotalForces(d->velocity_read(), d->q_read(), t);
        accFree += d->totalForces();

        if (d->LUMass()) {
          d->update_lu_mass();
          accFree = d->LUMass()->solve(accFree);
          // contains right (left limit) acceleration without contact force
        }
        DEBUG_PRINT(
            "accFree contains left limit acceleration at  t^-_{k+1} without contact force "
            ":\n");
        DEBUG_EXPR(accFree.display());
      } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
        auto& accFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];

        accFree.setZero();
        // get right state from memory
        // d->q()  contains position q_{k+1}
        // d->velocity() contains velocity v_{k+1}^- and not free velocity
        d->computeWrench(d->twist_read(), d->q_read(), t);
        accFree += d->wrench();

        if (d->LUMass()) {
          accFree = d->LUMass()->solve(accFree);
          // contains right (left limit) acceleration without contact force
        }

        DEBUG_PRINT(
            "accFree contains left limit acceleration at  t^-_{k+1} without contact force "
            ":\n");
        DEBUG_EXPR(accFree.display());
      } else
        THROW_EXCEPTION(
            "siconos::integrators::D1MinusLinearOSI::"
            "computeResiduHalfExplicitAccelerationLevel - only implemented for "
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

      _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(t, *indexSet2);

      if (_simulation->nonSmoothDynamicalSystem()->topology()->hasChanged()) {
        for (auto osns : *allOSNS) {
          osns->setHasBeenUpdated(false);
        }
      }

      if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]
               ->hasInteractions())) {
        (*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->compute(t);
        DEBUG_EXPR((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->display(););
        _simulation->nonSmoothDynamicalSystem()->updateInput(_simulation->nextTime(), 2);
      }
    }
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);
      auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
      // type of the current DS
      /* \warning the following conditional statement should be removed with a MechanicalDS
       * class */
      if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
        auto& residuFree = *workVectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];
        auto& accFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];

        // initialize *it->residuFree
        residuFree -= 0.5 * h * accFree;
        DEBUG_EXPR(accFree.display());

        if (d->p(2)) {
          // get right state from memory
          DEBUG_EXPR(d->LUMass()->display());
          DEBUG_EXPR(d->p(2)->display());
          siconos::algebra::SiconosVector dummy = d->p_read(2);  // value = contact force
          if (d->LUMass()) {
            d->update_lu_mass();
            dummy = d->LUMass()->solve(dummy);
          }

          residuFree -= 0.5 * h * dummy;
        }
        DEBUG_EXPR(residuFree.display());
      } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
        auto& residuFree = *workVectors[siconos::integrators::D1MinusLinearOSI::RESIDU_FREE];
        auto& accFree = *workVectors[siconos::integrators::D1MinusLinearOSI::FREE];
        // initialize *it->residuFree
        residuFree -= 0.5 * h * accFree;
        DEBUG_EXPR(accFree.display());

        if (d->p(2)) {
          // get right state from memory
          DEBUG_EXPR(d->LUMass()->display());
          DEBUG_EXPR(d->p(2)->display());
          siconos::algebra::SiconosVector dummy = d->p_read(2);  // value = contact force
          if (d->LUMass()) {
            dummy = d->LUMass()->solve(dummy);
          }
          residuFree -= 0.5 * h * dummy;
        }
      } else
        THROW_EXCEPTION(
            "siconos::integrators::D1MinusLinearOSI::"
            "computeResiduHalfExplicitAccelerationLevel - only implemented for "
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
      "siconos::integrators::D1MinusLinearOSI::computeResiduHalfExplicitAccelerationLevel()"
      "\n");

  return 0.;  // there is no Newton iteration and the residuum is assumed to vanish
}

void siconos::integrators::D1MinusLinearOSI::computeFreeOutputHalfExplicitAccelerationLevel(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) {
  DEBUG_BEGIN(
      "siconos::integrators::D1MinusLinearOSI::"
      "computeFreeOutputHalfExplicitAccelerationLevel\n");
  auto allOSNS = _simulation->oneStepNSProblems();  // all OSNSP
  auto indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  auto inter = indexSet->bundle(vertex_inter);
  auto& DSlink = inter->linkToDSVariables();
  auto& inter_work_block = *indexSet->properties(vertex_inter).workBlockVectors;

  // get relation and non smooth law information
  auto relationType = inter->relation()->getType();  // relation
  auto relationSubType = inter->relation()->getSubType();
  auto sizeY = inter->nonSmoothLaw()->size();  // related NSL

  // Jacobian of Relation with respect to degree of freedom
  std::shared_ptr<siconos::algebra::SiconosMatrix> C{nullptr};

  std::shared_ptr<siconos::algebra::BlockVector> Xfree{nullptr};  // free degree of freedom

  auto& osnsp_rhs = *(*indexSet->properties(vertex_inter)
                           .workVectors)[siconos::integrators::D1MinusLinearOSI::OSNSP_RHS];

  // define Xfree for velocity and acceleration level
  if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp) {
    // Xfree = inter->dataX();
    /* get the current velocity  of the aggregated ds */
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      Xfree = DSlink[tools::enum_to_index(modeling::LagrangianR::WorkDS::q1)];
      DEBUG_PRINT("Xfree = DSlink[siconos::modeling::LagrangianR::q1];\n");
    } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
      Xfree = DSlink[siconos::tools::enum_to_index(modeling::NewtonEulerR::WorkDS::velocity)];
      DEBUG_PRINT("Xfree = DSlink[NewtonEulerR::WorkDS::velocity];\n");
    } else
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::computeFreeOutput - unknown relation "
          "type.");

    DEBUG_EXPR(Xfree->display());
    assert(Xfree);
  } else if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp) {
    /* get the "free" acceleration of the aggregated ds */
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      Xfree = inter_work_block[siconos::integrators::D1MinusLinearOSI::xfree];
    } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
      Xfree = inter_work_block[siconos::integrators::D1MinusLinearOSI::xfree];
    } else
      THROW_EXCEPTION(
          "siconos::integrators::D1MinusLinearOSI::computeFreeOutput - unknown relation "
          "type.");
    DEBUG_PRINT("Xfree = DSlink[siconos::integrators::D1MinusLinearOSI::FREE];\n");
    DEBUG_EXPR(Xfree->display());
    assert(Xfree);
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
    // in osnsp_rhs the linear part of velocity or acceleration relation will be saved
    auto rel =
        std::static_pointer_cast<siconos::modeling::LagrangianR>(mainInteraction->relation());

    assert(Xfree);
    auto C =
        std::static_pointer_cast<siconos::modeling::LagrangianR>(mainInteraction->relation())
            ->jacobianhOver_q();
    siconos::algebra::matrixBlockVector_prod(C, *Xfree, osnsp_rhs, true);

    DEBUG_EXPR(osnsp_rhs.display(););

    // in osnsp_rhs corrections have to be added
    auto ID = std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeY);
    ID->setIdentity();

    if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp) {
      /*  explicit time dependence -> partial time derivative has to be added */

      // if (relationSubType == RheonomousR) // explicit time dependence -> partial time
      // derivative has to be added
      // {
      //   auto  q = *DSlink[tools::enum_to_index(modeling::LagrangianR::WorkDS::q0)];
      //   auto  z = *DSlink[siconos::modeling::LagrangianR::z];

      //   std::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->computehDot(simulation()->getTkp1(),
      //   q, z); *DSlink[siconos::modeling::LagrangianR::z] = z;
      // }

      if (relationSubType == siconos::modeling::RelationSubType::RheonomousR)
        THROW_EXCEPTION(
            "siconos::integrators::D1MinusLinearOSI::computeFreeOutput is not implemented  at "
            "velocity level for LagrangianRheonomousR.");
      auto nslEffectOnFreeOutput = std::make_shared<_NSLEffectOnFreeOutput>(
          osnsp, inter, indexSet->properties(vertex_inter));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }

    if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp) {
      if (relationSubType == siconos::modeling::RelationSubType::
                                 ScleronomousR)  // acceleration term involving Hessian matrix
                                                 // of Relation with respect to degree is added
      {
        DEBUG_PRINT(
            "siconos::integrators::D1MinusLinearOSI::computeFreeOutput. acceleration term "
            "involving Hessian matrix of Relation\n");
        DEBUG_EXPR(osnsp_rhs.display(););
        auto scler = std::static_pointer_cast<siconos::modeling::LagrangianScleronomousR>(
            inter->relation());
        scler->computeJacobianhOver_q_dot_X_qdot(simulation()->getTkp1(), *inter);
        osnsp_rhs += *ID * scler->jacobianhOver_q_dot_X_qdot();
      }
      DEBUG_EXPR(osnsp_rhs.display(););
    }
  } else if (relationType == siconos::modeling::RelationType::NewtonEuler) {
    auto ner =
        std::static_pointer_cast<siconos::modeling::NewtonEulerR>(mainInteraction->relation());
    auto CT = ner->jacobianhOver_q_prod_T();
    DEBUG_EXPR(CT->display());
    assert(Xfree);
    // creates a POINTER link between workX[ds] (xfree) and the
    // corresponding interactionBlock in each Interaction for each ds of the
    // current Interaction.
    // XXX Big quirks !!! -- xhub
    siconos::algebra::matrixBlockVector_prod(CT, *Xfree, osnsp_rhs, true);

    if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp) {
      // in osnsp_rhs corrections have to be added
      auto ID = std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeY);
      ID->setIdentity();

      DEBUG_PRINT(
          "siconos::integrators::D1MinusLinearOSI::computeFreeOutput.\n Adding the additional "
          "terms of the second order time derivative of constraints.\n");
      DEBUG_EXPR(osnsp_rhs.display(););

      /** Compute additional terms of the second order time derivative of constraints
       *
       * \f$ \nabla_q h(q) \dot T v + \frac{d}{dt}(\nabla_q h(q) ) T v \f$
       *
       */
      auto ds1 = indexSet->properties(vertex_inter).source;
      auto ds2 = indexSet->properties(vertex_inter).target;

      ner->computeSecondOrderTimeDerivativeTerms(simulation()->getTkp1(), *inter, ds1, ds2);

      DEBUG_EXPR(ner->secondOrderTimeDerivativeTerms()->display());
      osnsp_rhs += *ID * ner->secondOrderTimeDerivativeTerms();

      DEBUG_EXPR(osnsp_rhs.display(););
    }
    DEBUG_EXPR(osnsp_rhs.display(););

    if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() ==
        osnsp)  // impact terms are added
    {
      auto nslEffectOnFreeOutput = std::make_shared<_NSLEffectOnFreeOutput>(
          osnsp, inter, indexSet->properties(vertex_inter));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }
  } else
    THROW_EXCEPTION(
        "siconos::integrators::D1MinusLinearOSI::computeFreeOutput - not implemented for "
        "Relation of type " +
        siconos::tools::enum_to_string(relationType));

  DEBUG_EXPR(osnsp_rhs.display(););
  DEBUG_END(
      "siconos::integrators::D1MinusLinearOSI::computeFreeOutputHalfExplicitAccelerationLevel "
      "ends\n");
}

bool siconos::integrators::D1MinusLinearOSI::
    addInteractionInIndexSetHalfExplicitAccelerationLevel(
        std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  DEBUG_PRINT(
      "siconos::integrators::D1MinusLinearOSI::"
      "addInteractionInIndexSetHalfExplicitAccelerationLevel.\n");

  if (i == 1) {
    DEBUG_PRINT(" level == 1\n");
    if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
        siconos::modeling::Type::EqualityConditionNSL) {
      return true;
    }
    /* Active for impulses calculations? Contacts that are closed */
    double y = (*(inter->y(0)))(0);  // current position
    DEBUG_PRINTF("y= %24.16e\n", y);
    return (y <= DEFAULT_TOL_D1MINUS);
  } else if (i == 2) {
    DEBUG_PRINT(" level == 2\n");
    if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
        siconos::modeling::Type::EqualityConditionNSL) {
      return true;
    }
    /* Active for non-impulsive forces computation */
    auto y = (*(inter->y(0)))(0);     // current position
    auto yDot = (*(inter->y(1)))(0);  // current position
    DEBUG_PRINTF("y= %24.16e\n", y);
    DEBUG_PRINTF("yDot= %24.16e\n", yDot);
    DEBUG_EXPR(std::cout << std::boolalpha << (y <= DEFAULT_TOL_D1MINUS) << std::endl;);
    DEBUG_EXPR(std::cout << std::boolalpha << (yDot <= DEFAULT_TOL_D1MINUS) << std::endl;);
    return (y <= DEFAULT_TOL_D1MINUS) && (yDot <= DEFAULT_TOL_D1MINUS);
  } else if (i == 3) {
    if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
        siconos::modeling::Type::EqualityConditionNSL) {
      return false;
    }
    /*  Contacts which have been closing in the last time step */
    DEBUG_PRINT(" level == 3\n");
    auto y = (*(inter->y(0)))(0);   // current position
    auto y_k = (inter->y_k(0))(0);  // old position
    DEBUG_PRINTF("y= %24.16e\n", y);
    DEBUG_PRINTF("y_k= %24.16e\n", y_k);
    /* if Interaction has not been active in the previous calculation
       and now becomes active */
    return (y <= DEFAULT_TOL_D1MINUS && y_k > DEFAULT_TOL_D1MINUS);
  } else
    THROW_EXCEPTION(
        "siconos::integrators::D1MinusLinearOSI::"
        "addInteractionInIndexSetHalfExplicitAccelerationLevel, IndexSet[i > 3] does not "
        "exist.");
  return false;
}

bool siconos::integrators::D1MinusLinearOSI::
    removeInteractionFromIndexSetHalfExplicitAccelerationLevel(
        std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  DEBUG_PRINT(
      "siconos::integrators::D1MinusLinearOSI::"
      "removeInteractionFromIndexSetHalfExplicitAccelerationLevel.\n");

  return !(addInteractionInIndexSetHalfExplicitAccelerationLevel(inter, i));
}
