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

#include "D1MinusLinearOSI.hpp"
#include "SiconosAlgebraProd.hpp"
#include "Simulation.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianRheonomousR.hpp"
#include "LagrangianScleronomousR.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonImpactNSL.hpp"
#include "BlockVector.hpp"
#include "Topology.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OneStepNSProblem.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

/// @cond
using namespace RELATION;

/// for non-owned shared pointers (passing const siconos::algebra::SiconosVector into
/// functions that take std::shared_ptr<siconos::algebra::SiconosVector> without copy -- warning
/// const abuse!)
static void null_deleter(const siconos::algebra::SiconosVector *) {}
template <typename T> static std::shared_ptr<T> ptr(const T& a)
{
  return std::shared_ptr<siconos::algebra::SiconosVector>(&*(T*)&a, null_deleter);
}

double D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel()
{
  DEBUG_BEGIN("D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel(), starts\n");

  double t = _simulation->nextTime(); // end of the time step
  double told = _simulation->startingTime(); // beginning of the time step
  double h = _simulation->timeStep(); // time step length
  auto allOSNS  = _simulation->oneStepNSProblems(); // all OSNSP
  auto topo =  _simulation->nonSmoothDynamicalSystem()->topology();
  auto indexSet1 = topo->indexSet(1);

  /******************************************************************************************
   *  Step 1-  solve a LCP at velocity level for lambda^+_{k} for the last set indices
   *   if index1 is empty we should skip this step
   ******************************************************************************************/

  DEBUG_PRINT("\nEVALUATE LEFT HAND SIDE\n");

  // DEBUG_EXPR(std::cout<< "allOSNS->empty()   " << std::boolalpha << allOSNS->empty() << std::endl << std::endl);
  // DEBUG_EXPR(std::cout<< "allOSNS->size()   "  << allOSNS->size() << std::endl << std::endl);

// -- LEFT SIDE --

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    Type::Siconos dsType = Type::value(*ds);
    auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    std::shared_ptr<siconos::algebra::SiconosVector> work_tdg;

    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      auto d = std::static_pointer_cast<LagrangianDS> (ds);
      siconos::algebra::SiconosVector& vFree = *workVectors[D1MinusLinearOSI::FREE];
      /* POINTER CONSTRUCTOR : vFree will contain the velocity without contact force */
      vFree.zero();

      // get left state from memory
      const siconos::algebra::SiconosVector& qold = d->qMemory().getsiconos::algebra::SiconosVector(0);
      const siconos::algebra::SiconosVector& vold = d->velocityMemory().getsiconos::algebra::SiconosVector(0); // right limit

      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());

      /* compute the force and store in vFree */
      d->computeForces(told, ptr(qold), ptr(vold));
      DEBUG_EXPR(d->forces()->display());
      vFree += *(d->forces());

      /* Compute the acceleration due to the external force */
      /* vFree contains left (right limit) acceleration without contact force */
      if(d->inverseMass())
      {
        d->update_inverse_mass();
        d->inverseMass()->Solve(vFree);
      }

      work_tdg =  workVectors[D1MinusLinearOSI::FREE_TDG];
      work_tdg->zero();
      *work_tdg = vFree;


      /*Compute the right limit of the (free) velocity at  t^+_k with contact force :  */
      vFree  *= h ;
      vFree += vold;

      /* Compute a prediction for q  that will serve for computing new values
         of the Jacobian of the constraints */

      *d->q() = qold + h * vold;
      DEBUG_PRINT("vFree contains the right limit of the (free) velocity at  t^+_k with contact force :\n");
      DEBUG_EXPR(vFree.display());
      DEBUG_PRINT("work_tdg contains left (right limit) acceleration without contact forcework :\n");
      DEBUG_EXPR(work_tdg->display());
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      auto d = std::static_pointer_cast<NewtonEulerDS> (ds);
      siconos::algebra::SiconosVector& vFree = *workVectors[D1MinusLinearOSI::FREE];
      vFree.zero();

      // get left state from memory
      const siconos::algebra::SiconosVector& qold = d->qMemory().getsiconos::algebra::SiconosVector(0);
      const siconos::algebra::SiconosVector& vold = d->twistMemory().getsiconos::algebra::SiconosVector(0); // right limit
      DEBUG_EXPR(vFree.display());
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());

      d->computeForces(told, ptr(qold), ptr(vold));
      DEBUG_EXPR(d->forces()->display());
      vFree += *(d->forces());

      if(d->inverseMass())
      {
        d->update_inverse_mass();
        d->inverseMass()->Solve(vFree); // contains left (right limit) acceleration without contact force
      }

      work_tdg =  workVectors[D1MinusLinearOSI::FREE_TDG];;
      work_tdg->zero();
      DEBUG_EXPR(work_tdg->display());
      *work_tdg = vFree;

      /*Compute the right limit of the (free) velocity at  t^+_k with contact force :  */
      vFree *= h ;
      vFree += vold;

      // * q = * qold + h* *vold ; to be written consistently for Newton Euler DS
      DEBUG_PRINT("vFree contains the right limit of the (free) velocity at  t^+_k with contact force :\n");
      DEBUG_EXPR(vFree.display());
      DEBUG_PRINT("work_tdg contains left (right limit) acceleration without contact forcework :\n");
      DEBUG_EXPR(work_tdg->display());
    }
    else
    {
      THROW_EXCEPTION("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " +  std::to_string(dsType));
    }
  }


  if(!allOSNS->empty())
  {
    if(indexSet1->size() >0)
    {

      //_simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(t, *indexSet1);
      siconos::graphs::InteractionsGraph::VIterator ui, uiend;
      std::shared_ptr<siconos::modeling::Interaction> inter;
      for(std::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui)
      {
        inter = indexSet1->bundle(*ui);
        inter->relation()->computeJach(t, *inter);
        inter->relation()->computeJacg(told, *inter);
      }

      if(_simulation->nonSmoothDynamicalSystem()->topology()->hasChanged())
      {
	for(auto osns: *allOSNS)
        {
          osns->setHasBeenUpdated(false);
        }
      }
      assert((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]);

      if(((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->hasInteractions()))  /* it should be equivalent to indexSet1 */
      {
        DEBUG_PRINT("We compute lambda^+_{k} \n");
        (*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->compute(told);
        DEBUG_EXPR((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->display());
      }

    }
  }

  // Note Franck : at the time this results in a call to swapInMem of all Interactions of the NSDS
  // So let the simu do this.
  //(*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->saveInMemory(); // we push y and lambda in Memories
  _simulation->nonSmoothDynamicalSystem()->pushInteractionsInMemory();
  _simulation->nonSmoothDynamicalSystem()->updateInput(_simulation->nextTime(),2);

  /**************************************************************************************************************
   *  Step 2 -  compute v_{k,1}
   **************************************************************************************************************/


  DEBUG_PRINT("\n PREDICT RIGHT HAND SIDE\n");

  for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    // type of the current DS
    Type::Siconos dsType = Type::value(*ds);
    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      auto d = std::static_pointer_cast<LagrangianDS> (ds);
      // contains residu without nonsmooth effect
      siconos::algebra::SiconosVector& residuFree = *workVectors[D1MinusLinearOSI::RESIDU_FREE];
      siconos::algebra::SiconosVector& vFree = *workVectors[D1MinusLinearOSI::FREE];
      std::shared_ptr<siconos::algebra::SiconosVector> work_tdg = workVectors[D1MinusLinearOSI::FREE_TDG];



      // initialize *it->residuFree and predicted right velocity (left limit)
      siconos::algebra::SiconosVector& v = *d->velocity(); //contains velocity v_{k+1}^- and not free velocity
      v.zero();




      siconos::algebra::SiconosVector& p2 = *d->p(2);
      siconos::algebra::SiconosVector dummy(p2); // value = contact force
      DEBUG_EXPR(p2.display());
      /* we homogenize p(2) to a force for the user output   */
      p2 /= h;
      if(d->inverseMass())
      {
        d->inverseMass()->Solve(dummy);
        DEBUG_EXPR(d->inverseMass()->display(););
      }
      DEBUG_EXPR(vFree.display());
      DEBUG_EXPR(dummy.display());

      v = vFree + dummy;
      DEBUG_EXPR(v.display());


      // get left state from memory
      const siconos::algebra::SiconosVector& qold = d->qMemory().getsiconos::algebra::SiconosVector(0);
      const siconos::algebra::SiconosVector& vold = d->velocityMemory().getsiconos::algebra::SiconosVector(0);
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());

      siconos::algebra::SiconosVector& q = *d->q(); // POINTER CONSTRUCTOR : contains position q_{k+1}
      q = qold;

      siconos::algebra::scal(0.5 * h, vold + v, q, false);
      DEBUG_EXPR(q->display());


      residuFree.zero();
      residuFree -= 0.5 * (h * *work_tdg) + 0.5 * dummy;
      DEBUG_EXPR(residuFree.display());
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      auto d = std::static_pointer_cast<NewtonEulerDS> (ds);
      siconos::algebra::SiconosVector& vFree = *workVectors[D1MinusLinearOSI::FREE];
      siconos::algebra::SiconosVector& residuFree = *workVectors[D1MinusLinearOSI::RESIDU_FREE];// contains residu without nonsmooth effect
      std::shared_ptr<siconos::algebra::SiconosVector> work_tdg = workVectors[D1MinusLinearOSI::FREE_TDG];

      // get left state from memory
      const siconos::algebra::SiconosVector& qold = d->qMemory().getsiconos::algebra::SiconosVector(0);
      DEBUG_EXPR(const siconos::algebra::SiconosVector& vold = d->twistMemory().getsiconos::algebra::SiconosVector(0););

      // initialize *it->residuFree and predicted right velocity (left limit)
      siconos::algebra::SiconosVector& v = *d->twist(); //contains velocity v_{k+1}^- and not free velocity
      v.zero();


      siconos::algebra::SiconosVector& p2 = *d->p(2);
      siconos::algebra::SiconosVector dummy(p2); // value = contact force
      DEBUG_EXPR(p2.display());
      /* we homogenize p(2) to a force for the user output   */
      p2 /= h;

      if(d->inverseMass())
      {
        d->update_inverse_mass();
        d->inverseMass()->Solve(dummy);
      }

      DEBUG_EXPR(vFree.display());
      DEBUG_EXPR(qold.display());
      DEBUG_EXPR(vold.display());

      v = vFree +  dummy;


      DEBUG_EXPR(v.display());

      //first step consists in computing  \dot q.
      //second step consists in updating q.
      //
      auto T = d->T();
      siconos::algebra::SiconosVector& dotq = *d->dotq();
      siconos::algebra::prod(*T, v, dotq, true);

      const siconos::algebra::SiconosVector& dotqold = d->dotqMemory().getsiconos::algebra::SiconosVector(0);

      siconos::algebra::SiconosVector& q = *d->q(); // contains position q_{k+1}
      q = qold;

      siconos::algebra::scal(0.5 * h, dotqold + dotq, q, false);
      DEBUG_PRINT("new q before normalizing\n");
      DEBUG_EXPR(q.display());
      //q[3:6] must be normalized
      d->normalizeq();
      d->computeT();
      DEBUG_PRINT("new q after normalizing\n");
      DEBUG_EXPR(q.display());


      residuFree.zero();
      residuFree -= 0.5 * (h * *work_tdg) + 0.5 * dummy;
      DEBUG_EXPR(residuFree.display());



    }
    else
      THROW_EXCEPTION("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " +  std::to_string(dsType));


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
  if(!allOSNS->empty())
  {

    for(unsigned int level = levelMinForOutput(); level < levelMaxForOutput(); level++)
    {
      _simulation->nonSmoothDynamicalSystem()->updateOutput(_simulation->nextTime(),level);
    }
    _simulation->updateIndexSets();

    auto topo =  _simulation->nonSmoothDynamicalSystem()->topology();
    auto indexSet2 = topo->indexSet(2);

    if(indexSet2->size() > 0)
    {
      _isThereImpactInTheTimeStep = true;
    }
    else
    {
      _isThereImpactInTheTimeStep = false;
    }
  }

  /* If _isThereImpactInTheTimeStep = true;
   * we recompute residuFree by removing the contribution of the nonimpulsive contact forces.
   * We add the contribution of the external forces at the end
   * of the time--step
   * If _isThereImpactInTheTimeStep = false;
   * we recompute residuFree by adding   the contribution of the external forces at the end
   * and the contribution of the nonimpulsive contact forces that are computed by solving the osnsp.
   */
  if(_isThereImpactInTheTimeStep)
  {

    DEBUG_PRINT("There is an impact in the step. indexSet1->size() > 0.  _isThereImpactInTheTimeStep = true\n");
    for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);
      auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

      // type of the current DS
      Type::Siconos dsType = Type::value(*ds);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
      {
        auto d = std::static_pointer_cast<LagrangianDS> (ds);
        siconos::algebra::SiconosVector& residuFree = *workVectors[D1MinusLinearOSI::RESIDU_FREE];
        std::shared_ptr<siconos::algebra::SiconosVector> v = d->velocity();
        std::shared_ptr<siconos::algebra::SiconosVector> q = d->q();
        // const siconos::algebra::SiconosVector& qold = d->qMemory().getsiconos::algebra::SiconosVector(0);
        // const siconos::algebra::SiconosVector& vold = d->velocityMemory().getsiconos::algebra::SiconosVector(0); // right limit

        //residuFree->zero();
        //v->zero();
        std::shared_ptr<siconos::algebra::SiconosVector> work_tdg =  workVectors[D1MinusLinearOSI::FREE_TDG];
        assert(work_tdg);
        DEBUG_EXPR(work_tdg->display());
        residuFree =  - 0.5 * h * *work_tdg;

        DEBUG_EXPR(q->display(););
        DEBUG_EXPR(v->display(););
        d->computeForces(t, q, v);
        DEBUG_EXPR(d->forces()->display(););
        *work_tdg = *(d->forces());


        if(d->inverseMass())
        {
          d->update_inverse_mass();
          d->inverseMass()->Solve(*work_tdg);
          // contains right (left limit) acceleration without contact force
        }

        DEBUG_EXPR(work_tdg->display());
        residuFree -= 0.5 * h**work_tdg;
        DEBUG_EXPR(residuFree.display());
      }
      else if(dsType == Type::NewtonEulerDS)
      {
        auto d = std::static_pointer_cast<NewtonEulerDS> (ds);
        siconos::algebra::SiconosVector& residuFree = *workVectors[D1MinusLinearOSI::RESIDU_FREE];

        std::shared_ptr<siconos::algebra::SiconosVector> v = d->twist();
        std::shared_ptr<siconos::algebra::SiconosVector> q = d->q();
        // const siconos::algebra::SiconosVector& qold = d->qMemory().getsiconos::algebra::SiconosVector(0);
        // const siconos::algebra::SiconosVector& vold = d->twistMemory().getsiconos::algebra::SiconosVector(0); // right limit

        v->zero();
        std::shared_ptr<siconos::algebra::SiconosVector> work_tdg = workVectors[D1MinusLinearOSI::FREE_TDG];
        assert(work_tdg);
        residuFree = 0.5 * h * *work_tdg;
        work_tdg->zero();

        d->computeForces(t, q, v);
        *work_tdg += *(d->forces());

        if(d->inverseMass())
        {
          d->update_inverse_mass();
          d->inverseMass()->Solve(*work_tdg);
          // contains right (left limit) acceleration without contact force
        }

        residuFree -= 0.5 * h * *work_tdg;
        DEBUG_EXPR(residuFree.display());
      }
      else
        THROW_EXCEPTION("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " +  std::to_string(dsType));
    }
  }
  else
  {
    DEBUG_PRINT("There is no  impact in the step. indexSet1->size() = 0. _isThereImpactInTheTimeStep = false;\n");
    // -- RIGHT SIDE --
    // calculate acceleration without contact force
    for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);
      auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

      // type of the current DS
      Type::Siconos dsType = Type::value(*ds);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
      {

        auto d = std::static_pointer_cast<LagrangianDS> (ds);
        siconos::algebra::SiconosVector& vFree = *workVectors[D1MinusLinearOSI::FREE];
        vFree.zero();
        // get right state from memory
        std::shared_ptr<siconos::algebra::SiconosVector> q = d->q(); // contains position q_{k+1}
        std::shared_ptr<siconos::algebra::SiconosVector> v = d->velocity(); // contains velocity v_{k+1}^- and not free velocity

        DEBUG_EXPR(q->display());
        DEBUG_EXPR(v->display());
        // Lagrangian Nonlinear Systems
        if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
        {
          d->computeForces(t, q, v);
          DEBUG_EXPR(d->forces()->display());
          vFree += *(d->forces());

        }
        else
          THROW_EXCEPTION
          ("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " +  std::to_string(dsType));

        if(d->inverseMass())
        {
          d->update_inverse_mass();
          d->inverseMass()->Solve(vFree);
        }
        /* vFree contains right (left limit) acceleration without contact force */
        siconos::algebra::SiconosVector& residuFree = *workVectors[D1MinusLinearOSI::RESIDU_FREE];
        residuFree += -0.5 * h * vFree;
        DEBUG_EXPR(residuFree.display());
        /* Compute the right limit of the (free) velocity at  t^+_k with contact force : */
        const siconos::algebra::SiconosVector& vold = d->velocityMemory().getsiconos::algebra::SiconosVector(0);
        DEBUG_EXPR(vold.display());

        vFree = vold - residuFree;

        DEBUG_PRINT("vFree contains the right  limit of the (free) velocity at  t^-_{k+1} without contact force :\n");
        DEBUG_EXPR(vFree.display());

      }
      else if(dsType == Type::NewtonEulerDS)
      {
        auto d = std::static_pointer_cast<NewtonEulerDS> (ds);
        siconos::algebra::SiconosVector& vFree = *workVectors[D1MinusLinearOSI::FREE];
        vFree.zero();
        // get right state from memory
        std::shared_ptr<siconos::algebra::SiconosVector> q = d->q(); // contains position q_{k+1}
        std::shared_ptr<siconos::algebra::SiconosVector> v = d->twist(); // contains velocity v_{k+1}^- and not free velocity

        DEBUG_EXPR(q->display());
        DEBUG_EXPR(v->display());

        d->computeForces(t, q, v);
        vFree += *(d->forces());


        if(d->inverseMass())
        {
          d->update_inverse_mass();
          d->inverseMass()->Solve(vFree);
        }
        /* work_tdg contains right (left limit) acceleration without contact force */
        siconos::algebra::SiconosVector& residuFree = *workVectors[D1MinusLinearOSI::RESIDU_FREE];

        residuFree += -0.5 * h * vFree;

        /*  Compute the right limit of the (free) velocity at  t^+_k with contact force : */
        const siconos::algebra::SiconosVector& vold = d->twistMemory().getsiconos::algebra::SiconosVector(0);
        vFree = vold - residuFree;
        DEBUG_PRINT("vFree contains the  left limit of the (free) velocity at  t^-_{k+1} without contact force :\n");
        DEBUG_EXPR(vFree.display());
      }
      else
        THROW_EXCEPTION("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " +  std::to_string(dsType));

    }

    // solve a LCP at acceleration level only for contacts which have been active at the beginning of the time-step
    if(!allOSNS->empty())
    {
      // for (unsigned int level = _simulation->levelMinForOutput(); level < _simulation->levelMaxForOutput(); level++)
      // {
      //   _simulation->updateOutput(level);
      // }
      // _simulation->updateIndexSets();
      DEBUG_PRINT("We compute lambda^-_{k+1} \n");
      DEBUG_PRINTF("indexSet1->size() = %i\n",(int)indexSet1->size());

      _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(t, *indexSet1);

      if(_simulation->nonSmoothDynamicalSystem()->topology()->hasChanged())
      {
	for(auto osns: *allOSNS)
        {
          itOsns->setHasBeenUpdated(false);
        }
      }
      if(((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->hasInteractions()))
      {
        (*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->compute(t);
        DEBUG_EXPR((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]->display(););
      }
    }

    _simulation->nonSmoothDynamicalSystem()->updateInput(_simulation->nextTime(),2);
    for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      auto ds = _dynamicalSystemsGraph->bundle(*dsi);
      auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
      // type of the current DS
      Type::Siconos dsType = Type::value(*ds);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
      {
        auto d = std::static_pointer_cast<LagrangianDS> (ds);
        siconos::algebra::SiconosVector& residuFree = *workVectors[D1MinusLinearOSI::RESIDU_FREE];
        DEBUG_EXPR(residuFree.display(););

        if(d->p(2))
        {
          // get right state from memory
          DEBUG_EXPR(d->p(2)->display());
          siconos::algebra::SiconosVector& p2 = *d->p(2);
          siconos::algebra::SiconosVector dummy(*(d->p(2)));
          DEBUG_EXPR(p2.display());
          /* we homogenize p(2) to a force for the user output   */
          p2 *= 2.0/h;

          if(d->inverseMass())
          {
            d->update_inverse_mass();
            d->inverseMass()->Solve(dummy);
          }
          residuFree -=  dummy;

        }
        DEBUG_EXPR(residuFree.display());
      }
      else if(dsType == Type::NewtonEulerDS)
      {
        auto d = std::static_pointer_cast<NewtonEulerDS> (ds);
        siconos::algebra::SiconosVector& residuFree = *workVectors[D1MinusLinearOSI::RESIDU_FREE];

        if(d->p(2))
        {
          // get right state from memory
          siconos::algebra::SiconosVector& p2 = *d->p(2);
          siconos::algebra::SiconosVector dummy(*(d->p(2))); // value = contact force
          DEBUG_EXPR(p2.display());
          /* we homogenize p(2) to a force for the user output   */
          p2 *= 2.0/h;
          if(d->inverseMass())
          {
            d->update_inverse_mass();
            d->inverseMass()->Solve(dummy);
          }
          residuFree -=  dummy;

        }
        DEBUG_EXPR(residuFree.display());
      }
      else
        THROW_EXCEPTION("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " +  std::to_string(dsType));

      /**
       * \f[
       * \begin{cases}
       * v_{k,0} = \mbox{vold} \\
       * q_{k,0} = \mbox{qold} \\
       * F^+_{k} = \mbox{F(told,qold,vold)} \\
       * v_{k,1} = v_{k,0} + h M^{-1}_k (P^+_{2,k}+F^+_{k}) \\[2mm]
       * q_{k,1} = q_{k,0} + \frac{h}{2} (v_{k,0} + v_{k,1})  \\[2mm]
       * F^-_{k+1} = F(t^{-}_{k+1},q_{k,1},v_{k,1}) \\[2mm]
       * R_{free} = - \frac{h}{2}  M^{-1}_k (P^+_{2,k}+F^+_{k}) -  \frac{h}{2}  M^{-1}_{k+1} (P^-_{2,k+1}+F^-_{k+1}) \\[2mm]
       * \end{cases}
       * \f]
       **/

    }

  } // No impact


  DEBUG_END("D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel(), ends\n");

  return 0.; // there is no Newton iteration and the residuum is assumed to vanish
}

void D1MinusLinearOSI::computeFreeOutputHalfExplicitVelocityLevel(siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{

  DEBUG_PRINT("\n D1MinusLinearOSI::computeFreeOutputHalfExplicitVelocityLevel starts\n");
  auto allOSNS  = _simulation->oneStepNSProblems(); // all OSNSP
  auto indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  std::shared_ptr<siconos::modeling::Interaction> inter = indexSet->bundle(vertex_inter);
  auto& DSlink = inter->linkToDSVariables();
  auto& workBlockV = *indexSet->properties(vertex_inter).workBlockVectors;
  // get relation and non smooth law information
  auto relationType = inter->relation()->getType(); // relation
  auto relationSubType = inter->relation()->getSubType();
  unsigned int relativePosition = 0;
  unsigned int sizeY = inter->nonSmoothLaw()->size(); // related NSL

  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[3] = 0;
  coord[4] = 0;
  coord[5] = 0;
  coord[6] = 0;
  coord[7] = sizeY;
  std::shared_ptr<siconos::algebra::SiconosMatrix> C; // Jacobian of Relation with respect to degree of freedom
  std::shared_ptr<siconos::algebra::BlockVector> Xfree; // free degree of freedom



  siconos::algebra::SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[D1MinusLinearOSI::OSNSP_RHS];


  DEBUG_PRINT("osnsp_rhs before\n");
  DEBUG_EXPR(osnsp_rhs.display(););
  // define Xfree for velocity and acceleration level
  if(((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
  {
    /* get the current velocity  of the aggregated ds */
    if(relationType == siconos::modeling::RelationType::Lagrangian)
    {
      Xfree = DSlink[siconos::modeling::LagrangianR::q1];
      DEBUG_PRINT("Xfree = DSlink[siconos::modeling::LagrangianR::q1];\n");
    }
    else if(relationType == NewtonEuler)
    {
      Xfree = DSlink[NewtonEulerR::velocity];
      DEBUG_PRINT("Xfree = DSlink[NewtonEulerR::velocity];\n");
    }
    else
      THROW_EXCEPTION("D1MinusLinearOSI::computeFreeOutput - unknown relation type.");
    DEBUG_PRINT("Xfree contains the current velocity of the aggregated ds];\n");
    DEBUG_EXPR(Xfree->display());

  }
  else if(((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
  {
    /* get the free velocity of the aggregated ds */
    if(relationType == siconos::modeling::RelationType::Lagrangian)
    {
      Xfree = workBlockV[D1MinusLinearOSI::xfree];
    }
    else if(relationType == NewtonEuler)
    {
      Xfree = workBlockV[D1MinusLinearOSI::xfree];
    }
    else
      THROW_EXCEPTION("D1MinusLinearOSI::computeFreeOutput - unknown relation type.");
    assert(Xfree);
    DEBUG_PRINT("Xfree = DSlink[Lagrangian/NewtonEulerR::xfree];\n");
    DEBUG_PRINT("Xfree contains the free velocity of the aggregated ds];\n");
    DEBUG_EXPR(Xfree->display());

  }
  else
    THROW_EXCEPTION("D1MinusLinearOSI::computeFreeOutput - OSNSP neither on velocity nor on acceleration level.");

  // calculate data of interaction
  std::shared_ptr<siconos::modeling::Interaction> mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());

  // only Lagrangian Systems
  if(relationType == siconos::modeling::RelationType::Lagrangian)
  {
    // in osnsp_rhs the linear part of velocity or acceleration relation will be saved
    C = std::static_pointer_cast<LagrangianR>(mainInteraction->relation())->C();

    DEBUG_EXPR(C->display(););

    if(C)
    {
      assert(Xfree);
      coord[3] = C->size(1);
      coord[5] = C->size(1);
      siconos::algebra::subprod(*C, *Xfree, osnsp_rhs, coord, true);
    }
    DEBUG_EXPR(osnsp_rhs.display(););
    /*  explicit time dependence -> partial time derivative has to be added */
    if(relationSubType == RheonomousR)
    {
      THROW_EXCEPTION("D1MinusLinearOSI::computeFreeOutput is not implemented  at velocity level for LagrangianRheonomousR.");
    }

    /* add the contribution due to the coefficient of restitution*/
    if(((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
    {
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter, indexSet->properties(vertex_inter)));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }

    if(((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
    {
      /*Do nothing*/
    }
    DEBUG_EXPR(osnsp_rhs.display(););
  }
  /*Newton-Euler */
  else if(relationType == NewtonEuler)
  {
    auto CT =  std::static_pointer_cast<NewtonEulerR>(mainInteraction->relation())->jachqT();
    DEBUG_EXPR(CT->display());
    if(CT)
    {
      coord[3] = CT->size(1);
      coord[5] = CT->size(1);
      assert(Xfree);
      // creates a POINTER link between workX[ds] (xfree) and the
      // corresponding interactionBlock in each Interaction for each ds of the
      // current Interaction.
      // XXX Big quirks !!! -- xhub
      siconos::algebra::subprod(*CT, *Xfree, osnsp_rhs, coord, true);
    }

    /* add the contribution due to the coefficient of restitution*/
    if(((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
    {
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter, indexSet->properties(vertex_inter)));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }

    if(((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
    {
      /* Do nothing*/
    }
    DEBUG_EXPR(osnsp_rhs.display(););
  }



  else
    THROW_EXCEPTION("D1MinusLinearOSI::computeFreeOutput - not implemented for Relation of type " +  std::to_string(relationType));

  DEBUG_EXPR(osnsp_rhs.display(););



  DEBUG_PRINT("D1MinusLinearOSI::computeFreeOutputHalfExplicitVelocityLevel ends\n");

}


bool D1MinusLinearOSI::addInteractionInIndexSetHalfExplicitVelocityLevel(std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)
{
  DEBUG_PRINT("D1MinusLinearOSI::addInteractionInIndexSetHalfExplicitVelocityLevel.\n");

  if(i == 1)
  {
    DEBUG_PRINT(" level == 1\n");
    if(Type::value(*(inter->nonSmoothLaw())) == Type::EqualityConditionNSL)
    {
      return true;
    }
    /* Active for impulses calculations? Contacts that are closed */
    double y = (*(inter->y(0)))(0); // current position
    DEBUG_PRINTF("y= %18.14e\n", y);
    return (y <= DEFAULT_TOL_D1MINUS);

  }
  else if(i == 2)
  {
    if(Type::value(*(inter->nonSmoothLaw())) == Type::EqualityConditionNSL)
    {
      return false;
    }
    /*  Contacts which have been closing in the last time step */
    DEBUG_PRINT(" level == 2\n");
    double y = (*(inter->y(0)))(0); // current position
    double y_k = (inter->y_k(0))(0); // old position
    DEBUG_PRINTF("y= %18.14e\n", y);
    DEBUG_PRINTF("y_k= %18.14e\n", y_k);
    /* if Interaction has not been active in the previous calculation
       and now becomes active */
    return (y <= DEFAULT_TOL_D1MINUS && y_k > DEFAULT_TOL_D1MINUS);
  }
  else
    THROW_EXCEPTION("D1MinusLinearOSI::addInteractionInIndexSetHalfExplicitVelocityLevel, IndexSet[i > 2] does not exist.");
  return false;
}

bool D1MinusLinearOSI::removeInteractionFromIndexSetHalfExplicitVelocityLevel(std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i)
{
  DEBUG_PRINT("D1MinusLinearOSI::removeInteractionFromIndexSetHalfExplicitVelocityLevel.\n");

  return !(addInteractionInIndexSetHalfExplicitVelocityLevel(inter,i));
}
