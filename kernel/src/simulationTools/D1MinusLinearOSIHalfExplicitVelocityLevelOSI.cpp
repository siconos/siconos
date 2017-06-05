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

#include "D1MinusLinearOSI.hpp"
#include "Simulation.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianRheonomousR.hpp"
#include "LagrangianScleronomousR.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonImpactNSL.hpp"
#include "BlockVector.hpp"
#include "CxxStd.hpp"
#include "Topology.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OneStepNSProblem.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

/// @cond
using namespace RELATION;

double D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel()
{
  DEBUG_BEGIN("D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel(), starts\n");

  double t = _simulation->nextTime(); // end of the time step
  double told = _simulation->startingTime(); // beginning of the time step
  double h = _simulation->timeStep(); // time step length
  SP::OneStepNSProblems allOSNS  = _simulation->oneStepNSProblems(); // all OSNSP
  SP::Topology topo =  _simulation->nonSmoothDynamicalSystem()->topology();
  SP::InteractionsGraph indexSet1 = topo->indexSet(1);

  /******************************************************************************************
   *  Step 1-  solve a LCP at velocity level for lambda^+_{k} for the last set indices
   *   if index1 is empty we should skip this step
   ******************************************************************************************/

  DEBUG_PRINT("\nEVALUATE LEFT HAND SIDE\n");

  // DEBUG_EXPR(std::cout<< "allOSNS->empty()   " << std::boolalpha << allOSNS->empty() << std::endl << std::endl);
  // DEBUG_EXPR(std::cout<< "allOSNS->size()   "  << allOSNS->size() << std::endl << std::endl);

// -- LEFT SIDE --

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    Type::Siconos dsType = Type::value(*ds);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    SP::SiconosVector work_tdg;

    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
      SiconosVector& vFree = *workVectors[OneStepIntegrator::free];
      /* POINTER CONSTRUCTOR : vFree will contain the velocity without contact force */
      vFree.zero();

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit

      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());

      /* compute the force and store in vFree */
      d->computeForces(told, qold, vold);
      DEBUG_EXPR(d->forces()->display());
      vFree += *(d->forces());

      /* Compute the acceleration due to the external force */
      /* vFree contains left (right limit) acceleration without contact force */
      if(d->inverseMass())
      {
        d->update_inverse_mass();
        d->inverseMass()->PLUForwardBackwardInPlace(vFree);
      }

      /* Store the value of vFree in d->workspace(DynamicalSystem::free_tdg called work_tdg*/
      work_tdg =  workVectors[OneStepIntegrator::free_tdg];
      work_tdg->zero();
      *work_tdg = vFree;


      /*Compute the right limit of the (free) velocity at  t^+_k with contact force :  */
      vFree  *= h ;
      vFree += *vold;

      /* Compute a prediction for q  that will serve for computing new values
         of the Jacobian of the constraints */

      SP::SiconosVector q = d->q();

      * q = * qold + h* *vold ;
      DEBUG_PRINT("vFree contains the right limit of the (free) velocity at  t^+_k with contact force :\n");
      DEBUG_EXPR(vFree.display());
      DEBUG_PRINT("work_tdg contains left (right limit) acceleration without contact forcework :\n");
      DEBUG_EXPR(work_tdg->display());
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
      SiconosVector& vFree = *workVectors[OneStepIntegrator::free];
      vFree.zero();

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->twistMemory()->getSiconosVector(0); // right limit
      DEBUG_EXPR(vFree.display());
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());

      d->computeForces(told, qold, vold);
      DEBUG_EXPR(d->forces()->display());
      vFree += *(d->forces());

      if(d->inverseMass())
      {
        d->update_inverse_mass();
        d->inverseMass()->PLUForwardBackwardInPlace(vFree); // contains left (right limit) acceleration without contact force
      }

      work_tdg =  workVectors[OneStepIntegrator::free_tdg];;
      work_tdg->zero();
      DEBUG_EXPR(work_tdg->display());
      *work_tdg = vFree;

      /*Compute the right limit of the (free) velocity at  t^+_k with contact force :  */
      vFree *= h ;
      vFree += *vold;

      // * q = * qold + h* *vold ; to be written consistently for Newton Euler DS
      DEBUG_PRINT("vFree contains the right limit of the (free) velocity at  t^+_k with contact force :\n");
      DEBUG_EXPR(vFree.display());
      DEBUG_PRINT("work_tdg contains left (right limit) acceleration without contact forcework :\n");
      DEBUG_EXPR(work_tdg->display());
    }
    else
    {
      RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    }
  }


  if(!allOSNS->empty())
  {
    if(indexSet1->size() >0)
    {
      InteractionsGraph::VIterator ui, uiend;
      SP::Interaction inter;
      for(std11::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui)
      {
        inter = indexSet1->bundle(*ui);
        inter->relation()->computeJach(t, *inter, indexSet1->properties(*ui));
        inter->relation()->computeJacg(told, *inter, indexSet1->properties(*ui));
      }

      if(_simulation->nonSmoothDynamicalSystem()->topology()->hasChanged())
      {
        for(OSNSIterator itOsns = allOSNS->begin(); itOsns != allOSNS->end(); ++itOsns)
        {
          (*itOsns)->setHasBeenUpdated(false);
        }
      }
      assert((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]);

      if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->hasInteractions()))  /* it should be equivalent to indexSet1 */
      {
        DEBUG_PRINT("We compute lambda^+_{k} \n");
        (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->compute(told);
        DEBUG_EXPR((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->display());
      }

    }
  }

  // Note Franck : at the time this results in a call to swapInMem of all Interactions of the NSDS
  // So let the simu do this.
  //(*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->saveInMemory(); // we push y and lambda in Memories
  _simulation->nonSmoothDynamicalSystem()->pushInteractionsInMemory();
  _simulation->nonSmoothDynamicalSystem()->updateInput(_simulation->nextTime(),2);

  /**************************************************************************************************************
   *  Step 2 -  compute v_{k,1}
   **************************************************************************************************************/


  DEBUG_PRINT("\n PREDICT RIGHT HAND SIDE\n");

  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    // type of the current DS
    Type::Siconos dsType = Type::value(*ds);
    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
      // contains residu without nonsmooth effect
      SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
      SiconosVector& vFree = *workVectors[OneStepIntegrator::free];
      SP::SiconosVector work_tdg = workVectors[OneStepIntegrator::free_tdg];



      // initialize *it->residuFree and predicted right velocity (left limit)
      SP::SiconosVector v = d->velocity(); //contains velocity v_{k+1}^- and not free velocity
      v->zero();




      SP::SiconosVector p2 = d->p(2);
      SP::SiconosVector dummy(new SiconosVector(*p2)); // value = contact force
      DEBUG_EXPR(p2->display());
      /* we homogenize p(2) to a force for the user output   */
      *p2 /= h;
      if(d->inverseMass())
      {
        d->inverseMass()->PLUForwardBackwardInPlace(*dummy);
        DEBUG_EXPR(d->inverseMass()->display(););
      }
      DEBUG_EXPR(vFree.display());
      DEBUG_EXPR(dummy->display());

      *v = vFree + *dummy;
      DEBUG_EXPR(v->display());


      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());

      SP::SiconosVector q = d->q(); // POINTER CONSTRUCTOR : contains position q_{k+1}
      *q = *qold;

      scal(0.5 * h, *vold + *v, *q, false);
      DEBUG_EXPR(q->display());


      residuFree.zero();
      residuFree -= 0.5 * (h * *work_tdg) + 0.5* *dummy;
      DEBUG_EXPR(residuFree.display());
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
      SiconosVector& vFree = *workVectors[OneStepIntegrator::free];
      SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];// contains residu without nonsmooth effect
      SP::SiconosVector work_tdg = workVectors[OneStepIntegrator::free_tdg];

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->twistMemory()->getSiconosVector(0);

      // initialize *it->residuFree and predicted right velocity (left limit)
      SP::SiconosVector v = d->twist(); //contains velocity v_{k+1}^- and not free velocity

      v->zero();


      SP::SiconosVector p2 = d->p(2);
      SP::SiconosVector dummy(new SiconosVector(*p2)); // value = contact force
      DEBUG_EXPR(p2->display());
      /* we homogenize p(2) to a force for the user output   */
      *p2 /= h;

      if(d->inverseMass())
      {
        d->update_inverse_mass();
        d->inverseMass()->PLUForwardBackwardInPlace(*dummy);
      }

      DEBUG_EXPR(vFree.display());
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());

      *v = vFree +  *dummy;


      DEBUG_EXPR(v->display());

      //first step consists in computing  \dot q.
      //second step consists in updating q.
      //
      SP::SiconosMatrix T = d->T();
      SP::SiconosVector dotq = d->dotq();
      prod(*T, *v, *dotq, true);

      SP::SiconosVector dotqold = d->dotqMemory()->getSiconosVector(0);

      SP::SiconosVector q = d->q(); // POINTER CONSTRUCTOR : contains position q_{k+1}
      *q = *qold;

      scal(0.5 * h, *dotqold + *dotq, *q, false);
      DEBUG_PRINT("new q before normalizing\n");
      DEBUG_EXPR(q->display());
      //q[3:6] must be normalized
      d->normalizeq();
      d->computeT();
      DEBUG_PRINT("new q after normalizing\n");
      DEBUG_EXPR(q->display());


      residuFree.zero();
      residuFree -= 0.5 * (h* *work_tdg) + 0.5 * *dummy;
      DEBUG_EXPR(residuFree.display());



    }
    else
      RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);


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

    SP::Topology topo =  _simulation->nonSmoothDynamicalSystem()->topology();
    SP::InteractionsGraph indexSet2 = topo->indexSet(2);

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
    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
      VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

      // type of the current DS
      Type::Siconos dsType = Type::value(*ds);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
      {
        SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
        SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
        SP::SiconosVector v = d->velocity();
        SP::SiconosVector q = d->q();
        SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
        SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit

        //residuFree->zero();
        //v->zero();
        SP::SiconosVector work_tdg =  workVectors[OneStepIntegrator::free_tdg];
        assert(work_tdg);
        DEBUG_EXPR(work_tdg->display());
        residuFree =  - 0.5 * h**work_tdg;

        DEBUG_EXPR(q->display(););
        DEBUG_EXPR(v->display(););
        d->computeForces(t, q, v);
        DEBUG_EXPR(d->forces()->display(););
        *work_tdg = *(d->forces());


        if(d->inverseMass())
        {
          d->update_inverse_mass();
          d->inverseMass()->PLUForwardBackwardInPlace(*work_tdg);
          // contains right (left limit) acceleration without contact force
        }

        DEBUG_EXPR(work_tdg->display());
        residuFree -= 0.5 * h**work_tdg;
        DEBUG_EXPR(residuFree.display());
      }
      else if(dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
        SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];

        SP::SiconosVector v = d->twist();
        SP::SiconosVector q = d->q();
        SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
        SP::SiconosVector vold = d->twistMemory()->getSiconosVector(0); // right limit

        v->zero();
        SP::SiconosVector work_tdg = workVectors[OneStepIntegrator::free_tdg];
        assert(work_tdg);
        residuFree = 0.5 * h**work_tdg;
        work_tdg->zero();

        d->computeForces(t, q, v);
        *work_tdg += *(d->forces());

        if(d->inverseMass())
        {
          d->update_inverse_mass();
          d->inverseMass()->PLUForwardBackwardInPlace(*work_tdg);
          // contains right (left limit) acceleration without contact force
        }

        residuFree -= 0.5 * h**work_tdg;
        DEBUG_EXPR(residuFree.display());
      }
      else
        RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    }
  }
  else
  {
    DEBUG_PRINT("There is no  impact in the step. indexSet1->size() = 0. _isThereImpactInTheTimeStep = false;\n");
    // -- RIGHT SIDE --
    // calculate acceleration without contact force
    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
      VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

      // type of the current DS
      Type::Siconos dsType = Type::value(*ds);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
      {

        SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
        SiconosVector& vFree = *workVectors[OneStepIntegrator::free];
        vFree.zero();
        // get right state from memory
        SP::SiconosVector q = d->q(); // contains position q_{k+1}
        SP::SiconosVector v = d->velocity(); // contains velocity v_{k+1}^- and not free velocity

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
          RuntimeException::selfThrow
            ("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

        if(d->inverseMass())
        {
          d->update_inverse_mass();
          d->inverseMass()->PLUForwardBackwardInPlace(vFree);
        }
        /* vFree contains right (left limit) acceleration without contact force */
        SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
        residuFree += -0.5 * h * vFree;
        DEBUG_EXPR(residuFree.display());
        /* Compute the right limit of the (free) velocity at  t^+_k with contact force : */
        SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);
        DEBUG_EXPR(vold->display());

        vFree = *vold - residuFree;

        DEBUG_PRINT("vFree contains the right  limit of the (free) velocity at  t^-_{k+1} without contact force :\n");
        DEBUG_EXPR(vFree.display());

      }
      else if(dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
        SiconosVector& vFree = *workVectors[OneStepIntegrator::free];
        vFree.zero();
        // get right state from memory
        SP::SiconosVector q = d->q(); // contains position q_{k+1}
        SP::SiconosVector v = d->twist(); // contains velocity v_{k+1}^- and not free velocity

        DEBUG_EXPR(q->display());
        DEBUG_EXPR(v->display());

        d->computeForces(t, q, v);
        vFree += *(d->forces());


        if(d->inverseMass())
        {
          d->update_inverse_mass();
          d->inverseMass()->PLUForwardBackwardInPlace(vFree);
        }
        /* work_tdg contains right (left limit) acceleration without contact force */
        SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];

        residuFree += -0.5 * h * vFree;

        /*  Compute the right limit of the (free) velocity at  t^+_k with contact force : */
        SP::SiconosVector vold = d->twistMemory()->getSiconosVector(0);
        vFree = *vold - residuFree;
        DEBUG_PRINT("vFree contains the  left limit of the (free) velocity at  t^-_{k+1} without contact force :\n");
        DEBUG_EXPR(vFree.display());
      }
      else
        RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

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
      InteractionsGraph::VIterator ui, uiend;
      SP::Interaction inter;
      for(std11::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui)
      {
        inter = indexSet1->bundle(*ui);
        inter->relation()->computeJach(t, *inter, indexSet1->properties(*ui));
        inter->relation()->computeJacg(t, *inter, indexSet1->properties(*ui));
      }
      if(_simulation->nonSmoothDynamicalSystem()->topology()->hasChanged())
      {
        for(OSNSIterator itOsns = allOSNS->begin(); itOsns != allOSNS->end(); ++itOsns)
        {
          (*itOsns)->setHasBeenUpdated(false);
        }
      }
      if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->hasInteractions()))
      {
        (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->compute(t);
        DEBUG_EXPR((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->display(););
      }
    }

    _simulation->nonSmoothDynamicalSystem()->updateInput(_simulation->nextTime(),2);
    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
      VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
      // type of the current DS
      Type::Siconos dsType = Type::value(*ds);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
      {
        SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
        SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
        DEBUG_EXPR(residuFree.display(););

        if(d->p(2))
        {

          // get right state from memory
          DEBUG_EXPR(d->p(2)->display());
          SP::SiconosVector p2 = d->p(2);
          SiconosVector dummy(*(d->p(2)));
          DEBUG_EXPR(p2->display());
          /* we homogenize p(2) to a force for the user output   */
          *p2 *= 2.0/h;

          if(d->inverseMass())
          {
            d->update_inverse_mass();
            d->inverseMass()->PLUForwardBackwardInPlace(dummy);
          }
          residuFree -=  dummy;

        }
        DEBUG_EXPR(residuFree.display());
      }
      else if(dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
        SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];

        if(d->p(2))
        {
          // get right state from memory
          SP::SiconosVector p2 = d->p(2);
          SiconosVector dummy(*(d->p(2))); // value = contact force
          DEBUG_EXPR(p2->display());
          /* we homogenize p(2) to a force for the user output   */
          *p2 *= 2.0/h;
          if(d->inverseMass())
          {
            d->update_inverse_mass();
            d->inverseMass()->PLUForwardBackwardInPlace(dummy);
          }
          residuFree -=  dummy;

        }
        DEBUG_EXPR(residuFree.display());
      }
      else
        RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

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

void D1MinusLinearOSI::computeFreeOutputHalfExplicitVelocityLevel(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{

  DEBUG_PRINT("\n D1MinusLinearOSI::computeFreeOutputHalfExplicitVelocityLevel starts\n");
  SP::OneStepNSProblems allOSNS  = _simulation->oneStepNSProblems(); // all OSNSP
  SP::InteractionsGraph indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);
  VectorOfBlockVectors& DSlink = *indexSet->properties(vertex_inter).DSlink;
  // get relation and non smooth law information
  RELATION::TYPES relationType = inter->relation()->getType(); // relation
  RELATION::SUBTYPES relationSubType = inter->relation()->getSubType();
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
  SP::SiconosMatrix C; // Jacobian of Relation with respect to degree of freedom
  SP::BlockVector Xfree; // free degree of freedom



  SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[D1MinusLinearOSI::OSNSP_RHS];


  DEBUG_PRINT("osnsp_rhs before\n");
  DEBUG_EXPR(osnsp_rhs.display(););
  // define Xfree for velocity and acceleration level
  if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
  {
    /* get the current velocity  of the aggregated ds */
    if(relationType == Lagrangian)
    {
      Xfree = DSlink[LagrangianR::q1];
      DEBUG_PRINT("Xfree = DSlink[LagrangianR::q1];\n");
    }
    else if(relationType == NewtonEuler)
    {
      Xfree = DSlink[NewtonEulerR::velocity];
      DEBUG_PRINT("Xfree = DSlink[NewtonEulerR::velocity];\n");
    }
    else
      RuntimeException::selfThrow("D1MinusLinearOSI::computeFreeOutput - unknown relation type.");
    DEBUG_PRINT("Xfree contains the current velocity of the aggregated ds];\n");
    DEBUG_EXPR(Xfree->display());

  }
  else if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
  {
    /* get the free velocity of the aggregated ds */
    if(relationType == Lagrangian)
    {
      Xfree = DSlink[LagrangianR::xfree];
    }
    else if(relationType == NewtonEuler)
    {
      Xfree = DSlink[NewtonEulerR::xfree];
    }
    else
      RuntimeException::selfThrow("D1MinusLinearOSI::computeFreeOutput - unknown relation type.");
    assert(Xfree);
    DEBUG_PRINT("Xfree = DSlink[Lagrangian/NewtonEulerR::xfree];\n");
    DEBUG_PRINT("Xfree contains the free velocity of the aggregated ds];\n");
    DEBUG_EXPR(Xfree->display());

  }
  else
    RuntimeException::selfThrow("D1MinusLinearOSI::computeFreeOutput - OSNSP neither on velocity nor on acceleration level.");

  // calculate data of interaction
  SP::Interaction mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());

  // only Lagrangian Systems
  if(relationType == Lagrangian)
  {
    // in osnsp_rhs the linear part of velocity or acceleration relation will be saved
    C = std11::static_pointer_cast<LagrangianR>(mainInteraction->relation())->C();

    DEBUG_EXPR(C->display(););

    if(C)
    {
      assert(Xfree);
      coord[3] = C->size(1);
      coord[5] = C->size(1);
      subprod(*C, *Xfree, osnsp_rhs, coord, true);
    }
    DEBUG_EXPR(osnsp_rhs.display(););
    /*  explicit time dependence -> partial time derivative has to be added */
    if(relationSubType == RheonomousR)
    {
      RuntimeException::selfThrow("D1MinusLinearOSI::computeFreeOutput is not implemented  at velocity level for LagrangianRheonomousR.");
    }

    /* add the contribution due to the coefficient of restitution*/
    if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
    {
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter, indexSet->properties(vertex_inter)));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }

    if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
    {
      /*Do nothing*/
    }
    DEBUG_EXPR(osnsp_rhs.display(););
  }
  /*Newton-Euler */
  else if(relationType == NewtonEuler)
  {
    SP::SiconosMatrix CT =  std11::static_pointer_cast<NewtonEulerR>(mainInteraction->relation())->jachqT();
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
      subprod(*CT, *Xfree, osnsp_rhs, coord, true);
    }

    /* add the contribution due to the coefficient of restitution*/
    if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
    {
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter, indexSet->properties(vertex_inter)));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }

    if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
    {
      /* Do nothing*/
    }
    DEBUG_EXPR(osnsp_rhs.display(););
  }



  else
    RuntimeException::selfThrow("D1MinusLinearOSI::computeFreeOutput - not implemented for Relation of type " + relationType);

  DEBUG_EXPR(osnsp_rhs.display(););



  DEBUG_PRINT("D1MinusLinearOSI::computeFreeOutputHalfExplicitVelocityLevel ends\n");

}


bool D1MinusLinearOSI::addInteractionInIndexSetHalfExplicitVelocityLevel(SP::Interaction inter, unsigned int i)
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
    double yOld = (*(inter->yOld(0)))(0); // old position
    DEBUG_PRINTF("y= %18.14e\n", y);
    DEBUG_PRINTF("yOld= %18.14e\n", yOld);
    /* if Interaction has not been active in the previous calculation
       and now becomes active */
    return (y <= DEFAULT_TOL_D1MINUS && yOld > DEFAULT_TOL_D1MINUS);
  }
  else
    RuntimeException::selfThrow("D1MinusLinearOSI::addInteractionInIndexSetHalfExplicitVelocityLevel, IndexSet[i > 2] does not exist.");
  return false;
}

bool D1MinusLinearOSI::removeInteractionInIndexSetHalfExplicitVelocityLevel(SP::Interaction inter, unsigned int i)
{
  DEBUG_PRINT("D1MinusLinearOSI::removeInteractionInIndexSetHalfExplicitVelocityLevel.\n");

  return !(addInteractionInIndexSetHalfExplicitVelocityLevel(inter,i));
}
