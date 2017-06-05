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

double D1MinusLinearOSI::computeResiduHalfExplicitAccelerationLevel()
{
  DEBUG_BEGIN("\n D1MinusLinearOSI::computeResiduHalfExplicitAccelerationLevel()\n");

  double t = _simulation->nextTime(); // end of the time step
  double told = _simulation->startingTime(); // beginning of the time step
  double h = _simulation->timeStep(); // time step length
  SP::OneStepNSProblems allOSNS  = _simulation->oneStepNSProblems(); // all OSNSP
  SP::Topology topo =  _simulation->nonSmoothDynamicalSystem()->topology();
  SP::InteractionsGraph indexSet2 = topo->indexSet(2);

  /**************************************************************************************************************
   *  Step 1-  solve a LCP at acceleration level for lambda^+_{k} for the last set indices
   *   if index2 is empty we should skip this step
   **************************************************************************************************************/

  DEBUG_PRINT("\nEVALUATE LEFT HAND SIDE\n");

  DEBUG_EXPR(std::cout<< "allOSNS->empty()   " << std::boolalpha << allOSNS->empty() << std::endl << std::endl);
  DEBUG_EXPR(std::cout<< "allOSNS->size()   "  << allOSNS->size() << std::endl << std::endl);

// -- LEFT SIDE --

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    Type::Siconos dsType = Type::value(*ds);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    DEBUG_EXPR(ds->display());
    SP::SiconosVector work_tdg;

    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
      SiconosVector& accFree = *workVectors[OneStepIntegrator::free];
      /* POINTER CONSTRUCTOR : will contain the acceleration without contact force */
      accFree.zero();

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit

      DEBUG_EXPR(accFree.display());
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());

      /* compute the force and store in accFree */
      d->computeForces(told, qold, vold);
      DEBUG_EXPR(d->forces()->display());
      accFree += *(d->forces());

      /* Compute the acceleration due to the external force */
      /* accFree contains left (right limit) acceleration without contact force */
      if(d->inverseMass())
	{
	  d->update_inverse_mass();
	  d->inverseMass()->PLUForwardBackwardInPlace(accFree);
	}

      /* Store the value of accFree in d->workspace(DynamicalSystem::free_tdg called work_tdg*/
      work_tdg =  workVectors[OneStepIntegrator::free_tdg];
      work_tdg->zero();
      *work_tdg = accFree; // store the value in WorkFreeFree

      //d->addWorkVector(accFree,DynamicalSystem::free_tdg); // store the value in WorkFreeFree
      DEBUG_PRINT("accFree contains right limit acceleration at  t^+_k with contact force :\n");
      DEBUG_EXPR(accFree.display());
      DEBUG_PRINT("work_tdg contains right limit acceleration at t^+_k without contact force :\n");
      DEBUG_EXPR(work_tdg->display());
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
      SiconosVector& accFree = *workVectors[OneStepIntegrator::free];

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->twistMemory()->getSiconosVector(0); // right limit

      DEBUG_EXPR(accFree.display());
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());

      work_tdg =  workVectors[OneStepIntegrator::free_tdg];
      work_tdg->zero();
      DEBUG_EXPR(work_tdg->display());

      d->computeForces(told, qold, vold);
      DEBUG_EXPR(d->forces()->display());

      accFree += *(d->forces());

      if(d->inverseMass())
	{
	  d->update_inverse_mass();
	  d->inverseMass()->PLUForwardBackwardInPlace(accFree); // contains left (right limit) acceleration without contact force
	}
      *work_tdg = accFree; // store the value in WorkFreeFree

      DEBUG_PRINT("accFree contains right limit acceleration at  t^+_k with contact force :\n");
      DEBUG_EXPR(accFree.display());
      DEBUG_PRINT("work_tdg contains right limit acceleration at t^+_k without contact force :\n");
      DEBUG_EXPR(work_tdg->display());
    }
    else
    {
      RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    }
  }


  if(!allOSNS->empty())
  {
    if(indexSet2->size() >0)
    {
      InteractionsGraph::VIterator ui, uiend;
      SP::Interaction inter;
      for(std11::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui)
      {
        inter = indexSet2->bundle(*ui);
        inter->relation()->computeJach(t, *inter, indexSet2->properties(*ui));
        inter->relation()->computeJacg(told, *inter, indexSet2->properties(*ui));
      }

      if(_simulation->nonSmoothDynamicalSystem()->topology()->hasChanged())
      {
        for(OSNSIterator itOsns = allOSNS->begin(); itOsns != allOSNS->end(); ++itOsns)
        {
          (*itOsns)->setHasBeenUpdated(false);
        }
      }
      assert((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]);

      if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->hasInteractions()))  // it should be equivalent to indexSet2
      {
        DEBUG_PRINT("We compute lambda^+_{k} \n");
        (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->compute(told);
        DEBUG_EXPR((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->display());
      }


      // Note Franck : at the time this results in a call to swapInMem of all Interactions of the NSDS
      // So let the simu do this.
      //(*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->saveInMemory(); // we push y and lambda in Memories
      _simulation->nonSmoothDynamicalSystem()->pushInteractionsInMemory();
      _simulation->nonSmoothDynamicalSystem()->updateInput(_simulation->nextTime(),2);


      for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
      {
        if(!checkOSI(dsi)) continue;
        SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

        Type::Siconos dsType = Type::value(*ds);
        VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

        if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
        {
          SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
          SiconosVector& accFree = *workVectors[OneStepIntegrator::free];

          SP::SiconosVector dummy(new SiconosVector(*(d->p(2)))); // value = contact force
          if(d->inverseMass())
          {
            d->update_inverse_mass();
            d->inverseMass()->PLUForwardBackwardInPlace(*dummy);
          }
          accFree  += *(dummy);

          DEBUG_EXPR(d->p(2)->display());
        }
        else if(dsType == Type::NewtonEulerDS)
        {
          SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
          SiconosVector& accFree = *workVectors[OneStepIntegrator::free];


          SP::SiconosVector dummy(new SiconosVector(*(d->p(2)))); // value = contact force
          if(d->inverseMass())
          {
            d->update_inverse_mass();
            d->inverseMass()->PLUForwardBackwardInPlace(*dummy);
          }
          accFree  += *(dummy);

          DEBUG_EXPR(d->p(2)->display());

        }
        else
          RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

      }
    }
  }

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
      SiconosVector& accFree = *workVectors[OneStepIntegrator::free];

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      // initialize *it->residuFree and predicted right velocity (left limit)
      SP::SiconosVector v = d->velocity(); //contains velocity v_{k+1}^- and not free velocity
      residuFree.zero();
      v->zero();

      DEBUG_EXPR(accFree.display());
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());


      residuFree -= 0.5 * h* accFree;

      *v += h* accFree;
      *v += *vold;

      DEBUG_EXPR(residuFree.display());
      DEBUG_EXPR(v->display());

      SP::SiconosVector q = d->q(); // POINTER CONSTRUCTOR : contains position q_{k+1}
      *q = *qold;

      scal(0.5 * h, *vold + *v, *q, false);
      DEBUG_EXPR(q->display());
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
      SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];// contains residu without nonsmooth effect
      SiconosVector& accFree = *workVectors[OneStepIntegrator::free];
      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->twistMemory()->getSiconosVector(0);

      // initialize *it->residuFree and predicted right velocity (left limit)
      SP::SiconosVector v = d->twist(); //contains velocity v_{k+1}^- and not free velocity
      residuFree.zero();
      v->zero();

      DEBUG_EXPR(accFree.display());
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());


      residuFree -= 0.5 * h* accFree;

      *v += h* accFree;
      *v += *vold;

      DEBUG_EXPR(residuFree.display());
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
  if(!allOSNS->empty())
  {

    for(unsigned int level = levelMinForOutput();
        level < levelMaxForOutput(); level++)
    {
      _simulation->nonSmoothDynamicalSystem()->updateOutput(_simulation->nextTime(),level);
    }
    _simulation->updateIndexSets();

    SP::Topology topo =  _simulation->nonSmoothDynamicalSystem()->topology();
    SP::InteractionsGraph indexSet3 = topo->indexSet(3);

    if(indexSet3->size() > 0)
    {
      _isThereImpactInTheTimeStep = true;
      DEBUG_PRINT("There is an impact in the step. indexSet3->size() > 0. _isThereImpactInTheTimeStep = true;\n");
    }
    else
    {
      _isThereImpactInTheTimeStep = false;
      DEBUG_PRINT("There is no  impact in the step. indexSet3->size() = 0. _isThereImpactInTheTimeStep = false;\n");
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

    DEBUG_PRINT("There is an impact in the step. indexSet3->size() > 0.  _isThereImpactInTheTimeStep = true\n");
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
        //residuFree.zero();
        //v->zero();

        SP::SiconosVector work_tdg = workVectors[OneStepIntegrator::free_tdg];
        assert(work_tdg);
        residuFree =  - 0.5 * h* *work_tdg;


	d->computeForces(t, q, v);
	*work_tdg = *(d->forces());
	DEBUG_EXPR(d->forces()->display());

	if(d->inverseMass())
	  {
	    d->update_inverse_mass();
	    d->inverseMass()->PLUForwardBackwardInPlace(*work_tdg);
	  }
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

        //residuFree.zero();
        v->zero();
        SP::SiconosVector work_tdg = workVectors[OneStepIntegrator::free_tdg];
        assert(work_tdg);
        residuFree = 0.5 * h* *work_tdg;
        work_tdg->zero();

	d->computeForces(t, q, v);
	*work_tdg += *(d->forces());


	if(d->inverseMass())
	  {
	    d->update_inverse_mass();
	    d->inverseMass()->PLUForwardBackwardInPlace(*work_tdg);
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
    DEBUG_PRINT("There is no  impact in the step. indexSet3->size() = 0. _isThereImpactInTheTimeStep = false;\n");
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
        SiconosVector& accFree = *workVectors[OneStepIntegrator::free];
        accFree.zero();
        // get right state from memory
        SP::SiconosVector q = d->q(); // contains position q_{k+1}
        SP::SiconosVector v = d->velocity(); // contains velocity v_{k+1}^- and not free velocity
        DEBUG_EXPR(accFree.display());
        DEBUG_EXPR(q->display());
        DEBUG_EXPR(v->display());
        // Lagrangian Nonlinear Systems
        if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
        {
	  d->computeForces(t, q, v);
	  accFree += *(d->forces());
        }
        else
          RuntimeException::selfThrow
          ("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

	if(d->inverseMass())
	  {
	    d->update_inverse_mass();
	    d->inverseMass()->PLUForwardBackwardInPlace(accFree);// contains right (left limit) acceleration without contact force
	  }
	DEBUG_PRINT("accFree contains left limit acceleration at  t^-_{k+1} without contact force :\n");
	DEBUG_EXPR(accFree.display());
      }
      else if(dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
        SiconosVector& accFree = *workVectors[OneStepIntegrator::free];

        accFree.zero();
        // get right state from memory
        SP::SiconosVector q = d->q(); // contains position q_{k+1}
        SP::SiconosVector v = d->twist(); // contains velocity v_{k+1}^- and not free velocity
        DEBUG_EXPR(accFree.display());
        DEBUG_EXPR(q->display());
        DEBUG_EXPR(v->display());

	d->computeForces(t, q, v);
	accFree += *(d->forces());

	if(d->inverseMass())
	  {
	    d->update_inverse_mass();
	    d->inverseMass()->PLUForwardBackwardInPlace(accFree);// contains right (left limit) acceleration without contact force
	  }

        DEBUG_PRINT("accFree contains left limit acceleration at  t^-_{k+1} without contact force :\n");
        DEBUG_EXPR(accFree.display());
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
      InteractionsGraph::VIterator ui, uiend;
      SP::Interaction inter;
      for(std11::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui)
      {
        inter = indexSet2->bundle(*ui);
        inter->relation()->computeJach(t, *inter, indexSet2->properties(*ui));
        inter->relation()->computeJacg(t, *inter, indexSet2->properties(*ui));
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
        _simulation->nonSmoothDynamicalSystem()->updateInput(_simulation->nextTime(),2);
      }
    }
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
        SiconosVector& accFree = *workVectors[OneStepIntegrator::free];

        // initialize *it->residuFree
        residuFree -= 0.5 * h* accFree;
        DEBUG_EXPR(accFree.display());

        if(d->p(2))
        {

          // get right state from memory
          DEBUG_EXPR(d->inverseMass()->display());
          DEBUG_EXPR(d->p(2)->display());
          SiconosVector dummy(*(d->p(2))); // value = contact force
          if(d->inverseMass())
          {
            d->update_inverse_mass();
            d->inverseMass()->PLUForwardBackwardInPlace(dummy);
          }

          residuFree -= 0.5 * h*dummy;

        }
        DEBUG_EXPR(residuFree.display());
      }
      else if(dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
        SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
        SiconosVector& accFree = *workVectors[OneStepIntegrator::free];
        // initialize *it->residuFree
        residuFree -= 0.5 * h* accFree;
        DEBUG_EXPR(accFree.display());

        if(d->p(2))
        {
          // get right state from memory
          DEBUG_EXPR(d->inverseMass()->display());
          DEBUG_EXPR(d->p(2)->display());
          SiconosVector dummy(*(d->p(2))); // value = contact force
	  if(d->inverseMass())
	    {
	      d->update_inverse_mass();
	      d->inverseMass()->PLUForwardBackwardInPlace(dummy);
	    }
          residuFree -= 0.5 * h*dummy;

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


  DEBUG_END("D1MinusLinearOSI::computeResiduHalfExplicitAccelerationLevel()\n");

  return 0.; // there is no Newton iteration and the residuum is assumed to vanish
}


void D1MinusLinearOSI::computeFreeOutputHalfExplicitAccelerationLevel(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{

  DEBUG_BEGIN("D1MinusLinearOSI::computeFreeOutputHalfExplicitAccelerationLevel\n");
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


  // define Xfree for velocity and acceleration level
  if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
  {
    //Xfree = inter->dataX();
    /* get the current velocity  of the aggregated ds */
    if(relationType == Lagrangian)
    {
      Xfree = DSlink[LagrangianR::q1];
      DEBUG_PRINT("Xfree = DSlink[LagrangianR::q1];\n");
    }
    else if(relationType == NewtonEuler)
    {
      Xfree = DSlink[NewtonEulerR::velocity];
      DEBUG_PRINT("Xfree = DSlink[LagrangianR::veclocity];\n");
    }
    else
      RuntimeException::selfThrow("D1MinusLinearOSI::computeFreeOutput - unknown relation type.");

    DEBUG_EXPR(Xfree->display());
    assert(Xfree);

  }
  else if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
  {
    /* get the "free" acceleration of the aggregated ds */
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
    DEBUG_PRINT("Xfree = DSlink[LagrangianR::xfree];\n");
    DEBUG_EXPR(Xfree->display());
    assert(Xfree);
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

    if(C)
    {
      assert(Xfree);
      coord[3] = C->size(1);
      coord[5] = C->size(1);
      subprod(*C, *Xfree, osnsp_rhs, coord, true);
    }
    DEBUG_EXPR(osnsp_rhs.display(););

    // in osnsp_rhs corrections have to be added
    SP::SiconosMatrix ID(new SimpleMatrix(sizeY, sizeY));
    ID->eye();

    Index xcoord(8);
    xcoord[0] = 0;
    xcoord[1] = sizeY;
    xcoord[2] = 0;
    xcoord[3] = sizeY;
    xcoord[4] = 0;
    xcoord[5] = sizeY;
    xcoord[6] = 0;
    xcoord[7] = sizeY;


    if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
    {

      /*  explicit time dependence -> partial time derivative has to be added */

      // if (relationSubType == RheonomousR) // explicit time dependence -> partial time derivative has to be added
      // {
      //   SiconosVector q = *DSlink[LagrangianR::q0];
      //   SiconosVector z = *DSlink[LagrangianR::z];

      //   std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->computehDot(simulation()->getTkp1(), q, z);
      //   *DSlink[LagrangianR::z] = z;
      //   subprod(*ID, *(std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->hDot()), osnsp_rhs, xcoord, false);
      // }

      if(relationSubType == RheonomousR)
        RuntimeException::selfThrow("D1MinusLinearOSI::computeFreeOutput is not implemented  at velocity level for LagrangianRheonomousR.");
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter, indexSet->properties(vertex_inter)));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }


    if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
    {
      if(relationSubType == ScleronomousR)  // acceleration term involving Hessian matrix of Relation with respect to degree is added
      {
        DEBUG_PRINT("D1MinusLinearOSI::computeFreeOutput. acceleration term involving Hessian matrix of Relation\n");
        DEBUG_EXPR(osnsp_rhs.display(););
        std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->computedotjacqhXqdot(simulation()->getTkp1(), *inter, DSlink);
        subprod(*ID, *(std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->dotjacqhXqdot()), osnsp_rhs, xcoord, false);
      }
      DEBUG_EXPR(osnsp_rhs.display(););

    }


  }
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



    if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
    {
      // in osnsp_rhs corrections have to be added
      SP::SiconosMatrix ID(new SimpleMatrix(sizeY, sizeY));
      ID->eye();

      Index xcoord(8);
      xcoord[0] = 0;
      xcoord[1] = sizeY;
      xcoord[2] = 0;
      xcoord[3] = sizeY;
      xcoord[4] = 0;
      xcoord[5] = sizeY;
      xcoord[6] = 0;
      xcoord[7] = sizeY;

      DEBUG_PRINT("D1MinusLinearOSI::computeFreeOutput.\n Adding the additional terms of the second order time derivative of constraints.\n");
      DEBUG_EXPR(osnsp_rhs.display(););

      /** Compute additional terms of the second order time derivative of constraints
       *
       * \f$ \nabla_q h(q) \dot T v + \frac{d}{dt}(\nabla_q h(q) ) T v \f$
       *
       */
      SP::DynamicalSystem ds1 = indexSet->properties(vertex_inter).source;
      SP::DynamicalSystem ds2 = indexSet->properties(vertex_inter).target;

      std11::static_pointer_cast<NewtonEulerR>(inter->relation())->computeSecondOrderTimeDerivativeTerms(simulation()->getTkp1(), *inter, DSlink, ds1, ds2);

      DEBUG_EXPR((std11::static_pointer_cast<NewtonEulerR>(inter->relation())->secondOrderTimeDerivativeTerms())->display());

      subprod(*ID, *(std11::static_pointer_cast<NewtonEulerR>(inter->relation())->secondOrderTimeDerivativeTerms()), osnsp_rhs, xcoord, false);
      DEBUG_EXPR(osnsp_rhs.display(););


    }
    DEBUG_EXPR(osnsp_rhs.display(););

    if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)  // impact terms are added
    {
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter, indexSet->properties(vertex_inter)));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }
  }
  else
    RuntimeException::selfThrow("D1MinusLinearOSI::computeFreeOutput - not implemented for Relation of type " + relationType);

  DEBUG_EXPR(osnsp_rhs.display(););
  DEBUG_END("D1MinusLinearOSI::computeFreeOutputHalfExplicitAccelerationLevel ends\n");

}


bool D1MinusLinearOSI::addInteractionInIndexSetHalfExplicitAccelerationLevel(SP::Interaction inter, unsigned int i)
{
  DEBUG_PRINT("D1MinusLinearOSI::addInteractionInIndexSetHalfExplicitAccelerationLevel.\n");

  if(i == 1)
  {
    DEBUG_PRINT(" level == 1\n");
    if(Type::value(*(inter->nonSmoothLaw())) == Type::EqualityConditionNSL)
    {
      return true;
    }
    /* Active for impulses calculations? Contacts that are closed */
    double y = (*(inter->y(0)))(0); // current position
    DEBUG_PRINTF("y= %24.16e\n", y);
    return (y <= DEFAULT_TOL_D1MINUS);

  }
  else if(i == 2)
  {
    DEBUG_PRINT(" level == 2\n");
    if(Type::value(*(inter->nonSmoothLaw())) == Type::EqualityConditionNSL)
    {
      return true;
    }
    /* Active for non-impulsive forces computation */
    double y = (*(inter->y(0)))(0); // current position
    double yDot = (*(inter->y(1)))(0); // current position
    DEBUG_PRINTF("y= %24.16e\n", y);
    DEBUG_PRINTF("yDot= %24.16e\n", yDot);
    DEBUG_EXPR(std::cout << std::boolalpha << (y <= DEFAULT_TOL_D1MINUS) <<std::endl;);
    DEBUG_EXPR(std::cout << std::boolalpha << (yDot <= DEFAULT_TOL_D1MINUS) <<std::endl;);
    return (y <= DEFAULT_TOL_D1MINUS) && (yDot <= DEFAULT_TOL_D1MINUS);
  }
  else if(i == 3)
  {
    if(Type::value(*(inter->nonSmoothLaw())) == Type::EqualityConditionNSL)
    {
      return false;
    }
    /*  Contacts which have been closing in the last time step */
    DEBUG_PRINT(" level == 3\n");
    double y = (*(inter->y(0)))(0); // current position
    double yOld = (*(inter->yOld(0)))(0); // old position
    DEBUG_PRINTF("y= %24.16e\n", y);
    DEBUG_PRINTF("yOld= %24.16e\n", yOld);
    /* if Interaction has not been active in the previous calculation
       and now becomes active */
    return (y <= DEFAULT_TOL_D1MINUS && yOld > DEFAULT_TOL_D1MINUS);
  }
  else
    RuntimeException::selfThrow("D1MinusLinearOSI::addInteractionInIndexSetHalfExplicitAccelerationLevel, IndexSet[i > 3] does not exist.");
  return false;
}

bool D1MinusLinearOSI::removeInteractionInIndexSetHalfExplicitAccelerationLevel(SP::Interaction inter, unsigned int i)
{
  DEBUG_PRINT("D1MinusLinearOSI::removeInteractionInIndexSetHalfExplicitAccelerationLevel.\n");

  return !(addInteractionInIndexSetHalfExplicitAccelerationLevel(inter,i));
}
