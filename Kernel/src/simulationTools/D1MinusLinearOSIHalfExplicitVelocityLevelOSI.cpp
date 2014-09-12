/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

/// @cond
using namespace RELATION;



double D1MinusLinearOSI::computeResiduHalfExplicitVelocityLevel()
{
  double t = simulationLink->nextTime(); // end of the time step
  double told = simulationLink->startingTime(); // beginning of the time step
  double h = simulationLink->timeStep(); // time step length
  SP::OneStepNSProblems allOSNS  = simulationLink->oneStepNSProblems(); // all OSNSP
  SP::Topology topo =  simulationLink->model()->nonSmoothDynamicalSystem()->topology();
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);
  SP::InteractionsGraph indexSet1 = topo->indexSet(1);
  SP::InteractionsGraph indexSet2 = topo->indexSet(2);

  // Note FP : we must use graph rather than InteractionSet. I have update all the code
  // in this file by replacing allInteractions calls with proper call to indexSet0.
  // So it produces the same results as before my changes.
  // Anyway, I think that in some place it will be better to call indexSet2 rather than indexSet0 as it was done in
  // the original code?

  DEBUG_PRINTF("nextTime %f\n", t);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);

  DEBUG_PRINT("\nEVALUATE LEFT HAND SIDE\n");
  /** Step 1-  solve a LCP at acceleration level for lambda^+_{k} for the last set indices
   * if index2 is empty we skip this step
   *
   **/

  DEBUG_EXPR(std::cout<< "allOSNS->empty()   " << std::boolalpha << allOSNS->empty() << std::endl << std::endl);
  DEBUG_EXPR(std::cout<< "allOSNS->size()   "  << allOSNS->size() << std::endl << std::endl);

// -- LEFT SIDE --
  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // type of the current DS
    Type::Siconos dsType = Type::value(**it);
    SP::SiconosVector workFree;
    SP::SiconosVector workFreeFree;
    SP::SiconosMatrix Mold;


    if ((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
      workFree = d->workspace(DynamicalSystem::free); // POINTER CONSTRUCTOR : contains acceleration without contact force
      workFree->zero();

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit
      Mold = d->mass();

      DEBUG_EXPR(workFree->display());
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());
      DEBUG_EXPR(Mold->display());

      if (! d->workspace(DynamicalSystem::free_tdg))
      {
        d->allocateWorkVector(DynamicalSystem::free_tdg, d->getDim()) ;
      }
      workFreeFree = d->workspace(DynamicalSystem::free_tdg);
      workFreeFree->zero();
      DEBUG_EXPR(workFreeFree->display());

      if (d->forces())
      {
        d->computeForces(told, qold, vold);
        DEBUG_EXPR(d->forces()->display());

        *workFree += *(d->forces());
      }
      Mold->PLUForwardBackwardInPlace(*workFree); // contains left (right limit) acceleration without contact force
      d->addWorkVector(workFree,DynamicalSystem::free_tdg); // store the value in WorkFreeFree
    }

    else if(dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);
      workFree = d->workspace(DynamicalSystem::free); // POINTER CONSTRUCTOR : contains acceleration without contact force
      workFree->zero();

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit
      //Mold = d->mass();
      assert(!d->mass()->isPLUInversed());
      Mold.reset(new SimpleMatrix(*(d->mass()))); // we copy the mass matrix to avoid its factorization
      DEBUG_EXPR(workFree->display());
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());
      DEBUG_EXPR(Mold->display());

      if (! d->workspace(DynamicalSystem::free_tdg))
      {
        d->allocateWorkVector(DynamicalSystem::free_tdg, d->getDim()) ;
      }

      workFreeFree = d->workspace(DynamicalSystem::free_tdg);
      workFreeFree->zero();
      DEBUG_EXPR(workFreeFree->display());

      if (d->forces())
      {
        d->computeForces(told, qold, vold);
        DEBUG_EXPR(d->forces()->display());

        *workFree += *(d->forces());
      }
      Mold->PLUForwardBackwardInPlace(*workFree); // contains left (right limit) acceleration without contact force

      d->addWorkVector(workFree,DynamicalSystem::free_tdg); // store the value in WorkFreeFree

    }
    else
    {
      RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    }


    DEBUG_PRINT("workFree contains right limit acceleration at  t^+_k with contact force :\n");
    DEBUG_EXPR(workFree->display());
    DEBUG_PRINT("workFreeFree contains right limit acceleration at t^+_k without contact force :\n");
    DEBUG_EXPR(workFreeFree->display());

  }


  if (!allOSNS->empty())
  {
    if (indexSet2->size() >0)
    {
      InteractionsGraph::VIterator ui, uiend;
      SP::Interaction inter;
      for (std11::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui)
      {
        inter = indexSet2->bundle(*ui);
        inter->relation()->computeJach(t, *inter, indexSet2->properties(*ui));
        if (inter->relation()->getType() == NewtonEuler)
        {
          SP::DynamicalSystem ds1 = indexSet2->properties(*ui).source;
          SP::DynamicalSystem ds2 = indexSet2->properties(*ui).target;
          SP::NewtonEulerR ner = std11::static_pointer_cast<NewtonEulerR>(indexSet2->bundle(*ui)->relation());
          ner->computeJachqT(*inter, ds1, ds2);
        }
        inter->relation()->computeJacg(told, *inter, indexSet2->properties(*ui));
      }

      if (simulationLink->model()->nonSmoothDynamicalSystem()->topology()->hasChanged())
      {
        for (OSNSIterator itOsns = allOSNS->begin(); itOsns != allOSNS->end(); ++itOsns)
        {
          (*itOsns)->setHasBeenUpdated(false);
        }
      }
      assert((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]);

      if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->hasInteractions())) // it should be equivalent to indexSet2
      {
        DEBUG_PRINT("We compute lambda^+_{k} \n");
        (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->compute(told);
        DEBUG_EXPR((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->display());

        // Note Franck : at the time this results in a call to swapInMem of all Interactions of the NSDS
        // So let the simu do this.
        //(*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->saveInMemory(); // we push y and lambda in Memories
        simulationLink->pushInteractionsInMemory();
        simulationLink->updateInput(2);
      }

      for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
      {
        Type::Siconos dsType = Type::value(**it);
        if ((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
        {
          SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
          SP::SiconosVector workFree = d->workspace(DynamicalSystem::free); // POINTER CONSTRUCTOR : contains acceleration without contact force

          SP::SiconosVector dummy(new SiconosVector(*(d->p(2)))); // value = contact force
          SP::SiconosMatrix Mold = d->mass();
          Mold->PLUForwardBackwardInPlace(*dummy);
          *workFree  += *(dummy);

          DEBUG_EXPR(d->p(2)->display());
        }
        else if (dsType == Type::NewtonEulerDS)
        {
          SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);
          SP::SiconosVector workFree = d->workspace(DynamicalSystem::free); // POINTER CONSTRUCTOR : contains acceleration without contact force

          SP::SiconosVector dummy(new SiconosVector(*(d->p(2)))); // value = contact force
          SP::SiconosMatrix Mold(new SimpleMatrix(*(d->mass())));  // we copy the mass matrix to avoid its factorization
          DEBUG_EXPR(Mold->display());
          Mold->PLUForwardBackwardInPlace(*dummy);
          *workFree  += *(dummy);

          DEBUG_EXPR(d->p(2)->display());

        }
        else
          RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

      }
    }
  }

  /** Step 2 -  compute v_{k,1}
   *
   **/


  DEBUG_PRINT("\n PREDICT RIGHT HAND SIDE\n");



  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {

    // type of the current DS
    Type::Siconos dsType = Type::value(**it);
    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    if ((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
      SP::SiconosVector workFree = d->workspace(DynamicalSystem::free); // contains acceleration without contact force

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      // initialize *it->residuFree and predicted right velocity (left limit)
      SP::SiconosVector residuFree = (*it)->workspace(DynamicalSystem::freeresidu); // contains residu without nonsmooth effect
      SP::SiconosVector v = d->velocity(); //contains velocity v_{k+1}^- and not free velocity
      residuFree->zero();
      v->zero();

      DEBUG_EXPR(workFree->display());
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());


      *residuFree -= 0.5 * h**workFree;

      *v += h**workFree;
      *v += *vold;

      DEBUG_EXPR(residuFree->display());
      DEBUG_EXPR(v->display());

      SP::SiconosVector q = d->q(); // POINTER CONSTRUCTOR : contains position q_{k+1}
      *q = *qold;

      scal(0.5 * h, *vold + *v, *q, false);
      DEBUG_EXPR(q->display());
    }
    else if (dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);
      SP::SiconosVector workFree = d->workspace(DynamicalSystem::free);

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      // initialize *it->residuFree and predicted right velocity (left limit)
      SP::SiconosVector residuFree = (*it)->workspace(DynamicalSystem::freeresidu); // contains residu without nonsmooth effect
      SP::SiconosVector v = d->velocity(); //contains velocity v_{k+1}^- and not free velocity
      residuFree->zero();
      v->zero();

      DEBUG_EXPR(workFree->display());
      DEBUG_EXPR(qold->display());
      DEBUG_EXPR(vold->display());


      *residuFree -= 0.5 * h**workFree;

      *v += h**workFree;
      *v += *vold;

      DEBUG_EXPR(residuFree->display());
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
      d->updateT();
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
     * F_{k,+} = F(told,qol,vold) \\
     * Work_{freefree} =  M^{-1}_k (F^+_{k})  \mbox{stored in workFreeFree} \\
     * Work_{free} =  M^{-1}_k (P^+_{2,k}+F^+_{k})  \mbox{stored in workFree} \\
     * R_{free} = -h/2 * M^{-1}_k (P^+_{2,k}+F^+_{k})  \mbox{stored in ResiduFree} \\
     * v_{k,1} = v_{k,0} + h * M^{-1}_k (P^+_{2,k}+F^+_{k})  \mbox{stored in v} \\
     * q_{k,1} = q_{k,0} + \frac{h}{2} (v_{k,0} + v_{k,1}) \mbox{stored in q} \\
     * \end{cases}
     * \f]
     **/
  }

  DEBUG_PRINT("\n DECIDE STRATEGY\n");
  /** Decide of the strategy impact or smooth multiplier */
  _isThereImpactInTheTimeStep = false;
  if (!allOSNS->empty())
  {

    for (unsigned int level = simulationLink->levelMinForOutput(); level < simulationLink->levelMaxForOutput(); level++)
    {
      simulationLink->updateOutput(level);
    }
    simulationLink->updateIndexSets();

    SP::Topology topo =  simulationLink->model()->nonSmoothDynamicalSystem()->topology();
    SP::InteractionsGraph indexSet0 = topo->indexSet(0);
    SP::InteractionsGraph indexSet1 = topo->indexSet(1);
    SP::InteractionsGraph indexSet2 = topo->indexSet(2);

    if (indexSet1->size() > 0)
    {
      _isThereImpactInTheTimeStep = true;
      DEBUG_PRINT("There is an impact in the step. indexSet1->size() > 0. _isThereImpactInTheTimeStep = true;\n");
    }
    else
    {
      _isThereImpactInTheTimeStep = false;
      DEBUG_PRINT("There is no  impact in the step. indexSet1->size() = 0. _isThereImpactInTheTimeStep = false;\n");
    }
  }


  // We recompute the residu if _isThereImpactInTheTimeStep = true;
  if (_isThereImpactInTheTimeStep)
  {

    DEBUG_PRINT("There is an impact in the step. indexSet1->size() > 0.  _isThereImpactInTheTimeStep = true\n");

    for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    {
      // type of the current DS
      Type::Siconos dsType = Type::value(**it);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if ((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
      {
        SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
        SP::SiconosVector residuFree = d->workspace(DynamicalSystem::freeresidu);
        SP::SiconosVector v = d->velocity();
        SP::SiconosVector q = d->q();
        SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
        SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit

        SP::SiconosMatrix M = d->mass(); // POINTER CONSTRUCTOR : contains mass matrix

        //residuFree->zero();
        v->zero();
        SP::SiconosVector workFreeFree = d->workspace(DynamicalSystem::free_tdg);
        assert(workFreeFree);
        *residuFree = 0.5 * h**workFreeFree;
        workFreeFree->zero();

        d->computeMass();
        M->resetLU();
        DEBUG_EXPR(M->display());
        if (d->forces())
        {
          d->computeForces(t, q, v);
          *workFreeFree += *(d->forces());
        }

        M->PLUForwardBackwardInPlace(*workFreeFree); // contains right (left limit) acceleration without contact force
        *residuFree -= 0.5 * h**workFreeFree;
        DEBUG_EXPR(residuFree->display());
      }
      else if (dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);
        SP::SiconosVector residuFree = d->workspace(DynamicalSystem::freeresidu);
        SP::SiconosVector v = d->velocity();
        SP::SiconosVector q = d->q();
        SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
        SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit

        SP::SiconosMatrix M(new SimpleMatrix(*(d->mass()))); // we copy the mass matrix to avoid its factorization;
        DEBUG_EXPR(M->display());

        //residuFree->zero();
        v->zero();
        SP::SiconosVector workFreeFree = d->workspace(DynamicalSystem::free_tdg);
        assert(workFreeFree);
        *residuFree = 0.5 * h**workFreeFree;
        workFreeFree->zero();

        if (d->forces())
        {
          d->computeForces(t, q, v);
          *workFreeFree += *(d->forces());
        }

        M->PLUForwardBackwardInPlace(*workFreeFree); // contains right (left limit) acceleration without contact force
        *residuFree -= 0.5 * h**workFreeFree;
        DEBUG_EXPR(residuFree->display());
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
    for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    {
      // type of the current DS
      Type::Siconos dsType = Type::value(**it);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if ((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
      {

        SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
        SP::SiconosVector workFree = d->workspace(DynamicalSystem::free); // POINTER CONSTRUCTOR : contains acceleration without contact force
        workFree->zero();
        // get right state from memory
        SP::SiconosVector q = d->q(); // contains position q_{k+1}
        SP::SiconosVector v = d->velocity(); // contains velocity v_{k+1}^- and not free velocity
        SP::SiconosMatrix M = d->mass(); // POINTER CONSTRUCTOR : contains mass matrix

        DEBUG_EXPR(workFree->display());
        DEBUG_EXPR(q->display());
        DEBUG_EXPR(v->display());
        // Lagrangian Nonlinear Systems
        if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
        {
          d->computeMass();
          M->resetLU();
          DEBUG_EXPR(M->display());
          if (d->forces())
          {
            d->computeForces(t, q, v);
            *workFree += *(d->forces());
          }
        }
        else
          RuntimeException::selfThrow
          ("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

        M->PLUForwardBackwardInPlace(*workFree); // contains right (left limit) acceleration without contact force
        DEBUG_PRINT("workFree contains left limit acceleration at  t^-_{k+1} without contact force :\n");
        DEBUG_EXPR(workFree->display());
      }
      else if (dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);
        SP::SiconosVector workFree = d->workspace(DynamicalSystem::free); // POINTER CONSTRUCTOR : contains acceleration without contact force
        workFree->zero();
        // get right state from memory
        SP::SiconosVector q = d->q(); // contains position q_{k+1}
        SP::SiconosVector v = d->velocity(); // contains velocity v_{k+1}^- and not free velocity
        SP::SiconosMatrix M(new SimpleMatrix(*(d->mass()))); // we copy the mass matrix to avoid its factorization;

        DEBUG_EXPR(workFree->display());
        DEBUG_EXPR(q->display());
        DEBUG_EXPR(v->display());

        if (d->forces())
        {
          d->computeForces(t, q, v);
          *workFree += *(d->forces());
        }

        M->PLUForwardBackwardInPlace(*workFree); // contains right (left limit) acceleration without contact force
        DEBUG_PRINT("workFree contains left limit acceleration at  t^-_{k+1} without contact force :\n");
        DEBUG_EXPR(workFree->display());
      }
      else
        RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

    }

    // solve a LCP at acceleration level only for contacts which have been active at the beginning of the time-step
    if (!allOSNS->empty())
    {
      // for (unsigned int level = simulationLink->levelMinForOutput(); level < simulationLink->levelMaxForOutput(); level++)
      // {
      //   simulationLink->updateOutput(level);
      // }
      // simulationLink->updateIndexSets();
      DEBUG_PRINT("We compute lambda^-_{k+1} \n");
      InteractionsGraph::VIterator ui, uiend;
      SP::Interaction inter;
      for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
      {
        inter = indexSet0->bundle(*ui);
        inter->relation()->computeJach(t, *inter, indexSet0->properties(*ui));
        if (inter->relation()->getType() == NewtonEuler)
        {
          SP::DynamicalSystem ds1 = indexSet0->properties(*ui).source;
          SP::DynamicalSystem ds2 = indexSet0->properties(*ui).target;
          SP::NewtonEulerR ner = (std11::static_pointer_cast<NewtonEulerR>(inter->relation()));
          ner->computeJachqT(*inter, ds1, ds2);
        }
        inter->relation()->computeJacg(t, *inter, indexSet0->properties(*ui));
      }
      if (simulationLink->model()->nonSmoothDynamicalSystem()->topology()->hasChanged())
      {
        for (OSNSIterator itOsns = allOSNS->begin(); itOsns != allOSNS->end(); ++itOsns)
        {
          (*itOsns)->setHasBeenUpdated(false);
        }
      }

      if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->hasInteractions()))
      {
        (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->compute(t);
        DEBUG_EXPR((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->display(););
        simulationLink->updateInput(2);
      }
    }

    for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    {
      // type of the current DS
      Type::Siconos dsType = Type::value(**it);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if ((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
      {
        SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
        SP::SiconosVector workFree = d->workspace(DynamicalSystem::free); // contains acceleration without contact force

        // initialize *it->residuFree
        SP::SiconosVector residuFree = (*it)->workspace(DynamicalSystem::freeresidu); // POINTER CONSTRUCTOR : contains residu without nonsmooth effect
        *residuFree -= 0.5 * h**workFree;
        DEBUG_EXPR(workFree->display());

        if (d->p(2))
        {

          // get right state from memory
          SP::SiconosMatrix M = d->mass();
          DEBUG_EXPR(M->display());
          DEBUG_EXPR(d->p(2)->display());
          SP::SiconosVector dummy(new SiconosVector(*(d->p(2)))); // value = contact force

          M->PLUForwardBackwardInPlace(*dummy);
          *residuFree -= 0.5 * h**dummy;

        }
        DEBUG_EXPR(residuFree->display());
      }
      else if (dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);
        SP::SiconosVector workFree = d->workspace(DynamicalSystem::free); // contains acceleration without contact force

        // initialize *it->residuFree
        SP::SiconosVector residuFree = (*it)->workspace(DynamicalSystem::freeresidu); // POINTER CONSTRUCTOR : contains residu without nonsmooth effect
        *residuFree -= 0.5 * h**workFree;
        DEBUG_EXPR(workFree->display());

        if (d->p(2))
        {
          // get right state from memory
          SP::SiconosMatrix M(new SimpleMatrix(*(d->mass()))); // we copy the mass matrix to avoid its factorization;
          DEBUG_EXPR(M->display());
          DEBUG_EXPR(d->p(2)->display());
          SP::SiconosVector dummy(new SiconosVector(*(d->p(2)))); // value = contact force

          M->PLUForwardBackwardInPlace(*dummy);
          *residuFree -= 0.5 * h**dummy;

        }
        DEBUG_EXPR(residuFree->display());
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

  DEBUG_PRINT("D1MinusLinearOSI::computeResidu() ends\n");
  return 0.; // there is no Newton iteration and the residuum is assumed to vanish
}

