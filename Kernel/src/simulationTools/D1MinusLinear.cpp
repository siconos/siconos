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

#include "D1MinusLinear.hpp"
#include "Simulation.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianRheonomousR.hpp"
#include "LagrangianScleronomousR.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonImpactNSL.hpp"


#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include "debug.h"

using namespace std;
using namespace RELATION;

struct D1MinusLinear::_NSLEffectOnFreeOutput : public SiconosVisitor
{

  using SiconosVisitor::visit;

  OneStepNSProblem* _osnsp;
  SP::Interaction _inter;

  _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::Interaction inter) : _osnsp(p), _inter(inter) {};

  void visit(const NewtonImpactNSL& nslaw)
  {
    double e = nslaw.e();
    Index subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = _inter->getNonSmoothLawSize();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    subscal(e, *(_inter->y_k(_osnsp->levelMin())), *(_inter->yp()), subCoord, false);
  }
  void visit(const EqualityConditionNSL& nslaw)
  {
    ;
  }
};

D1MinusLinear::D1MinusLinear(SP::DynamicalSystem newDS) :
  OneStepIntegrator(OSI::D1MINUSLINEAR)
{
  OSIDynamicalSystems->insert(newDS);
}

D1MinusLinear::D1MinusLinear(DynamicalSystemsSet& newDS): OneStepIntegrator(OSI::D1MINUSLINEAR, newDS) {}

void D1MinusLinear::initialize()
{
  for(DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    Type::Siconos dsType = Type::value(**it);
    if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
      d->computeMass();
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      //SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);

    }
    else
      RuntimeException::selfThrow("D1MinusLinear::initialize - not implemented for Dynamical system type: " + dsType);
  }
}

double D1MinusLinear::computeResidu()
{

  DEBUG_PRINT("\n D1MinusLinear::computeResidu(), start\n");
  double t = simulationLink->nextTime(); // end of the time step
  double told = simulationLink->startingTime(); // beginning of the time step
  double h = simulationLink->timeStep(); // time step length
  SP::OneStepNSProblems allOSNS  = simulationLink->oneStepNSProblems(); // all OSNSP
  SP::Topology topo =  simulationLink->model()->nonSmoothDynamicalSystem()->topology();
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);
  SP::InteractionsGraph indexSet1 = topo->indexSet(1);
  SP::InteractionsGraph indexSet2 = topo->indexSet(2);




  SP::InteractionsSet allInteractions = simulationLink->model()->nonSmoothDynamicalSystem()->interactions(); // all Interactions

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
  for(DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // type of the current DS
    Type::Siconos dsType = Type::value(**it);
    SP::SiconosVector workFree;
    SP::SiconosVector workFreeFree;
    SP::SiconosMatrix Mold;


    if((dsType == Type::LagrangianDS) or (dsType == Type::LagrangianLinearTIDS))
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

      if(! d->workspace(DynamicalSystem::free_tdg))
      {
        d->allocateWorkVector(DynamicalSystem::free_tdg, d->getDim()) ;
      }
      workFreeFree = d->workspace(DynamicalSystem::free_tdg);
      workFreeFree->zero();
      DEBUG_EXPR(workFreeFree->display());

      if(d->forces())
      {
        d->computeForces(told, qold, vold);
        DEBUG_EXPR(d->forces()->display());

        *workFree += *(d->forces());
      }
      Mold->PLUForwardBackwardInPlace(*workFree); // contains left (right limit) acceleration without contact force
      d->addWorkVector(workFree,DynamicalSystem::free_tdg); // store the value in WorkFreeFree
    }
    else if((dsType == Type::NewtonEulerDS))
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

      if(! d->workspace(DynamicalSystem::free_tdg))
      {
        d->allocateWorkVector(DynamicalSystem::free_tdg, d->getDim()) ;
      }

      workFreeFree = d->workspace(DynamicalSystem::free_tdg);
      workFreeFree->zero();
      DEBUG_EXPR(workFreeFree->display());

      if(d->forces())
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
      RuntimeException::selfThrow("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    }


    DEBUG_PRINT("workFree contains right limit acceleration at  t^+_k with contact force :\n");
    DEBUG_EXPR(workFree->display());
    DEBUG_PRINT("workFreeFree contains right limit acceleration at t^+_k without contact force :\n");
    DEBUG_EXPR(workFreeFree->display());

  }


  if(!allOSNS->empty())
  {
    if(indexSet2->size() >0)
    {

      for(InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
      {
        (*it)->computeJach(told);
        (*it)->computeJacg(told);
      }

      if(simulationLink->model()->nonSmoothDynamicalSystem()->topology()->hasChanged())
      {
        for(OSNSIterator itOsns = allOSNS->begin(); itOsns != allOSNS->end(); ++itOsns)
        {
          (*itOsns)->setHasBeenUpdated(false);
        }
      }
      assert((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]);
      assert((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->interactions());

      if(!((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->interactions()->isEmpty()))  // it should be equivalent to indexSet2
      {
        DEBUG_PRINT("We compute lambda^+_{k} \n");
        (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->compute(told);
        DEBUG_EXPR((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->display());
        (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->saveInMemory(); // we push y and lambda in Memories
        simulationLink->updateInput(2);
      }

      for(DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
      {
        Type::Siconos dsType = Type::value(**it);
        if((dsType == Type::LagrangianDS) or (dsType == Type::LagrangianLinearTIDS))
        {
          SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
          SP::SiconosVector workFree = d->workspace(DynamicalSystem::free); // POINTER CONSTRUCTOR : contains acceleration without contact force

          SP::SiconosVector dummy(new SiconosVector(*(d->p(2)))); // value = contact force
          SP::SiconosMatrix Mold = d->mass();
          Mold->PLUForwardBackwardInPlace(*dummy);
          *workFree  += *(dummy);

          DEBUG_EXPR(d->p(2)->display());
        }
        else if(dsType == Type::NewtonEulerDS)
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
          RuntimeException::selfThrow("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);

      }
    }
  }

  /** Step 2 -  compute v_{k,1}
   *
   **/


  DEBUG_PRINT("\n PREDICT RIGHT HAND SIDE\n");



  for(DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {

    // type of the current DS
    Type::Siconos dsType = Type::value(**it);
    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    if((dsType == Type::LagrangianDS) or (dsType == Type::LagrangianLinearTIDS))
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
    else if(dsType == Type::NewtonEulerDS)
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
      DEBUG_EXPR(q->display());

    }
    else
      RuntimeException::selfThrow("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);


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
  if(!allOSNS->empty())
  {

    for(unsigned int level = simulationLink->levelMinForOutput(); level < simulationLink->levelMaxForOutput(); level++)
    {
      simulationLink->updateOutput(level);
    }
    simulationLink->updateIndexSets();

    SP::Topology topo =  simulationLink->model()->nonSmoothDynamicalSystem()->topology();
    SP::InteractionsGraph indexSet0 = topo->indexSet(0);
    SP::InteractionsGraph indexSet1 = topo->indexSet(1);
    SP::InteractionsGraph indexSet2 = topo->indexSet(2);

    if(indexSet1->size() > 0)
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
  if(_isThereImpactInTheTimeStep)
  {

    DEBUG_PRINT("There is an impact in the step. indexSet1->size() > 0.  _isThereImpactInTheTimeStep = true\n");
    
    for(DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    {
      // type of the current DS
      Type::Siconos dsType = Type::value(**it);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if((dsType == Type::LagrangianDS) or (dsType == Type::LagrangianLinearTIDS))
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
        if(d->forces())
        {
          d->computeForces(t, q, v);
          *workFreeFree += *(d->forces());
        }

        M->PLUForwardBackwardInPlace(*workFreeFree); // contains right (left limit) acceleration without contact force
        *residuFree -= 0.5 * h**workFreeFree;
        DEBUG_EXPR(residuFree->display());
      }
      else if(dsType == Type::NewtonEulerDS)
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

        if(d->forces())
        {
          d->computeForces(t, q, v);
          *workFreeFree += *(d->forces());
        }

        M->PLUForwardBackwardInPlace(*workFreeFree); // contains right (left limit) acceleration without contact force
        *residuFree -= 0.5 * h**workFreeFree;
        DEBUG_EXPR(residuFree->display());
      }
      else
        RuntimeException::selfThrow("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    }
  }
  else
  {
    DEBUG_PRINT("There is no  impact in the step. indexSet1->size() = 0. _isThereImpactInTheTimeStep = false;\n");
    // -- RIGHT SIDE --
    // calculate acceleration without contact force
    for(DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    {
      // type of the current DS
      Type::Siconos dsType = Type::value(**it);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if((dsType == Type::LagrangianDS) or (dsType == Type::LagrangianLinearTIDS))
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
        if(dsType == Type::LagrangianDS or dsType == Type::LagrangianLinearTIDS)
        {
          d->computeMass();
          M->resetLU();
          DEBUG_EXPR(M->display());
          if(d->forces())
          {
            d->computeForces(t, q, v);
            *workFree += *(d->forces());
          }
        }
        else
          RuntimeException::selfThrow
          ("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);

        M->PLUForwardBackwardInPlace(*workFree); // contains right (left limit) acceleration without contact force
        DEBUG_PRINT("workFree contains left limit acceleration at  t^-_{k+1} without contact force :\n");
        DEBUG_EXPR(workFree->display());
      }
      else if(dsType == Type::NewtonEulerDS)
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

        if(d->forces())
        {
          d->computeForces(t, q, v);
          *workFree += *(d->forces());
        }

        M->PLUForwardBackwardInPlace(*workFree); // contains right (left limit) acceleration without contact force
        DEBUG_PRINT("workFree contains left limit acceleration at  t^-_{k+1} without contact force :\n");
        DEBUG_EXPR(workFree->display());
      }
      else
        RuntimeException::selfThrow("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);

    }

    // solve a LCP at acceleration level only for contacts which have been active at the beginning of the time-step
    if(!allOSNS->empty())
    {
      // for (unsigned int level = simulationLink->levelMinForOutput(); level < simulationLink->levelMaxForOutput(); level++)
      // {
      //   simulationLink->updateOutput(level);
      // }
      // simulationLink->updateIndexSets();
      for(InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
      {
        (*it)->computeJach(t);
        (*it)->computeJacg(t);
      }

      if(simulationLink->model()->nonSmoothDynamicalSystem()->topology()->hasChanged())
      {
        for(OSNSIterator itOsns = allOSNS->begin(); itOsns != allOSNS->end(); ++itOsns)
        {
          (*itOsns)->setHasBeenUpdated(false);
        }
      }

      if(!((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->interactions())->isEmpty())
      {
        (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->compute(t);
        simulationLink->updateInput(2);
      }
    }

    for(DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    {
      // type of the current DS
      Type::Siconos dsType = Type::value(**it);
      /* \warning the following conditional statement should be removed with a MechanicalDS class */
      if((dsType == Type::LagrangianDS) or (dsType == Type::LagrangianLinearTIDS))
      {
        SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
        SP::SiconosVector workFree = d->workspace(DynamicalSystem::free); // contains acceleration without contact force

        // initialize *it->residuFree
        SP::SiconosVector residuFree = (*it)->workspace(DynamicalSystem::freeresidu); // POINTER CONSTRUCTOR : contains residu without nonsmooth effect
        *residuFree -= 0.5 * h**workFree;
        DEBUG_EXPR(workFree->display());

        if(d->p(2))
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
      else if(dsType == Type::NewtonEulerDS)
      {
        SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);
        SP::SiconosVector workFree = d->workspace(DynamicalSystem::free); // contains acceleration without contact force
        
        // initialize *it->residuFree
        SP::SiconosVector residuFree = (*it)->workspace(DynamicalSystem::freeresidu); // POINTER CONSTRUCTOR : contains residu without nonsmooth effect
        *residuFree -= 0.5 * h**workFree;
        DEBUG_EXPR(workFree->display());
        
        if(d->p(2))
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
        RuntimeException::selfThrow("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);

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

  DEBUG_PRINT("D1MinusLinear::computeResidu(), end\n");
  return 0.; // there is no Newton iteration and the residuum is assumed to vanish
}

void D1MinusLinear::computeFreeState()
{
  DEBUG_PRINT("\n D1MinusLinear::computeFreeState(), start\n");


  for(DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // type of the current DS
    Type::Siconos dsType = Type::value(**it);
    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    if((dsType == Type::LagrangianDS) or (dsType == Type::LagrangianLinearTIDS))
    {
      // Lagrangian Systems
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
      
      // get left state from memory
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit
      DEBUG_EXPR(vold->display());
      
      // get right information
      SP::SiconosMatrix M = d->mass();
      SP::SiconosVector vfree = d->velocity(); // POINTER CONSTRUCTOR : contains free velocity
      (*vfree) = *(d->workspace(DynamicalSystem::freeresidu));
      DEBUG_EXPR(d->workspace(DynamicalSystem::freeresidu)->display());
      // d->computeMass();
      // M->resetLU();
      // M->PLUForwardBackwardInPlace(*vfree);
      // DEBUG_EXPR(M->display());
      
      *vfree *= -1.;
      *vfree += *vold;
      DEBUG_EXPR(vfree->display());
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      // NewtonEuler Systems
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);
      
      // get left state from memory
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit
      DEBUG_EXPR(vold->display());
      
      // get right information
      SP::SiconosMatrix M(new SimpleMatrix(*(d->mass()))); // we copy the mass matrix to avoid its factorization;
      SP::SiconosVector vfree = d->velocity(); // POINTER CONSTRUCTOR : contains free velocity
      (*vfree) = *(d->workspace(DynamicalSystem::freeresidu));
      DEBUG_EXPR(d->workspace(DynamicalSystem::freeresidu)->display());

      *vfree *= -1.;
      *vfree += *vold;
      DEBUG_EXPR(vfree->display());
      
      
    }
    else
      RuntimeException::selfThrow("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);

  }


  DEBUG_PRINT("D1MinusLinear::computeFreeState(), end\n");


}

void D1MinusLinear::computeFreeOutput(SP::Interaction inter, OneStepNSProblem* osnsp)
{

  DEBUG_PRINT("D1MinusLinear::computeFreeOutput starts\n");


  SP::OneStepNSProblems allOSNS  = simulationLink->oneStepNSProblems(); // all OSNSP

  // get relation and non smooth law information
  RELATION::TYPES relationType = inter->getRelationType(); // relation
  RELATION::SUBTYPES relationSubType = inter->getRelationSubType();
  unsigned int relativePosition = 0;
  unsigned int sizeY = inter->getNonSmoothLawSize(); // related NSL

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
  SP::SiconosVector Yp = inter->yp(); // POINTER CONSTRUCTOR : contains output

  // define Xfree for velocity and acceleration level
  if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
  {
    //Xfree = inter->dataX();
    if(relationType == Lagrangian)
    {
      Xfree = inter->data(LagrangianR::q1);
    }
    else if(relationType == NewtonEuler)
    {
      Xfree = inter->data(NewtonEulerR::velocity);
    }
    else
      RuntimeException::selfThrow("D1MinusLinear::computeFreeOutput - unknown relation type.");

  }
  else if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
  {

    if(relationType == Lagrangian)
    {
      Xfree = inter->data(LagrangianR::free);
    }
    else if(relationType == NewtonEuler)
    {
      Xfree = inter->data(NewtonEulerR::free);
    }
    else
      RuntimeException::selfThrow("D1MinusLinear::computeFreeOutput - unknown relation type.");
    DEBUG_PRINT("Xfree = inter->data(LagrangianR::free);\n");
    DEBUG_EXPR(Xfree->display());
    assert(Xfree);
  }
  else
    RuntimeException::selfThrow("D1MinusLinear::computeFreeOutput - OSNSP neither on velocity nor on acceleration level.");

  // calculate data of interaction
  SP::Interaction mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());

  // only Lagrangian Systems
  if(relationType == Lagrangian)
  {
    // in Yp the linear part of velocity or acceleration relation will be saved
    C = std11::static_pointer_cast<LagrangianR>(mainInteraction->relation())->C();

    if(C)
    {
      assert(Xfree);
      assert(Yp);

      coord[3] = C->size(1);
      coord[5] = C->size(1);
      subprod(*C, *Xfree, *Yp, coord, true);
    }

    // in Yp corrections have to be added
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

    if(relationSubType == RheonomousR)  // explicit time dependence -> partial time derivative has to be added
    {
      if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
      {
        std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->computehDot(simulation()->getTkp1(), *inter);
        subprod(*ID, *(std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->hDot()), *Yp, xcoord, false);
      }
      else
        RuntimeException::selfThrow("D1MinusLinear::computeFreeOutput is only implemented  at velocity level for LagrangianRheonomousR.");
    }

    if(relationSubType == ScleronomousR)  // acceleration term involving Hessian matrix of Relation with respect to degree is added
    {
      DEBUG_PRINT("D1MinusLinear::computeFreeOutput. acceleration term involving Hessian matrix of Relation\n");
      DEBUG_EXPR(Yp->display(););

      if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
      {
        std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->computeNonLinearH2dot(simulation()->getTkp1(), *inter);
        subprod(*ID, *(std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->Nonlinearh2dot()), *Yp, xcoord, false);
      }
      DEBUG_EXPR(Yp->display(););
    }


  }

  else if(relationType == NewtonEuler)
  {
    SP::SiconosMatrix CT =  std11::static_pointer_cast<NewtonEulerR>(mainInteraction->relation())->jachqT();

    if(CT)
    {
      coord[3] = CT->size(1);
      coord[5] = CT->size(1);
      assert(Yp);
      assert(Xfree);
      // creates a POINTER link between workX[ds] (xfree) and the
      // corresponding interactionBlock in each Interaction for each ds of the
      // current Interaction.
      // XXX Big quirks !!! -- xhub
      subprod(*CT, *Xfree, *Yp, coord, true);
    }

    DEBUG_PRINT("D1MinusLinear::computeFreeOutput. acceleration term involving Hessian matrix of Relation\n");
    DEBUG_EXPR(Yp->display(););

    if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
    {
      // std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->computeNonLinearH2dot(simulation()->getTkp1(), *inter);
      // subprod(*ID, *(std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->Nonlinearh2dot()), *Yp, xcoord, false);
    }
    DEBUG_EXPR(Yp->display(););
    if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)  // impact terms are added
    {
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }
  }

  else
    RuntimeException::selfThrow("D1MinusLinear::computeFreeOutput - not implemented for Relation of type " + relationType);
  DEBUG_PRINT("D1MinusLinear::computeFreeOutput ends\n");

}

void D1MinusLinear::updateState(const unsigned int level)
{
  DEBUG_PRINTF("\n D1MinusLinear::updateState(const unsigned int level) start for level = %i\n",level);

  for(DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // type of the current DS
    Type::Siconos dsType = Type::value(**it);
    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    if((dsType == Type::LagrangianDS) or (dsType == Type::LagrangianLinearTIDS))
    {
      
      // Lagrangian Systems
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
      SP::SiconosMatrix M = d->mass();
      SP::SiconosVector v = d->velocity(); // POINTER CONSTRUCTOR : contains new velocity
      
      if(d->p(1))
      {
        SP::SiconosVector dummy(new SiconosVector(*(d->p(1)))); // value = nonsmooth impulse
        M->PLUForwardBackwardInPlace(*dummy); // solution for its velocity equivalent
        *v += *dummy; // add free velocity
        DEBUG_PRINT("\nRIGHT IMPULSE\n");
        DEBUG_EXPR(d->p(1)->display());
      }
      DEBUG_EXPR(d->q()->display());
      DEBUG_EXPR(d->velocity()->display());
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      // NewtonEuler Systems
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);
      SP::SiconosMatrix M(new SimpleMatrix(*(d->mass()))); // we copy the mass matrix to avoid its factorization;
      SP::SiconosVector v = d->velocity(); // POINTER CONSTRUCTOR : contains new velocity
      if(d->p(1))
      {
        SP::SiconosVector dummy(new SiconosVector(*(d->p(1)))); // value = nonsmooth impulse
        M->PLUForwardBackwardInPlace(*dummy); // solution for its velocity equivalent
        *v += *dummy; // add free velocity
        DEBUG_PRINT("\nRIGHT IMPULSE\n");
        DEBUG_EXPR(d->p(1)->display());
      }
      DEBUG_EXPR(d->q()->display());
      DEBUG_EXPR(d->velocity()->display());
    }
    else
      RuntimeException::selfThrow("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);

  }
    
    DEBUG_PRINT("\n D1MinusLinear::updateState(const unsigned int level) end\n");

}


bool D1MinusLinear::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  DEBUG_PRINT("D1MinusLinear::addInteractionInIndexSet.\n");
  assert((i == 1) || (i==2));
  // double h = simulationLink->timeStep();

  double y ;
  double yOld ;
  SP::Relation r = inter->relation();
  RELATION::TYPES relationType = r->getType();
  SP::LagrangianDS lds;
  if(relationType == Lagrangian)
  {

    // compute yold
    SP::BlockVector qold(new BlockVector());
    SP::BlockVector vold(new BlockVector());
    SP::BlockVector zold(new BlockVector());
    SP::BlockVector q(new BlockVector());
    SP::BlockVector v(new BlockVector());
    SP::BlockVector z(new BlockVector());

    for(DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
    {
      // check dynamical system type
      assert((Type::value(**it) == Type::LagrangianLinearTIDS ||
              Type::value(**it) == Type::LagrangianDS));

      // convert vDS systems into LagrangianDS and put them in vLDS
      lds = std11::static_pointer_cast<LagrangianDS> (*it);

      qold->insertPtr(lds->qMemory()->getSiconosVector(0));
      vold->insertPtr(lds->velocityMemory()->getSiconosVector(0));
      /** \warning Warning the value of z of not stored. */
      zold->insertPtr(lds->z());
      q->insertPtr(lds->q());
      v->insertPtr(lds->velocity());
      z->insertPtr(lds->z());
    }
    std11::static_pointer_cast<LagrangianScleronomousR>(r)->computeh(*inter,qold,zold);
    yOld = (inter->y(0))->getValue(0);
    // Compute current y (we assume that q stores q_{k,1} and v stores v_{k,1})
    // If not sure we have to store it into a new date in Interaction.
    std11::static_pointer_cast<LagrangianScleronomousR>(r)->computeh(*inter,q,z);
    y = (inter->y(0))->getValue(0);
  }

  DEBUG_PRINTF("D1MinusLinear::addInteractionInIndexSet of level = %i yOld=%e, y=%e \n", i,  yOld, y);

#if __cplusplus >= 201103L
  assert(!::isnan(y));
#else
  assert(!isnan(y));
#endif

  DEBUG_EXPR(
    if((yOld >0.0) && (y <= y))
    DEBUG_PRINT("D1MinusLinear::addInteractionInIndexSet contact are closing ((yOld >0.0) && (y <= y)).\n");
  );
  return ((yOld >0.0) && (y <= y));
}


bool D1MinusLinear::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(i == 1);
  double h = simulationLink->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
  double gamma = 1.0 / 2.0;
  // if (_useGamma)
  // {
  //   gamma = _gamma;
  // }
  DEBUG_PRINTF("D1MinusLinear::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
#if __cplusplus >= 201103L
  assert(!::isnan(y));
#else
  assert(!isnan(y));
#endif
  DEBUG_EXPR(
    if(y > 0)
    DEBUG_PRINT("D1MinusLinear::removeInteractionInIndexSet DEACTIVATE.\n");
  );

  return (y > 0.0);
}





void D1MinusLinear::insertDynamicalSystem(SP::DynamicalSystem ds)
{
  OSIDynamicalSystems->insert(ds);
}

D1MinusLinear* D1MinusLinear::convert(OneStepIntegrator* osi)
{
  return dynamic_cast<D1MinusLinear*>(osi);
}
