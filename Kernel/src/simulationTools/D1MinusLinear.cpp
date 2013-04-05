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
#include "LagrangianRheonomousR.hpp"
#include "LagrangianScleronomousR.hpp"
#include "NewtonImpactNSL.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
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
};

D1MinusLinear::D1MinusLinear(SP::DynamicalSystem newDS) :
  OneStepIntegrator(OSI::D1MINUSLINEAR)
{
  OSIDynamicalSystems->insert(newDS);
}

D1MinusLinear::D1MinusLinear(DynamicalSystemsSet& newDS): OneStepIntegrator(OSI::D1MINUSLINEAR, newDS) {}

void D1MinusLinear::initialize()
{
  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    Type::Siconos dsType = Type::value(**it);
    if (dsType != Type::LagrangianDS && dsType != Type::LagrangianLinearTIDS)
      RuntimeException::selfThrow("D1MinusLinear::initialize - not implemented for Dynamical system type: " + dsType);
    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
    d->computeMass();
  }
}

double D1MinusLinear::computeResidu()
{

  DEBUG_PRINT("\n D1MinusLinear::computeResidu(), start\n");
  double t = simulationLink->nextTime(); // end of the time step
  double told = simulationLink->startingTime(); // beginning of the time step
  double h = simulationLink->timeStep(); // time step length
  SP::OneStepNSProblems allOSNS  = simulationLink->oneStepNSProblems(); // all OSNSP
  SP::InteractionsSet allInteractions = simulationLink->model()->nonSmoothDynamicalSystem()->interactions(); // all Interactions

  DEBUG_PRINTF("nextTime %f\n", t);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);

  DEBUG_PRINT("\nLEFT SIDE\n");

  // -- LEFT SIDE --
  // calculate acceleration without contact force
  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // type of the current DS
    Type::Siconos dsType = Type::value(**it);
    if ((dsType != Type::LagrangianDS) and (dsType != Type::LagrangianLinearTIDS ))
    {
      RuntimeException::selfThrow("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    }


    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
    SP::SiconosVector workFree = d->workFree(); // POINTER CONSTRUCTOR : contains acceleration without contact force
    workFree->zero();

    // get left state from memory
    SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
    SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit
    SP::SiconosMatrix Mold = d->mass();

    DEBUG_EXPR(workFree->display());
    DEBUG_EXPR(qold->display());
    DEBUG_EXPR(vold->display());
    DEBUG_EXPR(Mold->display());

    // Lagrangian Nonlinear Systems
    if ((dsType == Type::LagrangianDS) or dsType == Type::LagrangianLinearTIDS)
    {
      if (d->forces())
      {
        d->computeForces(told, qold, vold);
        DEBUG_EXPR(d->forces()->display());


        *workFree += *(d->forces());
      }
    }
    // // Lagrangian Linear Systems
    // else if (dsType == Type::LagrangianLinearTIDS)
    // {
    //   SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS> (*it);

    //   SP::SiconosMatrix C = d->C(); // constant dissipation
    //   SP::SiconosMatrix K = d->K(); // constant stiffness
    //   SP::SiconosVector Fext = d->fExt(); // time dependent force

    //   if (K)
    //   {
    //     prod(*K, *qold, *workFree, false);
    //   }

    //   if (C)
    //   {
    //     prod(*C, *vold, *workFree, false);
    //   }

    //   *workFree *= -1.;

    //   if (Fext)
    //   {
    //     d->computeFExt(told);
    //     *workFree += *Fext;
    //   }
    // }
    else
    {
      RuntimeException::selfThrow("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    }

    Mold->PLUForwardBackwardInPlace(*workFree); // contains left (right limit) acceleration without contact force

    DEBUG_PRINT("workFree contains right limit acceleration at  t^+_k without contact force :\n")
    DEBUG_EXPR(workFree->display());
  }

  // solve a LCP at acceleration level

  // DEBUG_EXPR(std::cout<< "allOSNS->empty()   " << std::boolalpha << allOSNS->empty() << std::endl << std::endl);
  // DEBUG_EXPR(std::cout<< "allOSNS->size()   "  << allOSNS->size() << std::endl << std::endl);




  if (!allOSNS->empty())
  {
    for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
    {
      (*it)->computeJach(told);
      (*it)->computeJacg(told);
    }

    if (simulationLink->model()->nonSmoothDynamicalSystem()->topology()->hasChanged())
    {
      for (OSNSIterator itOsns = allOSNS->begin(); itOsns != allOSNS->end(); ++itOsns)
      {
        (*itOsns)->setHasBeenUpdated(false);
      }
    }
    assert((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]);
    assert((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->interactions());

    if (!((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->interactions()->isEmpty()))
    {
      (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->compute(told);
      simulationLink->updateInput(2);
    }
  }





  DEBUG_PRINT("\nADVANCE TO RIGHT SIDE\n");
  // ADVANCE TO RIGHT
  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // type of the current DS
    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
    SP::SiconosVector workFree = d->workFree(); // contains acceleration without contact force

    // get left state from memory
    SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
    SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit
    SP::SiconosMatrix Mold = d->mass();

    // initialize *it->residuFree and predicted right velocity (left limit)
    SP::SiconosVector residuFree = (*it)->residuFree(); // POINTER CONSTRUCTOR : contains residu without nonsmooth effect
    SP::SiconosVector v = d->velocity(); // POINTER CONSTRUCTOR : contains velocity v_{k+1}^- and not free velocity
    residuFree->zero();
    v->zero();

    DEBUG_EXPR(workFree->display());
    DEBUG_EXPR(qold->display());
    DEBUG_EXPR(vold->display());
    DEBUG_EXPR(Mold->display());


    if (d->p(2))
    {

      scal(-0.5 * h, *(d->p(2)), *residuFree, false);
      scal(h, *(d->p(2)), *v, false);

      Mold->PLUForwardBackwardInPlace(*residuFree);
      Mold->PLUForwardBackwardInPlace(*v);
      DEBUG_EXPR(d->p(2)->display());

      DEBUG_EXPR(residuFree->display());
      DEBUG_EXPR(v->display());

    }

    *residuFree -= 0.5 * h**workFree;

    *v += h**workFree;
    *v += *vold;

    DEBUG_EXPR(residuFree->display());
    DEBUG_EXPR(v->display());

    SP::SiconosVector q = d->q(); // POINTER CONSTRUCTOR : contains position q_{k+1}
    *q = *qold;

    DEBUG_EXPR(q->display());

    scal(0.5 * h, *vold + *v, *q, false);

    /** At this step, we obtain
     * \f[
     * \begin{cases}
     * v_{k,0} = \mbox{\tt vold} \\
     * q_{k,0} = qold \\
     * F_{k,+} = F(told,qol,vold) \\
     * R_{free} = -h/2 * M^{-1}_k (P^+_{2,k}+F^+_{k}) \\
     * v_{k,1} = v_{k,0} + h * M^{-1}_k (P^+_{2,k}+F^+_{k})  \mbox{stored in v} \\
     * q_{k,1} = q_{k,0} + \frac{h}{2} (v_{k,0} + v_{k,1}) \mbox{stored in q} \\
     * \end{cases}
     * \f]
     **/

    DEBUG_EXPR(Mold->display());
    DEBUG_EXPR(q->display());
  }

  DEBUG_PRINT("\nRIGHT SIDE\n");
  // -- RIGHT SIDE --
  // calculate acceleration without contact force
  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // type of the current DS
    Type::Siconos dsType = Type::value(**it);
    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
    SP::SiconosVector workFree = d->workFree(); // POINTER CONSTRUCTOR : contains acceleration without contact force
    workFree->zero();
    // get right state from memory
    SP::SiconosVector q = d->q(); // contains position q_{k+1}
    SP::SiconosVector v = d->velocity(); // contains velocity v_{k+1}^- and not free velocity
    SP::SiconosMatrix M = d->mass(); // POINTER CONSTRUCTOR : contains mass matrix

    DEBUG_EXPR(workFree->display());
    DEBUG_EXPR(q->display());
    DEBUG_EXPR(v->display());
    // Lagrangian Nonlinear Systems
    if (dsType == Type::LagrangianDS or dsType == Type::LagrangianLinearTIDS)
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
        ("D1MinusLinear::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    // // Lagrangian Linear Systems
    // else if (dsType == Type::LagrangianLinearTIDS)
    // {
    //   SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS> (*it);

    //   SP::SiconosMatrix C = d->C(); // constant dissipation
    //   SP::SiconosMatrix K = d->K(); // constant stiffness
    //   SP::SiconosVector Fext = d->fExt(); // time-dependent force

    //   if (K)
    //   {
    //     prod(*K, *q, *workFree, false);
    //   }

    //   if (C)
    //   {
    //     prod(*C, *v, *workFree, false);
    //   }

    //   *workFree *= -1.;

    //   if (Fext)
    //   {
    //     d->computeFExt(t);
    //     *workFree += *Fext;
    //   }
    // }

    //DEBUG_EXPR(M->display());

    M->PLUForwardBackwardInPlace(*workFree); // contains right (left limit) acceleration without contact force
    DEBUG_PRINT("workFree contains left limit acceleration at  t^-_{k+1} without contact force :\n")
    DEBUG_EXPR(workFree->display());
    DEBUG_EXPR(q->display());
    DEBUG_EXPR(v->display());
    //DEBUG_EXPR(M->display());
  }

  // solve a LCP at acceleration level only for contacts which have been active at the beginning of the time-step
  if (!allOSNS->empty())
  {
    for (unsigned int level = simulationLink->levelMinForOutput(); level < simulationLink->levelMaxForOutput(); level++)
    {
      simulationLink->updateOutput(level);
    }
    // special update to consider only contacts which have been active at the beginning of the time-step
    //
    // Maurice Bremond: indices must be recomputed
    // as we deal with dynamic graphs, vertices and edges are stored
    // in lists for fast add/remove during updateIndexSet(i)
    // we need indices of list elements to build the OSNS Matrix so we
    // need an update if graph has changed
    // this should be done in updateIndexSet(i) for all integrators only
    // if a graph has changed
    simulationLink->updateIndexSets();
    //simulationLink->model()->nonSmoothDynamicalSystem()->topology()->indexSet(1)->update_vertices_indices();
    //simulationLink->model()->nonSmoothDynamicalSystem()->topology()->indexSet(1)->update_edges_indices();
    //simulationLink->model()->nonSmoothDynamicalSystem()->topology()->indexSet(2)->update_vertices_indices();
    //simulationLink->model()->nonSmoothDynamicalSystem()->topology()->indexSet(2)->update_edges_indices();

    for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
    {
      (*it)->computeJach(t);
      (*it)->computeJacg(t);
    }

    if (simulationLink->model()->nonSmoothDynamicalSystem()->topology()->hasChanged())
    {
      for (OSNSIterator itOsns = allOSNS->begin(); itOsns != allOSNS->end(); ++itOsns)
      {
        (*itOsns)->setHasBeenUpdated(false);
      }
    }

    if (!((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->interactions())->isEmpty())
    {
      (*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]->compute(t);
      simulationLink->updateInput(2);
    }
  }

  DEBUG_PRINT("\nRESIDU\n");
  // -- RESIDU --
  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // type of the current DS
    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
    SP::SiconosVector workFree = d->workFree(); // contains acceleration without contact force

    // get right state from memory
    SP::SiconosMatrix M = d->mass();
    DEBUG_EXPR(M->display());

    // initialize *it->residuFree
    SP::SiconosVector residuFree = (*it)->residuFree(); // POINTER CONSTRUCTOR : contains residu without nonsmooth effect
    *residuFree -= 0.5 * h**workFree;

    DEBUG_EXPR(workFree->display());

    if (d->p(2))
    {
      SP::SiconosVector dummy(new SiconosVector(*(d->p(2)))); // value = contact force
      M->PLUForwardBackwardInPlace(*dummy);
      *residuFree -= 0.5 * h**dummy;

      DEBUG_EXPR(d->p(2)->display());
    }
    d->computeMass();
    M->resetLU();
    *residuFree = prod(*M, *residuFree);
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
    DEBUG_EXPR(residuFree->display());
  }
  DEBUG_PRINT("D1MinusLinear::computeResidu(), end\n");
  return 0.; // there is no Newton iteration and the residuum is assumed to vanish
}

void D1MinusLinear::computeFreeState()
{
  DEBUG_PRINT("\n D1MinusLinear::computeFreeState(), start\n");


  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // Lagrangian Systems
    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);

    // get left state from memory
    SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit
    DEBUG_EXPR(vold->display());

    // get right information
    SP::SiconosMatrix M = d->mass();
    SP::SiconosVector vfree = d->velocity(); // POINTER CONSTRUCTOR : contains free velocity
    (*vfree) = *(d->residuFree());
    DEBUG_EXPR(d->residuFree()->display());
    d->computeMass();
    M->resetLU();
    M->PLUForwardBackwardInPlace(*vfree);
    // DEBUG_EXPR(M->display());

    *vfree *= -1.;
    *vfree += *vold;
    DEBUG_EXPR(vfree->display());

  }


  DEBUG_PRINT("D1MinusLinear::computeFreeState(), end\n");


}

void D1MinusLinear::computeFreeOutput(SP::Interaction inter, OneStepNSProblem* osnsp)
{
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
  if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
  {
    //Xfree = inter->dataX();
    if  (relationType == Lagrangian)
    {
      Xfree = inter->data(LagrangianR::x);
    }
  }
  else if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
  {

    if  (relationType == Lagrangian)
    {
      Xfree = inter->data(LagrangianR::free);
    }
    // else if  (relationType == NewtonEuler)
    // {
    //   Xfree = inter->data(NewtonEulerR::free);
    // }

    assert(Xfree);
  }
  else
    RuntimeException::selfThrow("D1MinusLinear::computeFreeOutput - OSNSP neither on velocity nor on acceleration level.");

  // calculate data of interaction
  SP::Interaction mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());

  // only Lagrangian Systems
  if (relationType == Lagrangian)
  {
    // in Yp the linear part of velocity or acceleration relation will be saved
    C = std11::static_pointer_cast<LagrangianR>(mainInteraction->relation())->C();

    if (C)
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

    if (relationSubType == RheonomousR) // explicit time dependence -> partial time derivative has to be added
    {
      if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
      {
        std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->computehDot(simulation()->getTkp1(), *inter);
        subprod(*ID, *(std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->hDot()), *Yp, xcoord, false);
      }
      else
        RuntimeException::selfThrow("D1MinusLinear::computeFreeOutput is only implemented  at velocity level for LagrangianRheonomousR.");
    }
    if (relationSubType == ScleronomousR) // acceleration term involving Hesse matrix of Relation with respect to degree is added
    {
      if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
      {
        std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->computeNonLinearH2dot(simulation()->getTkp1(), *inter);
        subprod(*ID, *(std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->Nonlinearh2dot()), *Yp, xcoord, false);
      }
    }
    if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp) // impact terms are added
    {
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }
  }
  else
    RuntimeException::selfThrow("D1MinusLinear::computeFreeOutput - not implemented for Relation of type " + relationType);
}

void D1MinusLinear::updateState(const unsigned int level)
{
  DEBUG_PRINTF("\n D1MinusLinear::updateState(const unsigned int level) start for level = %i\n",level);

  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // Lagrangian Systems
    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
    SP::SiconosMatrix M = d->mass();
    SP::SiconosVector v = d->velocity(); // POINTER CONSTRUCTOR : contains new velocity

    if (d->p(1))
    {
      SP::SiconosVector dummy(new SiconosVector(*(d->p(1)))); // value = nonsmooth impulse
      M->PLUForwardBackwardInPlace(*dummy); // solution for its velocity equivalent
      *v += *dummy; // add free velocity
      DEBUG_PRINT("\nRIGHT IMPULSE\n");
      DEBUG_PRINTF("%f\n", (*(d->p(1)))(0));
      //DEBUG_EXPR(throw(1));
    }
    DEBUG_EXPR(d->q()->display());
    DEBUG_EXPR(d->velocity()->display());
  }

  DEBUG_PRINT("\n D1MinusLinear::updateState(const unsigned int level) end\n");

}


bool D1MinusLinear::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert((i == 1) || (i==2));
  double h = simulationLink->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity

  double gamma = 1.0 / 2.0;
  // if (_useGamma)
  // {
  //   gamma = _gamma;
  // }
  DEBUG_PRINTF("D1MinusLinear::addInteractionInIndexSet of level = %i yref=%e, yDot=%e, y_estimated=%e.\n", i,  y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;

#if __cplusplus >= 201103L
  assert(!::isnan(y));
#else
  assert(!isnan(y));
#endif

  DEBUG_EXPR(
    if (y <= 0)
      DEBUG_PRINT("D1MinusLinear::addInteractionInIndexSet ACTIVATE.\n");
    );
  return (y <= 0.0);
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
    if (y > 0)
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
