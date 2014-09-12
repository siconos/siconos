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


#define D1MINUSLINEAR_STD
/**#define D1MINUSLINEAR_FULL */ /** In this version, we keep the contribution of lambda(2)
                                  * even if there is an impact. The unique modification is in the computeResidu
                                  * method that is redevelopped.
                                  */

/// @cond
using namespace RELATION;

struct D1MinusLinearOSI::_NSLEffectOnFreeOutput : public SiconosVisitor
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
    subCoord[1] = _inter->nonSmoothLaw()->size();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    subscal(e, *(_inter->y_k(_osnsp->inputOutputLevel())), *(_inter->yForNSsolver()), subCoord, false);
  }
  void visit(const EqualityConditionNSL& nslaw)
  {
    ;
  }
};
/// @endcond

D1MinusLinearOSI::D1MinusLinearOSI() :
  OneStepIntegrator(OSI::D1MINUSLINEAROSI), _typeOfD1MinusLinearOSI(explicit_acceleration_level) {}


D1MinusLinearOSI::D1MinusLinearOSI(SP::DynamicalSystem newDS) :
  OneStepIntegrator(OSI::D1MINUSLINEAROSI)
{
  OSIDynamicalSystems->insert(newDS);
}


void D1MinusLinearOSI::initialize()
{
  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    Type::Siconos dsType = Type::value(**it);
    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*it);
      d->computeMass();
    }
    else if (dsType == Type::NewtonEulerDS)
    {
      //SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);

    }
    else
      RuntimeException::selfThrow("D1MinusLinearOSI::initialize - not implemented for Dynamical system type: " + dsType);
  }
}



double D1MinusLinearOSI::computeResidu()
{
  DEBUG_PRINT("\n D1MinusLinearOSI::computeResidu(), start\n");
  
  switch (_typeOfD1MinusLinearOSI)
  {
  case explicit_acceleration_level:
    return computeResiduExplicitAccelerationLevel();
  case explicit_acceleration_level_full:
    return computeResiduExplicitAccelerationLevelFull();
  case explicit_velocity_level:
    RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu() - not implemented for type of D1MinusLinearOSI: " + explicit_velocity_level);
    break;
  case halfexplicit_velocity_level:
    RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu() - not implemented for type of D1MinusLinearOSI: " + halfexplicit_velocity_level);
    break;
  }
  return 1;
}



void D1MinusLinearOSI::computeFreeState()
{
  DEBUG_PRINT("\n D1MinusLinearOSI::computeFreeState(), start\n");


  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // type of the current DS
    Type::Siconos dsType = Type::value(**it);
    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    if ((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
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
    else if (dsType == Type::NewtonEulerDS)
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
      RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

  }


  DEBUG_PRINT("D1MinusLinearOSI::computeFreeState(), end\n");


}

void D1MinusLinearOSI::updateState(const unsigned int level)
{
  DEBUG_PRINTF("\n D1MinusLinearOSI::updateState(const unsigned int level) start for level = %i\n",level);

  for (DSIterator it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    // type of the current DS
    Type::Siconos dsType = Type::value(**it);
    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    if ((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
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
        DEBUG_EXPR(d->p(1)->display());
      }
      DEBUG_EXPR(d->q()->display());
      DEBUG_EXPR(d->velocity()->display());
    }
    else if (dsType == Type::NewtonEulerDS)
    {
      // NewtonEuler Systems
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (*it);
      SP::SiconosMatrix M(new SimpleMatrix(*(d->mass()))); // we copy the mass matrix to avoid its factorization;
      SP::SiconosVector v = d->velocity(); // POINTER CONSTRUCTOR : contains new velocity
      if (d->p(1))
      {

        // Update the velocity
        SP::SiconosVector dummy(new SiconosVector(*(d->p(1)))); // value = nonsmooth impulse
        M->PLUForwardBackwardInPlace(*dummy); // solution for its velocity equivalent
        *v += *dummy; // add free velocity

        // update \f$ \dot q \f$
        SP::SiconosMatrix T = d->T();
        SP::SiconosVector dotq = d->dotq();
        prod(*T, *v, *dotq, true);

        DEBUG_PRINT("\nRIGHT IMPULSE\n");
        DEBUG_EXPR(d->p(1)->display());
      }
      DEBUG_EXPR(d->q()->display());
      DEBUG_EXPR(d->velocity()->display());
    }
    else
      RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

  }

  DEBUG_PRINT("\n D1MinusLinearOSI::updateState(const unsigned int level) end\n");

}

void D1MinusLinearOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{

  DEBUG_PRINT("D1MinusLinearOSI::computeFreeOutput starts\n");
  SP::OneStepNSProblems allOSNS  = simulationLink->oneStepNSProblems(); // all OSNSP
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
  SiconosVector& yForNSsolver = *inter->yForNSsolver();

  // define Xfree for velocity and acceleration level
  if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
  {
    //Xfree = inter->dataX();
    if (relationType == Lagrangian)
    {
      Xfree = DSlink[LagrangianR::q1];
    }
    else if (relationType == NewtonEuler)
    {
      Xfree = DSlink[NewtonEulerR::velocity];
    }
    else
      RuntimeException::selfThrow("D1MinusLinearOSI::computeFreeOutput - unknown relation type.");

  }
  else if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
  {

    if (relationType == Lagrangian)
    {
      Xfree = DSlink[LagrangianR::xfree];
    }
    else if (relationType == NewtonEuler)
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
  if (relationType == Lagrangian)
  {
    // in yForNSsolver the linear part of velocity or acceleration relation will be saved
    C = std11::static_pointer_cast<LagrangianR>(mainInteraction->relation())->C();

    if (C)
    {
      assert(Xfree);

      coord[3] = C->size(1);
      coord[5] = C->size(1);
      subprod(*C, *Xfree, yForNSsolver, coord, true);
    }

    // in yForNSsolver corrections have to be added
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
        SiconosVector q = *DSlink[LagrangianR::q0];
        SiconosVector z = *DSlink[LagrangianR::z];

        std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->computehDot(simulation()->getTkp1(), q, z);
        *DSlink[LagrangianR::z] = z;
        subprod(*ID, *(std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->hDot()), yForNSsolver, xcoord, false);
      }
      else
        RuntimeException::selfThrow("D1MinusLinearOSI::computeFreeOutput is only implemented  at velocity level for LagrangianRheonomousR.");
    }

    if (relationSubType == ScleronomousR) // acceleration term involving Hessian matrix of Relation with respect to degree is added
    {
      DEBUG_PRINT("D1MinusLinearOSI::computeFreeOutput. acceleration term involving Hessian matrix of Relation\n");
      DEBUG_EXPR(yForNSsolver.display(););

      if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
      {
        std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->computedotjacqhXqdot(simulation()->getTkp1(), *inter, DSlink);
        subprod(*ID, *(std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->dotjacqhXqdot()), yForNSsolver, xcoord, false);
      }
      DEBUG_EXPR(yForNSsolver.display(););
    }


  }

  else if (relationType == NewtonEuler)
  {
    SP::SiconosMatrix CT =  std11::static_pointer_cast<NewtonEulerR>(mainInteraction->relation())->jachqT();
    DEBUG_EXPR(CT->display());
    if (CT)
    {
      coord[3] = CT->size(1);
      coord[5] = CT->size(1);
      assert(Xfree);
      // creates a POINTER link between workX[ds] (xfree) and the
      // corresponding interactionBlock in each Interaction for each ds of the
      // current Interaction.
      // XXX Big quirks !!! -- xhub
      subprod(*CT, *Xfree, yForNSsolver, coord, true);
    }



    if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY + 1]).get() == osnsp)
    {
      // in yForNSsolver corrections have to be added
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
      DEBUG_EXPR(yForNSsolver.display(););

      /** Compute additional terms of the second order time derivative of constraints
       *
       * \f$ \nabla_q h(q) \dot T v + \frac{d}{dt}(\nabla_q h(q) ) T v \f$
       *
       */
      SP::DynamicalSystem ds1 = indexSet->properties(vertex_inter).source;
      SP::DynamicalSystem ds2 = indexSet->properties(vertex_inter).target;

      std11::static_pointer_cast<NewtonEulerR>(inter->relation())->computeSecondOrderTimeDerivativeTerms(simulation()->getTkp1(), *inter, DSlink, ds1, ds2);

      DEBUG_EXPR((std11::static_pointer_cast<NewtonEulerR>(inter->relation())->secondOrderTimeDerivativeTerms())->display());

      subprod(*ID, *(std11::static_pointer_cast<NewtonEulerR>(inter->relation())->secondOrderTimeDerivativeTerms()), yForNSsolver, xcoord, false);
      DEBUG_EXPR(yForNSsolver.display(););


    }
    DEBUG_EXPR(yForNSsolver.display(););
    if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp) // impact terms are added
    {
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }
  }

  else
    RuntimeException::selfThrow("D1MinusLinearOSI::computeFreeOutput - not implemented for Relation of type " + relationType);
  DEBUG_PRINT("D1MinusLinearOSI::computeFreeOutput ends\n");

}


bool D1MinusLinearOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  DEBUG_PRINT("D1MinusLinearOSI::addInteractionInIndexSet.\n");
  assert((i == 1) || (i==2));
  // double h = simulationLink->timeStep();

  double y = 0.0;
  double yOld =0.0;
  SP::Relation r = inter->relation();
  RELATION::TYPES relationType = r->getType();
  SP::LagrangianDS lds;
  if (relationType == Lagrangian)
  {

    // compute yold
    SP::BlockVector qoldB(new BlockVector());
    SP::BlockVector voldB(new BlockVector());
    SP::BlockVector zoldB(new BlockVector());
    SP::BlockVector qB(new BlockVector());
    SP::BlockVector vB(new BlockVector());
    SP::BlockVector zB(new BlockVector());

    for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
    {
      // check dynamical system type
      assert((Type::value(**it) == Type::LagrangianLinearTIDS ||
              Type::value(**it) == Type::LagrangianDS));

      // convert vDS systems into LagrangianDS and put them in vLDS
      lds = std11::static_pointer_cast<LagrangianDS> (*it);

      qoldB->insertPtr(lds->qMemory()->getSiconosVector(0));
      voldB->insertPtr(lds->velocityMemory()->getSiconosVector(0));
      /** \warning Warning the value of z of not stored. */
      zoldB->insertPtr(lds->z());
      qB->insertPtr(lds->q());
      vB->insertPtr(lds->velocity());
      zB->insertPtr(lds->z());
    }
    SiconosVector qold = *qoldB;
    SiconosVector zold = *zoldB;
    SiconosVector q = *qB;
    SiconosVector z = *zB;

    std11::static_pointer_cast<LagrangianScleronomousR>(r)->computeh(qold, zold, *inter->y(0));
    yOld = (inter->y(0))->getValue(0);
    // Compute current y (we assume that q stores q_{k,1} and v stores v_{k,1})
    // If not sure we have to store it into a new date in Interaction.
    std11::static_pointer_cast<LagrangianScleronomousR>(r)->computeh(q, z, *inter->y(0));
    y = (inter->y(0))->getValue(0);
  }

  DEBUG_PRINTF("D1MinusLinearOSI::addInteractionInIndexSet of level = %i yOld=%e, y=%e \n", i,  yOld, y);

  assert(!isnan(y));

  DEBUG_EXPR(
    if ((yOld >0.0) && (y <= y))
    DEBUG_PRINT("D1MinusLinearOSI::addInteractionInIndexSet contact are closing ((yOld >0.0) && (y <= y)).\n");
  );
  return ((yOld >0.0) && (y <= y));
}


bool D1MinusLinearOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
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
  DEBUG_PRINTF("D1MinusLinearOSI::removeInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
  assert(!isnan(y));
  DEBUG_EXPR(
    if (y > 0)
    DEBUG_PRINT("D1MinusLinearOSI::removeInteractionInIndexSet DEACTIVATE.\n");
  );

  return (y > 0.0);
}



double D1MinusLinearOSI::computeResiduExplicitAccelerationLevelFull()
{
RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu_explicit_acceleration_level has been removed due to obsolescence");
return 0.0;
}
