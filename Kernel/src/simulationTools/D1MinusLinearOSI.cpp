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

#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include "debug.h"


#define D1MINUSLINEAR_STD
/**#define D1MINUSLINEAR_FULL */ /** In this version, we keep the contribution of lambda(2)
                                  * even if there is an impact. The unique modification is in the computeResidu
                                  * method that is redevelopped.
                                  */

/// @cond
using namespace RELATION;

/// @endcond


void D1MinusLinearOSI::_NSLEffectOnFreeOutput::visit(const NewtonImpactNSL& nslaw)
{
  double e = nslaw.e();
  Index subCoord(4);
  subCoord[0] = 0;
  subCoord[1] = _inter->nonSmoothLaw()->size();
  subCoord[2] = 0;
  subCoord[3] = subCoord[1];
  subscal(e, *(_inter->y_k(_osnsp->inputOutputLevel())), *(_inter->yForNSsolver()), subCoord, false);
}


D1MinusLinearOSI::D1MinusLinearOSI() :
  OneStepIntegrator(OSI::D1MINUSLINEAROSI), _typeOfD1MinusLinearOSI(halfexplicit_acceleration_level) {}

D1MinusLinearOSI::D1MinusLinearOSI(unsigned int type) :
  OneStepIntegrator(OSI::D1MINUSLINEAROSI)
{
  setTypeOfD1MinusLinearOSI(type);
}


D1MinusLinearOSI::D1MinusLinearOSI(SP::DynamicalSystem newDS) :
  OneStepIntegrator(OSI::D1MINUSLINEAROSI), _typeOfD1MinusLinearOSI(halfexplicit_acceleration_level)
{
  OSIDynamicalSystems->insert(newDS);
}


void D1MinusLinearOSI::setTypeOfD1MinusLinearOSI(unsigned int type)
{
  if (type < numberOfTypeOfD1MinusLinearOSI )
  {
    _typeOfD1MinusLinearOSI = type;
  }
  else
  {
    RuntimeException::selfThrow("D1MinusLinearOSI::setTypeOfD1MinusLinearOSI");
  }
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
  case halfexplicit_acceleration_level:
    return computeResiduHalfExplicitAccelerationLevel();
  case halfexplicit_acceleration_level_full:
    return computeResiduHalfExplicitAccelerationLevelFull();
  case explicit_velocity_level:
    RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu() - not implemented for type of D1MinusLinearOSI: " + explicit_velocity_level);
    break;
  case halfexplicit_velocity_level:
    return computeResiduHalfExplicitVelocityLevel();
  }
  DEBUG_PRINT("D1MinusLinearOSI::computeResidu() ends\n");
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
  DEBUG_PRINT("\n D1MinusLinearOSI::computeFreeOutput(), start\n");
  switch (_typeOfD1MinusLinearOSI)
  {
  case halfexplicit_acceleration_level:
    return computeFreeOutputHalfExplicitAccelerationLevel(vertex_inter,osnsp);
  case halfexplicit_acceleration_level_full:
    return computeFreeOutputHalfExplicitAccelerationLevel(vertex_inter,osnsp);
  case explicit_velocity_level:
    RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu() - not implemented for type of D1MinusLinearOSI: " + explicit_velocity_level);
    break;
  case halfexplicit_velocity_level:
    RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu() - not implemented for type of D1MinusLinearOSI: " + halfexplicit_velocity_level);
    break;
  }
  DEBUG_PRINT("D1MinusLinearOSI::computeFreeOutput() ends\n");
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



double D1MinusLinearOSI::computeResiduHalfExplicitAccelerationLevelFull()
{
RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu_explicit_acceleration_level has been removed due to obsolescence");
return 0.0;
}
