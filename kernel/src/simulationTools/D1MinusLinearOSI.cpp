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

//#define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
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
//  subscal(e, *(_inter->y_k(_osnsp->inputOutputLevel())), *(_inter->yForNSsolver()), subCoord, false);
  subscal(e, *(_inter->yForNSsolver()), *(_inter->yForNSsolver()), subCoord, false);
}


D1MinusLinearOSI::D1MinusLinearOSI() :
  OneStepIntegrator(OSI::D1MINUSLINEAROSI), _typeOfD1MinusLinearOSI(halfexplicit_acceleration_level) {}

D1MinusLinearOSI::D1MinusLinearOSI(unsigned int type) :
  OneStepIntegrator(OSI::D1MINUSLINEAROSI)
{
  setTypeOfD1MinusLinearOSI(type);
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



unsigned int D1MinusLinearOSI::numberOfIndexSets() const
{
  switch (_typeOfD1MinusLinearOSI)
  {
  case halfexplicit_acceleration_level:
    return 4;
  case halfexplicit_acceleration_level_full:
    return 4;
  case halfexplicit_velocity_level:
    return 3;
  }
  RuntimeException::selfThrow("D1MinusLinearOSI::numberOfIndexSet - not implemented for D1minusLinear of type: " + _typeOfD1MinusLinearOSI);
  return 0;
}
void D1MinusLinearOSI::initialize(Model & m)
{
  DEBUG_BEGIN("D1MinusLinearOSI::initialize() \n");

  OneStepIntegrator::initialize(m);

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;

    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    Type::Siconos dsType = Type::value(*ds);
    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
      d->computeMass();
    }
    else if (dsType == Type::NewtonEulerDS)
    {
      //SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
    }
    else
      RuntimeException::selfThrow("D1MinusLinearOSI::initialize - not implemented for Dynamical system type: " + dsType);
  }

  DEBUG_PRINTF("D1MinusLinearOSI::initialize(). Type of OSI  %i ", _typeOfD1MinusLinearOSI );

  SP::OneStepNSProblems allOSNSP  = _simulation->oneStepNSProblems(); // all OSNSP

  bool isOSNSPinitialized = false ;
  switch (_typeOfD1MinusLinearOSI)
  {
  case halfexplicit_acceleration_level:
    // set evaluation levels (first is of velocity, second of acceleration type)
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY]->setIndexSetLevel(1);
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY]->setInputOutputLevel(1);
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY]->initialize(_simulation);

    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY + 1]->setIndexSetLevel(2);
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY + 1]->setInputOutputLevel(2);
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY + 1]->initialize(_simulation);
    isOSNSPinitialized = true ;
    DEBUG_EXPR((*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY + 1]->display());
    break;
  case halfexplicit_acceleration_level_full:
    // set evaluation levels (first is of velocity, second of acceleration type)
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY]->setIndexSetLevel(1);
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY]->setInputOutputLevel(1);
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY]->initialize(_simulation);

    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY + 1]->setIndexSetLevel(2);
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY + 1]->setInputOutputLevel(2);
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY + 1]->initialize(_simulation);
    isOSNSPinitialized = true ;
    break;
  case halfexplicit_velocity_level:
    // set evaluation levels (first is of velocity, second of acceleration type)
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY]->setIndexSetLevel(1);
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY]->setInputOutputLevel(1);
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY]->initialize(_simulation);

    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY + 1]->setIndexSetLevel(1); /** !!! */
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY + 1]->setInputOutputLevel(2);
    (*allOSNSP)[SICONOS_OSNSP_TS_VELOCITY + 1]->initialize(_simulation);
    isOSNSPinitialized = true ;
    break;
  }
  if (!isOSNSPinitialized)
  {
    RuntimeException::selfThrow("D1MinusLinearOSI::initialize() - not implemented for type of D1MinusLinearOSI: " + _typeOfD1MinusLinearOSI );
  }
  DEBUG_END("D1MinusLinearOSI::initialize() \n");
}



double D1MinusLinearOSI::computeResidu()
{

  DEBUG_PRINT("\n ******************************************************************\n");
  DEBUG_PRINT(" ******************************************************************\n");
  DEBUG_PRINT(" ******************************************************************\n");

  DEBUG_BEGIN("D1MinusLinearOSI::computeResidu()\n");


  DEBUG_PRINTF("nextTime %f\n", _simulation->nextTime());
  DEBUG_PRINTF("startingTime %f\n", _simulation->startingTime());
  DEBUG_PRINTF("time step size %f\n", _simulation->timeStep());

  switch (_typeOfD1MinusLinearOSI)
  {
  case halfexplicit_acceleration_level:
    return computeResiduHalfExplicitAccelerationLevel();
  case halfexplicit_acceleration_level_full:
    return computeResiduHalfExplicitAccelerationLevelFull();
  case halfexplicit_velocity_level:
    return computeResiduHalfExplicitVelocityLevel();
  }
  RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu() - not implemented for type of D1MinusLinearOSI: " + _typeOfD1MinusLinearOSI );
  DEBUG_END("D1MinusLinearOSI::computeResidu()\n");
  return 1;
}



void D1MinusLinearOSI::computeFreeState()
{
  DEBUG_BEGIN("D1MinusLinearOSI::computeFreeState()\n");
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

    Type::Siconos dsType = Type::value(*ds);
    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    if ((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      // Lagrangian Systems
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);

      // get left state from memory
      const SiconosVector& vold = d->velocityMemory().getSiconosVector(0); // right limit
      DEBUG_EXPR(vold->display());

      // get right information
      //SP::SiconosMatrix M = d->mass();
      SiconosVector& vfree = *d->velocity(); // POINTER CONSTRUCTOR : contains free velocity
      vfree = *(d->workspace(DynamicalSystem::freeresidu));
      DEBUG_EXPR(d->workspace(DynamicalSystem::freeresidu)->display());
      // d->computeMass();
      // M->resetLU();
      // M->PLUForwardBackwardInPlace(*vfree);
      // DEBUG_EXPR(M->display());

      vfree *= -1.;
      vfree += vold;
      DEBUG_EXPR(vfree.display());
    }
    else if (dsType == Type::NewtonEulerDS)
    {
      // NewtonEuler Systems
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);

      // get left state from memory
      const SiconosVector& vold = d->twistMemory().getSiconosVector(0); // right limit
      DEBUG_EXPR(vold.display());

      // get right information
      SP::SiconosMatrix M(new SimpleMatrix(*(d->mass()))); // we copy the mass matrix to avoid its factorization;
      SiconosVector& vfree = *d->twist(); // POINTER CONSTRUCTOR : contains free velocity
      vfree = *(d->workspace(DynamicalSystem::freeresidu));
      DEBUG_EXPR(d->workspace(DynamicalSystem::freeresidu)->display());

      vfree *= -1.;
      vfree += vold;
      DEBUG_EXPR(vfree.display());
    }
    else
      RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

  }


  DEBUG_END("D1MinusLinearOSI::computeFreeState()\n");


}

void D1MinusLinearOSI::updateState(const unsigned int level)
{
  DEBUG_BEGIN("D1MinusLinearOSI::updateState(const unsigned int level)\n");
  DEBUG_PRINTF("with level  = %i\n",level);

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;

    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

    Type::Siconos dsType = Type::value(*ds);

    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    /* Lagrangian DS*/
    if ((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
      SP::SiconosMatrix M = d->mass();
      SP::SiconosVector v = d->velocity();

      DEBUG_PRINT("Position and velocity before update\n");
      DEBUG_EXPR(d->q()->display());
      DEBUG_EXPR(d->velocity()->display());

      /* Add the contribution of the impulse if any */
      if (d->p(1))
      {
        DEBUG_EXPR(d->p(1)->display());
        /* copy the value of the impulse */
        SP::SiconosVector dummy(new SiconosVector(*(d->p(1))));
        /* Compute the velocity jump due to the impulse */
        M->PLUForwardBackwardInPlace(*dummy);
        /* Add the velocity jump to the free velocity */
        *v += *dummy;
      }

      DEBUG_PRINT("Position and velocity after update\n");
      DEBUG_EXPR(d->q()->display());
      DEBUG_EXPR(d->velocity()->display());
    }
    /*  NewtonEuler Systems */
    else if (dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
      SP::SiconosMatrix M(new SimpleMatrix(*(d->mass()))); // we copy the mass matrix to avoid its factorization;
      SP::SiconosVector v = d->twist(); // POINTER CONSTRUCTOR : contains new velocity
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

  DEBUG_END("\n D1MinusLinearOSI::updateState(const unsigned int level)\n");

}



void D1MinusLinearOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  DEBUG_PRINT("\n D1MinusLinearOSI::computeFreeOutput(), start\n");
  switch (_typeOfD1MinusLinearOSI)
  {
  case halfexplicit_acceleration_level:
    computeFreeOutputHalfExplicitAccelerationLevel(vertex_inter,osnsp);
    return;
  case halfexplicit_acceleration_level_full:
    computeFreeOutputHalfExplicitAccelerationLevel(vertex_inter,osnsp);
    return;
  case halfexplicit_velocity_level:
    computeFreeOutputHalfExplicitVelocityLevel(vertex_inter,osnsp);
    return;
  }
  RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu() - not implemented for type of D1MinusLinearOSI: " + _typeOfD1MinusLinearOSI);
  DEBUG_END("D1MinusLinearOSI::computeFreeOutput()\n");
}



bool D1MinusLinearOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  DEBUG_BEGIN("D1MinusLinearOSI::addInteractionInIndexSet.\n");
  DEBUG_END("D1MinusLinearOSI::addInteractionInIndexSet.\n");
  switch (_typeOfD1MinusLinearOSI)
  {
  case halfexplicit_acceleration_level:
    return addInteractionInIndexSetHalfExplicitAccelerationLevel(inter,i);
  case halfexplicit_velocity_level:
    return addInteractionInIndexSetHalfExplicitVelocityLevel(inter,i);
  }
  RuntimeException::selfThrow("D1MinusLinearOSI::addInteractionInIndexSet() - not implemented for type of D1MinusLinearOSI: " + _typeOfD1MinusLinearOSI);

  return 0;
}

bool D1MinusLinearOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  DEBUG_BEGIN("D1MinusLinearOSI::removeInteractionInIndexSet.\n");
  DEBUG_END("D1MinusLinearOSI::removeInteractionInIndexSet.\n");
  switch (_typeOfD1MinusLinearOSI)
  {
  case halfexplicit_acceleration_level:
    return removeInteractionInIndexSetHalfExplicitAccelerationLevel(inter,i);
  case halfexplicit_velocity_level:
    return removeInteractionInIndexSetHalfExplicitVelocityLevel(inter,i);
  }
  RuntimeException::selfThrow("D1MinusLinearOSI::removeInteractionInIndexSet() - not implemented for type of D1MinusLinearOSI: " + _typeOfD1MinusLinearOSI);
  return 0;
}


double D1MinusLinearOSI::computeResiduHalfExplicitAccelerationLevelFull()
{
RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu_explicit_acceleration_level has been removed due to obsolescence");
return 0.0;
}
