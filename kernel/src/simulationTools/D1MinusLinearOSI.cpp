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

// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
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
  SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[D1MinusLinearOSI::OSNSP_RHS];
  subscal(e, osnsp_rhs, osnsp_rhs, subCoord, false);
}


D1MinusLinearOSI::D1MinusLinearOSI() :
  OneStepIntegrator(OSI::D1MINUSLINEAROSI), _typeOfD1MinusLinearOSI(halfexplicit_acceleration_level)
{
  _steps =2; //Two evaluations of lambda(2) are made for each time--step
  _levelMinForOutput= 0;
  _levelMaxForOutput =2;
  _levelMinForInput =1;
  _levelMaxForInput =2;
}

D1MinusLinearOSI::D1MinusLinearOSI(unsigned int type) :
  OneStepIntegrator(OSI::D1MINUSLINEAROSI)
{
  _steps =2; //Two evaluations of lambda(2) are made for each time--step
  _levelMinForOutput= 0;
  _levelMaxForOutput =2;
  _levelMinForInput =1;
  _levelMaxForInput =2;
  setTypeOfD1MinusLinearOSI(type);
}

void D1MinusLinearOSI::setTypeOfD1MinusLinearOSI(unsigned int type)
{
  if(type < numberOfTypeOfD1MinusLinearOSI)
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
  switch(_typeOfD1MinusLinearOSI)
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
void D1MinusLinearOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)
{
  // Get work buffers from the graph
  VectorOfVectors& workVectors = *_initializeDSWorkVectors(ds);

  // Check dynamical system type
  Type::Siconos dsType = Type::value(*ds);
  assert(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS || Type::NewtonEulerDS);
  
  if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
  {
    SP::LagrangianDS lds = std11::static_pointer_cast<LagrangianDS> (ds);
    lds->init_generalized_coordinates(2); // acceleration is required for the ds
    lds->init_inverse_mass(); // invMass required to update post-impact velocity

    workVectors.resize(OneStepIntegrator::work_vector_of_vector_size);
    workVectors[OneStepIntegrator::residu_free].reset(new SiconosVector(lds->dimension()));
    workVectors[OneStepIntegrator::free].reset(new SiconosVector(lds->dimension()));
    workVectors[OneStepIntegrator::free_tdg].reset(new SiconosVector(lds->dimension()));
    lds->swapInMemory();
  }
  else if(dsType == Type::NewtonEulerDS)
  {
    SP::NewtonEulerDS neds = std11::static_pointer_cast<NewtonEulerDS> (ds);
    neds->init_inverse_mass(); // invMass required to update post-impact velocity
    workVectors.resize(OneStepIntegrator::work_vector_of_vector_size);
    workVectors[OneStepIntegrator::residu_free].reset(new SiconosVector(neds->dimension()));
    workVectors[OneStepIntegrator::free].reset(new SiconosVector(neds->dimension()));
    workVectors[OneStepIntegrator::free_tdg].reset(new SiconosVector(neds->dimension()));
    neds->swapInMemory();
  }
  else
    RuntimeException::selfThrow("D1MinusLinearOSI::initialize - not implemented for Dynamical system type: " + dsType);

  for (unsigned int k = _levelMinForInput ; k < _levelMaxForInput + 1; k++)
  {
    ds->initializeNonSmoothInput(k);
  }

}

void D1MinusLinearOSI::initialize_nonsmooth_problems()
{
  SP::OneStepNSProblems allOSNSP  = _simulation->oneStepNSProblems(); // all OSNSP

  bool isOSNSPinitialized = false ;
  switch(_typeOfD1MinusLinearOSI)
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

  if(!isOSNSPinitialized)
  {
    RuntimeException::selfThrow("D1MinusLinearOSI::initialize() - not implemented for type of D1MinusLinearOSI: " + _typeOfD1MinusLinearOSI);
  }
}

void D1MinusLinearOSI::fillDSLinks(Interaction &inter,
				     InteractionProperties& interProp,
				     DynamicalSystemsGraph & DSG)
{
  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;
  assert(ds1);
  assert(ds2);

  VectorOfVectors& workV = *interProp.workVectors;
  workV.resize(D1MinusLinearOSI::WORK_INTERACTION_LENGTH);
  workV[D1MinusLinearOSI::OSNSP_RHS].reset(new SiconosVector(inter.getSizeOfY()));

  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  assert(interProp.DSlink);

  Relation &relation =  *inter.relation();
  RELATION::TYPES relationType = relation.getType();

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  bool computeResidu = relation.requireResidu();
  inter.initializeMemory(computeResidu,_steps);

  if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2))))
  {
    std::cout << "checkOSI(DSG.descriptor(ds1)): "
              << std::boolalpha
              << checkOSI(DSG.descriptor(ds1)) << std::endl;
    std::cout << "checkOSI(DSG.descriptor(ds2)): "
              << std::boolalpha
              << checkOSI(DSG.descriptor(ds2)) << std::endl;


    RuntimeException::selfThrow("D1MinusLinearOSI::fillDSLinks. The implementation is not correct for two different OSI for one interaction");
  }


  /* allocate ant set work vectors for the osi */
  VectorOfVectors &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
  if (relationType == Lagrangian)
  {
    LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS> (ds1);
    DSlink[LagrangianR::p2].reset(new BlockVector());
    DSlink[LagrangianR::p2]->insertPtr(lds.p(2));

    DSlink[LagrangianR::xfree].reset(new BlockVector());
    DSlink[LagrangianR::xfree]->insertPtr(workVds1[OneStepIntegrator::free]);
  }
  else if (relationType == NewtonEuler)
  {
    DSlink[NewtonEulerR::xfree].reset(new BlockVector());
    DSlink[NewtonEulerR::xfree]->insertPtr(workVds1[OneStepIntegrator::free]);
  }

  if (ds1 != ds2)
  {
    VectorOfVectors &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
    if (relationType == Lagrangian)
    {
      DSlink[LagrangianR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);
      LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS> (ds2);
      DSlink[LagrangianR::p2]->insertPtr(lds.p(2));

    }
    else if (relationType == NewtonEuler)
    {
      DSlink[NewtonEulerR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);
    }
  }
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

  switch(_typeOfD1MinusLinearOSI)
  {
  case halfexplicit_acceleration_level:
    DEBUG_END("D1MinusLinearOSI::computeResidu()\n");
    return computeResiduHalfExplicitAccelerationLevel();
  case halfexplicit_acceleration_level_full:
    DEBUG_END("D1MinusLinearOSI::computeResidu()\n");
    return computeResiduHalfExplicitAccelerationLevelFull();
  case halfexplicit_velocity_level:
    DEBUG_END("D1MinusLinearOSI::computeResidu()\n");
    return computeResiduHalfExplicitVelocityLevel();
  }
  RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu() - not implemented for type of D1MinusLinearOSI: " + _typeOfD1MinusLinearOSI);
  DEBUG_END("D1MinusLinearOSI::computeResidu()\n");
  return 1;
}



void D1MinusLinearOSI::computeFreeState()
{
  DEBUG_BEGIN("D1MinusLinearOSI::computeFreeState()\n");
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    Type::Siconos dsType = Type::value(*ds);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      // Lagrangian Systems
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);


      // get left state from memory
      SiconosVector& vold = *d->velocityMemory()->getSiconosVector(0); // right limit
      DEBUG_EXPR(vold.display());
      SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
      SiconosVector &vfree  = *d->velocity(); // POINTER CONSTRUCTOR : contains free velocity

      // get right information
      vfree =  residuFree;
      DEBUG_EXPR(residuFree.display());
      // d->computeMass();
      // M->resetLU();
      // M->PLUForwardBackwardInPlace(vfree);
      // DEBUG_EXPR(M->display());

      vfree *= -1.;
      vfree += vold;
      DEBUG_EXPR(vfree.display());
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      // NewtonEuler Systems
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);

      // get left state from memory
      SiconosVector& vold = *d->twistMemory()->getSiconosVector(0); // right limit
      DEBUG_EXPR(vold.display());


      // get right information
      //SP::SiconosMatrix M(new SimpleMatrix(*(d->mass()))); // we copy the mass matrix to avoid its factorization;
      SiconosVector &vfree = *d->twist(); // POINTER CONSTRUCTOR : contains free velocity
      SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];

      vfree = residuFree;
      DEBUG_EXPR(residuFree.display());

      vfree *= -1.;
      vfree += vold;
      DEBUG_EXPR(vfree.display());


    }
    else
      RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

  }


  DEBUG_END("D1MinusLinearOSI::computeFreeState()\n");


}

void D1MinusLinearOSI::updateState(const unsigned int)
{
  DEBUG_BEGIN("D1MinusLinearOSI::updateState(const unsigned int)\n");

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;

    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

    Type::Siconos dsType = Type::value(*ds);

    /* \warning the following conditional statement should be removed with a MechanicalDS class */
    /* Lagrangian DS*/
    if((dsType == Type::LagrangianDS) || (dsType == Type::LagrangianLinearTIDS))
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
      SP::SiconosVector v = d->velocity();

      DEBUG_PRINT("Position and velocity before update\n");
      DEBUG_EXPR(d->q()->display());
      DEBUG_EXPR(d->velocity()->display());

      /* Add the contribution of the impulse if any */
      if(d->p(1))
      {
        DEBUG_EXPR(d->p(1)->display());
        /* copy the value of the impulse */
        SP::SiconosVector dummy(new SiconosVector(*(d->p(1))));
        /* Compute the velocity jump due to the impulse */
	if(d->inverseMass())
	  {
	    d->update_inverse_mass();
	    d->inverseMass()->PLUForwardBackwardInPlace(*dummy);
	  }
        /* Add the velocity jump to the free velocity */
        *v += *dummy;
      }

      DEBUG_PRINT("Position and velocity after update\n");
      DEBUG_EXPR(d->q()->display());
      DEBUG_EXPR(d->velocity()->display());
    }
    /*  NewtonEuler Systems */
    else if(dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
      SP::SiconosVector v = d->twist(); // POINTER CONSTRUCTOR : contains new velocity
      if(d->p(1))
      {

        // Update the velocity
        SP::SiconosVector dummy(new SiconosVector(*(d->p(1)))); // value = nonsmooth impulse
	if(d->inverseMass())
	  {
	    d->update_inverse_mass();
	    d->inverseMass()->PLUForwardBackwardInPlace(*dummy);
	  }
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

  DEBUG_END("\n D1MinusLinearOSI::updateState(const unsigned int )\n");

}



void D1MinusLinearOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  DEBUG_PRINT("D1MinusLinearOSI::computeFreeOutput(), start\n");
  switch(_typeOfD1MinusLinearOSI)
  {
  case halfexplicit_acceleration_level:
    computeFreeOutputHalfExplicitAccelerationLevel(vertex_inter,osnsp);
    DEBUG_END("D1MinusLinearOSI::computeFreeOutput()\n");
    return;
  case halfexplicit_acceleration_level_full:
    computeFreeOutputHalfExplicitAccelerationLevel(vertex_inter,osnsp);
    DEBUG_END("D1MinusLinearOSI::computeFreeOutput()\n");
    return;
  case halfexplicit_velocity_level:
    computeFreeOutputHalfExplicitVelocityLevel(vertex_inter,osnsp);
    DEBUG_END("D1MinusLinearOSI::computeFreeOutput()\n");
    return;
  }
  RuntimeException::selfThrow("D1MinusLinearOSI::computeResidu() - not implemented for type of D1MinusLinearOSI: " + _typeOfD1MinusLinearOSI);
  DEBUG_END("D1MinusLinearOSI::computeFreeOutput()\n");
}



bool D1MinusLinearOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  DEBUG_BEGIN("D1MinusLinearOSI::addInteractionInIndexSet.\n");
  DEBUG_END("D1MinusLinearOSI::addInteractionInIndexSet.\n");
  switch(_typeOfD1MinusLinearOSI)
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
  switch(_typeOfD1MinusLinearOSI)
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
