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
#include "MoreauJeanDirectProjectionOSI.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianDS.hpp"

#include "NewtonEulerR.hpp"
#include "LagrangianR.hpp"
#include "BlockVector.hpp"

#include "CxxStd.hpp"


#include "TypeName.hpp"
using namespace RELATION;
//#define STANDARD_ACTIVATION
#define FIRSTWAY_ACTIVATION
//#define SECONDWAY_ACTIVATION
//#define QFREE_ACTIVATION

// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>

#define SICONOS_MPC_DEFAULT_ACTIVATION_POS_THRESHOLD 1.e-7
#define SICONOS_MPC_DEFAULT_ACTIVATION_VEL_THRESHOLD 0.0
#define SICONOS_MPC_DEFAULT_DEACTIVATION_POS_THRESHOLD 1.e-7
#define SICONOS_MPC_DEFAULT_DEACTIVATION_VEL_THRESHOLD 0.0


MoreauJeanDirectProjectionOSI::MoreauJeanDirectProjectionOSI(double theta) : MoreauJeanOSI(theta)
{
  _levelMinForOutput= 0;
  _levelMaxForOutput =1;
  _levelMinForInput =0;
  _levelMaxForInput =1;
  _integratorType = OSI::MOREAUDIRECTPROJECTIONOSI;
  _deactivateYPosThreshold = SICONOS_MPC_DEFAULT_DEACTIVATION_POS_THRESHOLD;
  _deactivateYVelThreshold = SICONOS_MPC_DEFAULT_DEACTIVATION_VEL_THRESHOLD;
  _activateYPosThreshold =   SICONOS_MPC_DEFAULT_ACTIVATION_POS_THRESHOLD;
  _activateYVelThreshold =   SICONOS_MPC_DEFAULT_ACTIVATION_VEL_THRESHOLD;
}

MoreauJeanDirectProjectionOSI::MoreauJeanDirectProjectionOSI(double theta, double gamma) : MoreauJeanOSI(theta, gamma)
{
  _levelMinForOutput= 0;
  _levelMaxForOutput =1;
  _levelMinForInput =0;
  _levelMaxForInput =1;
  _integratorType = OSI::MOREAUDIRECTPROJECTIONOSI;
  _deactivateYPosThreshold = SICONOS_MPC_DEFAULT_DEACTIVATION_POS_THRESHOLD;
  _deactivateYVelThreshold = SICONOS_MPC_DEFAULT_DEACTIVATION_VEL_THRESHOLD;
  _activateYPosThreshold =   SICONOS_MPC_DEFAULT_ACTIVATION_POS_THRESHOLD;
  _activateYVelThreshold =   SICONOS_MPC_DEFAULT_ACTIVATION_VEL_THRESHOLD;
}

void MoreauJeanDirectProjectionOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)
{
  DEBUG_BEGIN("MoreauJeanDirectProjectionOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds) \n");
  MoreauJeanOSI::initializeDynamicalSystem(m, t, ds);
  const DynamicalSystemsGraph::VDescriptor& dsv = _dynamicalSystemsGraph->descriptor(ds);
  VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(dsv).workVectors;
  Type::Siconos dsType = Type::value(*ds);
  if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
  {
    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
    workVectors[OneStepIntegrator::qtmp].reset(new SiconosVector(d->dimension()));
  }
  else if(dsType == Type::NewtonEulerDS)
  {
    SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS>(ds);
    workVectors[OneStepIntegrator::qtmp].reset(new SiconosVector(d->getqDim()));
  }
  else
  {
    RuntimeException::selfThrow("MoreauJeanDirectProjectionOSI::initialize() - DS not of the right type");
  }
  for (unsigned int k = _levelMinForInput ; k < _levelMaxForInput + 1; k++)
  {
    DEBUG_PRINTF("ds->initializeNonSmoothInput(%i)\n", k);
    ds->initializeNonSmoothInput(k);
    DEBUG_EXPR_WE(
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
      if (d->p(k))
        std::cout << "d->p(" << k <<" ) exists" << std::endl;
      );

  }
  DEBUG_END("MoreauJeanDirectProjectionOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds) \n");

}

void MoreauJeanDirectProjectionOSI::fillDSLinks(Interaction &inter, InteractionProperties& interProp,
                                  DynamicalSystemsGraph & DSG)
{
  DEBUG_BEGIN("MoreauJeanDirectProjectionOSI::fillDSLinks(Interaction &inter, InteractionProperties& interProp, DynamicalSystemsGraph & DSG)\n");

  MoreauJeanOSI::fillDSLinks(inter, interProp,DSG);

  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;
  assert(ds1);
  assert(ds2);
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  Relation &relation =  *inter.relation();
  RELATION::TYPES relationType = relation.getType();


  unsigned int p0 =0;
  if (relationType == Lagrangian)
  {
    p0 = LagrangianR::p0;
  }
  else if (relationType == NewtonEuler)
  {
    p0 = NewtonEulerR::p0;
  }

  if (ds1 != ds2)
  {
    DEBUG_PRINT("ds1 != ds2\n");
    if ((!DSlink[p0]) || (DSlink[p0]->numberOfBlocks() !=2))
      DSlink[p0].reset(new BlockVector(2));
  }
  else
  {
    if ((!DSlink[p0]) || (DSlink[p0]->numberOfBlocks() !=1))
      DSlink[p0].reset(new BlockVector(1));
  }

  if(checkOSI(DSG.descriptor(ds1)))
  {
    DEBUG_PRINTF("ds1->number() %i is taken into account\n", ds1->number());
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    if (relationType == Lagrangian)
    {
      LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS> (ds1);
      DSlink[p0]->setVectorPtr(0,lds.p(0));
    }
    else if (relationType == NewtonEuler)
    {
      NewtonEulerDS& neds = *std11::static_pointer_cast<NewtonEulerDS> (ds1);
      DSlink[p0]->setVectorPtr(0,neds.p(0));
    }
  }
  DEBUG_PRINTF("ds1->number() %i\n",ds1->number());
  DEBUG_PRINTF("ds2->number() %i\n",ds2->number());

  if (ds1 != ds2)
  {
    DEBUG_PRINT("ds1 != ds2\n");
    if(checkOSI(DSG.descriptor(ds2)))
    {
      DEBUG_PRINTF("ds2->number() %i is taken into account\n",ds2->number());
      assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
      if (relationType == Lagrangian)
      {
        LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS> (ds2);
        DSlink[p0]->setVectorPtr(1,lds.p(0));
      }
      else if (relationType == NewtonEuler)
      {
        NewtonEulerDS& neds = *std11::static_pointer_cast<NewtonEulerDS> (ds2);
        DSlink[p0]->setVectorPtr(1,neds.p(0));
      }
    }
  }



  DEBUG_END("MoreauJeanDirectProjectionOSI::fillDSLinks(Interaction &inter, InteractionProperties& interProp, DynamicalSystemsGraph & DSG)\n");


}


void MoreauJeanDirectProjectionOSI::computeFreeState()
{
  MoreauJeanOSI::computeFreeState();
}

#ifdef STANDARD_ACTIVATION
bool MoreauJeanDirectProjectionOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  return MoreauJeanOSI::addInteractionInIndexSet(inter, i);
}

bool MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  return MoreauJeanOSI::removeInteractionInIndexSet(inter, i);
}
#endif


#ifdef FIRSTWAY_ACTIVATION
bool MoreauJeanDirectProjectionOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{

  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
  double gamma = 1.0 / 2.0;
  if(_useGamma)
  {
    gamma = _gamma;
  }
  DEBUG_PRINTF("MoreauJeanOSI::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;


  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::addInteractionInIndexSet yref=%e, yDot=%e.\n", y, yDot);

  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::addInteractionInIndexSet  _activateYPosThreshold =%e, _activateYVelThreshold=%e\n",
               _activateYPosThreshold ,
               _activateYVelThreshold);

  assert(!isnan(y));
#ifdef DEBUG_MESSAGES
  if(y <= _activateYPosThreshold)
    DEBUG_PRINT("MoreauJeanDirectProjectionOSI::addInteractionInIndexSet ACTIVATE.\n");
#endif
  return (y <= _activateYPosThreshold);
}

bool MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)

{
  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
  double gamma = 1.0 / 2.0;
  if(_useGamma)
  {
    gamma = _gamma;
  }
  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet yref=%e, yDot=%e .\n", y, yDot);
  y += gamma * h * yDot;

  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet  _deactivateYPosThreshold =%e, _deactivateYVelThreshold=%e\n",
               _deactivateYPosThreshold ,
               _deactivateYVelThreshold);

  assert(!isnan(y));
#ifdef DEBUG_MESSAGES
  if(y > _deactivateYPosThreshold && yDot >= _deactivateYVelThreshold)
    DEBUG_PRINT("MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet DEACTIVATE.\n");
#endif
  return (y > _deactivateYPosThreshold && yDot >= _deactivateYVelThreshold);
}

#endif



#ifdef SECONDWAY_ACTIVATION
bool MoreauJeanDirectProjectionOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{

  assert(i == 1);
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
#ifdef DEBUG_MESSAGES
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
#endif
  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::addInteractionInIndexSet yref=%e, yDot=%e.\n", y, yDot);

  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::addInteractionInIndexSet  _activateYPosThreshold =%e, _activateYVelThreshold=%e\n",
               _activateYPosThreshold ,
               _activateYVelThreshold);

  assert(!isnan(y));

  if(y <= _activateYPosThreshold)
    DEBUG_PRINT("MoreauJeanDirectProjectionOSI::addInteractionInIndexSet ACTIVATE.\n");
  return (y <= _activateYPosThreshold);
}

bool MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)

{
  assert(i == 1);
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
  double lambda = (inter->lambda(i))->getValue(0); // for i=1 y(i) is the velocity

  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet yref=%e, yDot=%e .\n", y, yDot);

  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet  _deactivateYPosThreshold =%e, _deactivateYVelThreshold=%e\n",
               _deactivateYPosThreshold ,
               _deactivateYVelThreshold);

  assert(!isnan(y));
  if(y > _deactivateYPosThreshold && lambda <= _deactivateYVelThreshold)
    DEBUG_PRINT("MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet DEACTIVATE.\n");
  return (y > _deactivateYPosThreshold && lambda <= _deactivateYVelThreshold);
}

#endif

#ifdef QFREE_ACTIVATION
bool MoreauJeanDirectProjectionOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{

  assert(i == 1);
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
#ifdef DEBUG_MESSAGES
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
#endif
  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::addInteractionInIndexSet yref=%e, yDot=%e.\n", y, yDot);

  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::addInteractionInIndexSet  _activateYPosThreshold =%e, _activateYVelThreshold=%e\n",
               _activateYPosThreshold ,
               _activateYVelThreshold);

  assert(!isnan(y));

  if(y <= _activateYPosThreshold)
    DEBUG_PRINT("MoreauJeanDirectProjectionOSI::addInteractionInIndexSet ACTIVATE.\n");
  return (y <= _activateYPosThreshold);
}

bool MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)

{
  assert(i == 1);
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity


  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet yref=%e, yDot=%e .\n", y, yDot);

  DEBUG_PRINTF("MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet  _deactivateYPosThreshold =%e, _deactivateYVelThreshold=%e\n",
               _deactivateYPosThreshold ,
               _deactivateYVelThreshold);

  assert(!isnan(y));
  if(y > _deactivateYPosThreshold)
    DEBUG_PRINT("MoreauJeanDirectProjectionOSI::removeInteractionInIndexSet DEACTIVATE.\n");
  return (y > _deactivateYPosThreshold);
}

#endif
