/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "MoreauJeanCombinedProjectionOSI.hpp"
#include "Simulation.hpp"
#include "LagrangianDS.hpp"
#include "NewtonEulerDS.hpp"

#include "NewtonEulerR.hpp"
#include "LagrangianR.hpp"
#include "BlockVector.hpp"

#include "TypeName.hpp"
using namespace RELATION;


//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>

MoreauJeanCombinedProjectionOSI::MoreauJeanCombinedProjectionOSI(double theta) : MoreauJeanOSI(theta)
{
  _levelMinForOutput= 0;
  _levelMaxForOutput =1;
  _levelMinForInput =0;
  _levelMaxForInput =1;
  //_integratorType = OSI::MOREAUDIRECTPROJECTIONOSI;
}


void MoreauJeanCombinedProjectionOSI::initializeWorkVectorsForDS(double t, SP::DynamicalSystem ds)
{
  DEBUG_BEGIN("MoreauJeanCombinedProjectionOSI::initializeWorkVectorsForDS( double t, SP::DynamicalSystem ds) \n");
 
  MoreauJeanOSI::initializeWorkVectorsForDS(t, ds);
  
  const DynamicalSystemsGraph::VDescriptor& dsv = _dynamicalSystemsGraph->descriptor(ds);
  VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(dsv).workVectors;
  Type::Siconos dsType = Type::value(*ds);
  if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
  {
    SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
    workVectors[MoreauJeanOSI::QTMP].reset(new SiconosVector(d->dimension()));
  }
  else if(dsType == Type::NewtonEulerDS)
  {
    SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS>(ds);
    workVectors[MoreauJeanOSI::QTMP].reset(new SiconosVector(d->getqDim()));
  }
  else
  {
    RuntimeException::selfThrow("MoreauJeanCombinedProjectionOSI::initialize() - DS not of the right type");
  }
  for (unsigned int k = _levelMinForInput ; k < _levelMaxForInput + 1; k++)
  {
    ds->initializeNonSmoothInput(k);
  }
  
  DEBUG_END("MoreauJeanCombinedProjectionOSI::initializeWorkVectorsForDS( double t, SP::DynamicalSystem ds) \n");
}

void MoreauJeanCombinedProjectionOSI::initializeWorkVectorsForInteraction(Interaction &inter, InteractionProperties& interProp,
                                  DynamicalSystemsGraph & DSG)
{
  DEBUG_BEGIN("MoreauJeanCombinedProjectionOSI::initializeWorkVectorsForInteraction(Interaction &inter, InteractionProperties& interProp, DynamicalSystemsGraph & DSG)\n");

  MoreauJeanOSI::initializeWorkVectorsForInteraction(inter, interProp,DSG);

  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;
  assert(ds1);
  assert(ds2);
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();

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



  DEBUG_END("MoreauJeanCombinedProjectionOSI::initializeWorkVectorsForInteraction(Interaction &inter, InteractionProperties& interProp, DynamicalSystemsGraph & DSG)\n");


}


bool MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(i == 1 || i == 2);
  //double h = _simulation->timeStep();
  if(i == 1)  // index set for resolution at the velocity
  {
    double y = (inter->y(0))->getValue(0); // y(0) is the position
    DEBUG_PRINTF("MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet yref=%e \n", y);
    DEBUG_EXPR(
      if(y <= 0)
      printf("MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet ACTIVATE in indexSet level = %i.\n", i);
    )
      return (y <= 0);
  }
  else if(i == 2)   //  special index for the projection
  {
    DEBUG_EXPR(
      double lambda = 0;
      lambda = (inter->lambda(1))->getValue(0); // lambda(1) is the contact impulse for MoreauJeanOSI scheme
      printf("MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet lambdaref=%e \n", lambda);
      if(lambda > 0)
      printf("MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet ACTIVATE in indexSet level = %i.\n", i);
    )
      //    return (lambda > 0);
      return true;
  }
  else
  {
    return false;
  }
}


bool MoreauJeanCombinedProjectionOSI::removeInteractionFromIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(0);
  return(0);
}

