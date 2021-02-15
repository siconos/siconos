/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "LagrangianScleronomousR.hpp"
#include "SiconosAlgebraProd.hpp"  // for matrix-vector prod
#include "Interaction.hpp"
#include "LagrangianDS.hpp"

#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
#include "siconos_debug.h"




using namespace RELATION;

// constructor from a set of data
LagrangianScleronomousR::LagrangianScleronomousR(const std::string& pluginh, const std::string& pluginJacobianhq):
  LagrangianR(ScleronomousR)
{
  _zeroPlugin();
  setComputehFunction(SSLH::getPluginName(pluginh), SSLH::getPluginFunctionName(pluginh));

  _pluginJachq->setComputeFunction(pluginJacobianhq);

// Warning: we cannot allocate memory for Jach[0] matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
}
// constructor from a data used for EventDriven scheme
LagrangianScleronomousR::LagrangianScleronomousR(const std::string& pluginh, const std::string& pluginJacobianhq, const std::string& pluginDotJacobianhq):
  LagrangianR(ScleronomousR)
{
  _zeroPlugin();
  setComputehFunction(SSLH::getPluginName(pluginh), SSLH::getPluginFunctionName(pluginh));

  _pluginJachq->setComputeFunction(pluginJacobianhq);

  _plugindotjacqh->setComputeFunction(pluginDotJacobianhq);
}



void LagrangianScleronomousR::_zeroPlugin()
{
  LagrangianR::_zeroPlugin();
  _pluginJachq.reset(new PluggedObject());
  _plugindotjacqh.reset(new PluggedObject());
}

void LagrangianScleronomousR::initialize(Interaction& inter)
{
  if(!_jachq)
  {
    unsigned int sizeY = inter.dimension();
    unsigned int sizeDS = inter.getSizeOfDS();
    _jachq.reset(new SimpleMatrix(sizeY, sizeDS));
  }
}


void LagrangianScleronomousR::checkSize(Interaction& inter)
{

}

void LagrangianScleronomousR::computeh(const BlockVector& q, BlockVector& z, SiconosVector& y)
{
  DEBUG_PRINT(" LagrangianScleronomousR::computeh(Interaction& inter, SP::BlockVector q, SP::BlockVector z)\n");
  if(_pluginh && _pluginh->fPtr)
  {
    auto qp = q.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FPtr3)(_pluginh->fPtr))(qp->size(), &(*qp)(0), y.size(), &(y(0)), zp->size(), &(*zp)(0));
    z = *zp;
    DEBUG_EXPR(y.display());

  }
  // else nothing
}

void LagrangianScleronomousR::computeJachq(const BlockVector& q, BlockVector& z)
{
  if(_jachq && _pluginJachq->fPtr)
  {
    auto qp = q.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    // get vector lambda of the current interaction
    ((FPtr3)(_pluginJachq->fPtr))(qp->size(), &(*qp)(0), _jachq->size(0), &(*_jachq)(0, 0), zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void LagrangianScleronomousR::computeDotJachq(const BlockVector& q, BlockVector& z, const BlockVector& qDot)
{
  if(_dotjachq && _plugindotjacqh->fPtr)
  {
    auto qp = q.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    auto qdotp = qDot.prepareVectorForPlugin();
    ((FPtr2)(_plugindotjacqh->fPtr))(qp->size(), &(*qp)(0), qdotp->size(), &(*qdotp)(0), &(*_dotjachq)(0, 0), zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void  LagrangianScleronomousR::computedotjacqhXqdot(double time, Interaction& inter, VectorOfBlockVectors& DSlink)
{
  DEBUG_PRINT("LagrangianScleronomousR::computeNonLinearH2dot starts");
  // Compute the H Jacobian dot
  LagrangianScleronomousR::computeDotJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z], *DSlink[LagrangianR::q1]);
  _dotjacqhXqdot.reset(new SiconosVector(_dotjachq->size(0)));
  DEBUG_EXPR(_dotjachq->display(););
  prod(*_dotjachq, *DSlink[LagrangianR::q1], *_dotjacqhXqdot);
  DEBUG_PRINT("LagrangianScleronomousR::computeNonLinearH2dot ends");
}

void LagrangianScleronomousR::computeOutput(double time, Interaction& inter,  unsigned int derivativeNumber)
{

  DEBUG_PRINTF("LagrangianScleronomousR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber) with time = %f and derivativeNumber = %i\n", time, derivativeNumber);
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  SiconosVector& y = *inter.y(derivativeNumber);
  if(derivativeNumber == 0)
  {
    computeh(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z], y);
  }
  else
  {
    computeJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);

    if(derivativeNumber == 1)
    {
      assert(_jachq);
      prod(*_jachq, *DSlink[LagrangianR::q1], y);
    }
    else if(derivativeNumber == 2)
    {

      assert(_jachq);
      prod(*_jachq, *DSlink[LagrangianR::q2], y);
      if(!_dotjachq)
      {
        unsigned int sizeY = inter.dimension();
        unsigned int sizeDS = inter.getSizeOfDS();
        _dotjachq.reset(new SimpleMatrix(sizeY, sizeDS));
      }
      computeDotJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z], *DSlink[LagrangianR::q1]);
      prod(*_dotjachq, *DSlink[LagrangianR::q1], y, false);
    }
    else
      THROW_EXCEPTION("LagrangianScleronomousR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber), index out of range");
  }
}


void LagrangianScleronomousR::computeInput(double time, Interaction& inter, unsigned int level)
{
  DEBUG_BEGIN("void LagrangianScleronomousR::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level) \n");

  DEBUG_PRINTF("level = %i\n", level);
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  computeJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);
  // get lambda of the concerned interaction
  SiconosVector& lambda = *inter.lambda(level);
  DEBUG_EXPR(lambda.display(););
  DEBUG_EXPR(_jachq->display(););
  // data[name] += trans(G) * lambda
  prod(lambda, *_jachq, *DSlink[LagrangianR::p0 + level], false);
  DEBUG_EXPR(DSlink[LagrangianR::p0 + level]->display(););
  DEBUG_END("void LagrangianScleronomousR::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level) \n");
}

void LagrangianScleronomousR::computeJach(double time, Interaction& inter)
{
  DEBUG_BEGIN("void LagrangianScleronomousR::computeJach(double time, Interaction& inter) \n");
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  DEBUG_EXPR(inter.display(););
  computeJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);
  // computeJachqDot(time, inter);
  if(!_dotjachq)
  {
    unsigned int sizeY = inter.dimension();
    unsigned int sizeDS = inter.getSizeOfDS();
    _dotjachq.reset(new SimpleMatrix(sizeY, sizeDS));
  }
  computeDotJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z], *DSlink[LagrangianR::q1]);
  // computeJachlambda(time, inter);
  // computehDot(time,inter);
  DEBUG_END("void LagrangianScleronomousR::computeJach(double time, Interaction& inter) \n");
}
