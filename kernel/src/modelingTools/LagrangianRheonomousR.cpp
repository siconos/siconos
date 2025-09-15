
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "LagrangianRheonomousR.hpp"
#include "SiconosAlgebraProd.hpp"  // for matrix-vector prod
#include "Interaction.hpp"
#include "LagrangianDS.hpp"

#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "siconos_debug.h"


using namespace RELATION;

// constructor from a set of data
LagrangianRheonomousR::LagrangianRheonomousR(const std::string& pluginh, const std::string& pluginJacobianhq, const std::string& pluginDoth):
  LagrangianR(RheonomousR)
{
  _zeroPlugin();
  // h
  setComputehFunction(SSLH::getPluginName(pluginh), SSLH::getPluginFunctionName(pluginh));

  _pluginJachq->setComputeFunction(pluginJacobianhq);

  // hDot
  setComputehDotFunction(SSLH::getPluginName(pluginDoth), SSLH::getPluginFunctionName(pluginDoth));
}

void LagrangianRheonomousR::initialize(Interaction& inter)
{
  if(!_jachq)
  {
    unsigned int sizeY = inter.dimension();
    unsigned int sizeDS = inter.getSizeOfDS();
    _jachq.reset(new SimpleMatrix(sizeY, sizeDS));
  }
}

void LagrangianRheonomousR::setComputehDotFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginhDot->setComputeFunction(pluginPath, functionName);
}

void LagrangianRheonomousR::_zeroPlugin()
{
  LagrangianR::_zeroPlugin();
  _pluginhDot.reset(new PluggedObject());
}

void LagrangianRheonomousR::computeh(double time, const BlockVector& q, BlockVector& z, SiconosVector& y)
{
  DEBUG_PRINT(" LagrangianRheonomousR::computeh(double time,Interaction& inter, SP::BlockVector q, SP::BlockVector z)");
  // arg= time. Unused in this function but required for interface.
  if(_pluginh->fPtr)
  {
    auto qp = q.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FPtr4)(_pluginh->fPtr))(qp->size(), &(*qp)(0), time, y.size(),  &(y)(0), zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void LagrangianRheonomousR::computehDot(double time, const BlockVector& q, BlockVector& z)
{
  if(_hDot && _pluginhDot->fPtr)
  {
    auto qp = q.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FPtr4)(_pluginhDot->fPtr))(qp->size(), &(*qp)(0), time, _hDot->size(),  &(*_hDot)(0), zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void LagrangianRheonomousR::computeJachq(double time,  const BlockVector& q, BlockVector& z)
{
  if(_jachq && _pluginJachq->fPtr)
  {
    auto qp = q.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FPtr4)(_pluginJachq->fPtr))(qp->size(), &(*qp)(0), time, _jachq->size(0), &(*_jachq)(0, 0), zp->size(), &(*zp)(0));
    z = *zp;
  }
}


void LagrangianRheonomousR::computeOutput(double time, Interaction& inter, unsigned int derivativeNumber)
{
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  SiconosVector& y = *inter.y(derivativeNumber);
  if(derivativeNumber == 0)
    computeh(time, *DSlink[LagrangianR::q0], *DSlink[LagrangianR::z], y);
  else
  {
    computeJachq(time, *DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);
    if(derivativeNumber == 1)
    {
      if(!_hDot)
      {
        unsigned int sizeY = inter.dimension();
        _hDot.reset(new SiconosVector(sizeY));
      }
      // Computation of the partial derivative w.r.t time of h(q,t)
      computehDot(time, *DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);
      assert(_jachq);
      // Computation of the partial derivative w.r.t q of h(q,t) : \nabla_q h(q,t) \dot q
      prod(*_jachq, *DSlink[LagrangianR::q1], y);
      // Sum of the terms
      y += *_hDot;
    }
    else if(derivativeNumber == 2)
    {
      assert(_jachq);
      prod(*_jachq, *DSlink[LagrangianR::q2], y); // Approx:,  ...
      // \warning : the computation of y[2] (in event-driven
      // simulation for instance) is approximated by y[2] =
      // Jach[0]q[2]. For the moment, other terms are neglected
      // (especially, partial derivatives with respect to time).
    }
    else
      THROW_EXCEPTION("LagrangianRheonomousR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber) index >2  not yet implemented.");
  }
}

void LagrangianRheonomousR::computeInput(double time, Interaction& inter,  unsigned int level)
{
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  computeJachq(time, *DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);
  // get lambda of the concerned interaction
  SiconosVector& lambda = *inter.lambda(level);
  // data[name] += trans(G) * lambda
  prod(lambda, *_jachq, *DSlink[LagrangianR::p0 + level], false);
}

void LagrangianRheonomousR::computeJach(double time, Interaction& inter)
{
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  computeJachq(time, *DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);
  // computeJachqDot(time, inter);
  //    computeDotJachq(time, q, z);
  // computeJachlambda(time, inter);
  computehDot(time, *DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);
}
