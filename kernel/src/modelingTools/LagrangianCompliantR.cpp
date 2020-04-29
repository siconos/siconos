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

#include "LagrangianCompliantR.hpp"
#include "SiconosAlgebraProd.hpp" // for matrix-vector prod
#include "Interaction.hpp"
#include "LagrangianDS.hpp"

#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

using namespace RELATION;

// constructor from a set of data
LagrangianCompliantR::LagrangianCompliantR(const std::string& pluginh, const std::string& pluginJacobianhq, const std::string& pluginJacobianhlambda) : LagrangianR(CompliantR)
{
  _zeroPlugin();
  setComputehFunction(SSLH::getPluginName(pluginh), SSLH::getPluginFunctionName(pluginh));
  _pluginJachq->setComputeFunction(SSLH::getPluginName(pluginJacobianhq), SSLH::getPluginFunctionName(pluginJacobianhq));
  _pluginJachlambda->setComputeFunction(SSLH::getPluginName(pluginJacobianhlambda), SSLH::getPluginFunctionName(pluginJacobianhlambda));
}

void LagrangianCompliantR::_zeroPlugin()
{
  _pluginJachq.reset(new PluggedObject());
  _pluginJachlambda.reset(new PluggedObject());
}


void LagrangianCompliantR::initialize(Interaction& inter)
{
  LagrangianR::initialize(inter);
  unsigned int sizeY = inter.dimension();

  if(! _jachlambda)
    _jachlambda.reset(new SimpleMatrix(sizeY, sizeY));
  else
    _jachlambda->resize(sizeY, sizeY);
}
void LagrangianCompliantR::checkSize(Interaction& inter)
{
}
void LagrangianCompliantR::computeh(double time, const BlockVector& q0, const SiconosVector& lambda, BlockVector& z, SiconosVector& y)
{

  if(_pluginh->fPtr)
  {
    auto qp = q0.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FPtr2)(_pluginh->fPtr))(qp->size(), &(*qp)(0), y.size(), lambda.getArray(), &(y)(0), zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void LagrangianCompliantR::computeJachq(double time, const BlockVector& q0, const SiconosVector& lambda, BlockVector& z)
{

  if(_jachq && _pluginJachq->fPtr)
  {
    auto qp = q0.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FPtr2)(_pluginJachq->fPtr))(qp->size(), &(*qp)(0), lambda.size(), lambda.getArray(), &(*_jachq)(0, 0), zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void LagrangianCompliantR::computeJachlambda(double time, const BlockVector& q0, const SiconosVector& lambda, BlockVector& z)
{

  if(_jachlambda && _pluginJachlambda->fPtr)
  {
    auto qp = q0.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FPtr2)_pluginJachlambda->fPtr)(qp->size(), &(*qp)(0), lambda.size(), lambda.getArray(), &(*_jachlambda)(0, 0), zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void LagrangianCompliantR::computeOutput(double time, Interaction& inter, unsigned int derivativeNumber)
{
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  if(derivativeNumber == 0)
  {
    SiconosVector& y = *inter.y(0);
    SiconosVector& lambda = *inter.lambda(0);
    computeh(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z], y);
  }
  else
  {
    SiconosVector& y = *inter.y(derivativeNumber);
    SiconosVector& lambda = *inter.lambda(derivativeNumber);
    computeJachq(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z]);
    computeJachlambda(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z]);
    if(derivativeNumber == 1)
    {
      // y = Jach[0] q1 + Jach[1] lambda
      prod(*_jachq, *DSlink[LagrangianR::q1], y);
      prod(*_jachlambda, lambda, y, false);
    }
    else if(derivativeNumber == 2)
      prod(*_jachq, *DSlink[LagrangianR::q2], y); // Approx: y[2] = Jach[0]q[2], other terms are neglected ...
    else
      RuntimeException::selfThrow("LagrangianCompliantR::computeOutput, index out of range or not yet implemented.");
  }
}

void LagrangianCompliantR::computeInput(double time, Interaction& inter, unsigned int level)
{
  // get lambda of the concerned interaction

  SiconosVector& lambda = *inter.lambda(level);
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  computeJachq(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z]);
  // data[name] += trans(G) * lambda
  prod(lambda, *_jachq, *DSlink[LagrangianR::p0 + level], false);
}

void LagrangianCompliantR::computeJach(double time, Interaction& inter)
{
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  SiconosVector& lambda = *inter.lambda(0);
  computeJachq(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z]);
  computeJachlambda(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z]);
}
