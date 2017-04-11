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

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "LagrangianCompliantR.hpp"
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

void LagrangianCompliantR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  LagrangianR::initComponents(inter, DSlink, workV, workM);
  unsigned int sizeY = inter.getSizeOfY();

  if (! _jachlambda)
    _jachlambda.reset(new SimpleMatrix(sizeY, sizeY));
  else
    _jachlambda->resize(sizeY, sizeY);
}

void LagrangianCompliantR::computeh(double time, SiconosVector& q0, SiconosVector& lambda, SiconosVector& z, SiconosVector& y)
{
  if (_pluginh->fPtr)
  {
    ((FPtr2)(_pluginh->fPtr))(q0.size(), &(q0)(0), y.size(), &(lambda)(0), &(y)(0), z.size(), &(z)(0));
  }
}

void LagrangianCompliantR::computeJachq(double time, SiconosVector& q0, SiconosVector& lambda, SiconosVector& z)
{

  if (_jachq && _pluginJachq->fPtr)
  {
    ((FPtr2)(_pluginJachq->fPtr))(q0.size(), &(q0)(0), lambda.size(), &(lambda)(0), &(*_jachq)(0, 0), z.size(), &(z)(0));
  }
}
void LagrangianCompliantR::computeJachlambda(double time, SiconosVector& q0, SiconosVector& lambda, SiconosVector& z)
{

  if (_jachlambda && _pluginJachlambda->fPtr)
  {
    ((FPtr2)_pluginJachlambda->fPtr)(q0.size(), &(q0)(0), lambda.size(), &(lambda)(0), &(*_jachlambda)(0, 0), z.size(), &(z)(0));
  }
}

void LagrangianCompliantR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber)
{
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  SiconosVector workZ = *DSlink[LagrangianR::z];
  if (derivativeNumber == 0)
  {
    SiconosVector& y = *inter.y(0);
    SiconosVector& lambda = *inter.lambda(0);
    SiconosVector workQ = *DSlink[LagrangianR::q0];

    computeh(time, workQ, lambda, workZ, y);
  }
  else
  {
    SiconosVector& y = *inter.y(derivativeNumber);
    SiconosVector& lambda = *inter.lambda(derivativeNumber);
    SiconosVector workQ = *DSlink[LagrangianR::q0];
    computeJachq(time, workQ, lambda, workZ);
    computeJachlambda(time, workQ, lambda, workZ);
    if (derivativeNumber == 1)
    {
      // y = Jach[0] q1 + Jach[1] lambda
      prod(*_jachq, *DSlink[LagrangianR::q1], y);
      prod(*_jachlambda, lambda, y, false);
    }
    else if (derivativeNumber == 2)
      prod(*_jachq, *DSlink[LagrangianR::q2], y); // Approx: y[2] = Jach[0]q[2], other terms are neglected ...
    else
      RuntimeException::selfThrow("LagrangianCompliantR::computeOutput, index out of range or not yet implemented.");
  }

  *DSlink[LagrangianR::z] = workZ;
}

void LagrangianCompliantR::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)
{
  // get lambda of the concerned interaction
  SiconosVector& lambda = *inter.lambda(level);
  VectorOfBlockVectors& DSlink = *interProp.DSlink;

  SiconosVector workQ = *DSlink[LagrangianR::q0];
  SiconosVector workZ = *DSlink[LagrangianR::z];
  computeJachq(time, workQ, lambda, workZ);
  // data[name] += trans(G) * lambda
  prod(lambda, *_jachq, *DSlink[LagrangianR::p0 + level], false);
  *DSlink[LagrangianR::z] = workZ;
}


void LagrangianCompliantR::computeJach(double time, Interaction& inter, InteractionProperties& interProp)
{
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  SiconosVector q = *DSlink[LagrangianR::q0];
  SiconosVector z = *DSlink[LagrangianR::z];
  SiconosVector& lambda = *inter.lambda(0);
  computeJachq(time, q, lambda, z);
  computeJachlambda(time, q, lambda, z);
}
