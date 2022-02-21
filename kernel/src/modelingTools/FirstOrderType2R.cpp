/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "FirstOrderType2R.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

FirstOrderType2R::FirstOrderType2R(const std::string& pluginh, const std::string& pluging):
  FirstOrderR(RELATION::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSLH::getPluginName(pluginh), SSLH::getPluginFunctionName(pluginh));
  setComputegFunction(SSLH::getPluginName(pluging), SSLH::getPluginFunctionName(pluging));
  // The jacobians are not set, and thus considered as null matrices at this point.
}

FirstOrderType2R::FirstOrderType2R(const std::string& pluginh, const std::string& pluging, const std::string& pluginJacobianhx, const std::string& pluginJacobianglambda):
  FirstOrderR(RELATION::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSLH::getPluginName(pluginh), SSLH::getPluginFunctionName(pluginh));
  setComputegFunction(SSLH::getPluginName(pluging), SSLH::getPluginFunctionName(pluging));

  setComputeJachxFunction(SSLH::getPluginName(pluginJacobianhx), SSLH::getPluginFunctionName(pluginJacobianhx));
  setComputeJacglambdaFunction(SSLH::getPluginName(pluginJacobianglambda), SSLH::getPluginFunctionName(pluginJacobianglambda));
}


void FirstOrderType2R::initialize(Interaction& inter)
{
  FirstOrderR::initialize(inter);

  unsigned int sizeY = inter.dimension();
  unsigned int sizeDS = inter.getSizeOfDS();
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  unsigned int sizeZ = DSlink[FirstOrderR::z]->size();
  VectorOfSMatrices& relationMat = inter.relationMatrices();


  if(!_C)
    relationMat[FirstOrderR::mat_C] = std::make_shared<SimpleMatrix>(sizeY, sizeDS);
  if(!_D)
    relationMat[FirstOrderR::mat_D] = std::make_shared<SimpleMatrix>(sizeY, sizeY);
  if(!_F)
    relationMat[FirstOrderR::mat_F] = std::make_shared<SimpleMatrix>(sizeY, sizeZ);
  if(!_B)
    relationMat[FirstOrderR::mat_B] = std::make_shared<SimpleMatrix>(sizeDS, sizeY);
  if(!_K)
    relationMat[FirstOrderR::mat_K] = std::make_shared<SimpleMatrix>(sizeDS, sizeDS);

//  if (!_jacgx)
//  {
//    relationMat[FirstOrderR::mat_K].reset(new SimpleMatrix(sizeDS, sizeDS));
  // TODO add this back to workV of the DS -> needed for X partial NS
//  }
}

void FirstOrderType2R::computeh(double time, const BlockVector& x, const SiconosVector& lambda, SiconosVector& y)
{
  auto xp = x.prepareVectorForPlugin();
  ((Type2PtrH)(_pluginh->fPtr))(xp->size(), xp->getArray(), lambda.size(), lambda.getArray(), y.size(), y.getArray());
}

void FirstOrderType2R::computeg(double time, const SiconosVector& lambda, BlockVector& r)
{
  auto rp = r.prepareVectorForPlugin();
  ((Type2PtrG)(_pluging->fPtr))(lambda.size(), lambda.getArray(), rp->size(), rp->getArray());
  r = *rp;
}

void FirstOrderType2R::computeOutput(double time, Interaction& inter, unsigned int level)
{
  DEBUG_BEGIN("FirstOrderType2R::computeOutput \n");
  auto& DSlink = inter.linkToDSVariables();
  BlockVector& x = *DSlink[FirstOrderR::x];
  // copy into Siconos continuous memory vector
  SiconosVector& y = *inter.y(level);
  SiconosVector& lambda = *inter.lambda(level);
  computeh(time, x, lambda, y);
  DEBUG_EXPR(y.display());
  DEBUG_END("FirstOrderType2R::computeOutput \n");
}

void FirstOrderType2R::computeInput(double time, Interaction& inter, unsigned int level)
{
  DEBUG_BEGIN("FirstOrderType2R::computeInput \n");
  auto& DSlink = inter.linkToDSVariables();
  // copy into Siconos continuous memory vector
  SiconosVector& lambda = *inter.lambda(level);
  computeg(time, lambda, *DSlink[FirstOrderR::r]);
  DEBUG_EXPR(DSlink[FirstOrderR::r]->display());
  DEBUG_END("FirstOrderType2R::computeInput \n");
}

void FirstOrderType2R::computeJachlambda(double time, const BlockVector& x, const SiconosVector& lambda, SimpleMatrix& D)
{
  THROW_EXCEPTION("FirstOrderType2R::computeJachlambda must be overload.");
}
void FirstOrderType2R::computeJachx(double time, const BlockVector& x, const SiconosVector& lambda, SimpleMatrix& C)
{
  THROW_EXCEPTION("FirstOrderType2R::computeJachx must be overload.");
  // Note FP: so this class should be virtual, isn't it?
}

void FirstOrderType2R::computeJach(double time, Interaction& inter)
{
  DEBUG_BEGIN("FirstOrderType2R::computeJach\n");
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  VectorOfSMatrices& relationMat = inter.relationMatrices();

  if(!_C)
  {
    computeJachx(time, *DSlink[FirstOrderR::x], *inter.lambda(0), *relationMat[FirstOrderR::mat_C]);
  }
  if(!_D)
  {
    computeJachlambda(time, *DSlink[FirstOrderR::x], *inter.lambda(0), *relationMat[FirstOrderR::mat_D]);
  }
  DEBUG_END("FirstOrderType2R::computeJach\n");
}

void FirstOrderType2R::computeJacglambda(double time, const SiconosVector& lambda, SimpleMatrix& B)
{
  THROW_EXCEPTION("FirstOrderType2R::computeJacglambda must be overload.");
}

void FirstOrderType2R::computeJacg(double time, Interaction& inter)
{
  DEBUG_BEGIN("FirstOrderType2R::computeJacg\n");
  if(!_B)
  {
    VectorOfSMatrices& relationMat = inter.relationMatrices();
    computeJacglambda(time, *inter.lambda(0), *relationMat[FirstOrderR::mat_B]);
  }
  DEBUG_END("FirstOrderType2R::computeJacg\n");
}
