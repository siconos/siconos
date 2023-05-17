/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "BlockVector.hpp"
#include "FirstOrderR.hpp"
#include "Interaction.hpp"
#include "PluggedObject.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::modeling::FirstOrderType2R::FirstOrderType2R(const std::string& pluginh,
                                                      const std::string& pluging)
    : FirstOrderR(RelationSubType::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(siconos::plugins::getPluginName(pluginh),
                      siconos::plugins::getPluginFunctionName(pluginh));
  setComputegFunction(siconos::plugins::getPluginName(pluging),
                      siconos::plugins::getPluginFunctionName(pluging));
  // The jacobians are not set, and thus considered as null matrices at this point.
}

siconos::modeling::FirstOrderType2R::FirstOrderType2R(const std::string& pluginh,
                                                      const std::string& pluging,
                                                      const std::string& pluginJacobianhx,
                                                      const std::string& pluginJacobianglambda)
    : FirstOrderR(RelationSubType::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(siconos::plugins::getPluginName(pluginh),
                      siconos::plugins::getPluginFunctionName(pluginh));
  setComputegFunction(siconos::plugins::getPluginName(pluging),
                      siconos::plugins::getPluginFunctionName(pluging));

  setComputeJachxFunction(siconos::plugins::getPluginName(pluginJacobianhx),
                          siconos::plugins::getPluginFunctionName(pluginJacobianhx));
  setComputeJacglambdaFunction(siconos::plugins::getPluginName(pluginJacobianglambda),
                               siconos::plugins::getPluginFunctionName(pluginJacobianglambda));
}

void siconos::modeling::FirstOrderType2R::initialize(Interaction& inter)
{
  FirstOrderR::initialize(inter);

  auto sizeY = inter.dimension();
  auto sizeDS = inter.getSizeOfDS();
  auto& DSlink = inter.linkToDSVariables();
  auto sizeZ = DSlink[FirstOrderR::z]->size();
  auto& relationMat = inter.relationMatrices();

  if (!_C)
    relationMat[FirstOrderR::mat_C] =
        std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeDS);
  if (!_D)
    relationMat[FirstOrderR::mat_D] =
        std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeY);
  if (!_F)
    relationMat[FirstOrderR::mat_F] =
        std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeZ);
  if (!_B)
    relationMat[FirstOrderR::mat_B] =
        std::make_shared<siconos::algebra::SimpleMatrix>(sizeDS, sizeY);
  if (!_K)
    relationMat[FirstOrderR::mat_K] =
        std::make_shared<siconos::algebra::SimpleMatrix>(sizeDS, sizeDS);

  //  if (!_jacgx)
  //  {
  //    relationMat[FirstOrderR::mat_K] =
  //    std::make_shared<siconos::algebra::siconos::algebra::SimpleMatrix>(sizeDS, sizeDS));
  // TODO add this back to workV of the DS -> needed for X partial NS
  //  }
}

void siconos::modeling::FirstOrderType2R::computeh(
    double time, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SiconosVector& y)
{
  auto xp = x.toSiconosVector();
  ((Type2PtrH)(_pluginh->fPtr))(xp->size(), xp->data(), lambda.size(), const_cast<double*>(lambda.data()),
                                y.size(), y.data());
}

void siconos::modeling::FirstOrderType2R::computeg(
    double time, const siconos::algebra::SiconosVector& lambda,
    siconos::algebra::BlockVector& r)
{
  auto rp = r.toSiconosVector();
  ((Type2PtrG)(_pluging->fPtr))(lambda.size(), const_cast<double*>(lambda.data()), rp->size(), rp->data());
  r = *rp;
}

void siconos::modeling::FirstOrderType2R::computeOutput(double time, Interaction& inter,
                                                        unsigned int level)
{
  DEBUG_BEGIN("siconos::modeling::FirstOrderType2R::computeOutput \n");
  auto& DSlink = inter.linkToDSVariables();
  siconos::algebra::BlockVector& x = *DSlink[FirstOrderR::x];
  // copy into Siconos continuous memory vector
  siconos::algebra::SiconosVector& y = *inter.y(level);
  siconos::algebra::SiconosVector& lambda = *inter.lambda(level);
  computeh(time, x, lambda, y);
  DEBUG_EXPR(y.display());
  DEBUG_END("siconos::modeling::FirstOrderType2R::computeOutput \n");
}

void siconos::modeling::FirstOrderType2R::computeInput(double time, Interaction& inter,
                                                       unsigned int level)
{
  DEBUG_BEGIN("siconos::modeling::FirstOrderType2R::computeInput \n");
  auto& DSlink = inter.linkToDSVariables();
  // copy into Siconos continuous memory vector
  siconos::algebra::SiconosVector& lambda = *inter.lambda(level);
  computeg(time, lambda, *DSlink[FirstOrderR::r]);
  DEBUG_EXPR(DSlink[FirstOrderR::r]->display());
  DEBUG_END("siconos::modeling::FirstOrderType2R::computeInput \n");
}

void siconos::modeling::FirstOrderType2R::computeJachlambda(
    double time, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SimpleMatrix& D)
{
  THROW_EXCEPTION("siconos::modeling::FirstOrderType2R::computeJachlambda must be overload.");
}
void siconos::modeling::FirstOrderType2R::computeJachx(
    double time, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::SimpleMatrix& C)
{
  THROW_EXCEPTION("siconos::modeling::FirstOrderType2R::computeJachx must be overload.");
  // Note FP: so this class should be virtual, isn't it?
}

void siconos::modeling::FirstOrderType2R::computeJach(double time, Interaction& inter)
{
  DEBUG_BEGIN("siconos::modeling::FirstOrderType2R::computeJach\n");
  auto& DSlink = inter.linkToDSVariables();
  auto& relationMat = inter.relationMatrices();

  if (!_C) {
    computeJachx(time, *DSlink[FirstOrderR::x], *inter.lambda(0),
                 *relationMat[FirstOrderR::mat_C]);
  }
  if (!_D) {
    computeJachlambda(time, *DSlink[FirstOrderR::x], *inter.lambda(0),
                      *relationMat[FirstOrderR::mat_D]);
  }
  DEBUG_END("siconos::modeling::FirstOrderType2R::computeJach\n");
}

void siconos::modeling::FirstOrderType2R::computeJacglambda(
    double time, const siconos::algebra::SiconosVector& lambda,
    siconos::algebra::SimpleMatrix& B)
{
  THROW_EXCEPTION("siconos::modeling::FirstOrderType2R::computeJacglambda must be overload.");
}

void siconos::modeling::FirstOrderType2R::computeJacg(double time, Interaction& inter)
{
  DEBUG_BEGIN("siconos::modeling::FirstOrderType2R::computeJacg\n");
  if (!_B) {
    auto& relationMat = inter.relationMatrices();
    computeJacglambda(time, *inter.lambda(0), *relationMat[FirstOrderR::mat_B]);
  }
  DEBUG_END("siconos::modeling::FirstOrderType2R::computeJacg\n");
}
