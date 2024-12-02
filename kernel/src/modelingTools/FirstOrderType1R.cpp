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
#include "FirstOrderType1R.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "PluggedObject.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

siconos::modeling::FirstOrderType1R::FirstOrderType1R(const std::string& pluginh,
                                                      const std::string& pluging)
    : FirstOrderR(RelationSubType::Type1R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(siconos::plugins::getPluginName(pluginh),
                      siconos::plugins::getPluginFunctionName(pluginh));
  setComputegFunction(siconos::plugins::getPluginName(pluging),
                      siconos::plugins::getPluginFunctionName(pluging));
  // The jacobians are not set, and thus considered as null matrices at this point.
}

siconos::modeling::FirstOrderType1R::FirstOrderType1R(const std::string& pluginh,
                                                      const std::string& pluging,
                                                      const std::string& pluginJachx,
                                                      const std::string& pluginJacglambda)
    : FirstOrderR(RelationSubType::Type1R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(siconos::plugins::getPluginName(pluginh),
                      siconos::plugins::getPluginFunctionName(pluginh));
  setComputegFunction(siconos::plugins::getPluginName(pluging),
                      siconos::plugins::getPluginFunctionName(pluging));
  setComputeJachxFunction(siconos::plugins::getPluginName(pluginJachx),
                          siconos::plugins::getPluginFunctionName(pluginJachx));
  setComputeJacglambdaFunction(siconos::plugins::getPluginName(pluginJacglambda),
                               siconos::plugins::getPluginFunctionName(pluginJacglambda));
}

void siconos::modeling::FirstOrderType1R::initialize(Interaction& inter)
{
  FirstOrderR::initialize(inter);

  // Check if an Interaction is connected to the Relation.
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
}

void siconos::modeling::FirstOrderType1R::computeh(double time,
                                                   const siconos::algebra::BlockVector& x,
                                                   siconos::algebra::BlockVector& z,
                                                   siconos::algebra::SiconosVector& y)
{
  assert(_pluginh &&
         "siconos::modeling::FirstOrderType1R::computeOutput() is not linked to a plugin "
         "function");
  auto xp = x.prepareVectorForPlugin();
  auto zp = z.prepareVectorForPlugin();
  ((Type1Ptr)(_pluginh->fPtr))(xp->size(), &(*xp)(0), y.size(), &(y)(0), zp->size(),
                               &(*zp)(0));
  z = *zp;
}

void siconos::modeling::FirstOrderType1R::computeg(
    double time, const siconos::algebra::SiconosVector& lambda,
    siconos::algebra::BlockVector& z, siconos::algebra::BlockVector& r)
{
  assert(_pluging &&
         "siconos::modeling::FirstOrderType1R::computeInput() is not linked to a plugin "
         "function");

  auto zp = z.prepareVectorForPlugin();
  auto rp = r.prepareVectorForPlugin();
  ((Type1Ptr)(_pluging->fPtr))(lambda.size(), lambda.getArray(), rp->size(), &(*rp)(0),
                               zp->size(), &(*zp)(0));
  z = *zp;
  r = *rp;
}
void siconos::modeling::FirstOrderType1R::computeOutput(double time, Interaction& inter,
                                                        unsigned int level)
{
  siconos::algebra::SiconosVector& y = *inter.y(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  auto& DSlink = inter.linkToDSVariables();
  // copy into Siconos continuous memory vector
  computeh(time, *DSlink[FirstOrderR::x], *DSlink[FirstOrderR::z], y);
}

void siconos::modeling::FirstOrderType1R::computeInput(double time, Interaction& inter,
                                                       unsigned int level)
{
  assert(_pluging &&
         "siconos::modeling::FirstOrderType1R::computeInput() is not linked to a plugin "
         "function");

  siconos::algebra::SiconosVector& lambda = *inter.lambda(level);

  auto& DSlink = inter.linkToDSVariables();
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  siconos::algebra::BlockVector& r = *DSlink[FirstOrderR::r];
  siconos::algebra::BlockVector& z = *DSlink[FirstOrderR::z];
  // copy into Siconos continuous memory vector
  computeg(time, lambda, z, r);
}

void siconos::modeling::FirstOrderType1R::computeJachx(double time,
                                                       const siconos::algebra::BlockVector& x,
                                                       siconos::algebra::BlockVector& z,
                                                       siconos::algebra::SimpleMatrix& C)
{
  //
  assert(_pluginJachx &&
         "siconos::modeling::FirstOrderType1R::computeJacobianH() failed; not linked to a "
         "plug-in function.");
  if (_pluginJachx) {
    auto xp = x.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((Type1Ptr)(_pluginJachx->fPtr))(xp->size(), &(*xp)(0), C.size(0), C.getArray(),
                                     zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::FirstOrderType1R::computeJachz(double time,
                                                       const siconos::algebra::BlockVector& x,
                                                       siconos::algebra::BlockVector& z,
                                                       siconos::algebra::SimpleMatrix& F)
{
  if (_pluginJachz && _pluginJachz->fPtr) {
    auto xp = x.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((Type1Ptr)(_pluginJachz->fPtr))(xp->size(), &(*xp)(0), F.size(0), F.getArray(),
                                     zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::FirstOrderType1R::computeJacglambda(
    double time, const siconos::algebra::SiconosVector& lambda,
    siconos::algebra::BlockVector& z, siconos::algebra::SimpleMatrix& B)
{
  assert(_pluginJacglambda &&
         "siconos::modeling::FirstOrderType1R::computeJacobiang() failed; not linked to a "
         "plug-in function.");
  if (_pluginJacglambda) {
    auto zp = z.prepareVectorForPlugin();
    ((Type1Ptr)(_pluginJacglambda->fPtr))(lambda.size(), lambda.getArray(), B.size(0),
                                          B.getArray(), zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::FirstOrderType1R::computeJach(double time, Interaction& inter)
{
  auto& DSlink = inter.linkToDSVariables();
  auto& relationMat = inter.relationMatrices();

  if (!_C) {
    computeJachx(time, *DSlink[FirstOrderR::x], *DSlink[FirstOrderR::z],
                 *relationMat[FirstOrderR::mat_C]);
  }
  if (!_F) {
    computeJachz(time, *DSlink[FirstOrderR::x], *DSlink[FirstOrderR::z],
                 *relationMat[FirstOrderR::mat_F]);
  }
}

void siconos::modeling::FirstOrderType1R::computeJacg(double time, Interaction& inter)
{
  auto& DSlink = inter.linkToDSVariables();
  auto& relationMat = inter.relationMatrices();
  if (!_B) {
    computeJacglambda(time, *inter.lambda(0), *DSlink[FirstOrderR::z],
                      *relationMat[FirstOrderR::mat_B]);
  }
}
