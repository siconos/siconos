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
#include "FirstOrderType2R.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <debug.h>


FirstOrderType2R::FirstOrderType2R():
  FirstOrderR(RELATION::Type2R)
{}

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

  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  unsigned int sizeZ = DSlink[FirstOrderR::z]->size();
  VectorOfSMatrices& relationMat = inter.relationMatrices();

  _vec_r.reset(new SiconosVector(sizeDS));
  _vec_x.reset(new SiconosVector(sizeDS));
  _vec_z.reset(new SiconosVector(sizeZ));



  if (!_C)
    relationMat[FirstOrderR::mat_C].reset(new SimpleMatrix(sizeY, sizeDS));
  if (!_D)
    relationMat[FirstOrderR::mat_D].reset(new SimpleMatrix(sizeY, sizeY));
  if (!_F)
    relationMat[FirstOrderR::mat_F].reset(new SimpleMatrix(sizeY, sizeZ));
  if (!_B)
    relationMat[FirstOrderR::mat_B].reset(new SimpleMatrix(sizeDS, sizeY));
  if (!_K)
    relationMat[FirstOrderR::mat_K].reset(new SimpleMatrix(sizeDS, sizeDS));

//  if (!_jacgx)
//  {
//    relationMat[FirstOrderR::mat_K].reset(new SimpleMatrix(sizeDS, sizeDS));
    // TODO add this back to workV of the DS -> needed for X partial NS
//  }
}



void FirstOrderType2R::checkSize(Interaction& inter){}

void FirstOrderType2R::computeh(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& y)
{
  ((Type2PtrH)(_pluginh->fPtr))(x.size(), x.getArray(), lambda.size(), lambda.getArray(), y.size(), y.getArray());
}

void FirstOrderType2R::computeg(double time, SiconosVector& lambda, SiconosVector& r)
{
  ((Type2PtrG)(_pluging->fPtr))(lambda.size(), lambda.getArray(), r.size(), r.getArray());
}

void FirstOrderType2R::computeOutput(double time, Interaction& inter, unsigned int level)
{
  DEBUG_BEGIN("FirstOrderType2R::computeOutput \n");
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  BlockVector& x = *DSlink[FirstOrderR::x];
  BlockVector& z = *DSlink[FirstOrderR::z];
  // copy into Siconos continuous memory vector
  SP::SiconosVector x_vec(new SiconosVector(x));
  SP::SiconosVector z_vec(new SiconosVector(z));
  SiconosVector& y = *inter.y(level);
  SiconosVector& lambda = *inter.lambda(level);
  computeh(time, *x_vec, lambda, y);
  DEBUG_EXPR(y.display());
  DEBUG_END("FirstOrderType2R::computeOutput \n");
}

void FirstOrderType2R::computeInput(double time, Interaction& inter, unsigned int level)
{
  DEBUG_BEGIN("FirstOrderType2R::computeInput \n");
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  // copy into Siconos continuous memory vector
  SP::SiconosVector r_vec(new SiconosVector(*DSlink[FirstOrderR::r]));
  SiconosVector& lambda = *inter.lambda(level);
  computeg(time, lambda,  *r_vec);
  *DSlink[FirstOrderR::r] = *r_vec;
  DEBUG_EXPR(DSlink[FirstOrderR::r]->display());
  DEBUG_END("FirstOrderType2R::computeInput \n");
}

void FirstOrderType2R::computeJachlambda(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJachlambda must be overload.");
}
void FirstOrderType2R::computeJachx(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJachx must be overload.");
}

void FirstOrderType2R::computeJach(double time, Interaction& inter)
{
  DEBUG_BEGIN("FirstOrderType2R::computeJach\n");
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  VectorOfSMatrices& relationMat = inter.relationMatrices();

  if (!_C)
  {
    *_vec_x = *DSlink[FirstOrderR::x];
    computeJachx(time, *_vec_x, *inter.lambda(0), *relationMat[FirstOrderR::mat_C]);
  }
  if (!_D)
  {
    *_vec_x = *DSlink[FirstOrderR::x];
    computeJachlambda(time, *_vec_x, *inter.lambda(0), *relationMat[FirstOrderR::mat_D]);
  }
  DEBUG_END("FirstOrderType2R::computeJach\n");
}

void FirstOrderType2R::computeJacglambda(double time, SiconosVector& lambda, SimpleMatrix& B)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJacglambda must be overload.");
}

void FirstOrderType2R::computeJacg(double time, Interaction& inter)
{
  DEBUG_BEGIN("FirstOrderType2R::computeJacg\n");
  if (!_B)
  {
    VectorOfSMatrices& relationMat = inter.relationMatrices();
    computeJacglambda(time, *inter.lambda(0), *relationMat[FirstOrderR::mat_B]);
  }
  DEBUG_END("FirstOrderType2R::computeJacg\n");
}
