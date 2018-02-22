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
#include "FirstOrderType1R.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"


FirstOrderType1R::FirstOrderType1R(const std::string& pluginh, const std::string& pluging):
  FirstOrderR(RELATION::Type1R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSLH::getPluginName(pluginh), SSLH::getPluginFunctionName(pluginh));
  setComputegFunction(SSLH::getPluginName(pluging), SSLH::getPluginFunctionName(pluging));
  // The jacobians are not set, and thus considered as null matrices at this point.
}

FirstOrderType1R::FirstOrderType1R(const std::string& pluginh, const std::string& pluging, const std::string& pluginJachx, const std::string& pluginJacglambda):
  FirstOrderR(RELATION::Type1R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSLH::getPluginName(pluginh), SSLH::getPluginFunctionName(pluginh));
  setComputegFunction(SSLH::getPluginName(pluging), SSLH::getPluginFunctionName(pluging));
  setComputeJachxFunction(SSLH::getPluginName(pluginJachx), SSLH::getPluginFunctionName(pluginJachx));
  setComputeJacglambdaFunction(SSLH::getPluginName(pluginJacglambda), SSLH::getPluginFunctionName(pluginJacglambda));
}


void FirstOrderType1R::initialize(Interaction& inter)
{

  FirstOrderR::initialize(inter);

  // Check if an Interaction is connected to the Relation.
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  unsigned int sizeZ = DSlink[FirstOrderR::z]->size();

  VectorOfSMatrices& relationMat = inter.relationMatrices();
  if (!_C)
    relationMat[FirstOrderR::mat_C].reset(new SimpleMatrix(sizeY, sizeDS));
  if (!_D)
    relationMat[FirstOrderR::mat_D].reset(new SimpleMatrix(sizeY, sizeY));
  if (!_F)
    relationMat[FirstOrderR::mat_F].reset(new SimpleMatrix(sizeY, sizeZ));
  if (!_B)
    relationMat[FirstOrderR::mat_B].reset(new SimpleMatrix(sizeDS, sizeY));
}


void FirstOrderType1R::checkSize(Interaction& inter)
{

}
void FirstOrderType1R::computeh(double time, SiconosVector& x, SiconosVector& z, SiconosVector& y)
{
  assert(_pluginh && "FirstOrderType1R::computeOutput() is not linked to a plugin function");
  
  ((Type1Ptr)(_pluginh->fPtr))(x.size(), &(x)(0), y.size(), &(y)(0), z.size(), &(z)(0));

}

void FirstOrderType1R::computeg(double time, SiconosVector& lambda, SiconosVector& z, SiconosVector& r)
{
  assert(_pluging && "FirstOrderType1R::computeInput() is not linked to a plugin function");

  ((Type1Ptr)(_pluging->fPtr))(lambda.size(), &(lambda)(0), r.size(), &(r)(0), z.size(), &(z)(0));

}
void FirstOrderType1R::computeOutput(double time, Interaction& inter, unsigned int level)
{
  SiconosVector& y = *inter.y(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  BlockVector& x = *DSlink[FirstOrderR::x];
  BlockVector& z = *DSlink[FirstOrderR::z];
  // copy into Siconos continuous memory vector
  SP::SiconosVector x_vec(new SiconosVector(x));
  SP::SiconosVector z_vec(new SiconosVector(z));
  
  computeh(time, *x_vec, *z_vec, y);

  *DSlink[FirstOrderR::z] = *z_vec;
}

void FirstOrderType1R::computeInput(double time, Interaction& inter, unsigned int level)
{
  assert(_pluging && "FirstOrderType1R::computeInput() is not linked to a plugin function");

  SiconosVector& lambda = *inter.lambda(level);

  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  BlockVector& r = *DSlink[FirstOrderR::r];
  BlockVector& z = *DSlink[FirstOrderR::z];
  // copy into Siconos continuous memory vector
  SP::SiconosVector r_vec(new SiconosVector(r));
  SP::SiconosVector z_vec(new SiconosVector(z));

  computeg(time, lambda, *z_vec, *r_vec);

  *DSlink[FirstOrderR::r] = *r_vec;
  *DSlink[FirstOrderR::z] = *z_vec;
}

void FirstOrderType1R::computeJachx(double time, SiconosVector& x, SiconosVector& z, SimpleMatrix& C)
{
  //
  assert(_pluginJachx && "FirstOrderType1R::computeJacobianH() failed; not linked to a plug-in function.");
  if(_C && _pluginJachx)
    ((Type1Ptr)(_pluginJachx->fPtr))(x.size(), &(x)(0), C.size(0), C.getArray(), z.size(), &(z)(0));

}

void FirstOrderType1R::computeJachz(double time, SiconosVector& x, SiconosVector& z, SimpleMatrix& F)
{
  if (_F && _pluginJachz && _pluginJachz->fPtr)
    ((Type1Ptr)(_pluginJachz->fPtr))(x.size(), &(x)(0), F.size(0), F.getArray(), z.size(), &(z)(0));

}

void FirstOrderType1R::computeJacglambda(double time, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& B)
{
  assert(_pluginJacglambda && "FirstOrderType1R::computeJacobiang() failed; not linked to a plug-in function.");
  if(_B && _pluginJacglambda)
    ((Type1Ptr)(_pluginJacglambda->fPtr))(lambda.size(), &(lambda)(0), B.size(0), B.getArray(), z.size(), &(z)(0));
}

void FirstOrderType1R::computeJach(double time, Interaction& inter)
{
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  VectorOfSMatrices& relationMat = inter.relationMatrices();

  *_vec_x = *DSlink[FirstOrderR::x];
  *_vec_z = *DSlink[FirstOrderR::z];
  if (!_C)
  {
    computeJachx(time, *_vec_x, *_vec_z, *relationMat[FirstOrderR::mat_C]);
  }
  if (!_F)
  {
    computeJachz(time, *_vec_x, *_vec_z, *relationMat[FirstOrderR::mat_F]);
  }
  *DSlink[FirstOrderR::z] = *_vec_z;
}

void FirstOrderType1R::computeJacg(double time, Interaction& inter)
{
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  VectorOfSMatrices& relationMat = inter.relationMatrices();

  *_vec_z = *DSlink[FirstOrderR::z];
  if (!_B)
  {
    computeJacglambda(time, *inter.lambda(0), *_vec_z, *relationMat[FirstOrderR::mat_B]);
  }
  *DSlink[FirstOrderR::z] = *_vec_z;
}
