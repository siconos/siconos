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
#include "FirstOrderNonLinearR.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

#include "debug.h"

typedef void (*FONLR_h)(double, unsigned, double*, unsigned, double*, double*, unsigned, double*);
typedef FONLR_h FONLR_g;
typedef FONLR_h FONLR_C;
typedef FONLR_h FONLR_B;
typedef FONLR_h FONLR_K;
typedef FONLR_h FONLR_D;


void FirstOrderNonLinearR::initialize(Interaction& inter)
{
  FirstOrderR::initialize(inter);

  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  unsigned int sizeZ = DSlink[FirstOrderR::z]->size();

  _vec_r.reset(new SiconosVector(sizeDS));
  _vec_x.reset(new SiconosVector(sizeDS));
  _vec_z.reset(new SiconosVector(sizeZ));



  VectorOfSMatrices& relationMat = inter.relationMatrices();

  relationMat[FirstOrderR::mat_C].reset(new SimpleMatrix(sizeY, sizeDS));
  relationMat[FirstOrderR::mat_D].reset(new SimpleMatrix(sizeY, sizeY));
  relationMat[FirstOrderR::mat_B].reset(new SimpleMatrix(sizeDS, sizeY));
  relationMat[FirstOrderR::mat_K].reset(new SimpleMatrix(sizeDS, sizeDS));

  // F ? e ?
}

void FirstOrderNonLinearR::checkSize(Interaction& inter)
{}

void FirstOrderNonLinearR::computeh(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& y)
{
  if (_pluginh)
  {
    ((FONLR_h)_pluginh->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), y.getArray(), z.size(), z.getArray());
  }
  else
  {
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeh - no plugin detected, you should provide one or derive this class and implement this function");
  }
}

void FirstOrderNonLinearR::computeg(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& r)
{
  if (_pluging)
  {
    ((FONLR_g)_pluging->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), r.getArray(), z.size(), z.getArray());
  }
  else
  {
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeg - no plugin detected, you should provide one or derive this class and implement this function");
  }
}

void FirstOrderNonLinearR::computeJachx(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& C)
{
  if (_pluginJachx)
  {
    ((FONLR_C)_pluginJachx->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), C.getArray(), z.size(), z.getArray());
  }
  else
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeJachx, you need to derive this function in order to use it");
}

void FirstOrderNonLinearR::computeJachlambda(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& D)
{
  if (_pluginJachlambda)
  {
    ((FONLR_D)_pluginJachlambda->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), D.getArray(), z.size(), z.getArray());
  }
  else
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeJachlambda, you need to either provide a matrix D or derive this function in order to use it");
}

void FirstOrderNonLinearR::computeJacglambda(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& B)
{
  if (_pluginJacglambda)
  {
    ((FONLR_B)_pluginJacglambda->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), B.getArray(), z.size(), z.getArray());
  }
  else
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeJacglambda, you need to either provide a matrix B or derive this function in order to use it");
}

void FirstOrderNonLinearR::computeJacgx(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& K)
{
  if (_pluginJacgx)
  {
    ((FONLR_K)_pluginJacgx->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), K.getArray(), z.size(), z.getArray());
  }
  else
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeJacgx, you need to either provide a matrix K or derive this function in order to use it");
}


void FirstOrderNonLinearR::computeOutput(double time, Interaction& inter, unsigned int level)
{
  DEBUG_PRINT("FirstOrderNonLinearR::computeOutput \n");
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  BlockVector& x = *DSlink[FirstOrderR::x];
  BlockVector& z = *DSlink[FirstOrderR::z];
  // copy into Siconos continuous memory vector
  *_vec_x = x ;
  *_vec_z = z ;
  SiconosVector& y = *inter.y(level);
  SiconosVector& lambda = *inter.lambda(level);
  computeh(time, *_vec_x, lambda, *_vec_z, y);
  DEBUG_END("FirstOrderNonLinearR::computeOutput \n");
}

void FirstOrderNonLinearR::computeInput(double time, Interaction& inter, unsigned int level)
{
  DEBUG_PRINT("FirstOrderNonLinearR::computeInput \n");
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  BlockVector& x = *DSlink[FirstOrderR::x];
  BlockVector& z = *DSlink[FirstOrderR::z];
  // copy into Siconos continuous memory vector
  *_vec_x = x ;
  *_vec_z = z ;
  SiconosVector& lambda = *inter.lambda(level);
  computeg(time, *_vec_x, lambda, *_vec_z, *_vec_r);
  *DSlink[FirstOrderR::r] = *_vec_r;
  DEBUG_END("FirstOrderNonLinearR::computeinput \n");
}

void FirstOrderNonLinearR::computeJach(double time, Interaction& inter)
{
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  VectorOfSMatrices& relationMat = inter.relationMatrices();

  *_vec_x = *DSlink[FirstOrderR::x];
  *_vec_z = *DSlink[FirstOrderR::z];
  SiconosVector& lambda = *inter.lambda(0);

  if (!_C)
  {
    computeJachx(time, *_vec_x, lambda, *_vec_z, *relationMat[FirstOrderR::mat_C]);
  }

  if (!_D)
  {
    computeJachlambda(time, *_vec_x, lambda, *_vec_z, *relationMat[FirstOrderR::mat_D]);
  }
  *DSlink[FirstOrderR::z] = *_vec_z;
}

void FirstOrderNonLinearR::computeJacg(double time, Interaction& inter)
{
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  VectorOfSMatrices& relationMat = inter.relationMatrices();

  *_vec_x = *DSlink[FirstOrderR::x];
  *_vec_z = *DSlink[FirstOrderR::z];
  SiconosVector& lambda = *inter.lambda(0);
  if (!_B)
  {
    computeJacglambda(time, *_vec_x, lambda, *_vec_z, *relationMat[FirstOrderR::mat_B]);
  }
  if (!_K)
  {
    computeJacgx(time, *_vec_x, lambda, *_vec_z, *relationMat[FirstOrderR::mat_K]);
  }
  *DSlink[FirstOrderR::z] = *_vec_z;
}
