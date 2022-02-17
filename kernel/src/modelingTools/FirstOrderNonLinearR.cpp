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
#include "FirstOrderNonLinearR.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

#include "siconos_debug.h"

typedef void (*FONLR_h)(double time, unsigned x_size, double *x, unsigned size_lambda, double* lambda, double *, unsigned z_size, double *z);
typedef FONLR_h FONLR_g;
typedef FONLR_h FONLR_C;
typedef FONLR_h FONLR_B;
typedef FONLR_h FONLR_K;
typedef FONLR_h FONLR_D;


void FirstOrderNonLinearR::initialize(Interaction& inter)
{
  FirstOrderR::initialize(inter);

  unsigned int sizeY = inter.dimension();
  unsigned int sizeDS = inter.getSizeOfDS();
  VectorOfSMatrices& relationMat = inter.relationMatrices();

  relationMat[FirstOrderR::mat_C] = std::make_shared<SimpleMatrix>(sizeY, sizeDS);
  relationMat[FirstOrderR::mat_D] = std::make_shared<SimpleMatrix>(sizeY, sizeY);
  relationMat[FirstOrderR::mat_B] = std::make_shared<SimpleMatrix>(sizeDS, sizeY);
  relationMat[FirstOrderR::mat_K] = std::make_shared<SimpleMatrix>(sizeDS, sizeDS);

  // F ?
}

void FirstOrderNonLinearR::checkSize(Interaction& inter)
{

}

void FirstOrderNonLinearR::computeh(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, SiconosVector& y)
{
  if(_pluginh)
  {
    auto xp = x.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FONLR_h)_pluginh->fPtr)(time, xp->size(), xp->getArray(), lambda.size(), lambda.getArray(), y.getArray(), zp->size(), zp->getArray());
    z = *zp;
  }
  else
  {
    THROW_EXCEPTION("FirstOrderNonLinearR::computeh - no plugin detected, you should provide one or derive this class and implement this function");
  }
}

void FirstOrderNonLinearR::computeg(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, BlockVector& r)
{
  if(_pluging)
  {
    auto xp = x.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    auto rp = r.prepareVectorForPlugin();
    ((FONLR_g)_pluging->fPtr)(time, xp->size(), xp->getArray(), lambda.size(), lambda.getArray(), rp->getArray(), zp->size(), zp->getArray());
    z = *zp;
    r = *rp;
  }
  else
  {
    THROW_EXCEPTION("FirstOrderNonLinearR::computeg - no plugin detected, you should provide one or derive this class and implement this function");
  }
}

void FirstOrderNonLinearR::computeJachx(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, SimpleMatrix& C)
{
  if(_pluginJachx)
  {
    auto xp = x.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FONLR_C)_pluginJachx->fPtr)(time, xp->size(), xp->getArray(), lambda.size(), lambda.getArray(), C.getArray(), zp->size(), zp->getArray());
    z = *zp;
  }
  else
    THROW_EXCEPTION("FirstOrderNonLinearR::computeJachx, you need to derive this function in order to use it");
}

void FirstOrderNonLinearR::computeJachlambda(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, SimpleMatrix& D)
{
  if(_pluginJachlambda)
  {
    auto xp = x.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FONLR_D)_pluginJachlambda->fPtr)(time, xp->size(), xp->getArray(), lambda.size(), lambda.getArray(), D.getArray(), zp->size(), zp->getArray());
    z = *zp;
  }
  else
    THROW_EXCEPTION("FirstOrderNonLinearR::computeJachlambda, you need to either provide a matrix D or derive this function in order to use it");
}

void FirstOrderNonLinearR::computeJacglambda(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, SimpleMatrix& B)
{
  if(_pluginJacglambda)
  {
    auto xp = x.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FONLR_B)_pluginJacglambda->fPtr)(time, xp->size(), xp->getArray(), lambda.size(), lambda.getArray(), B.getArray(), zp->size(), zp->getArray());
    z = *zp;
  }
  else
    THROW_EXCEPTION("FirstOrderNonLinearR::computeJacglambda, you need to either provide a matrix B or derive this function in order to use it");
}

void FirstOrderNonLinearR::computeJacgx(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, SimpleMatrix& K)
{
  if(_pluginJacgx)
  {
    auto xp = x.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((FONLR_K)_pluginJacgx->fPtr)(time, xp->size(), xp->getArray(), lambda.size(), lambda.getArray(), K.getArray(), zp->size(), zp->getArray());
    z = *zp;
  }
  else
    THROW_EXCEPTION("FirstOrderNonLinearR::computeJacgx, you need to either provide a matrix K or derive this function in order to use it");
}


void FirstOrderNonLinearR::computeOutput(double time, Interaction& inter, unsigned int level)
{
  DEBUG_PRINT("FirstOrderNonLinearR::computeOutput \n");
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  SiconosVector& y = *inter.y(level);
  SiconosVector& lambda = *inter.lambda(level);
  computeh(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z], y);
  DEBUG_END("FirstOrderNonLinearR::computeOutput \n");
}

void FirstOrderNonLinearR::computeInput(double time, Interaction& inter, unsigned int level)
{
  DEBUG_PRINT("FirstOrderNonLinearR::computeInput \n");
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  SiconosVector& lambda = *inter.lambda(level);
  computeg(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z], *DSlink[FirstOrderR::r]);
  DEBUG_END("FirstOrderNonLinearR::computeinput \n");
}

void FirstOrderNonLinearR::computeJach(double time, Interaction& inter)
{
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  VectorOfSMatrices& relationMat = inter.relationMatrices();

  SiconosVector& lambda = *inter.lambda(0);

  if(!_C)
  {
    computeJachx(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z], *relationMat[FirstOrderR::mat_C]);
  }

  if(!_D)
  {
    computeJachlambda(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z], *relationMat[FirstOrderR::mat_D]);
  }
}

void FirstOrderNonLinearR::computeJacg(double time, Interaction& inter)
{
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  VectorOfSMatrices& relationMat = inter.relationMatrices();

  SiconosVector& lambda = *inter.lambda(0);
  if(!_B)
  {
    computeJacglambda(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z], *relationMat[FirstOrderR::mat_B]);
  }
  if(!_K)
  {
    computeJacgx(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z], *relationMat[FirstOrderR::mat_K]);
  }
}
