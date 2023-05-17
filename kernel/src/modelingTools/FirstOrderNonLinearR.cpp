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

#include "FirstOrderNonLinearR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "PluggedObject.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
// #include "FirstOrderNonLinearDS.hpp"

// #include "siconos::algebra::BlockVector.hpp"
// #include "SimulationGraphs.hpp"

#include "siconos_debug.h"

typedef void (*FONLR_h)(double time, unsigned x_size, double* x, unsigned size_lambda,
                        double* lambda, double*, unsigned z_size, double* z);
using FONLR_g = FONLR_h;
using FONLR_C = FONLR_h;
using FONLR_B = FONLR_h;
using FONLR_K = FONLR_h;
using FONLR_D = FONLR_h;

void siconos::modeling::FirstOrderNonLinearR::initialize(Interaction& inter)
{
  FirstOrderR::initialize(inter);

  auto sizeY = inter.dimension();
  auto sizeDS = inter.getSizeOfDS();
  auto& relationMat = inter.relationMatrices();

  relationMat[FirstOrderR::mat_C] =
      std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeDS);
  relationMat[FirstOrderR::mat_D] =
      std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeY);
  relationMat[FirstOrderR::mat_B] =
      std::make_shared<siconos::algebra::SimpleMatrix>(sizeDS, sizeY);
  relationMat[FirstOrderR::mat_K] =
      std::make_shared<siconos::algebra::SimpleMatrix>(sizeDS, sizeDS);

  // F ?
}

void siconos::modeling::FirstOrderNonLinearR::checkSize(Interaction& inter) {}

void siconos::modeling::FirstOrderNonLinearR::computeh(
    double time, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::BlockVector& z,
    siconos::algebra::SiconosVector& y)
{
  if (_pluginh) {
    auto xp = x.toSiconosVector();
    auto zp = z.toSiconosVector();
    ((FONLR_h)_pluginh->fPtr)(time, xp->size(), xp->data(), lambda.size(),
                              const_cast<double*>(lambda.data()), y.data(), zp->size(), zp->data());
    z = *zp;
  }
  else {
    THROW_EXCEPTION(
        "siconos::modeling::FirstOrderNonLinearR::computeh - no plugin detected, you should "
        "provide one or derive this class and implement this function");
  }
}

void siconos::modeling::FirstOrderNonLinearR::computeg(
    double time, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::BlockVector& z,
    siconos::algebra::BlockVector& r)
{
  if (_pluging) {
    auto xp = x.toSiconosVector();
    auto zp = z.toSiconosVector();
    auto rp = r.toSiconosVector();
    ((FONLR_g)_pluging->fPtr)(time, xp->size(), xp->data(), lambda.size(),
                              const_cast<double*>(lambda.data()), rp->data(), zp->size(), zp->data());
    z = *zp;
    r = *rp;
  }
  else {
    THROW_EXCEPTION(
        "siconos::modeling::FirstOrderNonLinearR::computeg - no plugin detected, you should "
        "provide one or derive this class and implement this function");
  }
}

void siconos::modeling::FirstOrderNonLinearR::computeJachx(
    double time, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::BlockVector& z,
    siconos::algebra::SimpleMatrix& C)
{
  if (_pluginJachx) {
    auto xp = x.toSiconosVector();
    auto zp = z.toSiconosVector();
    ((FONLR_C)_pluginJachx->fPtr)(time, xp->size(), xp->data(), lambda.size(),
                                  const_cast<double*>(lambda.data()), C.data(), zp->size(), zp->data());
    z = *zp;
  }
  else
    THROW_EXCEPTION(
        "siconos::modeling::FirstOrderNonLinearR::computeJachx, you need to derive this "
        "function in order to use it");
}

void siconos::modeling::FirstOrderNonLinearR::computeJachlambda(
    double time, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::BlockVector& z,
    siconos::algebra::SimpleMatrix& D)
{
  if (_pluginJachlambda) {
    auto xp = x.toSiconosVector();
    auto zp = z.toSiconosVector();
    ((FONLR_D)_pluginJachlambda->fPtr)(time, xp->size(), xp->data(), lambda.size(),
                                       const_cast<double*>(lambda.data()), D.data(), zp->size(),
                                       zp->data());
    z = *zp;
  }
  else
    THROW_EXCEPTION(
        "siconos::modeling::FirstOrderNonLinearR::computeJachlambda, you need to either "
        "provide a matrix D or derive this function in order to use it");
}

void siconos::modeling::FirstOrderNonLinearR::computeJacglambda(
    double time, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::BlockVector& z,
    siconos::algebra::SimpleMatrix& B)
{
  if (_pluginJacglambda) {
    auto xp = x.toSiconosVector();
    auto zp = z.toSiconosVector();
    ((FONLR_B)_pluginJacglambda->fPtr)(time, xp->size(), xp->data(), lambda.size(),
                                       const_cast<double*>(lambda.data()), B.data(), zp->size(),
                                       zp->data());
    z = *zp;
  }
  else
    THROW_EXCEPTION(
        "siconos::modeling::FirstOrderNonLinearR::computeJacglambda, you need to either "
        "provide a matrix B or derive this function in order to use it");
}

void siconos::modeling::FirstOrderNonLinearR::computeJacgx(
    double time, const siconos::algebra::BlockVector& x,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::BlockVector& z,
    siconos::algebra::SimpleMatrix& K)
{
  if (_pluginJacgx) {
    auto xp = x.toSiconosVector();
    auto zp = z.toSiconosVector();
    ((FONLR_K)_pluginJacgx->fPtr)(time, xp->size(), xp->data(), lambda.size(),
                                  const_cast<double*>(lambda.data()), K.data(), zp->size(), zp->data());
    z = *zp;
  }
  else
    THROW_EXCEPTION(
        "siconos::modeling::FirstOrderNonLinearR::computeJacgx, you need to either provide a "
        "matrix K or derive this function in order to use it");
}

void siconos::modeling::FirstOrderNonLinearR::computeOutput(double time, Interaction& inter,
                                                            unsigned int level)
{
  DEBUG_PRINT("siconos::modeling::FirstOrderNonLinearR::computeOutput \n");
  auto& DSlink = inter.linkToDSVariables();
  auto& y = *inter.y(level);
  auto& lambda = *inter.lambda(level);
  computeh(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z], y);
  DEBUG_END("siconos::modeling::FirstOrderNonLinearR::computeOutput \n");
}

void siconos::modeling::FirstOrderNonLinearR::computeInput(double time, Interaction& inter,
                                                           unsigned int level)
{
  DEBUG_PRINT("siconos::modeling::FirstOrderNonLinearR::computeInput \n");
  auto& DSlink = inter.linkToDSVariables();
  auto& lambda = *inter.lambda(level);
  computeg(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z],
           *DSlink[FirstOrderR::r]);
  DEBUG_END("siconos::modeling::FirstOrderNonLinearR::computeinput \n");
}

void siconos::modeling::FirstOrderNonLinearR::computeJach(double time, Interaction& inter)
{
  auto& DSlink = inter.linkToDSVariables();
  auto& relationMat = inter.relationMatrices();

  auto& lambda = *inter.lambda(0);

  if (!_C) {
    computeJachx(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z],
                 *relationMat[FirstOrderR::mat_C]);
  }

  if (!_D) {
    computeJachlambda(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z],
                      *relationMat[FirstOrderR::mat_D]);
  }
}

void siconos::modeling::FirstOrderNonLinearR::computeJacg(double time, Interaction& inter)
{
  auto& DSlink = inter.linkToDSVariables();
  auto& relationMat = inter.relationMatrices();

  auto& lambda = *inter.lambda(0);
  if (!_B) {
    computeJacglambda(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z],
                      *relationMat[FirstOrderR::mat_B]);
  }
  if (!_K) {
    computeJacgx(time, *DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z],
                 *relationMat[FirstOrderR::mat_K]);
  }
}
