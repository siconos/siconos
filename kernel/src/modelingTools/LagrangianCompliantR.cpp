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

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "LagrangianCompliantR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "PluggedObject.hpp"
#include "PluginTypes.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

siconos::modeling::LagrangianCompliantR::LagrangianCompliantR(
    const std::string& pluginh, const std::string& pluginJacobianhq,
    const std::string& pluginJacobianhlambda)
    : LagrangianR(RelationSubType::CompliantR) {
  _zeroPlugin();
  setComputehFunction(siconos::plugins::getPluginName(pluginh),
                      siconos::plugins::getPluginFunctionName(pluginh));
  _pluginJachq->setComputeFunction(siconos::plugins::getPluginName(pluginJacobianhq),
                                   siconos::plugins::getPluginFunctionName(pluginJacobianhq));
  _pluginJachlambda->setComputeFunction(
      siconos::plugins::getPluginName(pluginJacobianhlambda),
      siconos::plugins::getPluginFunctionName(pluginJacobianhlambda));
}

void siconos::modeling::LagrangianCompliantR::_zeroPlugin() {
  _pluginJachq = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginJachlambda = std::make_shared<siconos::plugins::PluggedObject>();
}

void siconos::modeling::LagrangianCompliantR::initialize(Interaction& inter) {
  auto sizeY = inter.dimension();

  if (!_jachlambda)
    _jachlambda = std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeY);
  else
    _jachlambda->resize(sizeY, sizeY);
}
void siconos::modeling::LagrangianCompliantR::checkSize(Interaction& inter) {}
void siconos::modeling::LagrangianCompliantR::computeh(
    double time, const siconos::algebra::BlockVector& q0,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::BlockVector& z,
    siconos::algebra::SiconosVector& y) {
  if (_pluginh->fPtr) {
    auto qp = q0.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((siconos::plugins::FPtr2)(_pluginh->fPtr))(
        qp->size(), &(*qp)(0), y.size(), lambda.getArray(), &(y)(0), zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::LagrangianCompliantR::computeJachq(
    double time, const siconos::algebra::BlockVector& q0,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::BlockVector& z) {
  if (_jachq && _pluginJachq->fPtr) {
    auto qp = q0.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((siconos::plugins::FPtr2)(_pluginJachq->fPtr))(qp->size(), &(*qp)(0), lambda.size(),
                                                    lambda.getArray(), &(*_jachq)(0, 0),
                                                    zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::LagrangianCompliantR::computeJachlambda(
    double time, const siconos::algebra::BlockVector& q0,
    const siconos::algebra::SiconosVector& lambda, siconos::algebra::BlockVector& z) {
  if (_jachlambda && _pluginJachlambda->fPtr) {
    auto qp = q0.prepareVectorForPlugin();
    auto zp = z.prepareVectorForPlugin();
    ((siconos::plugins::FPtr2)_pluginJachlambda->fPtr)(
        qp->size(), &(*qp)(0), lambda.size(), lambda.getArray(), &(*_jachlambda)(0, 0),
        zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::LagrangianCompliantR::computeOutput(double time, Interaction& inter,
                                                            unsigned int derivativeNumber) {
  auto& DSlink = inter.linkToDSVariables();
  if (derivativeNumber == 0) {
    auto& y = *inter.y(0);
    auto& lambda = *inter.lambda(0);
    computeh(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z], y);
  } else {
    auto& y = *inter.y(derivativeNumber);
    auto& lambda = *inter.lambda(derivativeNumber);
    computeJachq(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z]);
    computeJachlambda(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z]);
    if (derivativeNumber == 1) {
      // y = Jach[0] q1 + Jach[1] lambda
      siconos::algebra::prod(*_jachq, *DSlink[LagrangianR::q1], y);
      siconos::algebra::prod(*_jachlambda, lambda, y, false);
    } else if (derivativeNumber == 2)
      siconos::algebra::prod(*_jachq, *DSlink[LagrangianR::q2],
                             y);  // Approx: y[2] = Jach[0]q[2], other terms are neglected ...
    else
      THROW_EXCEPTION(
          "siconos::modeling::LagrangianCompliantR::computeOutput, index out of range or not "
          "yet implemented.");
  }
}

void siconos::modeling::LagrangianCompliantR::computeInput(double time, Interaction& inter,
                                                           unsigned int level) {
  // get lambda of the concerned interaction

  auto& lambda = *inter.lambda(level);
  auto& DSlink = inter.linkToDSVariables();
  computeJachq(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z]);
  // data[name] += trans(G) * lambda
  siconos::algebra::prod(lambda, *_jachq, *DSlink[LagrangianR::p0 + level], false);
}

void siconos::modeling::LagrangianCompliantR::computeJach(double time, Interaction& inter) {
  auto& DSlink = inter.linkToDSVariables();
  auto& lambda = *inter.lambda(0);
  computeJachq(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z]);
  computeJachlambda(time, *DSlink[LagrangianR::q0], lambda, *DSlink[LagrangianR::z]);
}
