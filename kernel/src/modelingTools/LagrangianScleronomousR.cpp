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

#include "LagrangianScleronomousR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "PluggedObject.hpp"
#include "PluggedObject.hpp"          // getPluginFunctionname ...
#include "PluginTypes.hpp"            // FPtr2 ...
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
#include "SiconosException.hpp"
#include "siconos_debug.h"

// constructor from a set of data
siconos::modeling::LagrangianScleronomousR::LagrangianScleronomousR(
    const std::string& pluginh, const std::string& pluginJacobianhq)
    : LagrangianR(RelationSubType::ScleronomousR) {
  _zeroPlugin();
  setComputehFunction(siconos::plugins::getPluginName(pluginh),
                      siconos::plugins::getPluginFunctionName(pluginh));

  _pluginJachq->setComputeFunction(pluginJacobianhq);

  // Warning: we cannot allocate memory for Jach[0] matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
}
// constructor from a data used for EventDriven scheme
siconos::modeling::LagrangianScleronomousR::LagrangianScleronomousR(
    const std::string& pluginh, const std::string& pluginJacobianhq,
    const std::string& pluginDotJacobianhq)
    : LagrangianR(RelationSubType::ScleronomousR) {
  _zeroPlugin();
  setComputehFunction(siconos::plugins::getPluginName(pluginh),
                      siconos::plugins::getPluginFunctionName(pluginh));

  _pluginJachq->setComputeFunction(pluginJacobianhq);

  _plugindotjacqh->setComputeFunction(pluginDotJacobianhq);
}

void siconos::modeling::LagrangianScleronomousR::_zeroPlugin() {
  LagrangianR::_zeroPlugin();
  _pluginJachq = std::make_shared<siconos::plugins::PluggedObject>();
  _plugindotjacqh = std::make_shared<siconos::plugins::PluggedObject>();
}

void siconos::modeling::LagrangianScleronomousR::initialize(Interaction& inter) {
  if (!_jachq) {
    unsigned int sizeY = inter.dimension();
    unsigned int sizeDS = inter.getSizeOfDS();
    _jachq = std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeDS);
  }
}

void siconos::modeling::LagrangianScleronomousR::checkSize(Interaction& inter) {}

void siconos::modeling::LagrangianScleronomousR::computeh(
    const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z,
    siconos::algebra::SiconosVector& y) {
  DEBUG_PRINT(
      " siconos::modeling::LagrangianScleronomousR::computeh(Interaction& inter, "
      "siconos::algebra::BlockVector q, siconos::algebra::BlockVector z)\n");
  if (_pluginh && _pluginh->fPtr) {
    auto qp = q.toSiconosVector();
    auto zp = z.toSiconosVector();
    ((siconos::plugins::FPtr3)(_pluginh->fPtr))(qp->size(), &(*qp)(0), y.size(), &(y(0)),
                                                zp->size(), &(*zp)(0));
    z = *zp;
    DEBUG_EXPR(y.display());
  }
  // else nothing
}

void siconos::modeling::LagrangianScleronomousR::computeJachq(
    const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z) {
  if (_jachq && _pluginJachq->fPtr) {
    auto qp = q.toSiconosVector();
    auto zp = z.toSiconosVector();
    // get vector lambda of the current interaction
    ((siconos::plugins::FPtr3)(_pluginJachq->fPtr))(qp->size(), &(*qp)(0), _jachq->size(0),
                                                    &(*_jachq)(0, 0), zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::LagrangianScleronomousR::computeDotJachq(
    const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z,
    const siconos::algebra::BlockVector& qDot) {
  if (_dotjachq && _plugindotjacqh->fPtr) {
    auto qp = q.toSiconosVector();
    auto zp = z.toSiconosVector();
    auto qdotp = qDot.toSiconosVector();
    ((siconos::plugins::FPtr2)(_plugindotjacqh->fPtr))(qp->size(), &(*qp)(0), qdotp->size(),
                                                       &(*qdotp)(0), &(*_dotjachq)(0, 0),
                                                       zp->size(), &(*zp)(0));
    z = *zp;
  }
}

void siconos::modeling::LagrangianScleronomousR::computedotjacqhXqdot(double time,
                                                                      Interaction& inter) {
  auto& DSlink = inter.linkToDSVariables();
  DEBUG_PRINT("siconos::modeling::LagrangianScleronomousR::computeNonLinearH2dot starts");
  // Compute the H Jacobian dot
  computeDotJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z], *DSlink[LagrangianR::q1]);
  _dotjacqhXqdot = std::make_shared<siconos::algebra::SiconosVector>(_dotjachq->size(0));
  DEBUG_EXPR(_dotjachq->display(););
  siconos::algebra::prod(*_dotjachq, *DSlink[LagrangianR::q1], *_dotjacqhXqdot);
  DEBUG_PRINT("siconos::modeling::LagrangianScleronomousR::computeNonLinearH2dot ends");
}

void siconos::modeling::LagrangianScleronomousR::computeOutput(double time, Interaction& inter,
                                                               unsigned int derivativeNumber) {
  DEBUG_PRINTF(
      "siconos::modeling::LagrangianScleronomousR::computeOutput(double time, Interaction& "
      "inter, InteractionProperties& interProp, unsigned int derivativeNumber) with time = %f "
      "and derivativeNumber = %i\n",
      time, derivativeNumber);
  auto& DSlink = inter.linkToDSVariables();
  auto& y = *inter.y(derivativeNumber);
  if (derivativeNumber == 0) {
    computeh(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z], y);
  } else {
    computeJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);

    if (derivativeNumber == 1) {
      assert(_jachq);
      siconos::algebra::prod(*_jachq, *DSlink[LagrangianR::q1], y);
    } else if (derivativeNumber == 2) {
      assert(_jachq);
      siconos::algebra::prod(*_jachq, *DSlink[LagrangianR::q2], y);
      if (!_dotjachq) {
        auto sizeY = inter.dimension();
        auto sizeDS = inter.getSizeOfDS();
        _dotjachq = std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeDS);
      }
      computeDotJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z],
                      *DSlink[LagrangianR::q1]);
      siconos::algebra::prod(*_dotjachq, *DSlink[LagrangianR::q1], y, false);
    } else
      THROW_EXCEPTION(
          "siconos::modeling::LagrangianScleronomousR::computeOutput(double time, "
          "Interaction& inter, InteractionProperties& interProp, unsigned int "
          "derivativeNumber), index out of range");
  }
}

void siconos::modeling::LagrangianScleronomousR::computeInput(double time, Interaction& inter,
                                                              unsigned int level) {
  DEBUG_BEGIN(
      "void siconos::modeling::LagrangianScleronomousR::computeInput(double time, "
      "Interaction& inter, InteractionProperties& interProp, unsigned int level) \n");

  DEBUG_PRINTF("level = %i\n", level);
  auto& DSlink = inter.linkToDSVariables();
  computeJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  DEBUG_EXPR(lambda.display(););
  DEBUG_EXPR(_jachq->display(););
  // data[name] += trans(G) * lambda
  siconos::algebra::prod(lambda, *_jachq, *DSlink[LagrangianR::p0 + level], false);
  DEBUG_EXPR(DSlink[LagrangianR::p0 + level]->display(););
  DEBUG_END(
      "void siconos::modeling::LagrangianScleronomousR::computeInput(double time, "
      "Interaction& inter, InteractionProperties& interProp, unsigned int level) \n");
}

void siconos::modeling::LagrangianScleronomousR::computeJach(double time, Interaction& inter) {
  DEBUG_BEGIN(
      "void siconos::modeling::LagrangianScleronomousR::computeJach(double time, Interaction& "
      "inter) \n");
  auto& DSlink = inter.linkToDSVariables();
  DEBUG_EXPR(inter.display(););
  computeJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);
  // computeJachqDot(time, inter);
  if (!_dotjachq) {
    auto sizeY = inter.dimension();
    auto sizeDS = inter.getSizeOfDS();
    _dotjachq = std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, sizeDS);
  }
  computeDotJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z], *DSlink[LagrangianR::q1]);
  // computeJachlambda(time, inter);
  // computehDot(time,inter);
  DEBUG_END(
      "void siconos::modeling::LagrangianScleronomousR::computeJach(double time, Interaction& "
      "inter) \n");
}
