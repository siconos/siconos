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
#include "LagrangianLinearTIR.hpp"

#include <iostream>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

// Minimum data (C as pointer) constructor
siconos::modeling::LagrangianLinearTIR::LagrangianLinearTIR(
    std::shared_ptr<siconos::algebra::SimpleMatrix> C)
    : LagrangianR(RelationSubType::LinearTIR) {
  _jachq = C;
}

// Constructor from a complete set of data
siconos::modeling::LagrangianLinearTIR::LagrangianLinearTIR(
    std::shared_ptr<siconos::algebra::SimpleMatrix> C,
    std::shared_ptr<siconos::algebra::SimpleMatrix> F,
    std::shared_ptr<siconos::algebra::SiconosVector> e)
    : LagrangianR(RelationSubType::LinearTIR) {
  _jachq = C;
  _F = F;
  _e = e;
}

// Minimum data (C, e as pointers) constructor
siconos::modeling::LagrangianLinearTIR::LagrangianLinearTIR(
    std::shared_ptr<siconos::algebra::SimpleMatrix> C,
    std::shared_ptr<siconos::algebra::SiconosVector> e)
    : LagrangianR(RelationSubType::LinearTIR) {
  _jachq = C;
  _e = e;
}

void siconos::modeling::LagrangianLinearTIR::checkSize(Interaction& inter) {
  auto sizeY = inter.dimension();
  auto& DSlink = inter.linkToDSVariables();
  if (!(_jachq) || _jachq->size(1) != inter.getSizeOfDS() || _jachq->size(0) != sizeY)
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianLinearTIR::checkSize inconsistent sizes between H "
        "matrix and the interaction.");

  if ((_e) && _e->size() != sizeY)
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianLinearTIR::checkSize inconsistent sizes between e "
        "vector and the dimension of the interaction.");

  auto sizeZ = DSlink[LagrangianR::z]->size();
  if ((_F) && (_F->size(0) != sizeZ || _F->size(1) != sizeZ))
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianLinearTIR::checkSize inconsistent sizes between F "
        "matrix and the interaction.");
}
void siconos::modeling::LagrangianLinearTIR::computeOutput(double time, Interaction& inter,
                                                           unsigned int derivativeNumber) {
  DEBUG_BEGIN(
      "siconos::modeling::LagrangianLinearTIR::computeOutput(double time, Interaction& inter, "
      "unsigned int derivativeNumber)\n");
  // get y and lambda of the interaction
  auto& y = *inter.y(derivativeNumber);
  auto& DSlink = inter.linkToDSVariables();
  siconos::algebra::prod(*_jachq, *DSlink[LagrangianR::q0 + derivativeNumber], y);

  if (derivativeNumber == 0) {
    if (_e) y += *_e;
    if (_F) siconos::algebra::prod(*_F, *DSlink[LagrangianR::z], y, false);
  }

  if (_jachlambda) {
    auto& lambda = *inter.lambda(derivativeNumber);
    siconos::algebra::prod(*_jachlambda, lambda, y, false);
  }
  DEBUG_END(
      "siconos::modeling::LagrangianLinearTIR::computeOutput(double time, Interaction& inter, "
      "unsigned int derivativeNumber)\n");
}
void siconos::modeling::LagrangianLinearTIR::computeInput(double time, Interaction& inter,
                                                          unsigned int level) {
  DEBUG_BEGIN(
      "void siconos::modeling::LagrangianLinearTIR::computeInput(double time, Interaction& "
      "inter, unsigned int level)\n")
  // get lambda of the concerned interaction
  siconos::algebra::SiconosVector& lambda = *inter.lambda(level);
  auto& DSlink = inter.linkToDSVariables();
  // computation of p = Ht lambda
  DEBUG_EXPR(lambda.display(););
  DEBUG_EXPR(_jachq->display(););
  DEBUG_EXPR(DSlink[LagrangianR::p0 + level]->display(););
  siconos::algebra::prod(lambda, *_jachq, *DSlink[LagrangianR::p0 + level], false);
  DEBUG_END(
      "void siconos::modeling::LagrangianLinearTIR::computeInput(double time, Interaction& "
      "inter, unsigned int level)\n")
}

void siconos::modeling::LagrangianLinearTIR::display() const {
  LagrangianR::display();
  std::cout << "===== Lagrangian Linear Relation display ===== " << std::endl;
  std::cout << " C: " << std::endl;
  if (_jachq)
    _jachq->display();
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << " e: " << std::endl;
  if (_e)
    _e->display();
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << " F: " << std::endl;
  if (_F)
    _F->display();
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << "===================================== " << std::endl;
}
