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
#include "LagrangianLinearTIR.hpp"

#include <iostream>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "SiconosException.hpp"
#include "siconos_debug.h"

// Minimum data (C as pointer) constructor
siconos::modeling::LagrangianLinearTIR::LagrangianLinearTIR(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newC)
    : LagrangianR(RelationSubType::LinearTIR) {
  jacobianhOver_q_view_ =
      std::make_shared<siconos::algebra::MapType>(newC.data(), newC.rows(), newC.cols());
}

void siconos::modeling::LagrangianLinearTIR::checkSize(const Interaction& inter) const {
  auto sizeY = inter.dimension();

  if (!(jacobianhOver_q_view_) || jacobianhOver_q_view_->cols() != inter.getSizeOfDS() ||
      jacobianhOver_q_view_->rows() != sizeY)
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianLinearTIR::checkSize inconsistent sizes "
        "between H "
        "matrix and the interaction.");

  if ((eVector_view_) && eVector_view_->size() != sizeY)
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianCompliantLinearTIR::checkSize inconsistent sizes "
        "between e vector and the dimension of the interaction.");
}

void siconos::modeling::LagrangianLinearTIR::seteVector(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  eVector_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
}
void siconos::modeling::LagrangianLinearTIR::computeOutput(double time, Interaction& inter,
                                                           unsigned int derivativeNumber) {
  DEBUG_BEGIN(
      "siconos::modeling::LagrangianLinearTIR::computeOutput(double time, "
      "Interaction& inter, "
      "unsigned int derivativeNumber)\n");
  // get y and lambda of the interaction
  auto& y = *inter.y(derivativeNumber);
  auto& DSlink = inter.linkToDSVariables();
  siconos::algebra::matrixBlockVector_prod(*jacobianhOver_q_view_,
                                           *DSlink[LagrangianR::q0 + derivativeNumber], y);

  if (derivativeNumber == 0) {
    if (eVector_view_) y += *eVector_view_;
  }
  DEBUG_END(
      "siconos::modeling::LagrangianLinearTIR::computeOutput(double time, "
      "Interaction& inter, "
      "unsigned int derivativeNumber)\n");
}

void siconos::modeling::LagrangianLinearTIR::computeInput(double time, Interaction& inter,
                                                          unsigned int level) {
  DEBUG_BEGIN(
      "void siconos::modeling::LagrangianLinearTIR::computeInput(double time, "
      "Interaction& "
      "inter, unsigned int level)\n")
  // get lambda of the concerned interaction
  siconos::algebra::SiconosVector& lambda = *inter.lambda(level);
  auto& DSlink = inter.linkToDSVariables();
  // computation of p = Ht lambda
  DEBUG_EXPR(lambda.display(););
  DEBUG_EXPR(jacobianhOver_q_->display(););
  DEBUG_EXPR(DSlink[LagrangianR::p0 + level]->display(););
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *jacobianhOver_q_view_, *DSlink[LagrangianR::p0 + level], false);
  DEBUG_END(
      "void siconos::modeling::LagrangianLinearTIR::computeInput(double time, "
      "Interaction& "
      "inter, unsigned int level)\n")
}

void siconos::modeling::LagrangianLinearTIR::display() const {
  LagrangianR::display();
  std::cout << "===== Lagrangian Linear Relation display ===== " << std::endl;
  std::cout << " C:\n" << jacobianhOver_q_view_ << "\n";
  std::cout << " e: " << std::endl;
  if (eVector_view_)
    std::cout << eVector_view_ << "\n";
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << "===================================== " << std::endl;
}
