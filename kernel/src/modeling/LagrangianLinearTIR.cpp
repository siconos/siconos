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

#include "Interaction.hpp"
#include "SiconosAlgebraAddons.hpp"  // for matrix-vector prod
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "SiconosException.hpp"
#include "Tools.hpp"
#include "siconos_debug.h"

// Minimum data (C) constructor
siconos::modeling::LagrangianLinearTIR::LagrangianLinearTIR(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newC)
    : LagrangianR{RelationSubType::LinearTIR} {
  jacobianhOver_q_view_ =
      std::make_shared<siconos::algebra::MapType>(newC.data(), newC.rows(), newC.cols());
}

siconos::modeling::LagrangianLinearTIR::LagrangianLinearTIR(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newC,
    Eigen::Ref<siconos::algebra::SiconosVector> newe)
    : LagrangianR{RelationSubType::LinearTIR} {
  jacobianhOver_q_view_ =
      std::make_shared<siconos::algebra::MapType>(newC.data(), newC.rows(), newC.cols());

  eVector_view_ = std::make_shared<siconos::algebra::MapVectorType>(newe.data(), newe.size());
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
        "siconos::modeling::LagrangianLinearTIR::checkSize inconsistent sizes "
        "between e vector and the dimension of the interaction.");
}

void siconos::modeling::LagrangianLinearTIR::computeOutput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type derivativeNumber) {
  DEBUG_BEGIN("siconos::modeling::LagrangianLinearTIR::computeOutput()\n");
  // get y and lambda of the interaction
  auto& y = *inter.y(derivativeNumber);
  const auto& ds_vars = inter.read_dynamical_systems_variables();

  siconos::algebra::matrixBlockVector_prod(
      *jacobianhOver_q_view_,
      *ds_vars[tools::enum_to_index(LagrangianR::ds_var::q0) + derivativeNumber], y);

  if (derivativeNumber == 0) {
    if (eVector_view_) y += *eVector_view_;
  }

  DEBUG_END("siconos::modeling::LagrangianLinearTIR::computeOutput()\n");
}

void siconos::modeling::LagrangianLinearTIR::computeInput(double time, Interaction& inter,
                                                          siconos::algebra::blocks::size_type level) {
  DEBUG_BEGIN("void siconos::modeling::LagrangianLinearTIR::computeInput()\n")
  // get lambda of the concerned interaction
  siconos::algebra::SiconosVector& lambda = *inter.lambda(level);
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  // computation of p = Ht lambda
  DEBUG_EXPR(siconos::algebra::print(lambda););
  DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_););
  DEBUG_EXPR(siconos::algebra::print(
                 *ds_vars[tools::enum_to_index(LagrangianR::ds_var::p0) + level]););

  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *jacobianhOver_q_view_,
      *ds_vars[tools::enum_to_index(LagrangianR::ds_var::p0) + level], false);

  DEBUG_END("void siconos::modeling::LagrangianLinearTIR::computeInput()\n")
}

void siconos::modeling::LagrangianLinearTIR::display() const {
  LagrangianR::display();
  std::cout << "===== Lagrangian Linear Relation display ===== " << std::endl;
  std::cout << " C:\n";
  siconos::algebra::print(*jacobianhOver_q_view_);
  std::cout << " e: " << std::endl;
  if (eVector_view_)
    std::cout << eVector_view_ << "\n";
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << "===================================== " << std::endl;
}
