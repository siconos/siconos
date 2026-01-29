/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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
#include "StressLinearTIR.hpp"

#include <iostream>

#include "Interaction.hpp"
#include "Tools.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::mechanics::fem::StressLinearTIR::checkSize(
    const siconos::modeling::Interaction& inter) const {
  auto sizeY = inter.dimension();
  const auto& DSlink = inter.linkToDSVariables();

  DSlink[LagrangianR::DSlinkSize + SolidLinearTIDS::sigma]->size();
  if (!(jacobianhOver_q_view_) ||
      jacobianhOver_q_view_->cols() !=
          DSlink[LagrangianR::DSlinkSize + SolidLinearTIDS::sigma]->size() ||
      jacobianhOver_q_view_->rows() != sizeY)
    THROW_EXCEPTION(
        "siconos::mechanics::fem::StressLinearTIR::checkSize inconsistent sizes between H "
        "matrix and the interaction.");

  if ((eVector_view_) && eVector_view_->size() != sizeY)
    THROW_EXCEPTION(
        "siconos::mechanics::fem::StressLinearTIR::checkSize inconsistent sizes "
        "between e vector and the dimension of the interaction.");
}

void siconos::mechanics::fem::StressLinearTIR::computeOutput(
    double time, siconos::modeling::Interaction& inter, unsigned int derivativeNumber) {
  DEBUG_BEGIN("siconos::mechanics::fem::StressLinearTIR::computeOutput(...)\n");
  // get y and lambda of the interaction
  auto& y = *inter.y(derivativeNumber);
  auto& DSlink = inter.linkToDSVariables();
  siconos::algebra::matrixBlockVector_prod(
      *jacobianhOver_q_view_,
      *DSlink[tools::enum_to_index(LagrangianR::DSlinkSize + SolidLinearTIDS::sigma) +
              derivativeNumber],
      y);

  if (derivativeNumber == 0) {
    if (eVector_view_) y += *eVector_view_;
  }
  DEBUG_END("siconos::mechanics::fem::StressLinearTIR::computeOutput(...)\n");
}
void siconos::mechanics::fem::StressLinearTIR::computeInput(
    double time, siconos::modeling::Interaction& inter, unsigned int level) {
  // get lambda of the concerned interaction
  // Here lambda = plastic Rate

  siconos::algebra::SiconosVector& lambda = *inter.lambda(level);
  auto& DSlink = inter.linkToDSVariables();

  // computation of p = Ht lambda
  DEBUG_EXPR(lambda.display(););
  DEBUG_EXPR(_jachq->display(););
  DEBUG_EXPR(DSlink[SolidLinearTIDS::dotEpsilon]->display(););

  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *jacobianhOver_q_view_,
      *DSlink[tools::enum_to_index(LagrangianR::DSlinkSize + SolidLinearTIDS::epsilonp) + level],
      false);
}

void siconos::mechanics::fem::StressLinearTIR::display() const {
  LagrangianR::display();
  std::cout << "===== Lagrangian Linear Relation display ===== \n";
  std::cout << " C: " << std::endl;
  std::cout << jacobianhOver_q_view_ << "•n";

  std::cout << " e: \n";
  if (eVector_view_)
    std::cout << eVector_view_ << "\n";
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << "===================================== " << std::endl;
}

void siconos::mechanics::fem::StressLinearTIR::allocate_dslink_vectors(
    std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& DSlink) const {
  DSlink.resize(tools::enum_to_index(StressLinearTIR::WorkDS::DSlinkSize));

  // Default DSlink
  DSlink[tools::enum_to_index(StressLinearTIR::WorkDS::sigma)] =
      std::make_shared<siconos::algebra::BlockVector>();
  DSlink[tools::enum_to_index(StressLinearTIR::WorkDS::sigma1)] =
      std::make_shared<siconos::algebra::BlockVector>();

  DSlink[tools::enum_to_index(StressLinearTIR::WorkDS::epsilonp)] =
      std::make_shared<siconos::algebra::BlockVector>();
  DSlink[tools::enum_to_index(StressLinearTIR::WorkDS::epsilonp)] =
      std::make_shared<siconos::algebra::BlockVector>();
}