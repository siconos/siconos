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

#include <condition_variable>
#include <iostream>
#include <memory>

#include "BlockVector.hpp"
#include "DynamicalSystem.hpp"
#include "Interaction.hpp"
#include "SecondOrderDS.hpp"
#include "SiconosAlgebraAddons.hpp"  // matrixBlockVector_prod
#include "SolidLinearTIDS.hpp"
#include "Tools.hpp"  // enum_to_index
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::mechanics::fem::StressLinearTIR::checkSize(
    const siconos::modeling::Interaction& inter) const {
  auto sizeY = inter.dimension();
  const auto& ds_vars = inter.read_dynamical_systems_variables();

  if (!(jacobianhOver_q_view_) ||
      jacobianhOver_q_view_->cols() != ds_vars[tools::enum_to_index(ds_var::sigma)]->size() ||
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
    double time, siconos::modeling::Interaction& inter,
    siconos::algebra::blocks::size_type derivativeNumber) {
  DEBUG_BEGIN("siconos::mechanics::fem::StressLinearTIR::computeOutput(...)\n");
  // get y and lambda of the interaction
  auto& y = *inter.y(derivativeNumber);
  const auto& ds_vars = inter.read_dynamical_systems_variables();

  // y[i] = C * sigma_i (i=0 or 1)
  siconos::algebra::matrixBlockVector_prod(
      *jacobianhOver_q_view_, *ds_vars[tools::enum_to_index(ds_var::sigma) + derivativeNumber],
      y);

  // y[0] = C * sigma_0 + e
  if (derivativeNumber == 0) {
    if (eVector_view_) y += *eVector_view_;
  }
  DEBUG_END("siconos::mechanics::fem::StressLinearTIR::computeOutput(...)\n");
}
void siconos::mechanics::fem::StressLinearTIR::computeInput(
    double time, siconos::modeling::Interaction& inter,
    siconos::algebra::blocks::size_type level) {
  // get lambda of the concerned interaction
  // Here lambda = plastic Rate

  siconos::algebra::SiconosVector& lambda = *inter.lambda(level);
  const auto& ds_vars = inter.read_dynamical_systems_variables();

  // computation of p = Ht lambda
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *jacobianhOver_q_view_,
      *ds_vars[tools::enum_to_index(ds_var::epsilon_p) + level], false);
}

void siconos::mechanics::fem::StressLinearTIR::display() const {
  std::cout << "======= StressLinearTIR display ======= \n";
  LagrangianLinearTIR::display();
}

void siconos::mechanics::fem::StressLinearTIR::allocate_read_dynamical_systems_var_vectors(
    std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& ds_vars,
    modeling::DynamicalSystem& ds1, modeling::DynamicalSystem& ds2) const {
  // Put q, velocity of each DS into a block. (Pointers links, no copy!!)
  if (ds_vars.empty()) ds_vars.resize(tools::enum_to_index(StressLinearTIR::ds_var::size));

  // Default ds_vars.
  // We need:
  // - sigma, sigma1 to compute output (y)
  // - epsilon_p, epsilon_p1 to compute input (lambda)
  bool has_two_ds = &ds1 != &ds2;

  // if (!ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::sigma)])
  ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::sigma)] =
      std::make_shared<siconos::algebra::BlockVector>();
  // if (!ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::sigma1)])
  ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::sigma1)] =
      std::make_shared<siconos::algebra::BlockVector>();
  // if (!ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::epsilon_p)])
  ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::epsilon_p)] =
      std::make_shared<siconos::algebra::BlockVector>();
  // if (!ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::epsilon_p1)])
  ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::epsilon_p1)] =
      std::make_shared<siconos::algebra::BlockVector>();

  auto* solid1 = dynamic_cast<SolidLinearTIDS*>(&ds1);
  ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::sigma)]->insertPtr(solid1->stress());
  ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::sigma1)]->insertPtr(
      solid1->stressRate());
  ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::epsilon_p)]->insertPtr(
      solid1->plasticDeformation());
  ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::epsilon_p1)]->insertPtr(
      solid1->plasticRate());

  if (has_two_ds) {
    auto* solid2 = dynamic_cast<SolidLinearTIDS*>(&ds2);
    ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::sigma)]->insertPtr(solid2->stress());
    ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::sigma1)]->insertPtr(
        solid2->stressRate());
    ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::epsilon_p)]->insertPtr(
        solid2->plasticDeformation());
    ds_vars[tools::enum_to_index(StressLinearTIR::ds_var::epsilon_p1)]->insertPtr(
        solid2->plasticRate());
  }
}