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
#include "FirstOrderR.hpp"

#include "BlockVector.hpp"
#include "FirstOrderNonLinearDS.hpp"
#include "Interaction.hpp"

void siconos::modeling::FirstOrderR::initialize(Interaction& inter) {
  inter.relationVectors().resize(FirstOrderR::relationVectorsSize);
}

void siconos::modeling::FirstOrderR::display() const {
  std::cout << "=====> Relation of type "
            << static_cast<std::underlying_type<RelationType>::type>(_relationType)
            << " and subtype "
            << static_cast<std::underlying_type<RelationSubType>::type>(_subType) << "\n";

  if (jacobianhOver_state_view_) {
    std::cout << " jacobian h over state\n";
    if (hasConstantJacobianhOver_state_) std::cout << "(constant matrix) \n";
    std::cout << *jacobianhOver_state_view_ << "\n";
  }
  if (jacobianhOver_lambda_view_) {
    std::cout << " jacobian h over lambda\n";
    if (hasConstantJacobianhOver_lambda_) std::cout << "(constant matrix) \n";
    std::cout << *jacobianhOver_lambda_view_ << "\n";
  }
  if (jacobiangOver_state_view_) {
    std::cout << " jacobian g over state\n";
    if (hasConstantJacobiangOver_state_) std::cout << "(constant matrix) \n";
    std::cout << *jacobiangOver_state_view_ << "\n";
  }
  if (jacobiangOver_lambda_view_) {
    std::cout << " jacobian g over lambda\n";
    if (hasConstantJacobiangOver_lambda_) std::cout << "(constant matrix) \n";
    std::cout << *jacobiangOver_lambda_view_ << "\n";
  }
}

void siconos::modeling::FirstOrderR::allocate_read_dynamical_systems_var_vectors(
    std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& ds_vars,
    modeling::DynamicalSystem& ds1, modeling::DynamicalSystem& ds2) const {
  if (ds_vars.empty()) ds_vars.resize(FirstOrderR::size);

  bool has_two_ds = &ds1 != &ds2;

  ds_vars[FirstOrderR::Xxx] = std::make_shared<siconos::algebra::BlockVector>();
  ds_vars[FirstOrderR::Rrr] = std::make_shared<siconos::algebra::BlockVector>();
  auto* fods1 = dynamic_cast<FirstOrderNonLinearDS*>(&ds1);
  ds_vars[FirstOrderR::Xxx]->insertPtr(fods1->x());
  ds_vars[FirstOrderR::Rrr]->insertPtr(fods1->r());

  if (has_two_ds) {
    auto* fods2 = dynamic_cast<FirstOrderNonLinearDS*>(&ds2);
    ds_vars[FirstOrderR::Xxx]->insertPtr(fods2->x());
    ds_vars[FirstOrderR::Rrr]->insertPtr(fods2->r());
  }
}
