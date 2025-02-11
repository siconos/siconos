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

#include "Interaction.hpp"
#include "SiconosVisitor.hpp"

void siconos::modeling::FirstOrderR::initialize(Interaction& inter) {
  inter.relationVectors().resize(FirstOrderR::relationVectorsSize);
}

void siconos::modeling::FirstOrderR::display() const {
  std::cout << "=====> Relation of type "
            << static_cast<std::underlying_type<RelationType>::type>(_relationType)
            << " and subtype "
            << static_cast<std::underlying_type<RelationSubType>::type>(_subType) << "\n";

  if (jacobianhOver_state_view_) {
    std::cout << " jacobian h over state";
    if (hasConstantJacobianhOver_state_) std::cout << "(constant matrix) \n";
    std::cout << jacobianhOver_state_view_ << "\n";
  }
  if (jacobianhOver_lambda_view_) {
    std::cout << " jacobian h over lambda";
    if (hasConstantJacobianhOver_lambda_) std::cout << "(constant matrix) \n";
    std::cout << jacobianhOver_lambda_view_ << "\n";
  }
  if (jacobiangOver_state_view_) {
    std::cout << " jacobian g over state";
    if (hasConstantJacobiangOver_state_) std::cout << "(constant matrix) \n";
    std::cout << jacobiangOver_state_view_ << "\n";
  }
  if (jacobiangOver_lambda_view_) {
    std::cout << " jacobian g over lambda";
    if (hasConstantJacobiangOver_lambda_) std::cout << "(constant matrix) \n";
    std::cout << jacobiangOver_lambda_view_ << "\n";
  }
}