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
#include "LagrangianCompliantLinearTIR.hpp"

#include <iostream>

// #include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosAlgebraAddons.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Tools.hpp"

// Minimum data (C as pointer) constructor
siconos::modeling::LagrangianCompliantLinearTIR::LagrangianCompliantLinearTIR(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newC,
    Eigen::Ref<siconos::algebra::SiconosMatrix> newD)
    : LagrangianR{RelationSubType::CompliantLinearTIR} {
  jacobianhOver_q_view_ =  // C
      std::make_shared<siconos::algebra::MapType>(newC.data(), newC.rows(), newC.cols());

  DMatrix_view_ =
      std::make_shared<siconos::algebra::MapType>(newD.data(), newD.rows(), newD.cols());
}

siconos::modeling::LagrangianCompliantLinearTIR::LagrangianCompliantLinearTIR(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newC,
    Eigen::Ref<siconos::algebra::SiconosMatrix> newD,
    Eigen::Ref<siconos::algebra::SiconosVector> newe)
    : LagrangianR{RelationSubType::CompliantLinearTIR} {
  jacobianhOver_q_view_ =  // C
      std::make_shared<siconos::algebra::MapType>(newC.data(), newC.rows(), newC.cols());

  DMatrix_view_ =
      std::make_shared<siconos::algebra::MapType>(newD.data(), newD.rows(), newD.cols());

  eVector_view_ = std::make_shared<siconos::algebra::MapVectorType>(newe.data(), newe.size());
}

void siconos::modeling::LagrangianCompliantLinearTIR::initialize(Interaction& inter) {
  checkSize(inter);
}

void siconos::modeling::LagrangianCompliantLinearTIR::checkSize(
    const Interaction& inter) const {
  auto sizeY = inter.dimension();

  if (!(jacobianhOver_q_view_) || jacobianhOver_q_view_->cols() != inter.getSizeOfDS() ||
      jacobianhOver_q_view_->rows() != sizeY)
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianCompliantLinearTIR::checkSize inconsistent sizes "
        "between H matrix and the interaction.");

  if ((DMatrix_view_) && (DMatrix_view_->rows() != sizeY || DMatrix_view_->cols() != sizeY))
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianCompliantLinearTIR::checkSize inconsistent sizes "
        "between D matrix and the interaction.");

  if ((eVector_view_) && eVector_view_->size() != sizeY)
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianCompliantLinearTIR::checkSize inconsistent sizes "
        "between e vector and the dimension of the interaction.");
}

void siconos::modeling::LagrangianCompliantLinearTIR::computeInput(double time,
                                                                   Interaction& inter,
                                                                   unsigned int level) {
  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);
  auto& DSlink = inter.linkToDSVariables();
  // computation of p = Ht lambda
  siconos::algebra::transposeMatrixVector_prod_toBlock(
      lambda, *jacobianhOver_q_view_, *DSlink[tools::enum_to_index(WorkDS::p0) + level],
      false);
}
void siconos::modeling::LagrangianCompliantLinearTIR::computeOutput(
    double time, Interaction& inter, unsigned int derivativeNumber) {
  // get y and lambda of the interaction
  auto& y = *inter.y(derivativeNumber);
  auto& lambda = *inter.lambda(derivativeNumber);
  auto& DSlink = inter.linkToDSVariables();
  siconos::algebra::matrixBlockVector_prod(
      *jacobianhOver_q_view_, *DSlink[tools::enum_to_index(WorkDS::q0) + derivativeNumber], y);
  y += *DMatrix_view_ * lambda;
  if (derivativeNumber == 0) {
    if (eVector_view_) y += *eVector_view_;
  }
}

void siconos::modeling::LagrangianCompliantLinearTIR::display() const {
  LagrangianR::display();
  std::cout << "===== Lagrangian Linear Relation display ===== " << std::endl;
  std::cout << " C: \n" << *jacobianhOver_q_view_ << "\n";
  std::cout << " e: " << std::endl;
  if (eVector_view_)
    std::cout << *eVector_view_ << "\n";
  else
    std::cout << " -> nullptr " << std::endl;
  std::cout << " D: \n" << *DMatrix_view_ << "\n";
  std::cout << "===================================== " << std::endl;
}
