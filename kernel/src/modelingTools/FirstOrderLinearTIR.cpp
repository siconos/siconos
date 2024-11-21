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
#include "FirstOrderLinearTIR.hpp"

#include <iostream>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::modeling::FirstOrderLinearTIR::initialize(Interaction &inter) {
  DEBUG_PRINT("siconos::modeling::FirstOrderLinearTIR::initialize(Interaction & inter)\n");

  FirstOrderR::initialize(inter);  // ?

  // if C (or B, D, e)  is constant (following a call to setConstantC)
  // then
  // - no need to allocate internal storage
  // - map/view is already 'connected' to some external memory

  // Shouldn't we allow C==null or B==null?
  assert(jacobianhOver_state_view_);   // "C is null and is a required  ")
  assert(jacobiangOver_lambda_view_);  // "B is null and is a required  ")
  checkSize(inter);
}

void siconos::modeling::FirstOrderLinearTIR::checkSize(Interaction &inter) {
  // get inter and ds sizes
  auto sizeY = inter.dimension();
  auto sizeX = inter.getSizeOfDS();

  if (jacobianhOver_state_view_) {
    assert(jacobianhOver_state_view_->rows() == sizeY);
    assert(jacobianhOver_state_view_->cols() == sizeX);
  }
  if (jacobiangOver_lambda_view_) {
    assert(jacobiangOver_lambda_view_->rows() == sizeX);
    assert(jacobiangOver_lambda_view_->cols() == sizeY);
  }
  if (jacobianhOver_lambda_view_) {
    assert(jacobianhOver_lambda_view_->rows() == sizeY);
    assert(jacobianhOver_lambda_view_->cols() == sizeY);
  }
  if (eVector_view_) {
    assert(eVector_view_->size() == sizeY);
  }
}
void siconos::modeling::FirstOrderLinearTIR::setConstantB(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Warning: we can't check if newValue size is compliant with
   * the current relation/interaction. This will be done during
   * initialize / checkSize call!
   */

  jacobiangOver_lambda_internal_storage_ = nullptr;

  jacobiangOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantJacobiangOver_lambda_ = true;
}

void siconos::modeling::FirstOrderLinearTIR::setConstantC(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Warning: we can't check if newValue size is compliant with
   * the current relation/interaction. This will be done during
   * initialize / checkSize call!
   */

  jacobianhOver_state_internal_storage_ = nullptr;

  jacobianhOver_state_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantJacobianhOver_state_ = true;
}

void siconos::modeling::FirstOrderLinearTIR::setConstantD(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Warning: we can't check if newValue size is compliant with
   * the current relation/interaction. This will be done during
   * initialize / checkSize call!
   */

  jacobianhOver_lambda_internal_storage_ = nullptr;

  jacobianhOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantJacobianhOver_lambda_ = true;
}

void siconos::modeling::FirstOrderLinearTIR::setConstanteVector(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  eVector_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
  haseVector_ = true;
}

// void siconos::modeling::FirstOrderLinearTIR::computeh(
//     const siconos::algebra::BlockVector &x, const siconos::algebra::SiconosVector &lambda,
//     siconos::algebra::BlockVector &z, siconos::algebra::SiconosVector &y) {
//   // if (_C) C must be allocated. Checksize is there to ensure it.
//   siconos::algebra::matrixBlockVector_prod(*_C, x, y, true);
//   // else
//   //   y.setZero();

//   if (_D) siconos::algebra::prod(*_D, lambda, y, false);

//   if (_e) y += *_e;

//   if (_F) siconos::algebra::matrixBlockVector_prod(*_F, z, y, false);
// }

void siconos::modeling::FirstOrderLinearTIR::computeOutput(double time, Interaction &inter,
                                                           unsigned int level) {
  // We get y and lambda of the interaction (pointers)
  siconos::algebra::SiconosVector &y = *inter.y(level);
  siconos::algebra::SiconosVector &lambda = *inter.lambda(level);
  auto &DSlink = inter.linkToDSVariables();
  if (jacobianhOver_state_view_) {
    siconos::algebra::matrixBlockVector_prod(*jacobianhOver_state_view_, x, y, true);
  } else
    y.setZero();

  if (jacobianhOver_lambda_view_) {
    y += *jacobianhOver_lambda_view_ * lambda;
  }

  if (eVector_view_) y += *eVector_view_;
}

// void siconos::modeling::FirstOrderLinearTIR::computeg(
//     const siconos::algebra::SiconosVector &lambda, siconos::algebra::BlockVector &r) {
//   siconos::algebra::matrixVector_prod_toBlock(*_B, lambda, r, false);
// }

void siconos::modeling::FirstOrderLinearTIR::computeInput(double time, Interaction &inter,
                                                          unsigned int level) {
  DEBUG_BEGIN(
      "siconos::modeling::FirstOrderLinearTIR::computeInput(double time, Interaction& "
      "inter, unsigned int level)\n")
  auto &DSlink = inter.linkToDSVariables();
  DEBUG_EXPR(inter.lambda(level)->display(););
  DEBUG_EXPR(DSlink[FirstOrderR::r]->display(););
  if (jacobiangOver_lambda_view_) {
    auto &DSlink = inter.linkToDSVariables();
    *DSlink[FirstOrderR::r] += *jacobiangOver_lambda_view_ * *inter.lambda(0);
  }
  DEBUG_END(
      "siconos::modeling::FirstOrderLinearTIR::computeInput(double time, Interaction& "
      "inter, unsigned int level)\n")
}

void siconos::modeling::FirstOrderLinearTIR::display() const {
  std::cout << " ===== Linear Time Invariant relation display =====\n";
  std::cout << "| C \n";
  if (jacobianhOver_state_view_)
    std::cout << jacobianhOver_state_view_ << "\n";
  else
    std::cout << "->nullptr\n";

  std::cout << "| D\n";
  if (jacobianhOver_lambda_view_)
    std::cout << jacobianhOver_lambda_view_ << "\n";
  else
    std::cout << "->nullptr\n";

  std::cout << "| e\n";
  if (eVector_view_)
    std::cout << eVector_view_ << "\n";
  else
    std::cout << "->nullptr\n";

  std::cout << "| B \n";
  if (jacobiangOver_lambda_view_)
    std::cout << jacobiangOver_lambda_view_ << "\n";
  else
    std::cout << "->nullptr\n";
  std::cout << " ================================================== \n";
}
