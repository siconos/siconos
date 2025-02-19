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
#include "FirstOrderType1R.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

void siconos::modeling::FirstOrderType1R::initialize(Interaction& inter) {
  FirstOrderR::initialize(inter);

  auto sizeY = inter.dimension();
  auto sizeX = inter.getSizeOfDS();

  if (computejacobianhOver_state_) {
    if (!jacobianhOver_state_internal_storage_) {
      jacobianhOver_state_internal_storage_ =
          std::make_unique<std::vector<double>>(sizeY * sizeX);
    }
    jacobianhOver_state_view_ = std::make_shared<siconos::algebra::MapType>(
        jacobianhOver_state_internal_storage_->data(), sizeY, sizeX);
  }
  // if the jacobian is a constant matrix (following a call to setConstant...)
  // then
  // - no need to allocate internal storage
  // - map/view is already 'connected' to some external memory

  if (computejacobiangOver_lambda_) {
    if (!jacobiangOver_lambda_internal_storage_) {
      jacobiangOver_lambda_internal_storage_ =
          std::make_unique<std::vector<double>>(sizeX * sizeY);
    }
    jacobiangOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
        jacobiangOver_lambda_internal_storage_->data(), sizeX, sizeY);
  }
  checkSize(inter);
}

void siconos::modeling::FirstOrderType1R::checkSize(const Interaction& inter) const {
  // get inter and ds sizes
  auto sizeY = inter.dimension();
  auto sizeX = inter.getSizeOfDS();

  // Check if various operators sizes are consistent.
  // Reference: interaction.

  if (jacobianhOver_state_view_) {
    assert(jacobianhOver_state_view_->rows() == sizeY);
    assert(jacobianhOver_state_view_->cols() == sizeX);
  }
  if (jacobiangOver_lambda_view_) {
    assert(jacobiangOver_lambda_view_->rows() == sizeX);
    assert(jacobiangOver_lambda_view_->cols() == sizeY);
  }
}

void siconos::modeling::FirstOrderType1R::setComputehFunction(
    const siconos::modeling::func_prototypes::FunctionBV_V& h_func) {
  computeh_ = h_func;
}

void siconos::modeling::FirstOrderType1R::computeh(
    const siconos::algebra::BlockVector& state,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  if (computeh_) computeh_(state, y);
}

void siconos::modeling::FirstOrderType1R::setComputegFunction(
    const siconos::modeling::func_prototypes::FunctionV_BV& g_func) {
  computeg_ = g_func;
}

void siconos::modeling::FirstOrderType1R::computeg(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
    siconos::algebra::BlockVector& result) {
  if (computeg_) computeg_(lambda, result);
}

void siconos::modeling::FirstOrderType1R::setConstantJacobianhOver_state(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Warning: we can't check if newValue size is compliant with
   * the current relation/interaction. This will be done during
   * initialize / checkSize call!
   */

  jacobianhOver_state_internal_storage_ = nullptr;

  jacobianhOver_state_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantJacobianhOver_state_ = true;
  computejacobianhOver_state_ = nullptr;
}

void siconos::modeling::FirstOrderType1R::setComputeJacobianhOver_stateFunction(
    const siconos::modeling::func_prototypes::FunctionBV_M& func) {
  computejacobianhOver_state_ = func;
}

void siconos::modeling::FirstOrderType1R::computeJacobianhOver_state(
    const siconos::algebra::BlockVector& state) {
  if (computejacobianhOver_state_)
    computejacobianhOver_state_(state, *jacobianhOver_state_view_);
}

void siconos::modeling::FirstOrderType1R::setConstantJacobiangOver_lambda(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Warning: we can't check if newValue size is compliant with
   * the current relation/interaction. This will be done during
   * initialize / checkSize call!
   */

  jacobiangOver_lambda_internal_storage_ = nullptr;

  jacobiangOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantJacobiangOver_lambda_ = true;
  computejacobiangOver_lambda_ = nullptr;
}

void siconos::modeling::FirstOrderType1R::setComputeJacobiangOver_lambdaFunction(
    const siconos::modeling::func_prototypes::FunctionV_M& func) {
  hasConstantJacobiangOver_lambda_ = false;
  computejacobiangOver_lambda_ = func;
}

void siconos::modeling::FirstOrderType1R::computeJacobiangOver_lambda(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
  if (computejacobiangOver_lambda_)
    computejacobiangOver_lambda_(lambda, *jacobiangOver_lambda_view_);
}

void siconos::modeling::FirstOrderType1R::computeOutput(double time, Interaction& inter,
                                                        unsigned int level) {
  siconos::algebra::SiconosVector& y = *inter.y(0);
  auto& DSlink = inter.linkToDSVariables();
  computeh_(*DSlink[FirstOrderR::Xxx], y);
}

void siconos::modeling::FirstOrderType1R::computeInput(double time, Interaction& inter,
                                                       unsigned int level) {
  auto lambda = inter.lambda(level);

  auto& DSlink = inter.linkToDSVariables();
  computeg_(*lambda, *DSlink[FirstOrderR::Rrr]);
}

void siconos::modeling::FirstOrderType1R::computeJach(double time, Interaction& inter) {
  if (!hasConstantJacobianhOver_state_) {
    auto& DSlink = inter.linkToDSVariables();
    if (computejacobianhOver_state_)
      computejacobianhOver_state_(*DSlink[FirstOrderR::Xxx], *jacobianhOver_state_view_);
  }
}

void siconos::modeling::FirstOrderType1R::computeJacg(double time, Interaction& inter) {
  if (!hasConstantJacobiangOver_lambda_) {
    computejacobiangOver_lambda_(*inter.lambda(0), *jacobiangOver_lambda_view_);
  }
}
