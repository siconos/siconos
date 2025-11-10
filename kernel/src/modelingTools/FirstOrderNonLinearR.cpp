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

#include "FirstOrderNonLinearR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
//#include "siconos_debug.h"

void siconos::modeling::FirstOrderNonLinearR::initialize(Interaction& inter) {
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

  if (computejacobianhOver_lambda_) {
    if (!jacobianhOver_lambda_internal_storage_) {
      jacobianhOver_lambda_internal_storage_ =
          std::make_unique<std::vector<double>>(sizeY * sizeY);
    }
    jacobianhOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
        jacobianhOver_lambda_internal_storage_->data(), sizeY, sizeY);
  }

  if (computejacobiangOver_state_) {
    if (!jacobiangOver_state_internal_storage_) {
      jacobiangOver_state_internal_storage_ =
          std::make_unique<std::vector<double>>(sizeX * sizeX);
    }
    jacobiangOver_state_view_ = std::make_shared<siconos::algebra::MapType>(
        jacobiangOver_state_internal_storage_->data(), sizeX, sizeX);
  }

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

void siconos::modeling::FirstOrderNonLinearR::checkSize(const Interaction& inter) const {
  // get inter and ds sizes
  auto sizeY = inter.dimension();
  auto sizeX = inter.getSizeOfDS();

  // Check if various operators sizes are consistent.
  // Reference: interaction.

  assert(jacobianhOver_state_view_);   // We do not allow relations without this operator.
  assert(jacobiangOver_lambda_view_);  // We do not allow relations without this operator.

  if (jacobianhOver_state_view_) {
    assert(jacobianhOver_state_view_->rows() == sizeY);
    assert(jacobianhOver_state_view_->cols() == sizeX);
  }
  if (jacobianhOver_lambda_view_) {
    assert(jacobianhOver_lambda_view_->rows() == sizeY);
    assert(jacobianhOver_lambda_view_->cols() == sizeY);
  }

  if (jacobiangOver_state_view_) {
    assert(jacobiangOver_state_view_->rows() == sizeX);
    assert(jacobiangOver_state_view_->cols() == sizeX);
  }

  if (jacobiangOver_lambda_view_) {
    assert(jacobiangOver_lambda_view_->rows() == sizeX);
    assert(jacobiangOver_lambda_view_->cols() == sizeY);
  }
}

void siconos::modeling::FirstOrderNonLinearR::setComputehFunction(
    const siconos::modeling::func_prototypes::FunctionBVSV_V& h_func) {
  computeh_ = h_func;
}

void siconos::modeling::FirstOrderNonLinearR::computeh(
    const siconos::algebra::BlockVector& state, double time,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  if (computeh_) computeh_(state, time, lambda, y);
}

void siconos::modeling::FirstOrderNonLinearR::setComputegFunction(
    const siconos::modeling::func_prototypes::FunctionBVSV_BV& g_func) {
  computeg_ = g_func;
}

void siconos::modeling::FirstOrderNonLinearR::computeg(
    const siconos::algebra::BlockVector& state, double time,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
    siconos::algebra::BlockVector& result) {
  if (computeg_) computeg_(state, time, lambda, result);
}

void siconos::modeling::FirstOrderNonLinearR::setConstantJacobianhOver_state(
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue) {
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

void siconos::modeling::FirstOrderNonLinearR::setComputeJacobianhOver_stateFunction(
    const siconos::modeling::func_prototypes::FunctionBVSV_M& func) {
  computejacobianhOver_state_ = func;
}

void siconos::modeling::FirstOrderNonLinearR::computeJacobianhOver_state(
    const siconos::algebra::BlockVector& state, double time,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
  if (computejacobianhOver_state_)
    computejacobianhOver_state_(state, time, lambda, *jacobianhOver_state_view_);
}

void siconos::modeling::FirstOrderNonLinearR::setConstantJacobianhOver_lambda(
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue) {
  /**  Warning: we can't check if newValue size is compliant with
   * the current relation/interaction. This will be done during
   * initialize / checkSize call!
   */

  jacobianhOver_lambda_internal_storage_ = nullptr;

  jacobianhOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantJacobianhOver_lambda_ = true;
  computejacobianhOver_lambda_ = nullptr;
}

void siconos::modeling::FirstOrderNonLinearR::setComputeJacobianhOver_lambdaFunction(
    const siconos::modeling::func_prototypes::FunctionBVSV_M& func) {
  computejacobianhOver_lambda_ = func;
}

void siconos::modeling::FirstOrderNonLinearR::computeJacobianhOver_lambda(
    const siconos::algebra::BlockVector& state, double time,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
  if (computejacobianhOver_lambda_)
    computejacobianhOver_lambda_(state, time, lambda, *jacobianhOver_lambda_view_);
}

void siconos::modeling::FirstOrderNonLinearR::setConstantJacobiangOver_state(
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue) {
  /**  Warning: we can't check if newValue size is compliant with
   * the current relation/interaction. This will be done during
   * initialize / checkSize call!
   */

  jacobiangOver_state_internal_storage_ = nullptr;

  jacobiangOver_state_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantJacobiangOver_state_ = true;
  computejacobiangOver_state_ = nullptr;
}

void siconos::modeling::FirstOrderNonLinearR::setComputeJacobiangOver_stateFunction(
    const siconos::modeling::func_prototypes::FunctionBVSV_M& func) {
  computejacobiangOver_state_ = func;
}

void siconos::modeling::FirstOrderNonLinearR::computeJacobiangOver_state(
    const siconos::algebra::BlockVector& state, double time,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
  if (computejacobiangOver_state_)
    computejacobiangOver_state_(state, time, lambda, *jacobiangOver_state_view_);
}

void siconos::modeling::FirstOrderNonLinearR::setConstantJacobiangOver_lambda(
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue) {
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

void siconos::modeling::FirstOrderNonLinearR::setComputeJacobiangOver_lambdaFunction(
    const siconos::modeling::func_prototypes::FunctionBVSV_M& func) {
  computejacobiangOver_lambda_ = func;
}

void siconos::modeling::FirstOrderNonLinearR::computeJacobiangOver_lambda(
    const siconos::algebra::BlockVector& state, double time,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda) {
  if (computejacobiangOver_lambda_)
    computejacobiangOver_lambda_(state, time, lambda, *jacobiangOver_lambda_view_);
}

void siconos::modeling::FirstOrderNonLinearR::computeOutput(double time, Interaction& inter,
                                                            unsigned int level) {
  auto& DSlink = inter.linkToDSVariables();
  auto& y = *inter.y(level);
  auto& lambda = *inter.lambda(level);
  if (computeh_) computeh_(*DSlink[FirstOrderR::Xxx], time, lambda, y);
}

void siconos::modeling::FirstOrderNonLinearR::computeInput(double time, Interaction& inter,
                                                           unsigned int level) {
  auto& DSlink = inter.linkToDSVariables();
  auto& lambda = *inter.lambda(level);
  if (computeg_) computeg_(*DSlink[FirstOrderR::Xxx], time, lambda, *DSlink[FirstOrderR::Rrr]);
}

void siconos::modeling::FirstOrderNonLinearR::computeJach(double time, Interaction& inter) {
  auto& DSlink = inter.linkToDSVariables();
  auto& lambda = *inter.lambda(0);

  if (computejacobianhOver_state_) {
    computejacobianhOver_state_(*DSlink[FirstOrderR::Xxx], time, lambda,
                                *jacobianhOver_state_view_);
  }
  if (computejacobianhOver_lambda_) {
    computejacobianhOver_lambda_(*DSlink[FirstOrderR::Xxx], time, lambda,
                                 *jacobianhOver_lambda_view_);
  }
}

void siconos::modeling::FirstOrderNonLinearR::computeJacg(double time, Interaction& inter) {
  auto& DSlink = inter.linkToDSVariables();

  auto& lambda = *inter.lambda(0);
  if (computejacobiangOver_lambda_) {
    computejacobiangOver_lambda_(*DSlink[FirstOrderR::Xxx], time, lambda,
                                 *jacobiangOver_lambda_view_);
  }
  if (computejacobiangOver_state_) {
    computejacobiangOver_state_(*DSlink[FirstOrderR::Xxx], time, lambda,
                                *jacobiangOver_state_view_);
  }
}
