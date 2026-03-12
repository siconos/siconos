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
#include "FirstOrderLinearR.hpp"

#include <iostream>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "SiconosAlgebraAddons.hpp"
#include "SiconosMatrix.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <cassert>

#include "siconos_debug.h"

void siconos::modeling::FirstOrderLinearR::initialize(Interaction& inter) {
  FirstOrderR::initialize(inter);

  // get interesting size
  auto sizeY = inter.dimension();
  auto sizeX = inter.getSizeOfDS();

  if (computeC_) {
    if (!jacobianhOver_state_internal_storage_) {
      jacobianhOver_state_internal_storage_ =
          std::make_unique<std::vector<double>>(sizeY * sizeX);
    }
    jacobianhOver_state_view_ = std::make_shared<siconos::algebra::MapType>(
        jacobianhOver_state_internal_storage_->data(), sizeY, sizeX);
  }
  // if C is constant (following a call to setConstantC)
  // then
  // - no need to allocate internal storage
  // - map/view is already 'connected' to some external memory

  if (computeD_) {
    if (!jacobianhOver_lambda_internal_storage_) {
      jacobianhOver_lambda_internal_storage_ =
          std::make_unique<std::vector<double>>(sizeY * sizeY);
    }
    jacobianhOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
        jacobianhOver_lambda_internal_storage_->data(), sizeY, sizeY);
  }

  if (computeB_) {
    if (!jacobiangOver_lambda_internal_storage_) {
      jacobiangOver_lambda_internal_storage_ =
          std::make_unique<std::vector<double>>(sizeX * sizeY);
    }
    jacobiangOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
        jacobiangOver_lambda_internal_storage_->data(), sizeX, sizeY);
  }

  if (computeeVector_) {
    if (!eVector_internal_storage_) {
      eVector_internal_storage_ = std::make_unique<std::vector<double>>(sizeY);
    }
    eVector_view_ = std::make_shared<siconos::algebra::MapVectorType>(
        eVector_internal_storage_->data(), sizeY);
  }
  checkSize(inter);
}

void siconos::modeling::FirstOrderLinearR::checkSize(const Interaction& inter) const {
  // Check if various operators sizes are consistent.
  // Reference: interaction.

  if (jacobianhOver_state_view_) {
    assert(jacobianhOver_state_view_->rows() == inter.dimension());
    assert(jacobianhOver_state_view_->cols() == inter.getSizeOfDS());
  }
  if (jacobiangOver_lambda_view_) {
    assert(jacobiangOver_lambda_view_->rows() == inter.getSizeOfDS());
    assert(jacobiangOver_lambda_view_->cols() == inter.dimension());
  }
  if (jacobianhOver_lambda_view_) {
    assert(jacobianhOver_lambda_view_->rows() == inter.dimension());
    assert(jacobianhOver_lambda_view_->cols() == inter.dimension());
  }
  if (eVector_view_) {
    assert(eVector_view_->size() == inter.dimension());
  }
}

void siconos::modeling::FirstOrderLinearR::setConstantB(
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue) {
  /**  Warning: we can't check if newValue size is compliant with
   * the current relation/interaction. This will be done during
   * initialize / checkSize call!
   */

  jacobiangOver_lambda_internal_storage_ = nullptr;

  jacobiangOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantJacobiangOver_lambda_ = true;
  computeB_ = nullptr;
}

void siconos::modeling::FirstOrderLinearR::setComputeBFunction(
    const siconos::modeling::func_prototypes::FunctionS_M& func) {
  hasConstantJacobiangOver_lambda_ = false;
  computeB_ = func;
}

void siconos::modeling::FirstOrderLinearR::computeB(double time) {
  if (computeB_) computeB_(time, *jacobiangOver_lambda_view_);
}

void siconos::modeling::FirstOrderLinearR::setConstantC(
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue) {
  /**  Warning: we can't check if newValue size is compliant with
   * the current relation/interaction. This will be done during
   * initialize / checkSize call!
   */

  jacobianhOver_state_internal_storage_ = nullptr;

  jacobianhOver_state_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantJacobianhOver_state_ = true;
  computeC_ = nullptr;
}

void siconos::modeling::FirstOrderLinearR::setComputeCFunction(
    const siconos::modeling::func_prototypes::FunctionS_M& func) {
  hasConstantJacobianhOver_state_ = false;
  computeC_ = func;
}

void siconos::modeling::FirstOrderLinearR::computeC(double time) {
  if (computeC_) computeC_(time, *jacobianhOver_state_view_);
}

void siconos::modeling::FirstOrderLinearR::setConstantD(
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue) {
  /**  Warning: we can't check if newValue size is compliant with
   * the current relation/interaction. This will be done during
   * initialize / checkSize call!
   */

  jacobianhOver_lambda_internal_storage_ = nullptr;

  jacobianhOver_lambda_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantJacobianhOver_lambda_ = true;
  computeD_ = nullptr;
}

void siconos::modeling::FirstOrderLinearR::setComputeDFunction(
    const siconos::modeling::func_prototypes::FunctionS_M& func) {
  hasConstantJacobianhOver_lambda_ = false;
  computeD_ = func;
}

void siconos::modeling::FirstOrderLinearR::computeD(double time) {
  if (computeD_) computeD_(time, *jacobianhOver_lambda_view_);
}

void siconos::modeling::FirstOrderLinearR::setConstanteVector(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  eVector_internal_storage_ = nullptr;

  eVector_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
  haseVector_ = true;
  hasConstanteVector_ = true;
  computeeVector_ = nullptr;
}

void siconos::modeling::FirstOrderLinearR::setComputeeVectorFunction(
    const siconos::modeling::func_prototypes::FunctionS_V& fct) {
  haseVector_ = true;
  hasConstanteVector_ = false;
  computeeVector_ = fct;
}

void siconos::modeling::FirstOrderLinearR::computee(double time) {
  if (computeeVector_) computeeVector_(time, *eVector_view_);
}

// void siconos::modeling::FirstOrderLinearR::computeh(
//     double time, const siconos::algebra::BlockVector &x,
//     const siconos::algebra::SiconosVector &lambda, siconos::algebra::BlockVector &z,
//     siconos::algebra::SiconosVector &y) {
//   if (_C) {
//     computeC(time, z, *_C);
//     siconos::algebra::matrixBlockVector_prod(*_C, x, y, true);
//   } else
//     y.setZero();

//   if (_D) {
//     computeD(time, z, *_D);
//     y += *_D * lambda;
//   }
//   if (_F) {
//     computeF(time, z, *_F);
//     siconos::algebra::matrixBlockVector_prod(*_F, z, y, false);
//   }
//   if (_e) {
//     computee(time, z, *_e);
//     y += *_e;
//   }
// }

void siconos::modeling::FirstOrderLinearR::computeOutput(double time, Interaction& inter,
                                                         siconos::algebra::blocks::size_type) {
  DEBUG_BEGIN("siconos::modeling::FirstOrderLinearR::computeOutput \n");
  siconos::algebra::SiconosVector& y = *inter.y(0);
  siconos::algebra::SiconosVector& lambda = *inter.lambda(0);
  const auto& ds_vars = inter.read_dynamical_systems_variables();
  if (jacobianhOver_state_view_) {
    if (!hasConstantJacobianhOver_state_)  // C not constant
      computeC_(time, *jacobianhOver_state_view_);
    siconos::algebra::matrixBlockVector_prod(*jacobianhOver_state_view_,
                                             *ds_vars[FirstOrderR::Xxx], y, true);
  } else
    y.setZero();

  if (jacobianhOver_lambda_view_) {
    if (!hasConstantJacobianhOver_lambda_)  // D not constant
      computeD_(time, *jacobianhOver_lambda_view_);
    y += *jacobianhOver_lambda_view_ * lambda;
  }

  if (eVector_view_) {
    if (!hasConstanteVector_) computeeVector_(time, *eVector_view_);
    y += *eVector_view_;
  }
  DEBUG_END("siconos::modeling::FirstOrderLinearR::computeOutput \n");
}

// void siconos::modeling::FirstOrderLinearR::computeg(
//     double time, const siconos::algebra::SiconosVector &lambda,
//     siconos::algebra::BlockVector &z, siconos::algebra::BlockVector &r) {
//   computeB(time, z, *_B);
//   r += *_B * lambda;
// }

void siconos::modeling::FirstOrderLinearR::computeInput(double time, Interaction& inter,
                                                        siconos::algebra::blocks::size_type) {
  if (jacobiangOver_lambda_view_) {
    const auto& ds_vars = inter.read_dynamical_systems_variables();
    if (!hasConstantJacobiangOver_lambda_)  // B not constant
      computeB_(time, *jacobiangOver_lambda_view_);
    *ds_vars[FirstOrderR::Rrr] += *jacobiangOver_lambda_view_ * *inter.lambda(0);
  }
}

void siconos::modeling::FirstOrderLinearR::display() const {
  std::cout << " ===== Linear relation display ===== \n";
  std::cout << "| C \n";
  if (jacobianhOver_state_view_)
    std::cout << *jacobianhOver_state_view_ << "\n";
  else
    std::cout << "->nullptr\n";

  std::cout << "| D\n";
  if (jacobianhOver_lambda_view_)
    std::cout << *jacobianhOver_lambda_view_ << "\n";
  else
    std::cout << "->nullptr\n";

  std::cout << "| e\n";
  if (eVector_view_)
    std::cout << *eVector_view_ << "\n";
  else
    std::cout << "->nullptr\n";

  std::cout << "| B \n";
  if (jacobiangOver_lambda_view_)
    std::cout << *jacobiangOver_lambda_view_ << "\n";
  else
    std::cout << "->nullptr\n";
  std::cout << " ================================================== \n";
}
