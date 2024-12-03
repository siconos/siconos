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
#include "FirstOrderLinearDS.hpp"

#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // For prod matrix-vectors
#include "SiconosVector.hpp"
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include <iostream>

#include "siconos_debug.h"

// --- Constructors ---

siconos::modeling::FirstOrderLinearDS::FirstOrderLinearDS(
    Eigen::Ref<siconos::algebra::SiconosVector> x0,
    Eigen::Ref<siconos::algebra::SiconosMatrix> A,
    Eigen::Ref<siconos::algebra::SiconosVector> b)
    : FirstOrderNonLinearDS(x0) {
  setConstantA(A);
  setConstantbVector(b);
  isTimeInvariant_ = true;
}

// Copy constructor
siconos::modeling::FirstOrderLinearDS::FirstOrderLinearDS(const FirstOrderLinearDS &ds)
    : FirstOrderNonLinearDS(ds) {
  if (ds.bVector_view_) {
    bVector_internal_storage_ =
        std::make_unique<std::vector<double>>(ds.bVector_view_->size());
    bVector_view_ = std::make_shared<siconos::algebra::MapVectorType>(
        ds.bVector_view_->data(), ds.bVector_view_->size());
    if (ds.computebVector_) {
      computebVector_ = ds.computebVector_;
    } else
      hasConstantbVector_ = true;
  }
}

void siconos::modeling::FirstOrderLinearDS::setConstantA(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  setConstantJacobianfOver_x(newValue);
  computeA_ = nullptr;
}

void siconos::modeling::FirstOrderLinearDS::setComputeAFunction(
    const siconos::modeling::func_prototypes::FunctionS_M &func) {
  // Ensure that memory is properly allocated for jacobianfOver_x_
  if (!jacobianfOver_x_internal_storage_) {
    jacobianfOver_x_internal_storage_ =
        std::make_unique<std::vector<double>>(x_size_ * x_size_);
  }
  jacobianfOver_x_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianfOver_x_internal_storage_->data(), x_size_, x_size_);
  jacobianfOver_x_view_->setZero();

  hasJacobianfOver_x_ = true;
  hasConstantJacobianfOver_x_ = false;
  computejacobianfOver_x_ = nullptr;  // just in case ...
  computeA_ = func;
  isTimeInvariant_ = false;
}

void siconos::modeling::FirstOrderLinearDS::computeA(double time) {
  if (computeA_) {
    computeA_(time, *jacobianfOver_x_view_);
  }
}

void siconos::modeling::FirstOrderLinearDS::setConstantbVector(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  bVector_internal_storage_ = nullptr;
  assert(newValue.size() == x_size_);
  bVector_view_ = std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), x_size_);
  hasbVector_ = true;
  hasConstantbVector_ = true;
  computebVector_ = nullptr;
}

void siconos::modeling::FirstOrderLinearDS::setComputebVectorFunction(
    const siconos::modeling::func_prototypes::FunctionS_V &func) {
  // Ensure that memory is properly allocated for bVector
  if (!bVector_internal_storage_) {
    bVector_internal_storage_ = std::make_unique<std::vector<double>>(x_size_);
  }
  bVector_view_ = std::make_shared<siconos::algebra::MapVectorType>(
      bVector_internal_storage_->data(), x_size_);
  hasbVector_ = true;
  computebVector_ = func;
  isTimeInvariant_ = false;
}

void siconos::modeling::FirstOrderLinearDS::computeb(double time) {
  if (computebVector_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computebVector_(time, *bVector_view_);
}

void siconos::modeling::FirstOrderLinearDS::updatePlugins(double time) {
  if (isTimeInvariant_) return;
  if (computeMMatrix_) {
    computeMMatrix_(time, *MMatrix_view_);
    hasLU_M_ = false;  // LUM needs to be updated.
  }
  if (computebVector_) computebVector_(time, *bVector_view_);
  if (computeA_) computeA_(time, *jacobianfOver_x_view_);
}

/*This function is called only by LsodarOSI and eventDriven*/
void siconos::modeling::FirstOrderLinearDS::computeRhs(double time) {
  // second argument is useless at the time - Used in derived classes
  // compute A=jacobianfx

  *state_x_[1] = *rVector_;  // Warning: r update is done in Interactions/Relations

  if (jacobianfOver_x_view_) {  // if A
    if (computeA_) computeA_(time, *jacobianfOver_x_view_);
    *(state_x_[1]) += *jacobianfOver_x_view_ * *(state_x_[0]);
  }

  if (bVector_view_) {
    if (computebVector_) computebVector_(time, *bVector_view_);
    *(state_x_[1]) += *bVector_view_;
  }

  if (MMatrix_view_) {
    if (computeMMatrix_) {
      computeMMatrix_(time, *MMatrix_view_);
      hasLU_M_ = false;  // LUM needs to be updated.
    }
    // allocate invM if required
    if (!hasLU_M_) {
      LU_M_ =
          std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(*MMatrix_view_);
      hasLU_M_ = true;
    }
    *(state_x_[1]) = LU_M_->solve(*(state_x_[1]));
  }
}

void siconos::modeling::FirstOrderLinearDS::computeJacobianRhsOver_x(double time) {
  if (isTimeInvariant_ && !isFirstCall_) return;
  isFirstCall_ = false;
  if (jacobianfOver_x_view_) {
    if (computeA_) computeA_(time, *jacobianfOver_x_view_);
    if (MMatrix_view_) {
      if (hasConstantJacobianfOver_x_)  // else memory is shared between jacobianRhsOver_x_
                                        // and jacobianfOver_x_view_
        *jacobianRhsOver_x_->block(0, 0) = *jacobianfOver_x_view_;

      if (computeMMatrix_) {
        computeMMatrix_(time, *MMatrix_view_);
        hasLU_M_ = false;  // LUM needs to be updated.
      }
      if (!hasLU_M_) {
        LU_M_ = std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(
            *MMatrix_view_);
        hasLU_M_ = true;
      }
      *(jacobianRhsOver_x_->block(0, 0)) = LU_M_->solve(*(jacobianRhsOver_x_->block(0, 0)));
    }
  }
  // else 0
}

void siconos::modeling::FirstOrderLinearDS::display(bool brief) const {
  std::cout << "=== Linear system display, " << number_ << std::endl;
  std::cout << "- dimension : " << x_size_ << std::endl;
  std::cout << "- state :\n" << state_x_[0] << "\n";
  std::cout << "- initial state : \n" << *x0_view_ << "\n";
  std::cout << "- M matrix: \n";
  if (MMatrix_view_)
    std::cout << *MMatrix_view_ << "\n";
  else
    std::cout << "-> nullptr\n";
  std::cout << "- A matrix: \n";
  if (jacobianfOver_x_view_)
    std::cout << *jacobianfOver_x_view_ << "\n";
  else
    std::cout << "-> nullptr\n";
  std::cout << "- b vector: \n";
  if (bVector_view_)
    std::cout << *bVector_view_ << "\n";
  else
    std::cout << "-> nullptr\n";
  std::cout << "=============================" << std::endl;
}
