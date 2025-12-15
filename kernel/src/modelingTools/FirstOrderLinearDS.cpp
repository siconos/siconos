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

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include <iostream>

// #include "siconos_debug.h"

// --- Constructors ---

siconos::modeling::FirstOrderLinearDS::FirstOrderLinearDS(
    Eigen::Ref<siconos::algebra::SiconosVector> x0,
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> A,
    Eigen::Ref<siconos::algebra::SiconosVector> b, siconos::algebra::AliasTag)
    : FirstOrderNonLinearDS(x0, siconos::algebra::alias_t) {
  setConstantA(A, siconos::algebra::alias_t);
  setConstantbVector(b, siconos::algebra::alias_t);
  isTimeInvariant_ = true;
}

siconos::modeling::FirstOrderLinearDS::FirstOrderLinearDS(
    const siconos::algebra::SiconosVector& x0, const siconos::algebra::SiconosDenseMatrix& A,
    const siconos::algebra::SiconosVector& b, siconos::algebra::CopyTag)
    : FirstOrderNonLinearDS(x0, siconos::algebra::copy_t) {
  setConstantA(A, siconos::algebra::copy_t);
  setConstantbVector(b, siconos::algebra::copy_t);
  isTimeInvariant_ = true;
}

// Copy constructor
siconos::modeling::FirstOrderLinearDS::FirstOrderLinearDS(const FirstOrderLinearDS& ds)
    : FirstOrderNonLinearDS(ds) {
  if (ds.hasbVector()) {
    setConstantbVector(ds.bVector(), siconos::algebra::copy_t);
    if (ds.computebVector_) {
      setComputebVectorFunction(ds.computebVector_);
    }
  }

  if (ds.hasJacobianfOver_x()) {
    hasJacobianfOver_x_ = true;
    computejacobianfOver_x_ = nullptr;
    if (ds.computeA_) {
      setComputeAFunction(ds.computeA_);
    } else
      hasConstantJacobianfOver_x_ = true;
  }
}

void siconos::modeling::FirstOrderLinearDS::setConstantA(
    const siconos::algebra::SiconosDenseMatrix& newValue, siconos::algebra::CopyTag) {
  setConstantJacobianfOver_x(newValue, siconos::algebra::copy_t);
  computeA_ = nullptr;
}

void siconos::modeling::FirstOrderLinearDS::setConstantA(
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue, siconos::algebra::AliasTag) {
  setConstantJacobianfOver_x(newValue, siconos::algebra::alias_t);
  computeA_ = nullptr;
}

void siconos::modeling::FirstOrderLinearDS::setComputeAFunction(
    const siconos::modeling::func_prototypes::FunctionS_M& func) {
  // Ensure that memory is properly allocated for jacobianfOver_x_
  if (!std::holds_alternative<siconos::algebra::OwnedDense>(jacobianfOver_x_storage_)) {
    jacobianfOver_x_storage_ =
        std::make_unique<siconos::algebra::SiconosDenseMatrix>(x_size_, x_size_);
  }
  hasJacobianfOver_x_ = true;
  hasConstantJacobianfOver_x_ = false;
  computejacobianfOver_x_ = nullptr;  // just in case ...
  computeA_ = func;
  isTimeInvariant_ = false;
}

void siconos::modeling::FirstOrderLinearDS::computeA(double time) {
  if (computeA_) {
    use_jacobianfOver_x([&](auto& M) { computeA_(time, M); });
    is_jacobianRhsOver_x_uptodate_ = false;
  }
}

void siconos::modeling::FirstOrderLinearDS::setConstantbVector(
    const siconos::algebra::SiconosVector& newValue, siconos::algebra::CopyTag tag) {
  if (newValue.size() != x_size_)
    throw std::invalid_argument("setConstantbVector(copy): input vector has wrong size");

  // Deep copy into Owned storage
  bVector_storage_ = std::make_unique<siconos::algebra::SiconosVector>(newValue);
  hasbVector_ = true;
  hasConstantbVector_ = true;
  computebVector_ = nullptr;
}

void siconos::modeling::FirstOrderLinearDS::setConstantbVector(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue, siconos::algebra::AliasTag tag) {
  if (newValue.size() != x_size_)
    throw std::invalid_argument("setConstantbVector(alias): input vector has wrong size");

  // Aliasing: shared_ptr to a Map (view over external memory)
  bVector_storage_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
  hasbVector_ = true;
  hasConstantbVector_ = true;
  computebVector_ = nullptr;
}

void siconos::modeling::FirstOrderLinearDS::setComputebVectorFunction(
    const siconos::modeling::func_prototypes::FunctionS_V& func) {
  // Ensure that memory is properly allocated for bVector
  if (!std::holds_alternative<siconos::algebra::OwnedDenseVector>(bVector_storage_)) {
    bVector_storage_ = std::make_unique<siconos::algebra::SiconosVector>(x_size_);
  }
  hasbVector_ = true;
  computebVector_ = func;
  isTimeInvariant_ = false;
  hasConstantbVector_ = false;
}

void siconos::modeling::FirstOrderLinearDS::clearbVector() {
  computebVector_ = nullptr;  // Usefull for MatrixIntegrator and control stuff
  hasbVector_ = false;
  bVector_storage_ = std::monostate{};
  hasConstantbVector_ = false;
}

void siconos::modeling::FirstOrderLinearDS::computeb(double time) {
  if (computebVector_) use_bVector([&](auto& fv) { computebVector_(time, fv); });
}

void siconos::modeling::FirstOrderLinearDS::updatePlugins(double time) {
  if (isTimeInvariant_) return;
  if (computeMMatrix_) {
    use_MMatrix([&](auto& M) { computeMMatrix_(time, M); });
    hasLU_M_ = false;  // LUM needs to be updated.
  }
  if (computebVector_) use_bVector([&](auto& fv) { computebVector_(time, fv); });
  if (computeA_) use_jacobianfOver_x([&](auto& M) { computeA_(time, M); });
}

/*This function is called only by LsodarOSI and eventDriven*/
void siconos::modeling::FirstOrderLinearDS::computeRhs(double time) {
  // second argument is useless at the time - Used in derived classes
  // compute A=jacobianfx
  *state_x_[1] = *rVector_;  // Warning: r update is done in Interactions/Relations

  if (hasJacobianfOver_x()) {  // if A
    use_jacobianfOver_x([&](auto& mat) {
      if (computeA_) computeA_(time, mat);
      *(state_x_[1]) += mat * *(state_x_[0]);
    });
  }

  if (hasbVector()) {
    use_bVector([&](auto& bv) {
      if (computebVector_) computebVector_(time, bv);
      *(state_x_[1]) += bv;
    });
  }

  if (hasMMatrix()) {
    computeMMatrix(time);
    // allocate invM if required
    if (!hasLU_M_) {
      LU_M_ = std::make_shared<siconos::algebra::SiconosDenseLUMatrix>(MMatrix());
      hasLU_M_ = true;
    }

    *(state_x_[1]) = LU_M_->solve(*(state_x_[1]));
  }
}

void siconos::modeling::FirstOrderLinearDS::computeJacobianRhsOver_x(double time) {
  if (isTimeInvariant_ && !isFirstCall_) return;

  if (!hasJacobianfOver_x_) return;
  if (computeA_)  // if plugged, update
    use_jacobianfOver_x([&](auto& M) { computeA_(time, M); });

  // If M is defined ...
  if (hasMMatrix()) {
    use_MMatrix([&](auto& M) {
      if (computeMMatrix_) {
        // and plugged: update
        computeMMatrix_(time, M);
        hasLU_M_ = false;  // M has changed, LUM needs to be updated.
      }
      if (!hasLU_M_) {
        // When? If M has been updated or if it's the first call
        is_jacobianRhsOver_x_uptodate_ = false;
        LU_M_ = std::make_shared<siconos::algebra::SiconosDenseLUMatrix>(M);
        hasLU_M_ = true;
      }
    });
    // solve M*jacobianXRhS = jacobianfx
    if (!is_jacobianRhsOver_x_uptodate_) {
      // Solve M-1.jacobianfOver_x_view_ in temp result matrix
      use_jacobianfOver_x([&](const auto& jac) {
        siconos::algebra::SiconosDenseMatrix result = LU_M_->solve(jac);
        // and keep it as a flattened vector
        std::copy(result.data(), result.data() + x_size_ * x_size_,
                  jacobianRhsOver_x_.begin());
      });
      is_jacobianRhsOver_x_uptodate_ = true;
    }
  } else {  // No M. Just copy jacobianfOver_x_view_ into jacobianRhsOver_x_ (flat)
    if (!is_jacobianRhsOver_x_uptodate_) {
      use_jacobianfOver_x([&](const auto& jac) {
        std::copy(jac.data(), jac.data() + x_size_ * x_size_, jacobianRhsOver_x_.begin());
      });
      is_jacobianRhsOver_x_uptodate_ = true;
    }
  }
}

void siconos::modeling::FirstOrderLinearDS::display(bool brief) const {
  std::cout << "=== Linear system display, " << number_ << std::endl;
  std::cout << "- dimension : " << x_size_ << std::endl;
  std::cout << "- state :\n" << *state_x_[0] << "\n";
  std::cout << "- initial state : \n";
  use_x0([&](const auto& v) { siconos::algebra::print(v); });

  std::cout << "- M matrix: \n";
  if (hasMMatrix())
    use_MMatrix([&](const auto& M) { siconos::algebra::print(M); });
  else
    std::cout << "-> nullptr\n";
  std::cout << "- A matrix: \n";
  if (hasJacobianfOver_x())
    use_jacobianfOver_x([&](const auto& M) { siconos::algebra::print(M); });
  else
    std::cout << "-> nullptr\n";
  std::cout << "- b vector: \n";
  if (hasbVector())
    use_bVector([&](const auto& v) { siconos::algebra::print(v); });

  else
    std::cout << "-> nullptr\n";
  std::cout << "=============================" << std::endl;
}
