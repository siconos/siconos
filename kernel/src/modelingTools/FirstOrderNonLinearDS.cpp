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

#include "FirstOrderNonLinearDS.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include <iostream>
#include <memory>

#include "siconos_debug.h"

// From a minimum set of data
siconos::modeling::FirstOrderNonLinearDS::FirstOrderNonLinearDS(
    Eigen::Ref<siconos::algebra::SiconosVector> initial_state)
    : DynamicalSystem(initial_state.size()) {
  // Memory allocation only for required parts of the DS:
  // state (initial and current). All other operators are optional and
  // allocated with 'set'-like methods.
  assert(x_size_ > 0 && "dynamical system dimension should be greater than 0.");
  // Set initial conditions
  x0_view_ = std::make_shared<siconos::algebra::MapVectorType>(initial_state.data(),
                                                               x_size_);  // Shared memory!

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.
  // _x.resize(2); done in base class constructor.
  state_x_[0] = std::make_shared<siconos::algebra::SiconosVector>(*x0_view_);
  state_x_[1] = std::make_shared<siconos::algebra::SiconosVector>(x_size_);
  state_x_[1]->setZero();

  rVector_ = std::make_shared<siconos::algebra::SiconosVector>(x_size_);
  // FP: move this to initializeNonSmoothInput?

  rVector_->setZero();
  // dot x = r
}

// Copy constructor
siconos::modeling::FirstOrderNonLinearDS::FirstOrderNonLinearDS(
    const FirstOrderNonLinearDS &ds)
    : DynamicalSystem(ds) {
  if (ds.MMatrix_view_) {
    MMatrix_internal_storage_ =
        std::make_unique<std::vector<double>>(ds.MMatrix_view_->size());
    MMatrix_view_ = std::make_shared<siconos::algebra::MapType>(
        MMatrix_internal_storage_->data(), x_size_, x_size_);
    *MMatrix_view_ = ds.MMatrix();  // Copy
    hasMMatrix_ = true;
    if (ds.computeMMatrix_) {
      setComputeMMatrixFunction(ds.computeMMatrix_);
    } else
      hasConstantMMatrix_ = true;
  }

  if (ds.fVector_view_) {
    fVector_internal_storage_ =
        std::make_unique<std::vector<double>>(ds.fVector_view_->size());
    fVector_view_ = std::make_shared<siconos::algebra::MapVectorType>(
        fVector_internal_storage_->data(), fVector_internal_storage_->size());
    *fVector_view_ = ds.fVector();
    if (ds.computefVector_) {
      setComputefVectorFunction(ds.computefVector_);
    } else
      hasConstantfVector_ = true;
  }

  if (ds.jacobianfOver_x_view_) {
    jacobianfOver_x_internal_storage_ =
        std::make_unique<std::vector<double>>(ds.jacobianfOver_x_view_->size());
    jacobianfOver_x_view_ = std::make_shared<siconos::algebra::MapType>(
        jacobianfOver_x_internal_storage_->data(), x_size_, x_size_);
    *jacobianfOver_x_view_ = ds.jacobianfOver_x();  // Copy
    if (ds.computejacobianfOver_x_) {
      setComputeJacobianfOver_xFunction(ds.computejacobianfOver_x_);
    } else
      hasConstantJacobianfOver_x_ = true;
  }

  if (ds.LU_M()) {
    LU_M_ = std::make_shared<siconos::algebra::SiconosLUMatrix>(*MMatrix_view_);
    hasLU_M_ = true;
  }
  // Memory stuff to me moved to graph/osi
  if (ds.fbuffer_)
    fbuffer_ = std::make_shared<siconos::algebra::SiconosVector>(*(ds.fbuffer_));
  rMemory_ = ds.rMemory_;
}

void siconos::modeling::FirstOrderNonLinearDS::initRhs(double time) {
  computeRhs(time);
  jacobianRhsOver_x_.resize(x_size_ * x_size_);
  if (jacobianfOver_x_view_) {
    computeJacobianRhsOver_x(time);
  } else  // else no resize, jacobian of rhs is equal to 0.
    jacobianRhsOver_x_.setZero();

  is_jacobianRhsOver_x_uptodate_ = true;
}

void siconos::modeling::FirstOrderNonLinearDS::computeRhs(double time) {
  // second argument is useless at the time - Used in derived classes

  // compute rhs = M-1*( f + r ).

  *state_x_[1] = *rVector_;  // Warning: r update is done in Interactions/Relations

  if (fVector_view_) {
    if (computefVector_) computefVector_(*state_x_[0], time, *fVector_view_);
    *(state_x_[1]) += *fVector_view_;
  }

  if (MMatrix_view_) {
    computeMMatrix(time);
    if (!hasLU_M_) {
      // allocate invM at the first call of the present function
#ifdef SICONOS_SPARSE
      LU_M_ = std::make_shared<siconos::algebra::SparseLUMatrix>();
      LU_M_->analysePattern(*MMatrix_view_);
      LU_M_->factorize(*MMatrix_view_);
#else
      LU_M_ = std::make_shared<siconos::algebra::DenseLUMatrix>(*MMatrix_view_);
#endif
      *(state_x_[1]) = LU_M_->solve(*(state_x_[1]));
      hasLU_M_ = true;
    }
  }
}

void siconos::modeling::FirstOrderNonLinearDS::computeJacobianRhsOver_x(double time) {
  // second argument is useless here but used in derived classes

  // compute jacobian of rhs according to x, = M-1(jacobianfx + jacobianX(T.u))
  // At the time, second term is set to zero.

  // Cases:
  // 1- M not defined then jacobianRhsOver_x_ contains a flat version of jacobianfOver_x_view_
  // 2- M is defined then jacobianRhsOver_x_ contains a flat version of
  // M-1.jacobianfOver_x_view_

  if (!hasJacobianfOver_x_) return;

  if (computejacobianfOver_x_)  // if plugged, update
    computejacobianfOver_x_(*state_x_[0], time, *jacobianfOver_x_view_);

  // If M is defined ...
  if (MMatrix_view_) {
    computeMMatrix(time);
    if (!hasLU_M_) {
      // When? If M has been updated or if it's the first call
      is_jacobianRhsOver_x_uptodate_ = false;
#ifdef SICONOS_SPARSE
      LU_M_ = std::make_shared<siconos::algebra::SparseLUMatrix>();
      LU_M_->analysePattern(*MMatrix_view_);
      LU_M_->factorize(*MMatrix_view_);
#else
      LU_M_ = std::make_shared<siconos::algebra::DenseLUMatrix>(*MMatrix_view_);
#endif
      hasLU_M_ = true;
    }
    // solve M*jacobianXRhS = jacobianfx
    if (!is_jacobianRhsOver_x_uptodate_) {
      // Solve M-1.jacobianfOver_x_view_ in temp result matrix
      siconos::algebra::SiconosMatrix result = LU_M_->solve(*jacobianfOver_x_view_);
      // and keep it as a flattened vector
      jacobianRhsOver_x_.noalias() =
          Eigen::Map<const Eigen::VectorXd>(result.data(), jacobianRhsOver_x_.size());
      is_jacobianRhsOver_x_uptodate_ = true;
    }
  } else {
    if (!is_jacobianRhsOver_x_uptodate_) {
      // No M. Just copy jacobianfOver_x_view_ into jacobianRhsOver_x_ (flat)
      jacobianRhsOver_x_.noalias() = Eigen::Map<const Eigen::VectorXd>(
          jacobianfOver_x_view_->data(), jacobianRhsOver_x_.size());
      is_jacobianRhsOver_x_uptodate_ = true;
    }
  }
}

void siconos::modeling::FirstOrderNonLinearDS::setConstantMMatrix(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for MMatrix
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  MMatrix_internal_storage_ = nullptr;

  MMatrix_view_ =
      std::make_shared<siconos::algebra::MapType>(newValue.data(), x_size_, x_size_);
  hasMMatrix_ = true;
  hasConstantMMatrix_ = true;
  computeMMatrix_ = nullptr;
  is_jacobianRhsOver_x_uptodate_ = false;
}

void siconos::modeling::FirstOrderNonLinearDS::setComputeMMatrixFunction(
    const siconos::modeling::func_prototypes::FunctionS_M &new_func) {
  // Ensure that memory is properly allocated for MMatrix_
  if (!MMatrix_internal_storage_) {
    MMatrix_internal_storage_ = std::make_unique<std::vector<double>>(x_size_ * x_size_);
  }
  MMatrix_view_ = std::make_shared<siconos::algebra::MapType>(
      MMatrix_internal_storage_->data(), x_size_, x_size_);
  MMatrix_view_->setZero();
  hasMMatrix_ = true;
  hasConstantMMatrix_ = false;
  computeMMatrix_ = new_func;
  isTimeInvariant_ = false;
}

void siconos::modeling::FirstOrderNonLinearDS::computeMMatrix(double time) {
  if (computeMMatrix_) {
    computeMMatrix_(time, *MMatrix_view_);
    hasLU_M_ = false;  // M has changed, LUM needs to be updated.
    is_jacobianRhsOver_x_uptodate_ = false;
  }
}

void siconos::modeling::FirstOrderNonLinearDS::setConstantfVector(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  fVector_internal_storage_ = nullptr;
  assert(newValue.size() == x_size_);
  fVector_view_ = std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), x_size_);
  hasfVector_ = true;
  hasConstantfVector_ = true;
  computefVector_ = nullptr;
}

void siconos::modeling::FirstOrderNonLinearDS::setComputefVectorFunction(
    const siconos::modeling::func_prototypes::FunctionVS_V &func) {
  // Ensure that memory is properly allocated for fVector
  if (!fVector_internal_storage_) {
    fVector_internal_storage_ = std::make_unique<std::vector<double>>(x_size_);
  }
  fVector_view_ = std::make_shared<siconos::algebra::MapVectorType>(
      fVector_internal_storage_->data(), x_size_);
  hasfVector_ = true;
  computefVector_ = func;
  isTimeInvariant_ = false;
}

void siconos::modeling::FirstOrderNonLinearDS::computefVector(
    const Eigen::Ref<siconos::algebra::SiconosVector> &state, double time) {
  if (computefVector_) {
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computefVector_(state, time, *fVector_view_);
  }
}

void siconos::modeling::FirstOrderNonLinearDS::setConstantJacobianfOver_x(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for jacobianfOver_x
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  jacobianfOver_x_internal_storage_ = nullptr;

  jacobianfOver_x_view_ =
      std::make_shared<siconos::algebra::MapType>(newValue.data(), x_size_, x_size_);
  hasJacobianfOver_x_ = true;
  hasConstantJacobianfOver_x_ = true;
  computejacobianfOver_x_ = nullptr;
  is_jacobianRhsOver_x_uptodate_ = false;
}

void siconos::modeling::FirstOrderNonLinearDS::setComputeJacobianfOver_xFunction(
    const siconos::modeling::func_prototypes::FunctionVS_M &new_func) {
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
  computejacobianfOver_x_ = new_func;
  isTimeInvariant_ = false;
}

void siconos::modeling::FirstOrderNonLinearDS::computeJacobianfOver_x(
    const Eigen::Ref<siconos::algebra::SiconosVector> &state, double time) {
  if (computejacobianfOver_x_) {
    computejacobianfOver_x_(state, time, *jacobianfOver_x_view_);
    is_jacobianRhsOver_x_uptodate_ = false;
  }
}

void siconos::modeling::FirstOrderNonLinearDS::updatePlugins(double time) {
  computeMMatrix(time);
  if (computefVector_) computefVector_(*state_x_[0], time, *fVector_view_);
  if (computejacobianfOver_x_)
    computejacobianfOver_x_(*state_x_[0], time, *jacobianfOver_x_view_);
}

void siconos::modeling::FirstOrderNonLinearDS::initializeNonSmoothInput(unsigned int level) {
  /**\warning V.A. rVector_ should be initialized here and not in  the constructor
   * The level should also be used if we need more that one rVector_
   */
  if (!rVector_) rVector_ = std::make_shared<siconos::algebra::SiconosVector>(x_size_);
}

void siconos::modeling::FirstOrderNonLinearDS::resetToInitialState() {
  *(state_x_[0]) = *x0_view_;
}

void siconos::modeling::FirstOrderNonLinearDS::resetInitialStateToValue(double val) {
  x0_view_->setConstant(val);
  // FP: Should we reset state vector to x0 ? I'd say yes ...
  resetToInitialState();

  // FP: this function is only used in MatrixIntegrator for ZOH and friends ...
}

// ===== MEMORY MANAGEMENT FUNCTIONS =====

void siconos::modeling::FirstOrderNonLinearDS::initMemory(unsigned int steps) {
  DynamicalSystem::initMemory(steps);

  if (fVector_view_ && !fbuffer_)
    fbuffer_ = std::make_shared<siconos::algebra::SiconosVector>(x_size_);

  if (steps == 0)
    std::cout << "Warning : siconos::modeling::FirstOrderNonLinearDS::initMemory with size "
                 "equal to zero"
              << std::endl;
  else
    rMemory_.setMemorySize(steps, x_size_);
}

void siconos::modeling::FirstOrderNonLinearDS::swapInMemory() {
  DEBUG_BEGIN("void siconos::modeling::FirstOrderNonLinearDS::swapInMemory()\n");
  xMemory_.swap(*state_x_[0]);
  rMemory_.swap(*rVector_);
  if (fVector_view_) {
    assert(fbuffer_);
    *fbuffer_ = *fVector_view_;
  }
  DEBUG_EXPR(siconos::algebra::print(xMemory_));
  DEBUG_END("void siconos::modeling::FirstOrderNonLinearDS::swapInMemory()\n");
}

// ===== MISCELLANEOUS ====

void siconos::modeling::FirstOrderNonLinearDS::display(bool brief) const {
  std::cout << " =====> First Order Non Linear DS (number: " << number_ << ").\n";
  std::cout << "- dimension : " << x_size_ << std::endl;
  std::cout << "- state :\n" << *state_x_[0] << "\n";
  std::cout << "- initial state : \n" << *x0_view_ << "\n";
  std::cout << "- M matrix: \n";
  if (MMatrix_view_)
    std::cout << *MMatrix_view_ << "\n";
  else
    std::cout << "-> nullptr\n";
  std::cout << " ============================================\n";
}

void siconos::modeling::FirstOrderNonLinearDS::resetAllNonSmoothParts() {
  rVector_->setZero();
}

void siconos::modeling::FirstOrderNonLinearDS::resetNonSmoothPart(unsigned int level) {
  // V.A. 28/05/2012:  for the moment various level are not used for First Order systems
  // assert(0);
  rVector_->setZero();
}