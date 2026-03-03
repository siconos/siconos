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
#include "LagrangianSparseDS.hpp"

#include <iostream>
#include <memory>

#include "LagrangianR.hpp"
#include "SiconosConst.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"
#include "Tools.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES

#include "siconos_debug.h"

// Build from initial state only
siconos::modeling::LagrangianSparseDS::LagrangianSparseDS(
    Eigen::Ref<siconos::algebra::SiconosVector> qinit,
    Eigen::Ref<siconos::algebra::SiconosVector> vinit, siconos::algebra::AliasTag)
    : SecondOrderDS(2 * qinit.size(), vinit.size()) {
  assert(ndof_ > 0 && "lagrangian dynamical system dimension should be greater than 0.");

  // Set initial conditions

  q0_storage_ = std::make_shared<siconos::algebra::MapVectorType>(qinit.data(), ndof_);
  velocity0_storage_ = std::make_shared<siconos::algebra::MapVectorType>(vinit.data(), ndof_);

  // -- Memory allocation for vector and matrix members --
  state_q_[0] = std::make_shared<siconos::algebra::SiconosVector>(q0());
  state_q_[1] = std::make_shared<siconos::algebra::SiconosVector>(velocity0());

  /** \todo lazy Memory allocation */
  p_[1] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  p_[1]->setZero();
}

siconos::modeling::LagrangianSparseDS::LagrangianSparseDS(
    const siconos::algebra::SiconosVector& qinit, const siconos::algebra::SiconosVector& vinit,
    siconos::algebra::CopyTag)
    : SecondOrderDS(2 * qinit.size(), vinit.size()) {
  assert(ndof_ > 0 && "lagrangian dynamical system dimension should be greater than 0.");

  q0_storage_ = std::make_unique<siconos::algebra::SiconosVector>(qinit);
  velocity0_storage_ = std::make_unique<siconos::algebra::SiconosVector>(vinit);

  // -- Memory allocation for vector and matrix members --
  state_q_[0] = std::make_shared<siconos::algebra::SiconosVector>(q0());
  state_q_[1] = std::make_shared<siconos::algebra::SiconosVector>(velocity0());

  /** \todo lazy Memory allocation */
  p_[1] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  p_[1]->setZero();
}

void siconos::modeling::LagrangianSparseDS::initializeNonSmoothInput(
    siconos::algebra::blocks::size_type level) {
  if (!p_[level]) {
    p_[level] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
    p_[level]->setZero();
  }
}

void siconos::modeling::LagrangianSparseDS::resetToInitialState() {
  *(state_q_[0]) = q0();
  *(state_q_[1]) = velocity0();
}

void siconos::modeling::LagrangianSparseDS::initMemoryForGeneralizedCoordinates(
    siconos::algebra::blocks::size_type level) {
  assert(level > 1);
  if (!state_q_[level])
    state_q_[level] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
}

void siconos::modeling::LagrangianSparseDS::init_lu_mass() {
  //  if (mass_mat_ && !hasConstantMass_) {
  if (!hasMass_) return;
  if (!hasConstantMass_) computeMass(*state_q_[0]);
  if (!hasLUMass_) {
    // LU factorization
    use_mass([&](const auto& M) {
      LUMass_ = std::make_shared<siconos::algebra::SiconosSparseLUMatrix>(M);
    });
    hasLUMass_ = true;
  }
}

void siconos::modeling::LagrangianSparseDS::initRhs(double time) {
  DEBUG_BEGIN("siconos::modeling::LagrangianSparseDS::initRhs(double time)\n");
  // dim
  x_size_ = 2 * ndof_;

  // All links between DS and LagrangianSparseDS class members are pointer links,
  // which means that no useless memory is allocated when connection is
  // established. One exception: zero and identity matrices, used to filled in M
  // and jacobianfx.

  // Initial conditions and state

  // WARNING : this function is supposed to be called
  // by the OneStepIntegrator, and maybe several times for the same DS
  // if the system is involved in more than one interaction. So, we must check
  // if p2 and q2 already exist to be sure that ds_vars won't be lost.
  x0_storage_ = std::make_unique<siconos::algebra::SiconosVector>(x_size_);

  use_x0([&](auto& xinit) {
    xinit.head(ndof_) = q0();  // COPY !
    xinit.tail(ndof_) = velocity0();
  });

  state_x_[0] = std::make_shared<siconos::algebra::SiconosVector>(x_size_);
  *(state_x_[0]) << *state_q_[0], *state_q_[1];

  if (!state_q_[2]) {
    state_q_[2] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
    state_q_[2]->setZero();
  }
  state_x_[1] = std::make_shared<siconos::algebra::SiconosVector>(x_size_);
  *(state_x_[1]) << *state_q_[1], *state_q_[2];

  // Everything concerning rhs and its jacobian is handled in initRhs and
  // computeXXX related functions.

  if (!p_[2]) {
    p_[2] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
    p_[2]->setZero();
  }

  if (fint_ || hasFext_ || fgyr_) {
    if (!totalForces_) totalForces_ = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  }

  // Compute mass and LU factorization
  if (hasMass_) {
    computeMass(*state_q_[0]);
    // LU factorization
    use_mass([&](const auto& M) {
      LUMass_ = std::make_shared<siconos::algebra::SiconosSparseLUMatrix>(M);
    });
    hasLUMass_ = true;
  }

  computeRhs(time);
  // -- Jacobian Rhs over x --

  // The jacobian is saved in a flattened version, as a vector
  // which contains the concatenation of the columns of the real matrix
  jacobianRhsOver_x_.resize(x_size_ * x_size_);
  jacobianRhsOver_x_.setZero();

  // A lambda to compute the position in jacobian vector from the
  // coordinates in the "real" jacobian matrix.
  auto idx = [&](int row, int col) { return col * (2 * ndof_) + row; };
  // Upper-left block of the jacobian=0. Nothing to be done
  // Upper-right block = identity
  for (Eigen::Index i = 0; i < ndof_; ++i) jacobianRhsOver_x_(idx(i, ndof_ + i)) = 1.0;

  if (!jacobianTotalForcesOver_q_ && !jacobianTotalForcesOver_velocity_) {
    is_jacobianRhsOver_x_uptodate_ = true;

    return;  // The jacobian is required only if at least one of the two jacobians above
             // exists.
  }

  // - Fill parts corresponding to the jacobians of total forces -
  // mass and lu_mass are up to date since we have already called init_lu_mass
  // in computeRhs()

  if (hasMass()) {
    // In that case, we need to compute inv(Mass).jacobian_v or _q Forces
    siconos::algebra::SiconosSparseMatrix jacq;  // to save M^-1 . jac totalforces/q
    siconos::algebra::SiconosSparseMatrix jacv;  // to save M^-1 . jac totalforces/v

    if (hasJacobianTotalForcesOver_q()) {
      // Update if required
      computeJacobianTotalForcesOver_q(*state_q_[1], *state_q_[0], time);
      // Solve MjacobianX(1,0) =
      // jacobianFL[0]
      jacq = LUMass_->solve(*jacobianTotalForcesOver_q_);
    }
    if (hasJacobianTotalForcesOver_velocity()) {
      // Update if required
      computeJacobianTotalForcesOver_velocity(*state_q_[1], *state_q_[0], time);
      // Solve MjacobianX(1,1) =
      // jacobianFL[1]
      jacv = LUMass_->solve(*jacobianTotalForcesOver_velocity_);
    }
    // Now fill in jacobianRhsOver_x_
    // Bottom-left block (jacobian / q)
    if (hasJacobianTotalForcesOver_q()) {
      for (Eigen::Index k = 0; k < jacq.outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(jacq, k); it; ++it) {
          double value = it.value();
          jacobianRhsOver_x_(idx(ndof_ + it.row(), it.col())) = value;
        }
      }
    }
    // Bottom-right block (jacobian / vel)
    if (hasJacobianTotalForcesOver_velocity()) {
      for (Eigen::Index k = 0; k < jacv.outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(jacv, k); it; ++it) {
          double value = it.value();
          jacobianRhsOver_x_(idx(ndof_ + it.row(), ndof_ + it.col())) = value;
        }
      }
    }
  } else  // No mass
  {
    if (hasJacobianTotalForcesOver_q()) {
      // Update if required
      computeJacobianTotalForcesOver_q(*state_q_[1], *state_q_[0], time);
    }
    if (hasJacobianTotalForcesOver_velocity()) {
      // Update if required
      computeJacobianTotalForcesOver_velocity(*state_q_[1], *state_q_[0], time);
    }
    // Now fill in jacobianRhsOver_x_

    // Bottom-left block (jacobian / q)
    if (hasJacobianTotalForcesOver_q()) {
      for (Eigen::Index k = 0; k < jacobianTotalForcesOver_q_->outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(
                 *jacobianTotalForcesOver_q_, k);
             it; ++it) {
          double value = it.value();
          jacobianRhsOver_x_(idx(ndof_ + it.row(), it.col())) = value;
        }
      }
    }
    // Bottom-right block (jacobian / vel)
    if (hasJacobianTotalForcesOver_velocity()) {
      for (Eigen::Index k = 0; k < jacobianTotalForcesOver_velocity_->outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(
                 *jacobianTotalForcesOver_velocity_, k);
             it; ++it) {
          double value = it.value();
          jacobianRhsOver_x_(idx(ndof_ + it.row(), ndof_ + it.col())) = value;
        }
      }
    }
  }
  is_jacobianRhsOver_x_uptodate_ = true;
  DEBUG_EXPR(display(););
  DEBUG_END("siconos::modeling::LagrangianSparseDS::initRhs(double time)\n");
}

// ////  MASS ////

void siconos::modeling::LagrangianSparseDS::setConstantMass(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue, siconos::algebra::AliasTag) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for mass
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  mass_storage_ =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);
  hasMass_ = true;
  hasConstantMass_ = true;
  computemass_ = nullptr;
  computemass_python_ = nullptr;
}

void siconos::modeling::LagrangianSparseDS::computeMass(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position) {
  auto isComputed = siconos::algebra::computeSparse(
      mass_storage_, computemass_, computemass_python_, "mass_storage", position);
  if (isComputed) {
    hasLUMass_ = false;  // LU factorisation needs update
    is_jacobianRhsOver_x_uptodate_ = false;
  }
}

void siconos::modeling::LagrangianSparseDS::setComputeFintFunction(
    const siconos::modeling::func_prototypes::FunctionVVS_V& fint_func) {
  // Ensure that memory is properly allocated for fint_
  if (!fint_) {
    fint_ = std::make_unique<siconos::algebra::SiconosVector>(ndof_);
  }
  hasFint_ = true;
  computefint_ = fint_func;
}

void siconos::modeling::LagrangianSparseDS::computeFint(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position, double time) {
  if (computefint_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computefint_(velocity, position, time, *fint_);
}

void siconos::modeling::LagrangianSparseDS::setConstantJacobianFintOver_q(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue, siconos::algebra::AliasTag) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for the jacobian
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  jacobianFintOver_q_storage_ =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);
  hasJacobianFintOver_q_ = true;
  hasConstantJacobianFintOver_q_ = true;
  computejacobianFintOver_q_ = nullptr;
  computejacobianFintOver_q_python_ = nullptr;
  is_jacobianRhsOver_x_uptodate_ = false;
}

void siconos::modeling::LagrangianSparseDS::computeJacobianFintOver_q(
    const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<siconos::algebra::SiconosVector>& position, double time) {
  auto isComputed =
      siconos::algebra::computeSparse(jacobianFintOver_q_storage_, computejacobianFintOver_q_,
                                      computejacobianFintOver_q_python_,
                                      "jacobianFintOver_q_storage_", velocity, position, time);

  if (isComputed) is_jacobianRhsOver_x_uptodate_ = false;
}

void siconos::modeling::LagrangianSparseDS::setConstantJacobianFintOver_velocity(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue, siconos::algebra::AliasTag) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for the jacobian
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  jacobianFintOver_velocity_storage_ =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);
  hasJacobianFintOver_velocity_ = true;
  hasConstantJacobianFintOver_velocity_ = true;
  computejacobianFintOver_velocity_ = nullptr;
  computejacobianFintOver_velocity_python_ = nullptr;
  is_jacobianRhsOver_x_uptodate_ = false;
}

void siconos::modeling::LagrangianSparseDS::computeJacobianFintOver_velocity(
    const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<siconos::algebra::SiconosVector>& position, double time) {
  auto isComputed = siconos::algebra::computeSparse(
      jacobianFintOver_velocity_storage_, computejacobianFintOver_velocity_,
      computejacobianFintOver_velocity_python_, "jacobianFintOver_velocity_storage_", velocity,
      position, time);

  if (isComputed) is_jacobianRhsOver_x_uptodate_ = false;
}

void siconos::modeling::LagrangianSparseDS::setComputeFgyrFunction(
    const siconos::modeling::func_prototypes::FunctionVV_V& fgyr_func) {
  // Ensure that memory is properly allocated for fgyr_
  if (!fgyr_) {
    fgyr_ = std::make_unique<siconos::algebra::SiconosVector>(ndof_);
  }

  hasFgyr_ = true;
  computefgyr_ = fgyr_func;
}

void siconos::modeling::LagrangianSparseDS::computeFgyr(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position) {
  if (computefgyr_) computefgyr_(velocity, position, *fgyr_);
}

void siconos::modeling::LagrangianSparseDS::setConstantJacobianFgyrOver_q(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue, siconos::algebra::AliasTag) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for the jacobian
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  jacobianFgyrOver_q_storage_ =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);
  hasJacobianFgyrOver_q_ = true;
  hasConstantJacobianFgyrOver_q_ = true;
  computejacobianFgyrOver_q_ = nullptr;
  computejacobianFgyrOver_q_python_ = nullptr;
  is_jacobianRhsOver_x_uptodate_ = false;
}

void siconos::modeling::LagrangianSparseDS::computeJacobianFgyrOver_q(
    const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<siconos::algebra::SiconosVector>& position) {
  auto isComputed = siconos::algebra::computeSparse(
      jacobianFgyrOver_q_storage_, computejacobianFgyrOver_q_,
      computejacobianFgyrOver_q_python_, "jacobianFgyrOver_q_storage_", velocity, position);

  if (isComputed) is_jacobianRhsOver_x_uptodate_ = false;
}

void siconos::modeling::LagrangianSparseDS::setConstantJacobianFgyrOver_velocity(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue, siconos::algebra::AliasTag) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for the jacobian
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  jacobianFgyrOver_velocity_storage_ =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);
  hasJacobianFgyrOver_velocity_ = true;
  hasConstantJacobianFgyrOver_velocity_ = true;
  computejacobianFgyrOver_velocity_ = nullptr;
  computejacobianFgyrOver_velocity_python_ = nullptr;
  is_jacobianRhsOver_x_uptodate_ = false;
}

void siconos::modeling::LagrangianSparseDS::computeJacobianFgyrOver_velocity(
    const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<siconos::algebra::SiconosVector>& position) {
  auto isComputed = siconos::algebra::computeSparse(
      jacobianFgyrOver_velocity_storage_, computejacobianFgyrOver_velocity_,
      computejacobianFgyrOver_velocity_python_, "jacobianFgyrOver_velocity_storage_", velocity,
      position);

  if (isComputed) is_jacobianRhsOver_x_uptodate_ = false;
}

void siconos::modeling::LagrangianSparseDS::setConstantFext(
    const siconos::algebra::SiconosVector& newValue, siconos::algebra::CopyTag tag) {
  if (newValue.size() != ndof_)
    throw std::invalid_argument("setConstantFext(copy): input vector has wrong size");

  // Deep copy into Owned storage
  fext_storage_ = std::make_unique<siconos::algebra::SiconosVector>(newValue);
  hasFext_ = true;
  hasConstantFext_ = true;
  computefext_ = nullptr;
}

void siconos::modeling::LagrangianSparseDS::setConstantFext(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue, siconos::algebra::AliasTag tag) {
  if (newValue.size() != ndof_)
    throw std::invalid_argument("setConstantFext(alias): input vector has wrong size");

  // Aliasing: shared_ptr to a Map (view over external memory)
  fext_storage_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
  hasFext_ = true;
  hasConstantFext_ = true;
  computefext_ = nullptr;
}

void siconos::modeling::LagrangianSparseDS::setComputeFextFunction(
    const siconos::modeling::func_prototypes::FunctionS_V& fext_func) {
  // Ensure that memory is properly allocated for fext_
  if (!std::holds_alternative<siconos::algebra::OwnedDenseVector>(fext_storage_)) {
    fext_storage_ = std::make_unique<siconos::algebra::SiconosVector>(ndof_);
  }
  hasFext_ = true;
  hasConstantFext_ = false;
  computefext_ = fext_func;
}

void siconos::modeling::LagrangianSparseDS::computeFext(double time) {
  if (computefext_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    // siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
    //     fext_storage_, [&](auto& fext) { computefext_(time, fext); }, "fext_storage_");

    use_fext([&](auto& fext) { computefext_(time, fext); });
}

void siconos::modeling::LagrangianSparseDS::computeRhs(double time) {
  DEBUG_BEGIN("siconos::modeling::LagrangianSparseDS::computeRhs(double time)");
  *state_q_[2] = *(p_[2]);  // Warning: r/p update is done in Interactions/Relations

  computeTotalForces(*state_q_[1], *state_q_[0], time);
  if (totalForces_) {
    *state_q_[2] += *totalForces_;
    DEBUG_EXPR(std::cout << *totalForces_ << "\n";);
  }

  // Computes q[2] = inv(mass)*(fL+p) by solving Mq[2]=fL+p.
  // mass and inv(mass) updates
  init_lu_mass();

  //  if(mass->isPlugged()) : mass may be not plugged in LagrangianSparseDS children
  if (hasLUMass_) *state_q_[2] = LUMass_->solve(*state_q_[2]);

  state_x_[1]->head(ndof_) = *state_q_[1];
  state_x_[1]->tail(ndof_) = *state_q_[2];
  DEBUG_END("siconos::modeling::LagrangianSparseDS::computeRhs(double time)");
}

void siconos::modeling::LagrangianSparseDS::computeJacobianRhsOver_x(double time) {
  if (is_jacobianRhsOver_x_uptodate_ && hasConstantMass_ &&
      hasConstantJacobianTotalForcesOver_q() && hasConstantJacobianTotalForcesOver_velocity())
    return;  // Mass and all jacobian constants and jacobian rhs already up to date
  // --> this has importance for the LagrangianSparseLinearTIDS case, where this method is not
  // reimplemented

  is_jacobianRhsOver_x_uptodate_ = true;
  // Any compute(...) below if active, will turn this to false

  // We only need to update the lower part of the jacobian matrix

  if (!hasConstantJacobianTotalForcesOver_q())
    computeJacobianTotalForcesOver_q(*state_q_[1], *state_q_[0], time);

  if (!hasConstantJacobianTotalForcesOver_velocity())
    computeJacobianTotalForcesOver_velocity(*state_q_[1], *state_q_[0], time);

  if (hasMass()) init_lu_mass();  // Update mass value if required

  // If none of the updates above has changed anything, then just return
  if (is_jacobianRhsOver_x_uptodate_) return;

  // At this point, all operators are up to date

  // A lambda to compute the position in jacobian vector from the
  // coordinates in the "real" jacobian matrix.
  // Remind that jacobianRhsOver_x_ is a vector which contains
  // the concatenation of the columns of the real matrix
  auto idx = [&](int row, int col) { return col * (2 * ndof_) + row; };

  // - Fill parts corresponding to the jacobians of total forces -
  if (hasMass()) {
    // In that case, we need to compute inv(Mass).jacobian_v or _q Forces
    siconos::algebra::SiconosSparseMatrix jacq;  // to save M^-1 . jac totalforces/q
    siconos::algebra::SiconosSparseMatrix jacv;  // to save M^-1 . jac totalforces/v
    // Now fill in jacobianRhsOver_x_
    // Bottom-left block (jacobian / q)
    if (hasJacobianTotalForcesOver_q()) {
      // Solve MjacobianX(1,0) =
      // jacobianFL[0]
      jacq = LUMass_->solve(*jacobianTotalForcesOver_q_);
      for (Eigen::Index k = 0; k < jacq.outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(jacq, k); it; ++it) {
          double value = it.value();
          jacobianRhsOver_x_(idx(ndof_ + it.row(), it.col())) = value;
        }
      }
    }
    // Bottom-right block (jacobian / vel)
    if (hasJacobianTotalForcesOver_velocity()) {
      // Solve MjacobianX(1,1) =
      // jacobianFL[1]
      jacv = LUMass_->solve(*jacobianTotalForcesOver_velocity_);
      for (Eigen::Index k = 0; k < jacv.outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(jacv, k); it; ++it) {
          double value = it.value();
          jacobianRhsOver_x_(ndof_ + it.row(), ndof_ + it.col()) = value;
        }
      }
    }
  } else  // No mass
  {
    // Bottom-left block (jacobian / q)
    if (hasJacobianTotalForcesOver_q()) {
      for (Eigen::Index k = 0; k < jacobianTotalForcesOver_q_->outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(
                 *jacobianTotalForcesOver_q_, k);
             it; ++it) {
          double value = it.value();
          jacobianRhsOver_x_(idx(ndof_ + it.row(), it.col())) = value;
        }
      }
    }
    // Bottom-right block (jacobian / vel)
    if (hasJacobianTotalForcesOver_velocity()) {
      for (Eigen::Index k = 0; k < jacobianTotalForcesOver_velocity_->outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(
                 *jacobianTotalForcesOver_velocity_, k);
             it; ++it) {
          double value = it.value();
          jacobianRhsOver_x_(ndof_ + it.row(), ndof_ + it.col()) = value;
        }
      }
    }
  }
  is_jacobianRhsOver_x_uptodate_ = true;
}

void siconos::modeling::LagrangianSparseDS::computeTotalForces(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position, double time) {
  if (hasFint_ or hasFext_ or hasFgyr_) {
    if (!totalForces_) {
      totalForces_ = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
    }
    totalForces_->setZero();
  } else
    return;

  // 1 - Computes the required function
  computeFint(velocity, position, time);
  computeFext(time);
  computeFgyr(velocity, position);

  if (fint_) *totalForces_ -= *fint_;

  if (hasFext_) use_fext([&](auto const& fext) { *totalForces_ += fext; });

  if (fgyr_) *totalForces_ -= *fgyr_;
}

void siconos::modeling::LagrangianSparseDS::computeJacobianTotalForcesOver_q(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position, double time) {
  if (hasJacobianFintOver_q_ or hasJacobianFgyrOver_q_) {
    if (!jacobianTotalForcesOver_q_) {
      jacobianTotalForcesOver_q_ =
          std::make_unique<siconos::algebra::SiconosSparseMatrix>(ndof_, ndof_);
      jacobianTotalForcesOver_q_->setZero();
    }
  } else {
    return;
  }

  computeJacobianFintOver_q(*state_q_[1], *state_q_[0], time);
  computeJacobianFgyrOver_q(*state_q_[1], *state_q_[0]);
  jacobianTotalForcesOver_q_->setZero();
  if (hasJacobianFintOver_q_)
    useJacobianFintOver_q([&](const auto& M) { *jacobianTotalForcesOver_q_ -= M; });

  if (hasJacobianFgyrOver_q_)
    useJacobianFgyrOver_q([&](const auto& M) { *jacobianTotalForcesOver_q_ -= M; });
}

void siconos::modeling::LagrangianSparseDS::computeJacobianTotalForcesOver_velocity(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position, double time) {
  if (hasJacobianFintOver_velocity_ or hasJacobianFgyrOver_velocity_) {
    if (!jacobianTotalForcesOver_velocity_) {
      jacobianTotalForcesOver_velocity_ =
          std::make_unique<siconos::algebra::SiconosSparseMatrix>(ndof_, ndof_);
      jacobianTotalForcesOver_velocity_->setZero();
    }
  } else
    return;

  computeJacobianFintOver_velocity(*state_q_[1], *state_q_[0], time);
  computeJacobianFgyrOver_velocity(*state_q_[1], *state_q_[0]);
  jacobianTotalForcesOver_velocity_->setZero();
  if (hasJacobianFintOver_velocity_)
    useJacobianFintOver_velocity(
        [&](const auto& M) { *jacobianTotalForcesOver_velocity_ -= M; });
  if (hasJacobianFgyrOver_q_)
    useJacobianFgyrOver_velocity(
        [&](const auto& M) { *jacobianTotalForcesOver_velocity_ -= M; });
}

void siconos::modeling::LagrangianSparseDS::display(bool brief) const {
  brief = false;  // Full display always.
  std::cout << "=====> LagrangianSparseDS  System display (number: " << number_ << ").\n";
  std::cout << "- ndof_ : " << ndof_ << "\n";
  std::cout << "- q \n";
  siconos::algebra::print(*state_q_[0]);
  std::cout << "- q0 \n";
  use_q0([&](const auto& v) { siconos::algebra::print(v); });

  std::cout << "- velocity\n ";
  siconos::algebra::print(*state_q_[1]);
  std::cout << "- acceleration \n";
  if (state_q_[2])
    siconos::algebra::print(*state_q_[2]);
  else
    std::cout << "-> nullptr\n";

  std::cout << "- v0 \n";
  use_velocity0([&](const auto& v) { siconos::algebra::print(v); });

  std::cout << "- p[0] \n";
  if (p_[0])
    siconos::algebra::print(*p_[0]);
  else
    std::cout << "-> nullptr\n";
  std::cout << "- p[1]\n";
  siconos::algebra::print(*p_[1]);
  std::cout << "- p[2]\n";
  if (p_[2])
    siconos::algebra::print(*p_[2]);
  else
    std::cout << "-> nullptr\n";

  if (!brief) {
    if (hasMass_) {
      std::cout << "- Mass\n ";
      use_mass([&](const auto& M) { siconos::algebra::print(M); });
      std::cout << "\n";
    }

    std::cout << "- Forces \n";
    if (totalForces_)
      siconos::algebra::print(*totalForces_);
    else
      std::cout << "-> nullptr\n";

    std::cout << "- jacobian of TotalForces over q \n";
    if (jacobianTotalForcesOver_q_)
      siconos::algebra::print(*jacobianTotalForcesOver_q_);
    else
      std::cout << "-> nullptr\n";

    std::cout << "- jacobian of TotalForces over velocity \n";
    if (jacobianTotalForcesOver_velocity_)
      siconos::algebra::print(*jacobianTotalForcesOver_velocity_);
    else
      std::cout << "-> nullptr\n";

    if (fint_) std::cout << "- FInt \n" << *fint_ << "\n";
    if (hasJacobianFintOver_q_) {
      std::cout << "- jacobian of fint over q\n";
      useJacobianFintOver_q([&](const auto& M) { siconos::algebra::print(M); });
      std::cout << "\n";
    }
    if (hasJacobianFintOver_velocity_) {
      std::cout << "- jacobian of fint over velocity\n";
      useJacobianFintOver_velocity([&](const auto& M) { siconos::algebra::print(M); });
      std::cout << "\n";
    }

    if (fgyr_) std::cout << "- FGyr \n" << *fgyr_ << "\n";
    if (hasJacobianFgyrOver_q_) {
      std::cout << "- jacobian of fgyr over q\n";
      useJacobianFgyrOver_q([&](const auto& M) { siconos::algebra::print(M); });
      std::cout << "\n";
    }
    if (hasJacobianFgyrOver_velocity_) {
      std::cout << "- jacobian of fgyr over velocity\n";
      useJacobianFgyrOver_velocity([&](const auto& M) { siconos::algebra::print(M); });
      std::cout << "\n";
    }
  }

  std::cout << "===================================== \n";
}

// --- Functions for memory handling ---
void siconos::modeling::LagrangianSparseDS::initMemory(
    siconos::algebra::blocks::size_type steps) {
  if (steps == 0) {
    std::cout << "Warning : LagragianDS::initMemory with size equal to zero\n";
    return;
  }
  qMemory_.setMemorySize(steps, ndof_);
  velocityMemory_.setMemorySize(steps, ndof_);
  totalForcesMemory_.setMemorySize(steps, ndof_);
  pMemory_.resize(3);

  // TODO : initMemory in graph + f(OSI/level)
  for (siconos::algebra::blocks::size_type level = 0; level < 3; ++level) {
    if (pMemory_[level].size() == 0) pMemory_[level].setMemorySize(steps, ndof_);
  }
}

void siconos::modeling::LagrangianSparseDS::swapInMemory() {
  qMemory_.swap(*state_q_[0]);
  velocityMemory_.swap(*state_q_[1]);
  if (totalForces_) totalForcesMemory_.swap(*totalForces_);

  // initialization of the reaction force due to the non smooth law
  // note: these are a no-op if either memory or vector is null
  pMemory_[0].swap(p_[0]);
  pMemory_[1].swap(p_[1]);
  pMemory_[2].swap(p_[2]);
  xMemory_.swap(state_x_[0]);
}

void siconos::modeling::LagrangianSparseDS::resetAllNonSmoothParts() {
  if (p_[0]) p_[0]->setZero();
  if (p_[1]) p_[1]->setZero();
  if (p_[2]) p_[2]->setZero();
}

void siconos::modeling::LagrangianSparseDS::resetNonSmoothPart(
    siconos::algebra::blocks::size_type level) {
  if (level < siconos::internal::LEVELMAX)
    if (p_[level]) p_[level]->setZero();
}

void siconos::modeling::LagrangianSparseDS::computePostImpactVelocity() {
  // When this function is call, q[1] is supposed to be pre-impact velocity.
  // We solve M(v+ - v-) = p - The result is saved in(place of) p[1].
  DEBUG_BEGIN("siconos::modeling::LagrangianSparseDS::computePostImpactV()\n");
  siconos::algebra::SiconosVector tmp(*p_[1]);
  if (LUMass_) tmp = LUMass_->solve(tmp);  // we assume mass/lumass are uptodate?
  *state_q_[1] += tmp;                     // v+ = v- + p
  DEBUG_BEGIN("siconos::modeling::LagrangianSparseDS::computePostImpactV() END \n");
}

double siconos::modeling::LagrangianSparseDS::computeKineticEnergy() {
  DEBUG_BEGIN("NewtonEulerDS::computeKineticEnergy()\n");
  double K;
  if (hasMass_)
    use_mass([&](const auto& M) { K = 0.5 * (M * *state_q_[1]).dot(*state_q_[1]); });
  // warning: sparse*dense --> dense, a temp. vector is allocated.
  else
    K = 0.5 * state_q_[1]->dot(*state_q_[1]);

  DEBUG_PRINTF("Kinetic Energy = %e\n", K);
  DEBUG_END("siconos::modeling::LagrangianSparseDS::computeKineticEnergy()\n");
  return K;
}
