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
#include "NewtonEulerDS.hpp"

#include <boost/math/quaternion.hpp>
#include <iostream>

#include "NewtonEulerR.hpp"
#include "RotationQuaternion.hpp"
#include "SiconosAlgebraAddons.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"
// #include "Tools.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::modeling::NewtonEulerDS::NewtonEulerDS(
    Eigen::Ref<siconos::algebra::SiconosVector7> initial_position,
    Eigen::Ref<siconos::algebra::SiconosVector6> initial_twist, double mass,
    Eigen::Ref<siconos::algebra::SiconosMatrix33> inertia, siconos::algebra::AliasTag)
    : SecondOrderDS(13, 6), scalarMass_{mass} {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::NewtonEulerDS(...)\n");

  // --- Initial conditions ---
  // q0 : contains the center of mass coordinate, and the quaternion initial. (dim(Q0)=7)
  // twist0 : contains the initial velocity of center of mass and the omega initial.
  // (dim(VTwist0)=6)

  q0_storage_ = std::make_shared<siconos::algebra::MapVector7Type>(initial_position.data(),
                                                                   initial_position.size());
  twist0_storage_ = std::make_shared<siconos::algebra::MapVector6Type>(initial_twist.data(),
                                                                       initial_twist.size());

  // --- Current state (q and twist variables) ---
  state_q_ = std::make_shared<siconos::algebra::SiconosVector>(q0());
  twist_ = std::make_shared<siconos::algebra::SiconosVector>(twist0());

  dotq_ = std::make_shared<siconos::algebra::SiconosVector>();
  dotq_->setZero();

  /** \todo lazy Memory allocation */
  p_.resize(3);
  p_[1] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);  // Needed in NewtonEulerR
  p_[1]->setZero();
  //  -- Total Inertia Matrix --
  //   Remind that inertial matrix is accessed with inertialMatrix() method, a view on
  //   the block
  totalInertiaMatrix_ = std::make_shared<siconos::algebra::SiconosMatrix66>();
  totalInertiaMatrix_->setZero();
  (*totalInertiaMatrix_)(0, 0) = scalarMass_;
  (*totalInertiaMatrix_)(1, 1) = scalarMass_;
  (*totalInertiaMatrix_)(2, 2) = scalarMass_;
  totalInertiaMatrix_->block<3, 3>(3, 3) = inertia;  // COPY
  hasLUMass_ = false;

  // --- T(q) matrix ---

  T_ = std::make_unique<siconos::algebra::SiconosMatrix76>();  // qDim_, ndof_);
  T_->setZero();
  (*T_)(0, 0) = 1.0;
  (*T_)(1, 1) = 1.0;
  (*T_)(2, 2) = 1.0;
  siconos::modeling::newton_euler::computeT(*state_q_, *T_);

  // -- Tdot --
  // lazy init, at first call of computeTdot()

  // --- Wrench ---
  wrench_ = std::make_shared<siconos::algebra::SiconosVector6>();
  wrench_->setZero();
  /** The follwing jacobian are always allocated since we have always
   * Gyroscopical forces that has non linear forces
   * This should be remove if the integration is explicit or _nullifyMGyr(false) is set to true
   */

  jacobianWrenchOver_twist_ =
      std::make_shared<siconos::algebra::SiconosMatrix66>();  // ndof_, ndof_);
  jacobianWrenchOver_twist_->setZero();
  // jacobianWrenchOver_q_ will be allocated only if required
  // (if fint and/or mint are defined or if mext is expressed in the inertial frame.)

  DEBUG_END("siconos::modeling::NewtonEulerDS::NewtonEulerDS(...)\n");
}

siconos::modeling::NewtonEulerDS::NewtonEulerDS(
    const siconos::algebra::SiconosVector7& initial_position,
    const siconos::algebra::SiconosVector6& initial_twist, double mass,
    const siconos::algebra::SiconosMatrix33& inertia, siconos::algebra::CopyTag)
    : SecondOrderDS(13, 6), scalarMass_{mass} {  // --- Initial conditions ---
  // q0 : contains the center of mass coordinate, and the quaternion initial. (dim(Q0)=7)
  // twist0 : contains the initial velocity of center of mass and the omega initial.
  // (dim(VTwist0)=6)

  q0_storage_ = std::make_unique<siconos::algebra::SiconosVector7>(initial_position);
  twist0_storage_ = std::make_unique<siconos::algebra::SiconosVector6>(initial_twist);

  // --- Current state (q and twist variables) ---
  state_q_ = std::make_shared<siconos::algebra::SiconosVector>(q0());
  twist_ = std::make_shared<siconos::algebra::SiconosVector>(twist0());

  dotq_ = std::make_shared<siconos::algebra::SiconosVector>();
  dotq_->setZero();

  /** \todo lazy Memory allocation */
  p_.resize(3);
  p_[1] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);  // Needed in NewtonEulerR
  p_[1]->setZero();
  //  -- Total Inertia Matrix --
  //   Remind that inertial matrix is accessed with inertialMatrix() method, a view on
  //   the bloxk
  totalInertiaMatrix_ = std::make_shared<siconos::algebra::SiconosMatrix66>();
  totalInertiaMatrix_->setZero();
  (*totalInertiaMatrix_)(0, 0) = scalarMass_;
  (*totalInertiaMatrix_)(1, 1) = scalarMass_;
  (*totalInertiaMatrix_)(2, 2) = scalarMass_;
  totalInertiaMatrix_->block<3, 3>(3, 3) = inertia;  // COPY
  hasLUMass_ = false;

  // --- T(q) matrix ---

  T_ = std::make_unique<siconos::algebra::SiconosMatrix76>();  // qDim_, ndof_);
  T_->setZero();
  (*T_)(0, 0) = 1.0;
  (*T_)(1, 1) = 1.0;
  (*T_)(2, 2) = 1.0;
  siconos::modeling::newton_euler::computeT(*state_q_, *T_);

  // -- Tdot --
  // lazy init, at first call of computeTdot()

  // --- Wrench ---
  wrench_ = std::make_shared<siconos::algebra::SiconosVector6>();
  wrench_->setZero();
  /** The follwing jacobian are always allocated since we have always
   * Gyroscopical forces that has non linear forces
   * This should be remove if the integration is explicit or _nullifyMGyr(false) is set to true
   */

  jacobianWrenchOver_twist_ =
      std::make_shared<siconos::algebra::SiconosMatrix66>();  // ndof_, ndof_);
  jacobianWrenchOver_twist_->setZero();
  // jacobianWrenchOver_q_ will be allocated only if required
  // (if fint and/or mint are defined or if mext is expressed in the inertial frame.)
}

void siconos::modeling::NewtonEulerDS::resetToInitialState() {
  // set q and q[1] to q0 and Twist0
  *state_q_ = q0();
  *twist_ = twist0();
}

void siconos::modeling::NewtonEulerDS::initRhs(double time) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::initRhs(double time)\n");
  // dim
  x_size_ = qDim_ + 6;
  x0_storage_ = std::make_unique<siconos::algebra::SiconosVector>(x_size_);

  use_x0([&](auto& xinit) {
    xinit.head(ndof_) = q0();  // COPY !
    xinit.tail(ndof_) = twist0();
  });

  state_x_[0] = std::make_shared<siconos::algebra::SiconosVector>(x_size_);
  *(state_x_[0]) << *state_q_, *twist_;

  if (!acceleration_) {
    acceleration_ = std::make_shared<siconos::algebra::SiconosVector>(6);
    acceleration_->setZero();
  }

  // Compute dotq_
  siconos::modeling::newton_euler::computeT(*state_q_, *T_);
  *dotq_ = *T_ * *state_q_;

  state_x_[1] = std::make_shared<siconos::algebra::SiconosVector>(x_size_);
  *(state_x_[1]) << *dotq_, *acceleration_;

  // Nothing to do for the initialization of the wrench

  // Everything concerning rhs and its jacobian is handled in initRhs and computeXXX
  // related functions.

  if (!p_[2]) {
    p_[2] = std::make_shared<siconos::algebra::SiconosVector>(6);
    p_[2]->setZero();
  }

  init_lu_mass();

  computeRhs(time);

  // -- Jacobian Rhs over x --
  /** \warning the derivative of T w.r.t to q is neglected */

  // The jacobian is saved in a flattened version, as a vector
  jacobianRhsOver_x_.resize(x_size_ * x_size_);
  // initialize to zero
  jacobianRhsOver_x_.setZero();

  // Set top-right block to T
  for (siconos::algebra::Index j = 0; j < 6; ++j) {
    for (siconos::algebra::Index i = 0; i < qDim_; ++i) {
      jacobianRhsOver_x_((j + qDim_) * x_size_ + i) = (*T_)(i, j);
    }
  }

  // - Fill parts corresponding to the jacobians of total forces -
  // mass and lu_mass are up to date since we have already called init_lu_mass
  // In that case, we'll need a buffer to save inv(Mass).jacobian_qForces and
  // inv(Mass).jacobian_v Forces
  if (hasJacobianWrenchOver_q() || hasJacobianWrenchOver_twist())
    buffer_.resize(6 * qDim_ + 36);

  // bottom-right block
  if (hasJacobianWrenchOver_q()) {
    // Update if required
    computeJacobianWrenchOver_q(*twist_, *state_q_, time);
    // View onto left part of buffer_
    Eigen::Map<siconos::algebra::SiconosMatrix67> jacq(buffer_.data(), ndof_, qDim_);
    // Solve MjacobianX(1,0) = jacobianFL[0]
    siconos::algebra::solve(*LUMass_, *totalInertiaMatrix_, *jacobianWrenchOver_q_, jacq);
    for (siconos::algebra::Index j = 0; j < qDim_; ++j) {
      for (siconos::algebra::Index i = 0; i < 6; ++i) {
        jacobianRhsOver_x_(j * x_size_ + i + qDim_) = jacq(i, j);
      }
    }
  }

  // bottom-left block
  if (hasJacobianWrenchOver_twist()) {
    // Update if required
    computeJacobianWrenchOver_twist(*twist_, *state_q_, time);
    // View onto right part of buffer_
    Eigen::Map<siconos::algebra::SiconosMatrix66> jacv(buffer_.data() + 6 * qDim_, ndof_,
                                                       ndof_);
    // Solve MjacobianX(1,1) = jacobianFL[1]
    siconos::algebra::solve(*LUMass_, *totalInertiaMatrix_, *jacobianWrenchOver_twist_, jacv);
    for (siconos::algebra::Index j = 0; j < 6; ++j) {
      for (siconos::algebra::Index i = 0; i < 6; ++i) {
        jacobianRhsOver_x_((j + qDim_) * x_size_ + i + qDim_) = jacv(i, j);
      }
    }
  }
  is_jacobianRhsOver_x_uptodate_ = true;

  DEBUG_EXPR(display(););
  DEBUG_END("siconos::modeling::NewtonEulerDS::initRhs(double time)\n");
}

void siconos::modeling::NewtonEulerDS::initializeNonSmoothInput(
    siconos::algebra::blocks::size_type level) {
  DEBUG_PRINTF(
      "siconos::modeling::NewtonEulerDS::initializeNonSmoothInput(...) for "
      "level = %i\n",
      level);

  if (!p_[level]) {
    if (level == 0) {
      p_[level] = std::make_shared<siconos::algebra::SiconosVector>(qDim_);
    } else
      p_[level] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  }
  p_[level]->setZero();
}

void siconos::modeling::NewtonEulerDS::computeRhs(double time) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeRhs(double time)");
  *acceleration_ = *(p_[2]);  // Warning: r/p update is done in Interactions/Relations

  computeWrench(*twist_, *state_q_, time);
  *acceleration_ += *wrench_;
  DEBUG_EXPR(siconos::algebra::print(*wrench_););

  if (LUMass_)
    siconos::algebra::solve(*LUMass_, *totalInertiaMatrix_, *acceleration_, *acceleration_);

  // Compute dotq_
  siconos::modeling::newton_euler::computeT(*state_q_, *T_);
  *dotq_ = *T_ * *twist_;

  state_x_[1]->head(qDim_) = *dotq_;
  state_x_[1]->tail(qDim_) = *acceleration_;
}

void siconos::modeling::NewtonEulerDS::computeJacobianRhsOver_x(double time) {
  siconos::modeling::newton_euler::computeT(*state_q_, *T_);
  for (siconos::algebra::Index j = 0; j < 6; ++j) {
    for (siconos::algebra::Index i = 0; i < qDim_; ++i) {
      jacobianRhsOver_x_((j + qDim_) * x_size_ + i) = (*T_)(i, j);
    }
  }
  // bottom-right block
  if (hasJacobianWrenchOver_q()) {
    // Update if required
    computeJacobianWrenchOver_q(*twist_, *state_q_, time);
    // View onto left part of buffer_
    Eigen::Map<siconos::algebra::SiconosMatrix67> jacq(buffer_.data(), ndof_, ndof_);
    // Solve MjacobianX(1,0) = jacobianFL[0]
    siconos::algebra::solve(*LUMass_, *totalInertiaMatrix_, *jacobianWrenchOver_q_, jacq);

    for (siconos::algebra::Index j = 0; j < qDim_; ++j) {
      for (siconos::algebra::Index i = 0; i < 6; ++i) {
        jacobianRhsOver_x_(j * x_size_ + i + qDim_) = jacq(i, j);
      }
    }
  }
  // bottom-left block
  if (hasJacobianWrenchOver_twist()) {
    // Update if required
    computeJacobianWrenchOver_twist(*twist_, *state_q_, time);
    // View onto right part of buffer_
    Eigen::Map<siconos::algebra::SiconosMatrix66> jacv(buffer_.data(), ndof_, ndof_);
    // Solve MjacobianX(1,1) = jacobianFL[1]
    siconos::algebra::solve(*LUMass_, *totalInertiaMatrix_, *jacobianWrenchOver_twist_, jacv);

    for (siconos::algebra::Index j = 0; j < 6; ++j) {
      for (siconos::algebra::Index i = 0; i < 6; ++i) {
        jacobianRhsOver_x_((j + qDim_) * x_size_ + i + qDim_) = jacv(i, j);
      }
    }
  }
}

void siconos::modeling::NewtonEulerDS::setInertia(double ix, double iy, double iz) {
  (*totalInertiaMatrix_)(3, 3) = ix;
  (*totalInertiaMatrix_)(4, 4) = iy;
  (*totalInertiaMatrix_)(5, 5) = iz;
  hasLUMass_ = false;
}

void siconos::modeling::NewtonEulerDS::init_lu_mass() {
  if (totalInertiaMatrix_ && !LUMass_) {
    LUMass_ = std::make_shared<siconos::algebra::SiconosDenseLUMatrix>(*totalInertiaMatrix_);
    hasLUMass_ = true;
  }
}

/////////////////// FEXT  ////////////////////

void siconos::modeling::NewtonEulerDS::setConstantFext(
    const siconos::algebra::SiconosVector3& newValue, siconos::algebra::CopyTag tag) {
  // Deep copy into Owned storage
  fext_storage_ = std::make_unique<siconos::algebra::SiconosVector3>(newValue);
  hasFext_ = true;
  hasConstantFext_ = true;
  computefext_ = nullptr;
}

void siconos::modeling::NewtonEulerDS::setConstantFext(
    Eigen::Ref<siconos::algebra::SiconosVector3> newValue, siconos::algebra::AliasTag tag) {
  fext_storage_ =
      std::make_shared<siconos::algebra::MapVector3Type>(newValue.data(), newValue.size());

  hasFext_ = true;
  hasConstantFext_ = true;
  computefext_ = nullptr;
}

void siconos::modeling::NewtonEulerDS::setComputeFextFunction(
    const func_prototypes::FunctionS_V& fext_func) {
  if (!std::holds_alternative<siconos::algebra::OwnedDenseVector3>(fext_storage_)) {
    fext_storage_ = std::make_unique<siconos::algebra::SiconosVector3>();
    use_fext([&](auto& vec) { vec.setZero(); });
  }
  hasFext_ = true;
  hasConstantFext_ = false;
  computefext_ = fext_func;
}

/////////////////// MEXT  ////////////////////

void siconos::modeling::NewtonEulerDS::setConstantMext(
    const siconos::algebra::SiconosVector3& newValue, siconos::algebra::CopyTag tag) {
  // Deep copy into Owned storage
  mext_storage_ = std::make_unique<siconos::algebra::SiconosVector3>(newValue);

  hasMext_ = true;
  hasConstantMext_ = true;
  computemext_ = nullptr;
}

void siconos::modeling::NewtonEulerDS::setConstantMext(
    Eigen::Ref<siconos::algebra::SiconosVector3> newValue, siconos::algebra::AliasTag tag) {
  mext_storage_ =
      std::make_shared<siconos::algebra::MapVector3Type>(newValue.data(), newValue.size());

  hasMext_ = true;
  hasConstantMext_ = true;
  computemext_ = nullptr;
}

void siconos::modeling::NewtonEulerDS::setComputeMextFunction(
    const func_prototypes::FunctionS_V& mext_func) {
  if (!std::holds_alternative<siconos::algebra::OwnedDenseVector3>(mext_storage_)) {
    mext_storage_ = std::make_unique<siconos::algebra::SiconosVector3>();
    use_mext([&](auto& vec) { vec.setZero(); });
  }
  hasMext_ = true;
  hasConstantMext_ = false;
  computemext_ = mext_func;
}

/////////////////// FINT and its Jacobians  ////////////////////

void siconos::modeling::NewtonEulerDS::setComputeFintFunction(
    const func_prototypes::FunctionVVS_V& fint_func) {
  // No need to ensure memory alloc -> computed directly into wrench_
  hasFint_ = true;
  computefint_ = fint_func;
  if (!jacobianWrenchOver_q_) {
    jacobianWrenchOver_q_ =
        std::make_shared<siconos::algebra::SiconosMatrix67>();  // ndof_, qDim_);
    jacobianWrenchOver_q_->setZero();
  }
}

void siconos::modeling::NewtonEulerDS::setComputeJacobianFintOver_qFunction(
    const func_prototypes::FunctionVVS_M& new_func) {
  // No need to ensure memory alloc -> computed directly into jacobianWrenchOver_q_
  hasJacobianFintOver_q_ = true;
  computejacobianFintOver_q_ = new_func;
  computeJacobianFintOver_q_byFD_ = false;
  if (!jacobianWrenchOver_q_) {
    jacobianWrenchOver_q_ =
        std::make_shared<siconos::algebra::SiconosMatrix67>();  // ndof_, qDim_);
    jacobianWrenchOver_q_->setZero();
  }
}

void siconos::modeling::NewtonEulerDS::setComputeJacobianFintOver_twistFunction(
    const func_prototypes::FunctionVVS_M& new_func) {
  // No need to ensure memory alloc -> computed directly into jacobianWrenchOver_q_
  hasJacobianFintOver_twist_ = true;
  computejacobianFintOver_twist_ = new_func;
  computeJacobianFintOver_twist_byFD_ = false;
}

/////////////////// MINT and its Jacobians  ////////////////////

void siconos::modeling::NewtonEulerDS::setComputeMintFunction(
    const func_prototypes::FunctionVVS_V& mint_func) {
  // No need to ensure memory alloc -> computed directly into wrench_
  hasMint_ = true;
  computemint_ = mint_func;
  if (!jacobianWrenchOver_q_) {
    jacobianWrenchOver_q_ = std::make_shared<siconos::algebra::SiconosMatrix67>();
    jacobianWrenchOver_q_->setZero();
  }
}

void siconos::modeling::NewtonEulerDS::setComputeJacobianMintOver_qFunction(
    const func_prototypes::FunctionVVS_M& new_func) {
  // No need to ensure memory alloc -> computed directly into jacobianWrenchOver_q_
  computejacobianMintOver_q_ = new_func;
  computeJacobianMintOver_q_byFD_ = false;
  if (!jacobianWrenchOver_q_) {
    jacobianWrenchOver_q_ = std::make_shared<siconos::algebra::SiconosMatrix67>();
    jacobianWrenchOver_q_->setZero();
  }
}

void siconos::modeling::NewtonEulerDS::setComputeJacobianMintOver_twistFunction(
    const func_prototypes::FunctionVVS_M& new_func) {
  // No need to ensure memory alloc -> computed directly into jacobianWrenchOver_q_
  computejacobianMintOver_twist_ = new_func;
  computeJacobianMintOver_twist_byFD_ = false;
}

void siconos::modeling::NewtonEulerDS::setComputeJacobianFintOver_q_byFD(bool value) {
  computeJacobianFintOver_q_byFD_ = value;
}

void siconos::modeling::NewtonEulerDS::setComputeJacobianFintOver_twist_byFD(bool value) {
  computeJacobianFintOver_twist_byFD_ = value;
}

void siconos::modeling::NewtonEulerDS::setComputeJacobianMintOver_q_byFD(bool value) {
  computeJacobianMintOver_q_byFD_ = value;
}

void siconos::modeling::NewtonEulerDS::setComputeJacobianMintOver_twist_byFD(bool value) {
  computeJacobianMintOver_twist_byFD_ = value;
}

void siconos::modeling::NewtonEulerDS::computeWrench(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q, double time) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeWrench(...)\n")

  assert(wrench_);
  wrench_->setZero();

  // wrench = [fext(t)  - fint(t) , mext(t) - mgyr(twist) - mint(twist, q, t)]

  // A buffer (fixed size, on the stack)
  // used to compute temp. values for fext, fint and so on
  // Warning: if parallel computation is in action, mind
  // competitive access to the memory
  siconos::algebra::SiconosVector3 buffer = siconos::algebra::SiconosVector3::Zero();
  if (hasFext_) {
    if (computefext_) {
      use_fext([&](auto& fext) {
        computefext_(time, fext);
        wrench_->head(3) += fext;
      });
    } else  // if (hasConstantFext_)
      use_fext([&](auto& fext) { wrench_->head(3) += fext; });
  }

  if (hasFint_) {
    computefint_(twist, q, time, buffer);
    // wrench[0:2] -= fint
    wrench_->head(3) -= buffer;
  }

  if (hasMext_) {
    if (computemext_)  // mext not constant --> we need to use the buffer
      use_mext([&](auto& mext) { computemext_(time, mext); });
    if (isMextExpressedInInertialFrame_) {
      use_mext([&](auto& mext) {
        siconos::geometry::rewriteVectorFromAbsoluteToBodyFrame(q, mext);
      });
      use_mext([&](auto& mext) { wrench_->tail(3) += mext; });
      // wrench[3:6] += mext
    }
  }

  if (hasMgyr_) {
    siconos::modeling::newton_euler::computeMgyr(twist, *totalInertiaMatrix_, buffer);
    wrench_->tail(3) -= buffer;
  }
  if (hasMint_) {
    computemint_(twist, q, time, buffer);
    wrench_->tail(3) -= buffer;
  }

  DEBUG_EXPR(siconos::algebra::print(*wrench_));
  DEBUG_END("siconos::modeling::NewtonEulerDS::computeWrench(...) \n");
}

void siconos::modeling::NewtonEulerDS::computeJacobianWrenchOver_q(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q, double time) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeJacobianWrenchOver_q(double time) \n");

  if (!jacobianWrenchOver_q_) return;

  jacobianWrenchOver_q_->setZero();

  siconos::algebra::SiconosMatrix37 matrix_buffer;  //{ 3, 7};
  // TMP. TODO: external setup for this memory
  matrix_buffer.setZero();

  // Jacobian fint ?
  if (hasFint_) {
    if (computejacobianFintOver_q_) {  // if user-defined function
      computejacobianFintOver_q_(twist, q, time, matrix_buffer);
    } else {
      siconos::modeling::newton_euler::computeJacobianFOver_q_byFD(
          twist, q, time, std::sqrt(std::numeric_limits<double>::epsilon()), computefint_,
          matrix_buffer);
    }
    jacobianWrenchOver_q_->topRows(3) -= matrix_buffer;
  }
  // Jacobian mint ?
  if (hasMint_) {
    if (computejacobianMintOver_q_)  // if user-defined function
      computejacobianMintOver_q_(twist, q, time, matrix_buffer);

    else  // FD computation
      siconos::modeling::newton_euler::computeJacobianFOver_q_byFD(
          twist, q, time, std::sqrt(std::numeric_limits<double>::epsilon()), computemint_,
          matrix_buffer);

    jacobianWrenchOver_q_->bottomRows(3) -= matrix_buffer;
  }

  if (hasMext_) {
    if (isMextExpressedInInertialFrame_) {
      // tmp value to compute mext. We do not update mext_ of the present object
      siconos::algebra::SiconosVector3 mext = siconos::algebra::SiconosVector3::Zero();
      if (computemext_)
        computemext_(time, mext);
      else
        use_mext([&](auto& val) { mext = val; });
      newton_euler::computeJacobianMExtqExpressedInInertialFrame(*state_q_, time, mext,
                                                                 matrix_buffer);
      jacobianWrenchOver_q_->bottomRows(3) += matrix_buffer;
    }
    DEBUG_EXPR(siconos::algebra::print(*_jacobianWrenchq););
  }
}

void siconos::modeling::NewtonEulerDS::computeJacobianWrenchOver_twist(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q, double time) {
  if (!jacobianWrenchOver_twist_) return;

  jacobianWrenchOver_twist_->setZero();

  siconos::algebra::SiconosMatrix36 matrix_buffer;  //{3, 6};
  // TMP. TODO: external setup for this memory
  matrix_buffer.setZero();

  if (hasFint_) {
    if (computejacobianFintOver_twist_) {  // if user-defined function
      computejacobianFintOver_twist_(twist, q, time, matrix_buffer);
    } else {
      siconos::modeling::newton_euler::computeJacobianFOver_twist_byFD(
          twist, q, time, std::sqrt(std::numeric_limits<double>::epsilon()), computefint_,
          matrix_buffer);
    }
    jacobianWrenchOver_twist_->topRows(3) -= matrix_buffer;
  }

  // Jacobian mint ?
  if (hasMint_) {
    if (computejacobianMintOver_twist_)  // if user-defined function
      computejacobianMintOver_twist_(twist, q, time, matrix_buffer);

    else  // FD computation
      siconos::modeling::newton_euler::computeJacobianFOver_twist_byFD(
          twist, q, time, std::sqrt(std::numeric_limits<double>::epsilon()), computemint_,
          matrix_buffer);

    jacobianWrenchOver_twist_->bottomRows(3) -= matrix_buffer;
  }

  if (hasMgyr_) {
    // siconos::modeling::newton_euler::computeJacobianMGyrOver_twist_byFD(...);
    siconos::modeling::newton_euler::computeJacobianMGyrOver_twist(twist, totalInertiaMatrix(),
                                                                   matrix_buffer);
    jacobianWrenchOver_twist_->bottomRows(3) -= matrix_buffer;
  }
}

void siconos::modeling::NewtonEulerDS::display(bool brief) const {
  std::cout << "=====> NewtonEuler System display (number: " << number_ << ").\n";
  std::cout << "- q \n";
  siconos::algebra::print(*state_q_);

  std::cout << "- initial state: \n";
  use_q0([&](const auto& v) { siconos::algebra::print(v); });

  std::cout << "- twist \n";
  siconos::algebra::print(*twist_);

  std::cout << "- twist0 \n ";
  use_twist0([&](const auto& v) { siconos::algebra::print(v); });

  std::cout << "- dotq \n";
  if (dotq_)
    siconos::algebra::print(*dotq_);
  else
    std::cout << "-> nullptr\n";

  std::cout << "- p[0] \n";
  if (p_[0])
    siconos::algebra::print(*p_[0]);
  else
    std::cout << "-> nullptr\n";

  std::cout << "- p[1] \n";
  if (p_[1])
    siconos::algebra::print(*p_[1]);
  else
    std::cout << "-> nullptr\n";

  std::cout << "- p[2] \n";
  if (p_[2])
    siconos::algebra::print(*p_[2]);
  else
    std::cout << "-> nullptr\n";

  std::cout << "total inertia matrix :\n" << *totalInertiaMatrix_ << std::endl;

  std::cout << "===================================== \n";
}

void siconos::modeling::NewtonEulerDS::setIsMextExpressedInInertialFrame(bool value) {
  isMextExpressedInInertialFrame_ = value;
  if (!jacobianWrenchOver_q_) {
    jacobianWrenchOver_q_ =
        std::make_shared<siconos::algebra::SiconosMatrix67>();  // ndof_, qDim_);
    jacobianWrenchOver_q_->setZero();
  }
}

// --- Functions for memory handling ---
void siconos::modeling::NewtonEulerDS::initMemory(siconos::algebra::blocks::size_type steps) {
  DynamicalSystem::initMemory(steps);

  if (steps == 0)
    std::cout << "Warning : siconos::modeling::NewtonEulerDS::initMemory with size equal to "
                 "zero\n";
  else {
    qMemory_.setMemorySize(steps, qDim_);
    twistMemory_.setMemorySize(steps, ndof_);
    wrenchMemory_.setMemorySize(steps, ndof_);
    dotqMemory_.setMemorySize(steps, qDim_);
    //    swapInMemory(); Useless, done in osi->initializeWorkVectorsForDS
  }
}

void siconos::modeling::NewtonEulerDS::swapInMemory() {
  //  xMemory_->swap(_x[0]);
  qMemory_.swap(*state_q_);
  twistMemory_.swap(*twist_);
  dotqMemory_.swap(*dotq_);
  wrenchMemory_.swap(*wrench_);
}

void siconos::modeling::NewtonEulerDS::resetAllNonSmoothParts() {
  if (!p_[1]) p_[1] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  p_[1]->setZero();
}

void siconos::modeling::NewtonEulerDS::resetNonSmoothPart(
    siconos::algebra::blocks::size_type level) {
  if (p_[level]) p_[level]->setZero();
}

void siconos::modeling::NewtonEulerDS::computeT(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q) {
  siconos::modeling::newton_euler::computeT(q, *T_);
}

void siconos::modeling::NewtonEulerDS::computeTdot() {
  if (!Tdot_) {
    Tdot_ = std::make_unique<siconos::algebra::SiconosMatrix76>();  // qDim_, ndof_);
    Tdot_->setZero();
  }
  // Update Tdot[3:6,3:5]
  siconos::modeling::newton_euler::computeT(*dotq_, *Tdot_);
}

void siconos::modeling::NewtonEulerDS::normalizeq() {
  siconos::geometry::normalizeq(*state_q_);
}

double siconos::modeling::NewtonEulerDS::computeKineticEnergy() {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeKineticEnergy()\n");
  assert(twist_);
  assert(totalInertiaMatrix_);
  DEBUG_EXPR(siconos::algebra::print(*twist_));

  double K = 0.5 * (*totalInertiaMatrix_ * *twist_).dot(*twist_);

  DEBUG_PRINTF("Kinetic Energy = %e\n", K);
  DEBUG_END("siconos::modeling::NewtonEulerDS::computeKineticEnergy()\n");
  return K;
}

siconos::algebra::SiconosVector3 siconos::modeling::NewtonEulerDS::linearVelocityInBodyFrame()
    const {
  siconos::algebra::SiconosVector3 v = twist_->head(3);  // copy
  siconos::geometry::rewriteVectorFromAbsoluteToBodyFrame(*state_q_, v);
  return v;  // RVO, no copy!
}

siconos::algebra::SiconosVector3 siconos::modeling::NewtonEulerDS::angularVelocityInBodyFrame()
    const {
  siconos::algebra::SiconosVector3 w = twist_->tail(3);
  siconos::geometry::rewriteVectorFromBodyToAbsoluteFrame(*state_q_, w);
  return w;  // RVO, no copy!
}

void siconos::modeling::NewtonEulerDS::setScalarMass(double mass) {
  scalarMass_ = mass;
  (*totalInertiaMatrix_)(0, 0) = scalarMass_;
  (*totalInertiaMatrix_)(1, 1) = scalarMass_;
  (*totalInertiaMatrix_)(2, 2) = scalarMass_;
  hasLUMass_ = false;
};

void siconos::modeling::NewtonEulerDS::addExtForceAtPos(
    const Eigen::Ref<siconos::algebra::SiconosVector3>& force, bool forceAbsRef,
    const Eigen::Ref<siconos::algebra::SiconosVector3>& pos, bool posAbsRef) {
  use_fext([&](auto& fext) {
    newton_euler::computeFextForceAtPos(*state_q_, force, forceAbsRef, fext, true);
  });
  use_mext([&](auto& mext) {
    newton_euler::computeMextForceAtPos(*state_q_, isMextExpressedInInertialFrame_, force,
                                        forceAbsRef, pos, posAbsRef, mext, true);
  });
}

///////////////////////////////////////////////////////////////////////////////////////
/////////// Free functions in the namespace siconos::modeling::newton_euler ///////////
///////////////////////////////////////////////////////////////////////////////////////

void siconos::modeling::newton_euler::computeT(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q,
    Eigen::Ref<siconos::algebra::SiconosMatrix76> T) {
  DEBUG_BEGIN("computeT(...)\n")

  // Warning : we assume that
  // - T is properly allocated and of size 7x6
  // - T[0:3,0:3] is already initialized to I_{3x3} (identity)

  double q0 = q(3) / 2.0;
  double q1 = q(4) / 2.0;
  double q2 = q(5) / 2.0;
  double q3 = q(6) / 2.0;
  T(3, 3) = -q1;
  T(3, 4) = -q2;
  T(3, 5) = -q3;
  T(4, 3) = q0;
  T(4, 4) = -q3;
  T(4, 5) = q2;
  T(5, 3) = q3;
  T(5, 4) = q0;
  T(5, 5) = -q1;
  T(6, 3) = -q2;
  T(6, 4) = q1;
  T(6, 5) = q0;
  DEBUG_END("computeT(...)\n")
}

void siconos::modeling::newton_euler::computeMextForceAtPos(
    const Eigen::Ref<siconos::algebra::SiconosVector7>& q, bool isMextExpressedInInertialFrame,
    const Eigen::Ref<siconos::algebra::SiconosVector3>& force, bool forceAbsRef,
    const Eigen::Ref<siconos::algebra::SiconosVector3>& pos, bool posAbsRef,
    Eigen::Ref<siconos::algebra::SiconosVector3> mExt, bool accumulate) {
  assert(force.size() == 3);
  assert(mExt.size() == 3);

  siconos::algebra::SiconosVector3 local_frc(force);  // copy

  if (forceAbsRef) {
    siconos::geometry::rewriteVectorFromAbsoluteToBodyFrame(q, local_frc);
  }

  siconos::algebra::SiconosVector3 moment = siconos::algebra::SiconosVector3::Zero();
  if (posAbsRef) {
    siconos::algebra::SiconosVector3 local_pos = pos - q.head<3>();
    siconos::geometry::rewriteVectorFromAbsoluteToBodyFrame(q, local_pos);
    moment = local_pos.cross(local_frc);
  } else {
    moment = pos.cross(local_frc);
  }

  if (isMextExpressedInInertialFrame)
    siconos::geometry::rewriteVectorFromBodyToAbsoluteFrame(q, moment);

  if (accumulate)
    mExt += moment;
  else
    mExt = moment;
}

void siconos::modeling::newton_euler::computeFextForceAtPos(
    const Eigen::Ref<siconos::algebra::SiconosVector7>& q,
    const Eigen::Ref<siconos::algebra::SiconosVector3>& force, bool forceAbsRef,
    Eigen::Ref<siconos::algebra::MapVector3Type> fext, bool accumulate) {
  siconos::algebra::SiconosVector3 abs_frc(force);

  if (!forceAbsRef) siconos::geometry::rewriteVectorFromBodyToAbsoluteFrame(q, abs_frc);
  if (accumulate)
    fext += abs_frc;
  else
    fext = abs_frc;
}

void siconos::modeling::newton_euler::computeMgyr(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
    const Eigen::Ref<const siconos::algebra::SiconosMatrix66>& inertiaMatrix,
    Eigen::Ref<siconos::algebra::SiconosVector3> result) {
  auto omega = twist.tail<3>();
  auto inertia = inertiaMatrix.block<3, 3>(3, 3);
  result = omega.cross(inertia * omega);
}

void siconos::modeling::newton_euler::computeJacobianMGyrOver_twist(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
    const Eigen::Ref<const siconos::algebra::SiconosMatrix66>& inertia,
    Eigen::Ref<siconos::algebra::SiconosMatrix36> result) {
  result.setZero();
  // input is assumed to be a 3x3 matrix
  // result(i) = ei x Inertia.Omega  + Omega x I ei

  auto omega = twist.tail<3>();
  auto iomega = inertia.block<3, 3>(3, 3) * omega;

  siconos::algebra::SiconosVector3 ei;
  /*See equation of DevNotes.pdf, equation with label eq:NE_nablaFL1*/

  for (int i = 0; i < 3; i++) {
    ei.setZero();
    ei(i) = 1.0;
    result.col(i + 3) = ei.cross(iomega) + omega.cross(inertia.block<3, 3>(3, 3) * ei);
  }
  // Check if Jacobian is valid. Warning to the transpose operation in
  // _jacobianMGyrtwist->setValue(3 + j, 3 + i, ei_Iomega(j) +
  // omega_Iei(j));
}

void siconos::modeling::newton_euler::computeJacobianMGyrOver_twist_byFD(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist, double epsilonFD,
    const Eigen::Ref<const siconos::algebra::SiconosMatrix66>& inertia,
    const siconos::modeling::func_prototypes::FunctionMV_V& mgyr_func,
    Eigen::Ref<siconos::algebra::SiconosMatrix36> result) {
  siconos::algebra::SiconosVector3 mgyr0 = siconos::algebra::SiconosVector3::Zero();
  siconos::algebra::SiconosVector3 mgyr = siconos::algebra::SiconosVector3::Zero();
  mgyr_func(inertia, twist, mgyr0);

  result.setZero();
  siconos::algebra::SiconosVector6 veps = twist;  // copy
  for (int j = 0; j < 3; j++) {
    veps(0) += epsilonFD;
    mgyr_func(inertia, veps, mgyr);
    result.col(j + 3) = (mgyr - mgyr0) / epsilonFD;
    veps(j) -= epsilonFD;
  }
}

void siconos::modeling::newton_euler::computeJacobianMExtqExpressedInInertialFrame(
    const Eigen::Ref<siconos::algebra::SiconosVector7>& q, double time,
    const siconos::algebra::SiconosVector3& mext,
    Eigen::Ref<siconos::algebra::SiconosMatrix37> result) {
  DEBUG_EXPR(siconos::algebra::print(q));
  DEBUG_EXPR(siconos::algebra::print(mext));

  double q0 = q(3);
  double q1 = q(4);
  double q2 = q(5);
  double q3 = q(6);
  double v0 = mext(0);
  double v1 = mext(1);
  double v2 = mext(2);
  // This routine compute the jacobian with respect to p of R^T(p)v
  // See 11.225 ... in devNotes.pdf
  result.setZero();
  result(0, 3) = q0 * v0 + q3 * v1 - q2 * v2;
  result(0, 4) = q1 * v0 + q2 * v1 + q3 * v2;
  result(0, 5) = -q2 * v0 + q1 * v1 - q0 * v2;
  result(0, 6) = -q3 * v0 + q0 * v1 + q1 * v2;

  result(1, 3) = -q3 * v0 + q0 * v1 + q1 * v2;
  result(1, 4) = q2 * v0 - q1 * v1 + q0 * v2;
  result(1, 5) = q1 * v0 + q2 * v1 + q3 * v2;
  result(1, 6) = -q0 * v0 - q3 * v1 + q2 * v2;

  result(2, 3) = q2 * v0 - q1 * v1 + q0 * v2;
  result(2, 4) = q3 * v0 - q0 * v1 - q1 * v2;
  result(2, 5) = q0 * v0 + q3 * v1 - q2 * v2;
  result(2, 6) = q1 * v0 + q2 * v1 + q3 * v2;

  result *= 2.0;
}

void siconos::modeling::newton_euler::computeJacobianMExtqExpressedInInertialFrameByFD(
    const Eigen::Ref<siconos::algebra::SiconosVector7>& q, double time,
    const siconos::modeling::func_prototypes::FunctionS_V& mext_func,
    bool isMextExpressedInInertialFrame, double epsilonFD,
    Eigen::Ref<siconos::algebra::SiconosMatrix33> result)

{
  /* The computation of Jacobian of R^T mExt is somehow very rough since the pertubation
   * that we apply to q  that gives qeps does not provide a unit quaternion. The rotation
   * is computed assuming that the quaternion is unit (see quaternionRotate...()
   */

  siconos::algebra::SiconosVector3 mext = siconos::algebra::SiconosVector3::Zero();
  siconos::algebra::SiconosVector3 mext0 = siconos::algebra::SiconosVector3::Zero();
  mext_func(time, mext0);
  if (isMextExpressedInInertialFrame)
    siconos::geometry::rewriteVectorFromAbsoluteToBodyFrame(q, mext0);

  auto qeps = q;  // copy
  result.setZero();

  for (int j = 3; j < 7; j++) {
    qeps(j) += epsilonFD;
    mext_func(time, mext);
    if (isMextExpressedInInertialFrame)
      siconos::geometry::rewriteVectorFromAbsoluteToBodyFrame(qeps, mext);
    result.col(j) = (mext - mext0) / epsilonFD;
    qeps(j) -= epsilonFD;
  }
}

void siconos::modeling::newton_euler::computeJacobianFOver_twist_byFD(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q, double time, double epsilonFD,
    const siconos::modeling::func_prototypes::FunctionVVS_V& f_func,
    Eigen::Ref<siconos::algebra::SiconosMatrix36> result) {
  siconos::algebra::SiconosVector3 mint0 = siconos::algebra::SiconosVector3::Zero();
  siconos::algebra::SiconosVector3 mint = siconos::algebra::SiconosVector3::Zero();
  f_func(twist, q, time, mint0);
  result.setZero();
  siconos::algebra::SiconosVector6 veps = twist;  // copy
  auto inveps = 1. / epsilonFD;
  for (int j = 0; j < veps.size(); j++) {
    veps(j) += epsilonFD;
    f_func(veps, q, time, mint);
    result.col(j) = (mint - mint0) * inveps;
    veps(j) -= epsilonFD;
  }
}

void siconos::modeling::newton_euler::computeJacobianFOver_q_byFD(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q, double time, double epsilonFD,
    const siconos::modeling::func_prototypes::FunctionVVS_V& f_func,
    Eigen::Ref<siconos::algebra::SiconosMatrix37> result) {
  siconos::algebra::SiconosVector3 mint0 = siconos::algebra::SiconosVector3::Zero();
  siconos::algebra::SiconosVector3 mint = siconos::algebra::SiconosVector3::Zero();
  f_func(twist, q, time, mint0);
  result.setZero();
  siconos::algebra::SiconosVector qeps = q;  // copy
  auto inveps = 1. / epsilonFD;
  for (int j = 0; j < qeps.size(); j++) {
    qeps(j) += epsilonFD;
    f_func(twist, qeps, time, mint);
    result.col(j) = (mint - mint0) * inveps;
    qeps(j) -= epsilonFD;
  }
}
