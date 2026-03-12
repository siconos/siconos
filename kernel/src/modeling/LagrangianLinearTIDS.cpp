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
#include "LagrangianLinearTIDS.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <iostream>

#include "StorageTools.hpp"
#include "siconos_debug.h"

siconos::modeling::LagrangianLinearTIDS::LagrangianLinearTIDS(
    Eigen::Ref<siconos::algebra::SiconosVector> q0,
    Eigen::Ref<siconos::algebra::SiconosVector> v0,
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newmass, siconos::algebra::AliasTag)
    : LagrangianDS(q0, v0, siconos::algebra::alias_t) {
  setConstantMass(newmass, siconos::algebra::alias_t);
}

siconos::modeling::LagrangianLinearTIDS::LagrangianLinearTIDS(
    const siconos::algebra::SiconosVector& q0, const siconos::algebra::SiconosVector& v0,
    const siconos::algebra::SiconosDenseMatrix& newmass, siconos::algebra::CopyTag tag)
    : LagrangianDS(q0, v0, siconos::algebra::copy_t) {
  setConstantMass(newmass, siconos::algebra::copy_t);
}

void siconos::modeling::LagrangianLinearTIDS::initRhs(double time) {
  // dim
  x_size_ = 2 * ndof_;
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

  // Compute mass and LU factorization
  init_lu_mass();

  computeRhs(time);

  // The jacobian is saved in a flattened version, as a vector
  jacobianRhsOver_x_.resize(x_size_ * x_size_);

  // Fill null and identity part
  jacobianRhsOver_x_.setZero();
  for (siconos::algebra::Index j = 0; j < ndof_; ++j) {
    jacobianRhsOver_x_((ndof_ + j) * x_size_ + j) = 1.0;
  }
  // - Fill parts corresponding to the jacobians of total forces -
  // mass and lu_mass are (must be) up to date
  if (hasMass()) {
    // In that case, we'll need a buffer to save inv(Mass).jacobian_qForces and
    // inv(Mass).jacobian_v Forces
    if (hasStiffnessMatrix() || hasDampingMatrix()) buffer_.resize(ndof_ * ndof_ * 2);
    // View onto left part of buffer_
    Eigen::Map<siconos::algebra::SiconosDenseMatrix> jacq(buffer_.data(), ndof_, ndof_);
    // View onto right part of buffer_
    Eigen::Map<siconos::algebra::SiconosDenseMatrix> jacv(buffer_.data() + ndof_ * ndof_,
                                                          ndof_, ndof_);

    if (hasStiffnessMatrix()) {
      // Solve MjacobianX(1,0) = jacobianFL[0]
      useStiffness([&](const auto& K) { jacq = LUMass_->solve(-1. * K); });
    }
    if (hasDampingMatrix()) {
      // Solve MjacobianX(1,1) = jacobianFL[1]
      useDamping([&](const auto& C) { jacv = LUMass_->solve(-1. * C); });
    }
    // Now fill in jacobianRhsOver_x_
    for (siconos::algebra::Index j = 0; j < ndof_; ++j) {
      // Bottom-left block (jacobian / q)
      if (hasStiffnessMatrix()) {
        for (siconos::algebra::Index i = 0; i < ndof_; ++i)
          jacobianRhsOver_x_(j * x_size_ + i + ndof_) = jacq(i, j);
      }
      // Bottom-right block (jacobian / vel)
      if (hasDampingMatrix()) {
        for (siconos::algebra::Index i = 0; i < ndof_; ++i)
          jacobianRhsOver_x_((j + ndof_) * x_size_ + i + ndof_) = jacv(i, j);
      }
    }
  } else  // No mass
  {       // ==> no buffer
    //  fill in jacobianRhsOver_x_
    for (siconos::algebra::Index j = 0; j < ndof_; ++j) {
      // Bottom-left block (jacobian / q)
      if (hasStiffnessMatrix()) {
        useStiffness([&](const auto& K) {
          for (siconos::algebra::Index i = 0; i < ndof_; ++i)

            jacobianRhsOver_x_(j * x_size_ + i + ndof_) = -K(i, j);
        });
      }
      // Bottom-right block (jacobian / vel)
      if (hasDampingMatrix()) {
        useDamping([&](const auto& C) {
          for (siconos::algebra::Index i = 0; i < ndof_; ++i)
            jacobianRhsOver_x_((j + ndof_) * x_size_ + i + ndof_) = -C(i, j);
        });
      }
    }
  }
  is_jacobianRhsOver_x_uptodate_ = true;
}

void siconos::modeling::LagrangianLinearTIDS::setStiffnessMatrix(
    const siconos::algebra::SiconosDenseMatrix& newValue, siconos::algebra::CopyTag tag) {
  assert(newValue.rows() == newValue.cols());
  assert(newValue.rows() == ndof_);

  stiffnessMatrix_storage_ = std::make_unique<siconos::algebra::SiconosDenseMatrix>(newValue);
  hasFint_ = true;
}

void siconos::modeling::LagrangianLinearTIDS::setStiffnessMatrix(
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue,
    siconos::algebra::AliasTag tag) {
  assert(newValue.rows() == newValue.cols());
  assert(newValue.rows() == ndof_);

  stiffnessMatrix_storage_ =
      std::make_shared<siconos::algebra::MapType>(newValue.data(), ndof_, ndof_);
  hasFint_ = true;
}

void siconos::modeling::LagrangianLinearTIDS::setDampingMatrix(
    const siconos::algebra::SiconosDenseMatrix& newValue, siconos::algebra::CopyTag tag) {
  assert(newValue.rows() == newValue.cols());
  assert(newValue.rows() == ndof_);

  dampingMatrix_storage_ = std::make_unique<siconos::algebra::SiconosDenseMatrix>(newValue);
  hasFint_ = true;
}

void siconos::modeling::LagrangianLinearTIDS::setDampingMatrix(
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue,
    siconos::algebra::AliasTag tag) {
  assert(newValue.rows() == newValue.cols());
  assert(newValue.rows() == ndof_);

  dampingMatrix_storage_ =
      std::make_shared<siconos::algebra::MapType>(newValue.data(), ndof_, ndof_);
  hasFint_ = true;
}

void siconos::modeling::LagrangianLinearTIDS::display(bool brief) const {
  LagrangianDS::display(brief);
  std::cout << "===== Lagrangian Linear Time Invariant System display ===== \n";
  if (hasStiffnessMatrix()) {
    std::cout << "- Stiffness matrix\n ";
    useStiffness([&](const auto& M) { siconos::algebra::print(M); });
    std::cout << "\n";
  }
  if (hasDampingMatrix()) {
    std::cout << "- Damping matrix\n ";
    useDamping([&](const auto& M) { siconos::algebra::print(M); });
    std::cout << "\n";
  }
  std::cout << "=========================================================== \n";
}

void siconos::modeling::LagrangianLinearTIDS::computeTotalForces(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time) {
  DEBUG_PRINT("siconos::modeling::LagrangianLinearTIDS::computeTotalForces(v,q,t) \n");
  if (hasFint_ or hasFext_) {  // hasFint_ = true if K or C
    if (!totalForces_) {
      totalForces_ = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
      totalForces_->setZero();
    } else
      totalForces_->setZero();
  } else
    return;

  if (hasFext_) {
    computeFext(time);
    use_fext([&](auto const& fext) { *totalForces_ += fext; });
  }

  if (hasStiffnessMatrix()) useStiffness([&](auto const& K) { *totalForces_ -= K * q; });
  if (hasDampingMatrix()) useDamping([&](auto const& C) { *totalForces_ -= C * velocity; });
}
