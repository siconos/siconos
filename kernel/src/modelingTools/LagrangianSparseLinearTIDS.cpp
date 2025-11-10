/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
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
#include "LagrangianSparseLinearTIDS.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <iostream>

#include "siconos_debug.h"

siconos::modeling::LagrangianSparseLinearTIDS::LagrangianSparseLinearTIDS(
    Eigen::Ref<siconos::algebra::SiconosVector> q0,
    Eigen::Ref<siconos::algebra::SiconosVector> v0,
    const siconos::algebra::SiconosSparseMatrix& newmass)
    : LagrangianSparseDS(q0, v0) {
  hasConstantMass_ = true;
  hasMass_ = true;
  computemass_ = nullptr;
  setConstantMassCopy(newmass);
};

siconos::modeling::LagrangianSparseLinearTIDS::LagrangianSparseLinearTIDS(
    Eigen::Ref<siconos::algebra::SiconosVector> q0,
    Eigen::Ref<siconos::algebra::SiconosVector> v0,
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newmass)
    : LagrangianSparseDS(q0, v0) {
  hasConstantMass_ = true;
  hasMass_ = true;
  computemass_ = nullptr;
  setConstantMassAlias(newmass);
};

void siconos::modeling::LagrangianSparseLinearTIDS::setStiffnessMatrixAlias(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for mass
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  stiffnessMatrix_storage =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);
}

void siconos::modeling::LagrangianSparseLinearTIDS::setDampingMatrixAlias(
    Eigen::Map<siconos::algebra::SiconosSparseMatrix>& newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for mass
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  dampingMatrix_storage =
      std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(newValue);
}

// Note FP: if required, add another constructor with shared memory between input M and
// internal mass matrix See for example LagrangianSparseDS::setConstantMassAlias.

void siconos::modeling::LagrangianSparseLinearTIDS::initRhs(double time) {
  // dim
  x_size_ = 2 * ndof_;
  // WARNING : this function is supposed to be called
  // by the OneStepIntegrator, and maybe several times for the same DS
  // if the system is involved in more than one interaction. So, we must check
  // if p2 and q2 already exist to be sure that DSlink won't be lost.

  x0_internal_storage_ = std::make_unique<std::vector<double>>(x_size_);
  x0_view_ = std::make_shared<siconos::algebra::MapVectorType>(x0_internal_storage_->data(),
                                                               x0_internal_storage_->size());
  x0_view_->head(ndof_) = *q0_view_;  // COPY !
  x0_view_->tail(ndof_) = *velocity0_view_;

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
  if (hasMass_) {
    // LU factorization
    LUMass_ = std::make_shared<siconos::algebra::SiconosSparseLUMatrix>(mass());
    hasLUMass_ = true;
  }

  computeRhs(time);

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

  // - Fill parts corresponding to the jacobians of total forces -
  // mass and lu_mass are up to date since we have already called init_lu_mass
  // in computeRhs()

  if (hasMass()) {
    // In that case, we'll need a buffer to save inv(Mass).jacobian_qForces and
    // inv(Mass).jacobian_v Forces

    // In that case, we need to compute inv(Mass).jacobian_v or _q Forces
    siconos::algebra::SiconosSparseMatrix jacq;  // to save -M^-1 . stiffness
    siconos::algebra::SiconosSparseMatrix jacv;  // to save -M^-1 . damping

    if (hasStiffnessMatrix()) {
      // Solve MjacobianX(1,0) = jacobianFL[0]
      useStiffness([&](const auto& K) { jacq = LUMass_->solve(-1. * K); });
    }
    if (hasDampingMatrix()) {
      // Solve MjacobianX(1,1) = jacobianFL[1]
      useDamping([&](const auto& C) { jacv = LUMass_->solve(-1. * C); });
    }

    // Now fill in jacobianRhsOver_x_
    if (hasStiffnessMatrix()) {
      for (Eigen::Index k = 0; k < jacq.outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(jacq, k); it; ++it) {
          double value = it.value();
          jacobianRhsOver_x_(idx(ndof_ + it.row(), it.col())) = value;
        }
      }
    }
    if (hasDampingMatrix()) {
      for (Eigen::Index k = 0; k < jacv.outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(jacv, k); it; ++it) {
          double value = it.value();
          jacobianRhsOver_x_(idx(ndof_ + it.row(), ndof_ + it.col())) = value;
        }
      }
    }

  } else  // No mass
  {       //  fill in jacobianRhsOver_x_
    if (hasStiffnessMatrix()) {
      useStiffness([&](const auto& K) {
        for (Eigen::Index k = 0; k < K.outerSize(); ++k) {
          for (Eigen::InnerIterator it(K, k); it; ++it) {
            double value = -it.value();
            jacobianRhsOver_x_(idx(ndof_ + it.row(), it.col())) = value;
          }
        }
      });
    }
    if (hasDampingMatrix()) {
      useStiffness([&](const auto& C) {
        for (Eigen::Index k = 0; k < C.outerSize(); ++k) {
          for (Eigen::InnerIterator it(C, k); it; ++it) {
            double value = -it.value();
            jacobianRhsOver_x_(idx(ndof_ + it.row(), ndof_ + it.col())) = value;
          }
        }
      });
    }
  }
  is_jacobianRhsOver_x_uptodate_ = true;
}

void siconos::modeling::LagrangianSparseLinearTIDS::display(bool brief) const {
  LagrangianSparseDS::display(brief);
  std::cout << "===== Lagrangian Linear Time Invariant System display ===== \n";

  if (hasStiffnessMatrix()) {
    std::cout << "- Stiffness Matrix K:\n";
    useStiffness([&](const auto& K) { siconos::algebra::print(K); });
    std::cout << "\n";
  }

  if (hasDampingMatrix()) {
    std::cout << "- Viscosity Matrix C:\n";
    useDamping([&](const auto& C) { siconos::algebra::print(C); });
    std::cout << "\n";
  }
  std::cout << "=========================================================== \n";
}

void siconos::modeling::LagrangianSparseLinearTIDS::computeTotalForces(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time) {
  DEBUG_PRINT("siconos::modeling::LagrangianSparseLinearTIDS::computeTotalForces(v,q,t) \n");
  if (hasFint_ or hasFext_) {  // hasFint_ = true if K or C
    if (!totalForces_) {
      totalForces_ = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
      totalForces_->setZero();
    } else
      totalForces_->setZero();
  } else
    return;

  if (fext_view_) {
    computeFext(time);
    *totalForces_ += *fext_view_;
  }

  if (hasStiffnessMatrix()) useStiffness([&](const auto& K) { *totalForces_ -= K * q; });

  if (hasDampingMatrix()) useDamping([&](const auto& C) { *totalForces_ -= C * velocity; });
}
