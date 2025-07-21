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
  mass_mat_ = std::make_shared<siconos::algebra::SiconosSparseMatrix>(newmass);
};

// Note FP: if required, add another constructor with shared memory between input M and
// internal mass matrix See for example LagrangianSparseDS::setConstantMass.

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
  if (mass_mat_) {
    // LU factorization
    LUMass_ = std::make_shared<siconos::algebra::SiconosSparseLUMatrix>(*mass_mat_);
    hasLUMass_ = true;
  }

  computeRhs(time);

  // The jacobian is saved in a flattened version, as a vector
  jacobianRhsOver_x_.resize(x_size_ * x_size_);

  // Fill null and identity part
  jacobianRhsOver_x_.setZero();
  for (unsigned int j = 0; j < ndof_; ++j) {
    jacobianRhsOver_x_((ndof_ + j) * x_size_ + j) = 1.0;
  }
  // - Fill parts corresponding to the jacobians of total forces -
  // mass and lu_mass are up to date since we have already called init_lu_mass
  if (hasMass()) {
    // In that case, we'll need a buffer to save inv(Mass).jacobian_qForces and
    // inv(Mass).jacobian_v Forces

    siconos::algebra::SiconosSparseMatrix jacq;
    siconos::algebra::SiconosSparseMatrix jacv;

    if (hasStiffnessMatrix()) {
      // Solve MjacobianX(1,0) = jacobianFL[0]
      jacq = LUMass_->solve(-1. * *stiffnessMatrix_);
    }
    if (hasDampingMatrix()) {
      // Solve MjacobianX(1,1) = jacobianFL[1]
      jacv = LUMass_->solve(-1. * *dampingMatrix_);
    }
    // Now fill in jacobianRhsOver_x_
    if (hasStiffnessMatrix()) {
      for (int k = 0; k < jacq.outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(jacq, k); it; ++it) {
          int i = it.row();
          int j = it.col();
          double value = it.value();
          int idx = j * x_size_ + i + ndof_;
          jacobianRhsOver_x_(idx) = value;
        }
      }
    }
    if (hasDampingMatrix()) {
      for (int k = 0; k < jacv.outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(jacv, k); it; ++it) {
          int i = it.row();
          int j = it.col();
          double value = it.value();
          int idx = (j + ndof_) * x_size_ + i + ndof_;
          jacobianRhsOver_x_(idx) = value;
        }
      }
    }

  } else  // No mass
  {       //  fill in jacobianRhsOver_x_
    if (hasStiffnessMatrix()) {
      for (int k = 0; k < stiffnessMatrix_->outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(*stiffnessMatrix_, k); it;
             ++it) {
          int i = it.row();
          int j = it.col();
          double value = -it.value();
          int idx = j * x_size_ + i + ndof_;
          jacobianRhsOver_x_(idx) = value;
        }
      }
    }
    if (hasDampingMatrix()) {
      for (int k = 0; k < dampingMatrix_->outerSize(); ++k) {
        for (siconos::algebra::SiconosSparseMatrix::InnerIterator it(*dampingMatrix_, k); it;
             ++it) {
          int i = it.row();
          int j = it.col();
          double value = -it.value();
          int idx = (j + ndof_) * x_size_ + i + ndof_;
          jacobianRhsOver_x_(idx) = value;
        }
      }
    }
  }
  is_jacobianRhsOver_x_uptodate_ = true;
}

void siconos::modeling::LagrangianSparseLinearTIDS::setStiffnessMatrix(
    siconos::algebra::SiconosSparseMatrix& newValue) {
  assert(newValue.rows() == newValue.cols());
  assert(newValue.rows() == ndof_);
  stiffnessMatrix_.reset(&newValue, [](siconos::algebra::SiconosSparseMatrix*) {
    // No-op deleter: the shared ptr does not own the matrix memory
    // Be cautious !!!
  });

  hasFint_ = true;
}
void siconos::modeling::LagrangianSparseLinearTIDS::setStiffnessMatrixWithCopy(
    const siconos::algebra::SiconosSparseMatrix& newValue) {
  assert(newValue.rows() == ndof_);
  assert(newValue.cols() == ndof_);
  stiffnessMatrix_ =
      std::make_shared<siconos::algebra::SiconosSparseMatrix>(newValue);  // copy
  hasFint_ = true;
}
void siconos::modeling::LagrangianSparseLinearTIDS::setDampingMatrix(
    siconos::algebra::SiconosSparseMatrix& newValue) {
  assert(newValue.rows() == newValue.cols());
  assert(newValue.rows() == ndof_);

  dampingMatrix_.reset(&newValue, [](siconos::algebra::SiconosSparseMatrix*) {
    // No-op deleter: the shared ptr does not own the matrix memory
    // Be cautious !!!
  });

  hasFint_ = true;
}
void siconos::modeling::LagrangianSparseLinearTIDS::setDampingMatrixWithCopy(
    const siconos::algebra::SiconosSparseMatrix& newValue) {
  assert(newValue.rows() == ndof_);
  assert(newValue.cols() == ndof_);
  dampingMatrix_ = std::make_shared<siconos::algebra::SiconosSparseMatrix>(newValue);  // copy
  hasFint_ = true;
}

void siconos::modeling::LagrangianSparseLinearTIDS::display(bool brief) const {
  LagrangianSparseDS::display(brief);
  std::cout << "===== Lagrangian Linear Time Invariant System display ===== \n";

  if (stiffnessMatrix_) {
    std::cout << "- Stiffness Matrix K:\n";
    siconos::algebra::print(*stiffnessMatrix_);
    std::cout << "\n";
  }

  if (dampingMatrix_) {
    std::cout << "- Viscosity Matrix C:\n";
    siconos::algebra::print(*dampingMatrix_);
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

  if (stiffnessMatrix_) *totalForces_ -= *stiffnessMatrix_ * q;
  if (dampingMatrix_) *totalForces_ -= *dampingMatrix_ * velocity;
}
