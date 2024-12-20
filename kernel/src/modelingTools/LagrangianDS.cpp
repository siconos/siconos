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
#include "LagrangianDS.hpp"

#include <iostream>
#include <memory>

#include "BlockMatrix.hpp"
#include "BlockVector.hpp"
#include "SiconosConst.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "SiconosVisitor.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES

#include "siconos_debug.h"

// Build from initial state only
siconos::modeling::LagrangianDS::LagrangianDS(Eigen::Ref<siconos::algebra::SiconosVector> q0,
                                              Eigen::Ref<siconos::algebra::SiconosVector> v0)
    : SecondOrderDS(2 * q0.size(), v0.size()) {
  assert(ndof_ > 0 && "lagrangian dynamical system dimension should be greater than 0.");

  // Set initial conditions
  q0_view_ = std::make_shared<siconos::algebra::MapVectorType>(q0.data(), ndof_);
  velocity0_view_ = std::make_shared<siconos::algebra::MapVectorType>(v0.data(), ndof_);

  // -- Memory allocation for vector and matrix members --
  state_q_[0] = std::make_shared<siconos::algebra::SiconosVector>(*q0_view_);
  state_q_[1] = std::make_shared<siconos::algebra::SiconosVector>(*velocity0_view_);

  /** \todo lazy Memory allocation */
  p_[1] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  p_[1]->setZero();
}

void siconos::modeling::LagrangianDS::initializeNonSmoothInput(unsigned int level) {
  if (!p_[level]) {
    p_[level] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
    p_[level]->setZero();
  }
}

void siconos::modeling::LagrangianDS::resetToInitialState() {
  if (q0_view_) {
    *(state_q_[0]) = *q0_view_;
  } else
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianDS::resetToInitialState - initial "
        "position q0_view_ is null");
  if (velocity0_view_) {
    *(state_q_[1]) = *velocity0_view_;
  } else
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianDS::resetToInitialState - initial "
        "velocity velocity0_view_ "
        "is null");
}

void siconos::modeling::LagrangianDS::initMemoryForGeneralizedCoordinates(unsigned int level) {
  assert(level > 1);
  if (!state_q_[level])
    state_q_[level] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
}

void siconos::modeling::LagrangianDS::init_lu_mass() {
  if (mass_view_) {
    computeMass(*state_q_[0]);
    // LU factorization
    LUMass_ = std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(*mass_view_);
    hasLUMass_ = true;
  }
}

void siconos::modeling::LagrangianDS::update_lu_mass() {
  if (mass_view_ && !hasConstantMass_) {
    computeMass(*state_q_[0]);
    LUMass_ = std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(*mass_view_);
  }
}

void siconos::modeling::LagrangianDS::initRhs(double time) {
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::initRhs(double time)\n");
  // dim
  x_size_ = 2 * ndof_;

  // All links between DS and LagrangianDS class members are pointer links,
  // which means that no useless memory is allocated when connection is
  // established. One exception: zero and identity matrices, used to filled in M
  // and jacobianfx.

  // Initial conditions and state

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

  if (fint_ || fext_view_ || fgyr_) {
    if (!totalForces_) totalForces_ = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  }

  if (fint_ || fgyr_) {
    hasJacobianTotalForces_ = true;
    if (!jacobianTotalForcesOver_q_)
      jacobianTotalForcesOver_q_ =
          std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_);
    if (!jacobianTotalForcesOver_velocity_)
      jacobianTotalForcesOver_velocity_ =
          std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_);
  }

  init_lu_mass();

  computeRhs(time);

  bool flag1 = false, flag2 = false;
  if (jacobianTotalForcesOver_q_) {
    // Solve MjacobianX(1,0) = jacobianFL[0]
    computeJacobianTotalForcesOver_q(*state_q_[1], *state_q_[0], time);

    rhsMatrices_[jacobianXBloc10_] =
        std::make_shared<siconos::algebra::SiconosMatrix>(*jacobianTotalForcesOver_q_);
    *rhsMatrices_[jacobianXBloc10_] = LUMass_->solve(*rhsMatrices_[jacobianXBloc10_]);
    flag1 = true;
  }

  if (jacobianTotalForcesOver_velocity_) {
    // Solve MjacobianX(1,1) = jacobianFL[1]
    computeJacobianTotalForcesOver_velocity(*state_q_[1], *state_q_[0], time);
    rhsMatrices_[jacobianXBloc11_] =
        std::make_shared<siconos::algebra::SiconosMatrix>(*jacobianTotalForcesOver_velocity_);
    *rhsMatrices_[jacobianXBloc11_] = LUMass_->solve(*rhsMatrices_[jacobianXBloc11_]);
    flag2 = true;
  }

  if (!rhsMatrices_[zeroMatrix_]) {
    rhsMatrices_[zeroMatrix_] =
        std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_);
    rhsMatrices_[zeroMatrix_]->setZero();
  }
  if (!rhsMatrices_[idMatrix_]) {
    rhsMatrices_[idMatrix_] = std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_);
    rhsMatrices_[idMatrix_]->setIdentity();
  }
  if (flag1 && flag2)
    jacobianRhsOver_x_ = std::make_shared<siconos::algebra::BlockMatrix>(
        rhsMatrices_[zeroMatrix_], rhsMatrices_[idMatrix_], rhsMatrices_[jacobianXBloc10_],
        rhsMatrices_[jacobianXBloc11_]);
  else if (flag1)  // flag2 = false
    jacobianRhsOver_x_ = std::make_shared<siconos::algebra::BlockMatrix>(
        rhsMatrices_[zeroMatrix_], rhsMatrices_[idMatrix_], rhsMatrices_[jacobianXBloc10_],
        rhsMatrices_[zeroMatrix_]);
  else if (flag2)  // flag1 = false
    jacobianRhsOver_x_ = std::make_shared<siconos::algebra::BlockMatrix>(
        rhsMatrices_[zeroMatrix_], rhsMatrices_[idMatrix_], rhsMatrices_[zeroMatrix_],
        rhsMatrices_[jacobianXBloc11_]);
  else
    jacobianRhsOver_x_ = std::make_shared<siconos::algebra::BlockMatrix>(
        rhsMatrices_[zeroMatrix_], rhsMatrices_[idMatrix_], rhsMatrices_[zeroMatrix_],
        rhsMatrices_[zeroMatrix_]);
  DEBUG_EXPR(display(););
  DEBUG_END("siconos::modeling::LagrangianDS::initRhs(double time)\n");
}

////  MASS ////

void siconos::modeling::LagrangianDS::setConstantMass(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for mass
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  mass_internal_storage_ = nullptr;

  mass_view_ = std::make_shared<siconos::algebra::MapType>(newValue.data(), ndof_, ndof_);
  hasMass_ = true;
  hasConstantMass_ = true;
  computemass_ = nullptr;
}
void siconos::modeling::LagrangianDS::setComputeMassFunction(
    const siconos::modeling::func_prototypes::FunctionV_M& new_func) {
  // Ensure that memory is properly allocated for mass_
  if (!mass_internal_storage_) {
    mass_internal_storage_ = std::make_unique<std::vector<double>>(ndof_ * ndof_);
  }
  mass_view_ = std::make_shared<siconos::algebra::MapType>(mass_internal_storage_->data(),
                                                           ndof_, ndof_);
  mass_view_->setZero();
  hasMass_ = true;
  hasConstantMass_ = false;
  computemass_ = new_func;
}

void siconos::modeling::LagrangianDS::computeMass(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position) {
  if (computemass_) {
    computemass_(position, *mass_view_);
  }
}

////////////////

void siconos::modeling::LagrangianDS::setComputeFintFunction(
    const siconos::modeling::func_prototypes::FunctionVVS_V& fint_func) {
  // Ensure that memory is properly allocated for fint_
  if (!fint_) {
    fint_ = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  }
  hasFint_ = true;
  computefint_ = fint_func;
}

void siconos::modeling::LagrangianDS::computeFint(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position, double time) {
  if (computefint_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computefint_(velocity, position, time, *fint_);
}

void siconos::modeling::LagrangianDS::setConstantJacobianFintOver_q(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for jacobianFintOver_q
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  jacobianFintOver_q_internal_storage_ = nullptr;

  jacobianFintOver_q_view_ =
      std::make_shared<siconos::algebra::MapType>(newValue.data(), ndof_, ndof_);
  hasJacobianFintOver_q_ = true;
  hasConstantJacobianFintOver_q_ = true;
  computejacobianFintOver_q_ = nullptr;
}

void siconos::modeling::LagrangianDS::setComputeJacobianFintOver_qFunction(
    const siconos::modeling::func_prototypes::FunctionVVS_M& new_func) {
  // Ensure that memory is properly allocated for jacobianFintOver_q_
  if (!jacobianFintOver_q_internal_storage_) {
    jacobianFintOver_q_internal_storage_ =
        std::make_unique<std::vector<double>>(ndof_ * ndof_);
  }
  jacobianFintOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianFintOver_q_internal_storage_->data(), ndof_, ndof_);
  jacobianFintOver_q_view_->setZero();
  hasJacobianFintOver_q_ = true;
  hasConstantJacobianFintOver_q_ = false;
  computejacobianFintOver_q_ = new_func;
}

void siconos::modeling::LagrangianDS::computeJacobianFintOver_q(
    const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<siconos::algebra::SiconosVector>& position, double time) {
  if (computejacobianFintOver_q_) {
    computejacobianFintOver_q_(velocity, position, time, *jacobianFintOver_q_view_);
  }
}

void siconos::modeling::LagrangianDS::setConstantJacobianFintOver_velocity(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for jacobianFintOver_velocity
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  jacobianFintOver_velocity_internal_storage_ = nullptr;

  jacobianFintOver_velocity_view_ =
      std::make_shared<siconos::algebra::MapType>(newValue.data(), ndof_, ndof_);
  hasJacobianFintOver_velocity_ = true;
  hasConstantJacobianFintOver_velocity_ = true;
  computejacobianFintOver_velocity_ = nullptr;
}

void siconos::modeling::LagrangianDS::setComputeJacobianFintOver_velocityFunction(
    const siconos::modeling::func_prototypes::FunctionVVS_M& new_func) {
  // Ensure that memory is properly allocated for jacobianFintOver_velocity_
  if (!jacobianFintOver_velocity_internal_storage_) {
    jacobianFintOver_velocity_internal_storage_ =
        std::make_unique<std::vector<double>>(ndof_ * ndof_);
  }
  jacobianFintOver_velocity_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianFintOver_velocity_internal_storage_->data(), ndof_, ndof_);
  jacobianFintOver_velocity_view_->setZero();
  hasJacobianFintOver_velocity_ = true;
  hasConstantJacobianFintOver_velocity_ = false;
  computejacobianFintOver_velocity_ = new_func;
}

void siconos::modeling::LagrangianDS::computeJacobianFintOver_velocity(
    const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<siconos::algebra::SiconosVector>& position, double time) {
  if (computejacobianFintOver_velocity_) {
    computejacobianFintOver_velocity_(velocity, position, time,
                                      *jacobianFintOver_velocity_view_);
  }
}

void siconos::modeling::LagrangianDS::setComputeFgyrFunction(
    const siconos::modeling::func_prototypes::FunctionVV_V& fgyr_func) {
  // Ensure that memory is properly allocated for fgyr_
  if (!fgyr_) {
    fgyr_ = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  }

  hasFgyr_ = true;
  computefgyr_ = fgyr_func;
}

void siconos::modeling::LagrangianDS::computeFgyr(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position) {
  if (computefgyr_) computefgyr_(velocity, position, *fgyr_);
}

void siconos::modeling::LagrangianDS::setConstantJacobianFgyrOver_q(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for jacobianFgyrOver_q
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  jacobianFgyrOver_q_internal_storage_ = nullptr;

  jacobianFgyrOver_q_view_ =
      std::make_shared<siconos::algebra::MapType>(newValue.data(), ndof_, ndof_);
  hasJacobianFgyrOver_q_ = true;
  hasConstantJacobianFgyrOver_q_ = true;
  computejacobianFgyrOver_q_ = nullptr;
}

void siconos::modeling::LagrangianDS::setComputeJacobianFgyrOver_qFunction(
    const siconos::modeling::func_prototypes::FunctionVV_M& new_func) {
  // Ensure that memory is properly allocated for jacobianFgyrOver_q_
  if (!jacobianFgyrOver_q_internal_storage_) {
    jacobianFgyrOver_q_internal_storage_ =
        std::make_unique<std::vector<double>>(ndof_ * ndof_);
  }
  jacobianFgyrOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianFgyrOver_q_internal_storage_->data(), ndof_, ndof_);
  jacobianFgyrOver_q_view_->setZero();
  hasJacobianFgyrOver_q_ = true;
  hasConstantJacobianFgyrOver_q_ = false;
  computejacobianFgyrOver_q_ = new_func;
}

void siconos::modeling::LagrangianDS::computeJacobianFgyrOver_q(
    const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<siconos::algebra::SiconosVector>& position) {
  if (computejacobianFgyrOver_q_) {
    computejacobianFgyrOver_q_(velocity, position, *jacobianFgyrOver_q_view_);
  }
}

void siconos::modeling::LagrangianDS::setConstantJacobianFgyrOver_velocity(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for jacobianFgyrOver_velocity
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  jacobianFgyrOver_velocity_internal_storage_ = nullptr;

  jacobianFgyrOver_velocity_view_ =
      std::make_shared<siconos::algebra::MapType>(newValue.data(), ndof_, ndof_);
  hasJacobianFgyrOver_velocity_ = true;
  hasConstantJacobianFgyrOver_velocity_ = true;
  computejacobianFgyrOver_velocity_ = nullptr;
}

void siconos::modeling::LagrangianDS::setComputeJacobianFgyrOver_velocityFunction(
    const siconos::modeling::func_prototypes::FunctionVV_M& new_func) {
  // Ensure that memory is properly allocated for jacobianFgyrOver_velocity_
  if (!jacobianFgyrOver_velocity_internal_storage_) {
    jacobianFgyrOver_velocity_internal_storage_ =
        std::make_unique<std::vector<double>>(ndof_ * ndof_);
  }
  jacobianFgyrOver_velocity_view_ = std::make_shared<siconos::algebra::MapType>(
      jacobianFgyrOver_velocity_internal_storage_->data(), ndof_, ndof_);
  jacobianFgyrOver_velocity_view_->setZero();
  hasJacobianFgyrOver_velocity_ = true;
  hasConstantJacobianFgyrOver_velocity_ = false;
  computejacobianFgyrOver_velocity_ = new_func;
}

void siconos::modeling::LagrangianDS::computeJacobianFgyrOver_velocity(
    const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<siconos::algebra::SiconosVector>& position) {
  if (computejacobianFgyrOver_velocity_) {
    computejacobianFgyrOver_velocity_(velocity, position, *jacobianFgyrOver_velocity_view_);
  }
}

void siconos::modeling::LagrangianDS::setConstantFext(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for Fext_
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  fext_internal_storage_ = nullptr;

  fext_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
  hasFext_ = true;
  hasConstantFext_ = true;
  computefext_ = nullptr;
}

void siconos::modeling::LagrangianDS::setComputeFextFunction(
    const siconos::modeling::func_prototypes::FunctionS_V& fext_func) {
  // Ensure that memory is properly allocated for fext_
  if (!fext_internal_storage_) {
    fext_internal_storage_ = std::make_unique<std::vector<double>>(ndof_);
  }
  fext_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(fext_internal_storage_->data(), ndof_);

  hasFext_ = true;
  hasConstantFext_ = false;
  computefext_ = fext_func;
}

void siconos::modeling::LagrangianDS::computeFext(double time) {
  if (computefext_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computefext_(time, *fext_view_);
}

void siconos::modeling::LagrangianDS::computeRhs(double time) {
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::computeRhs(double time)");
  *state_q_[2] = *(p_[2]);  // Warning: r/p update is done in Interactions/Relations

  // if(totalForces_)
  //   {
  computeTotalForces(*state_q_[1], *state_q_[0], time);
  if (totalForces_) {
    *state_q_[2] += *totalForces_;
    DEBUG_EXPR(std::cout << *totalForces_ << "\n";);
  }

  // Computes q[2] = inv(mass)*(fL+p) by solving Mq[2]=fL+p.
  // -- Case 1: if mass is constant, then a copy of imass is LU-factorized
  // during initialization and saved into LUMass_
  // -- Case 2: mass is not constant, we copy it into LUMass_
  // Then we proceed with PLUForwardBackward.
  // mass and inv(mass) computation
  update_lu_mass();

  //  if(mass->isPlugged()) : mass may be not plugged in LagrangianDS children
  if (LUMass_) *state_q_[2] = LUMass_->solve(*state_q_[2]);

  state_x_[1]->head(ndof_) = *state_q_[1];
  state_x_[1]->tail(ndof_) = *state_q_[2];
  DEBUG_END("siconos::modeling::LagrangianDS::computeRhs(double time)");
}

void siconos::modeling::LagrangianDS::computeJacobianRhsOver_x(double time) {
  computeMass(*state_q_[0]);

  if (jacobianTotalForcesOver_q_ || jacobianTotalForcesOver_velocity_) {
    update_lu_mass();
  }

  if (jacobianTotalForcesOver_q_) {
    /** \warning the Jacobian of the inverse of the mass matrix
     * w.r.t q is not taken into account */

    std::shared_ptr<siconos::algebra::SiconosMatrix> bloc10 = jacobianRhsOver_x_->block(1, 0);
    computeJacobianTotalForcesOver_q(*state_q_[1], *state_q_[0], time);
    *bloc10 = *jacobianTotalForcesOver_q_;
    *bloc10 = LUMass_->solve(*bloc10);
  }

  if (jacobianTotalForcesOver_velocity_) {
    std::shared_ptr<siconos::algebra::SiconosMatrix> bloc11 = jacobianRhsOver_x_->block(1, 1);
    computeJacobianTotalForcesOver_velocity(*state_q_[1], *state_q_[0], time);
    *bloc11 = *jacobianTotalForcesOver_velocity_;
    *bloc11 = LUMass_->solve(*bloc11);
  }
}

void siconos::modeling::LagrangianDS::computeTotalForces(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position, double time) {
  if (hasFint_ or hasFext_ or hasFgyr_) {
    if (!totalForces_) {
      totalForces_ = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
      totalForces_->setZero();
    } else
      totalForces_->setZero();
  } else
    return;

  // 1 - Computes the required function
  computeFint(velocity, position, time);
  computeFext(time);
  computeFgyr(position, velocity);

  if (fint_) *totalForces_ -= *fint_;

  if (fext_view_) *totalForces_ += *fext_view_;

  if (fgyr_) *totalForces_ -= *fgyr_;
}

void siconos::modeling::LagrangianDS::computeJacobianTotalForcesOver_q(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position, double time) {
  if (hasJacobianFintOver_q_ or hasJacobianFgyrOver_q_) {
    if (!jacobianTotalForcesOver_q_) {
      jacobianTotalForcesOver_q_ =
          std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_);
      jacobianTotalForcesOver_q_->setZero();
    }
  } else
    return;

  computeJacobianFintOver_q(*state_q_[1], *state_q_[0], time);
  computeJacobianFgyrOver_q(*state_q_[1], *state_q_[0]);
  jacobianTotalForcesOver_q_->setZero();
  if (hasJacobianFintOver_q_) *jacobianTotalForcesOver_q_ -= *jacobianFintOver_q_view_;
  if (hasJacobianFgyrOver_q_) *jacobianTotalForcesOver_q_ -= *jacobianFgyrOver_q_view_;
}

void siconos::modeling::LagrangianDS::computeJacobianTotalForcesOver_velocity(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& position, double time) {
  if (hasJacobianFintOver_velocity_ or hasJacobianFgyrOver_velocity_) {
    if (!jacobianTotalForcesOver_velocity_) {
      jacobianTotalForcesOver_velocity_ =
          std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_);
      jacobianTotalForcesOver_velocity_->setZero();
    }
  } else
    return;

  computeJacobianFintOver_velocity(*state_q_[1], *state_q_[0], time);
  computeJacobianFgyrOver_velocity(*state_q_[1], *state_q_[0]);
  jacobianTotalForcesOver_velocity_->setZero();
  if (hasJacobianFintOver_velocity_)
    *jacobianTotalForcesOver_velocity_ -= *jacobianFintOver_velocity_view_;
  if (hasJacobianFgyrOver_q_)
    *jacobianTotalForcesOver_velocity_ -= *jacobianFgyrOver_velocity_view_;
}

void siconos::modeling::LagrangianDS::display(bool brief) const {
  std::cout << "=====> Lagrangian System display (number: " << number_ << ").\n";
  std::cout << "- ndof_ : " << ndof_ << "\n";
  std::cout << "- q " << *state_q_[0] << "\n";
  std::cout << "- q0 " << *q0_view_ << "\n";
  std::cout << "- velocity " << *state_q_[1] << "\n";
  std::cout << "- acceleration \n";
  if (state_q_[2])
    std::cout << *state_q_[2] << "\n";
  else
    std::cout << "-> nullptr\n";

  std::cout << "- v0 " << *velocity0_view_ << "\n";

  std::cout << "- p[0] \n";
  if (p_[0])
    std::cout << *p_[0] << "\n";
  else
    std::cout << "-> nullptr\n";
  std::cout << "- p[1] " << *p_[1] << "\n";
  if (p_[2])
    std::cout << *p_[2] << "\n";
  else
    std::cout << "-> nullptr\n";

  if (!brief) {
    if (mass_view_) std::cout << "- Mass " << mass_view_ << "\n";

    std::cout << "- Forces \n";
    if (totalForces_)
      std::cout << *totalForces_ << "\n";
    else
      std::cout << "-> nullptr\n";

    if (fint_) std::cout << "- FInt " << *fint_ << "\n";

    std::cout << "- jacobian of TotalForces over q \n";
    if (jacobianTotalForcesOver_q_)
      std::cout << *jacobianTotalForcesOver_q_ << "\n";
    else
      std::cout << "-> nullptr\n";

    if (hasJacobianFintOver_q_)
      std::cout << "- jacobian of fint over q " << jacobianFintOver_q_view_ << "\n";

    std::cout << "- jacobian of TotalForces over velocity \n";
    if (jacobianTotalForcesOver_velocity_)
      std::cout << *jacobianTotalForcesOver_velocity_ << "\n";
    else
      std::cout << "-> nullptr\n";
  }

  std::cout << "===================================== \n";
}

// --- Functions for memory handling ---
void siconos::modeling::LagrangianDS::initMemory(unsigned int steps) {
  DEBUG_PRINTF(
      "siconos::modeling::LagrangianDS::initMemory(unsigned int steps) with "
      "steps = %i\n",
      steps);
  if (steps == 0)
    std::cout << "Warning : LagragianDS::initMemory with size equal to zero\n";
  else {
    qMemory_.setMemorySize(steps, ndof_);
    velocityMemory_.setMemorySize(steps, ndof_);
    totalForcesMemory_.setMemorySize(steps, ndof_);
    pMemory_.resize(3);

    // TODO : initMemory in graph + f(OSI/level)
    for (unsigned int level = 0; level < 3; ++level) {
      if (pMemory_[level].size() == 0) pMemory_[level].setMemorySize(steps, ndof_);
    }

    // swapInMemory();
  }
}

void siconos::modeling::LagrangianDS::swapInMemory() {
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

void siconos::modeling::LagrangianDS::resetAllNonSmoothParts() {
  if (p_[0]) p_[0]->setZero();
  if (p_[1]) p_[1]->setZero();
  if (p_[2]) p_[2]->setZero();
}

void siconos::modeling::LagrangianDS::resetNonSmoothPart(unsigned int level) {
  if (level < siconos::internal::LEVELMAX)
    if (p_[level]) p_[level]->setZero();
}

void siconos::modeling::LagrangianDS::computePostImpactVelocity() {
  // When this function is call, q[1] is supposed to be pre-impact velocity.
  // We solve M(v+ - v-) = p - The result is saved in(place of) p[1].
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::computePostImpactV()\n");
  siconos::algebra::SiconosVector tmp(*p_[1]);
  if (LUMass_) tmp = LUMass_->solve(tmp);  // we assume mass/lumass are uptodate?
  *state_q_[1] += tmp;                     // v+ = v- + p
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::computePostImpactV() END \n");
}

double siconos::modeling::LagrangianDS::computeKineticEnergy() {
  DEBUG_BEGIN("NewtonEulerDS::computeKineticEnergy()\n");
  double K;
  if (mass_view_)
    K = 0.5 * state_q_[1]->dot(*mass_view_ * *state_q_[1]);
  else
    K = 0.5 * state_q_[1]->dot(*state_q_[1]);

  DEBUG_PRINTF("Kinetic Energy = %e\n", K);
  DEBUG_END("siconos::modeling::LagrangianDS::computeKineticEnergy()\n");
  return K;
}

void siconos::modeling::LagrangianDS::acceptSP(
    std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const {
  tourist->visit(*this);
}
