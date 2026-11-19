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

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "NewtonEulerR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "NewtonEulerDS.hpp"  // computeT ...
#include "SiconosAlgebraAddons.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Tools.hpp"
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::modeling::NewtonEulerR::initialize(Interaction& inter) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerR::initialize(Interaction& inter)\n");

  auto ySize = inter.dimension();
  auto Hcols = 7;  // 7 * number of DS in the interaction
  if (inter.has2Bodies()) Hcols *= 2;
  // setHMatrix is the only way to set H. So if it's required, storage
  // is already allocated at this point.

  // No relation without H
  // - either a shared view with an external object (if set with setConstantConstraintJacobian)
  // - or memory allocated in derived class (constructor or in initialize)
  if (!constraintJacobian_view_) {
    // It means that setConstantConstraintJacobian has not been call and that internal storage
    // is required.
    constraintJacobian_storage_.resize(ySize, Hcols);
    constraintJacobian_view_ = std::make_shared<siconos::algebra::MapType>(
        constraintJacobian_storage_.data(), ySize, Hcols);
    constraintJacobian_view_->setZero();
  }

  // JacobianhOver_q_dot is optional. Allocation only if a user-defined function has been set.
  if (computeConstraintJacobianDot_) {
    // True only if setComputeContactMatrixDotFunction has been called
    // Ensure that memory is properly allocated for H_dot
    if (!constraintJacobianDot_) {
      constraintJacobianDot_ = std::make_shared<siconos::algebra::SiconosMatrix>(ySize, Hcols);
      constraintJacobianDot_->setZero();
    }
  }
  auto qSize = inter.getSizeOfDS();  // sum of considered DS sizes
  if (!contactMatrix_)
    contactMatrix_ = std::make_shared<siconos::algebra::SiconosMatrix>(ySize, qSize);
  contactMatrix_->setZero();

  // Allocate internal buffer, used to save T(q)
  T_buffer_.setZero();
  T_buffer_(0, 0) = 1.;
  T_buffer_(1, 1) = 1.0;
  T_buffer_(2, 2) = 1.0;

  const auto& ds_vars = inter.read_dynamical_systems_variables();
  contactForce_.resize(ds_vars[siconos::tools::enum_to_index(ds_var::p1)]->size());
  contactForce_.setZero();

  checkSize(inter);

  DEBUG_END("siconos::modeling::NewtonEulerR::initialize(Interaction& inter)\n");
}

void siconos::modeling::NewtonEulerR::checkSize(const Interaction& inter) const {
  [[maybe_unused]] auto ySize = inter.dimension();
  [[maybe_unused]] auto qSize = inter.getSizeOfDS();  // sum of considered DS sizes
  [[maybe_unused]] auto Hcols = 7 * (qSize / 6);      // 7 * number of DS in the interaction

  if (constraintJacobian_view_) {
    assert(constraintJacobian_view_->rows() == ySize);
    assert(constraintJacobian_view_->cols() == Hcols);
  }

  if (constraintJacobianDot_) {
    assert(constraintJacobianDot_->rows() == ySize);
    assert(constraintJacobianDot_->cols() == Hcols);
  }

  if (haseVector_) assert(eVector_view_->size() == ySize);

  if (contactMatrix_) {
    assert(contactMatrix_->rows() == ySize);
    assert(contactMatrix_->cols() == qSize);
  }
}

void siconos::modeling::NewtonEulerR::setConstantConstraintJacobian(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  // Warning: mat. dims are not checked since we have no clue here regarding the size of the
  // Interaction.

  constraintJacobian_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
  hasConstantConstraintJacobian_ = true;
}

void siconos::modeling::NewtonEulerR::setConstanteVector(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  // Warning: mat. dims are not checked since we have no clue here regarding the size of the
  // Interaction.

  eVector_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
  haseVector_ = true;
}

void siconos::modeling::NewtonEulerR::setComputeConstraintJacobianDotFunction(
    const siconos::modeling::func_prototypes::FunctionBVBV_M& func) {
  // constraintJacobianDot_ is allocated in initialize, or at first use if the function is set
  // later.
  computeConstraintJacobianDot_ = func;
}

void siconos::modeling::NewtonEulerR::allocateConstraintJacobianDot_() {
  if (!constraintJacobianDot_) {
    assert(constraintJacobian_view_ &&
           "NewtonEulerR: G must be allocated before constraintJacobianDot");
    constraintJacobianDot_ = std::make_shared<siconos::algebra::SiconosMatrix>(
        constraintJacobian_view_->rows(), constraintJacobian_view_->cols());
    constraintJacobianDot_->setZero();
  }
}

void siconos::modeling::NewtonEulerR::computeConstraintJacobianDot(
    const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot) {
  if (!computeConstraintJacobianDotLegacy_) {  // dispatch to a class overriding the deprecated
                                               // name
    detail::LegacyCallGuard guard(computeConstraintJacobianDotLegacy_);
    SICONOS_IGNORE_DEPRECATED_BEGIN
    computeH_NE_dot(q, qdot);
    SICONOS_IGNORE_DEPRECATED_END
    return;
  }
  if (computeConstraintJacobianDot_) {
    allocateConstraintJacobianDot_();
    computeConstraintJacobianDot_(q, qdot, *constraintJacobianDot_);
  }
}

void siconos::modeling::NewtonEulerR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  // Default implementation: convert to BlockVector and call existing method
  y.noalias() = constraintJacobian_view_->leftCols(7) * q1;
  if (q2) {
    y.noalias() += constraintJacobian_view_->rightCols(7) * *q2;
  }
  if (haseVector_) y += *eVector_view_;
}

void siconos::modeling::NewtonEulerR::computeOutput(
    double time, Interaction& inter, siconos::algebra::blocks::size_type derivativeNumber) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerR::computeOutput(...)\n");
  DEBUG_PRINTF("with time = %f and derivativeNumber = %i starts\n", time, derivativeNumber);

  const auto& ds_vars = inter.read_dynamical_systems_variables();
  auto& y = *inter.y(derivativeNumber);
  auto q = ds_vars[siconos::tools::enum_to_index(ds_var::q0)];
  if (derivativeNumber == 0) {
    if (q->numberOfBlocks() == 2)
      computeh(*q->vector(0), *q->vector(1), y);
    else
      computeh(*q->vector(0), std::nullopt, y);

  } else {
    if (derivativeNumber == 1) {
      assert(contactMatrix_);
      assert(ds_vars[siconos::tools::enum_to_index(ds_var::velocity)]);
      DEBUG_EXPR(siconos::algebra::print(*contactMatrix_););
      DEBUG_EXPR(siconos::algebra::print(
          (*ds_vars[siconos::tools::enum_to_index(ds_var::velocity)])));
      ;
      /* \warning  V.A. 15/04/2016
       * We decide finally not to update the Jacobian there. To be discussed
       */
      //      if (!hasConstantConstraintJacobian_) computeConstraintJacobian(time, inter, q);
      //      computeContactMatrix(inter, q);
      auto& vel = *ds_vars[siconos::tools::enum_to_index(ds_var::velocity)];
      siconos::algebra::matrixBlockVector_prod(*contactMatrix_, vel, y);
      DEBUG_EXPR(siconos::algebra::print(y););
    } else if (derivativeNumber == 2) {
      THROW_EXCEPTION(
          "Call siconos::modeling::NewtonEulerR::computeOutput(..., derivativeNumber=2) is "
          "not allowed.");
    } else
      THROW_EXCEPTION(
          "siconos::modeling::NewtonEulerR::computeOutput(..., derivativeNumber) "
          "derivativeNumber out of range or not yet implemented.");
  }
  DEBUG_END("siconos::modeling::NewtonEulerR::computeOutput(...)\n");
}

void siconos::modeling::NewtonEulerR::computeInput(double time, Interaction& inter,
                                                   siconos::algebra::blocks::size_type level) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerR::computeInput(...)\n")
  DEBUG_PRINTF("with time = %f and level = %i starts\n", time, level);
  DEBUG_EXPR(printf("interaction %p\n", &inter););
  DEBUG_EXPR(siconos::algebra::print(inter););
  const auto& ds_vars = inter.read_dynamical_systems_variables();

  // get lambda of the concerned interaction
  const auto& lambda = *inter.lambda(level);

  DEBUG_EXPR(siconos::algebra::print(lambda););
  DEBUG_EXPR(
      siconos::algebra::print(*ds_vars[siconos::tools::enum_to_index(ds_var::p0) + level]););
  if (level == 1 ||
      level == 2) /* \warning : we assume that ContactForce is given by lambda[level] */
  {
    contactForce_.noalias() = contactMatrix_->transpose() * lambda;
    DEBUG_EXPR(siconos::algebra::print(contactForce_););

    siconos::algebra::transposeMatrixVector_prod_toBlock(
        lambda, *contactMatrix_, *ds_vars[siconos::tools::enum_to_index(ds_var::p0) + level],
        false);

    DEBUG_EXPR(siconos::algebra::print(*contactMatrix_););
    DEBUG_EXPR(
        siconos::algebra::print(*ds_vars[siconos::tools::enum_to_index(ds_var::p0) + level]););
  } else if (level == 0) {
    siconos::algebra::transposeMatrixVector_prod_toBlock(
        lambda, *constraintJacobian_view_,
        *ds_vars[siconos::tools::enum_to_index(ds_var::p0) + level], false);
  } else
    THROW_EXCEPTION(
        "siconos::modeling::NewtonEulerR::computeInput(double time, Interaction& inter, "
        "InteractionProperties& interProp, siconos::algebra::blocks::size_type level)  not "
        "yet implemented for "
        "level "
        "> 1");
  DEBUG_END("siconos::modeling::NewtonEulerR::computeInput(...)\n");
}

void siconos::modeling::NewtonEulerR::computeConstraintJacobian(
    double time, Interaction& inter, const siconos::algebra::BlockVector& q0) {
  // Default: G is constant. Dispatch to a class overriding the deprecated name.
  if (!computeConstraintJacobianLegacy_) {
    detail::LegacyCallGuard guard(computeConstraintJacobianLegacy_);
    SICONOS_IGNORE_DEPRECATED_BEGIN
    computeH_NE_(time, inter, q0);
    SICONOS_IGNORE_DEPRECATED_END
  }
}

void siconos::modeling::NewtonEulerR::computeH_NE_(double time, Interaction& inter,
                                                   const siconos::algebra::BlockVector& q0) {
  computeConstraintJacobian(time, inter, q0);
}

void siconos::modeling::NewtonEulerR::computeContactMatrix(
    const Interaction& inter, const siconos::algebra::BlockVector& q0) {
  if (!computeContactMatrixLegacy_) {  // dispatch to a class overriding the deprecated name
    detail::LegacyCallGuard guard(computeContactMatrixLegacy_);
    SICONOS_IGNORE_DEPRECATED_BEGIN
    computeH_NE_prod_T(inter, q0);
    SICONOS_IGNORE_DEPRECATED_END
    return;
  }
  computeContactMatrixFromConstraintJacobian_(q0);
}

void siconos::modeling::NewtonEulerR::computeH_NE_prod_T(
    const Interaction& inter, const siconos::algebra::BlockVector& q0) {
  computeContactMatrix(inter, q0);
}

void siconos::modeling::NewtonEulerR::computeH_NE_dot(
    const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot) {
  computeConstraintJacobianDot(q, qdot);
}

void siconos::modeling::NewtonEulerR::computeContactMatrixFromConstraintJacobian_(
    const siconos::algebra::BlockVector& q0) {
  // H = G T(q), body by body
  for (Eigen::Index body = 0; body < static_cast<Eigen::Index>(q0.numberOfBlocks()); ++body) {
    siconos::modeling::newton_euler::computeT(*q0.vector(body), T_buffer_);
    contactMatrix_->middleCols(6 * body, 6) =
        constraintJacobian_view_->middleCols(7 * body, 7) * T_buffer_;
  }
  DEBUG_EXPR(siconos::algebra::print(*contactMatrix_););
}

void siconos::modeling::NewtonEulerR::computeJach(double time, Interaction& inter) {
  DEBUG_BEGIN("NewtonEulerR::computeJachq(double time, Interaction& inter, ...) \n");
  DEBUG_PRINTF("with time =  %f\n", time);
  DEBUG_PRINTF("with inter =  %p\n", &inter);

  const auto& ds_vars = inter.read_dynamical_systems_variables();

  if (!hasConstantConstraintJacobian_)
    computeConstraintJacobian(time, inter,
                              *ds_vars[siconos::tools::enum_to_index(ds_var::q0)]);
  computeContactMatrix(inter, *ds_vars[siconos::tools::enum_to_index(ds_var::q0)]);

  // computeJachqDot(time, inter); // This is not needed here
  // computeDotJachq(time, inter);
  DEBUG_END("NewtonEulerR::computeJachq(double time, Interaction& inter, ...) \n");
}

void siconos::modeling::NewtonEulerR::computeContactMatrixDotV(
    double time, Interaction& inter, std::shared_ptr<DynamicalSystem> ds1,
    std::shared_ptr<DynamicalSystem> ds2) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerR::computeContactMatrixDotV\n");

  if (!computeConstraintJacobianDot_)
    THROW_EXCEPTION(
        "NewtonEulerR::computeContactMatrixDotV: no function to compute "
        "H_NE_dot, call setComputeConstraintJacobianDotFunction.");

  const auto& ds_vars = inter.read_dynamical_systems_variables();
  const auto& q = *ds_vars[siconos::tools::enum_to_index(ds_var::q0)];
  const auto& dotq = *ds_vars[siconos::tools::enum_to_index(ds_var::dotq)];
  const auto& v = *ds_vars[siconos::tools::enum_to_index(ds_var::velocity)];

  // Jdot * dotq (dotq is assumed to be up to date)
  allocateConstraintJacobianDot_();
  computeConstraintJacobianDot_(q, dotq, *constraintJacobianDot_);
  DEBUG_EXPR(siconos::algebra::print(*constraintJacobianDot_););
  contactMatrixDotV_.resize(constraintJacobianDot_->rows());
  siconos::algebra::matrixBlockVector_prod(*constraintJacobianDot_, dotq, contactMatrixDotV_,
                                           true);

  // + J * Tdot * v, body by body
  auto addJTdotv = [&](DynamicalSystem& ds, Eigen::Index body) {
    auto& neds = dynamic_cast<NewtonEulerDS&>(ds);
    neds.computeTdot();
    contactMatrixDotV_.noalias() +=
        constraintJacobian_view_->middleCols(7 * body, 7) * (neds.Tdot() * *v.vector(body));
  };
  addJTdotv(*ds1, 0);
  if (inter.has2Bodies()) addJTdotv(*ds2, 1);

  DEBUG_EXPR(siconos::algebra::print(contactMatrixDotV_););
  DEBUG_END("siconos::modeling::NewtonEulerR::computeContactMatrixDotV\n");
}

void siconos::modeling::NewtonEulerR::allocate_read_dynamical_systems_var_vectors(
    std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& ds_vars,
    modeling::DynamicalSystem& ds1, modeling::DynamicalSystem& ds2) const {
  ds_vars.resize(siconos::tools::enum_to_index(NewtonEulerR::ds_var::size));
  ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::q0)] =
      std::make_shared<siconos::algebra::BlockVector>();  // displacement
  ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::velocity)] =
      std::make_shared<siconos::algebra::BlockVector>();  // velocity
  ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::dotq)] =
      std::make_shared<siconos::algebra::BlockVector>();  // qdot
  ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::p0)] =
      std::make_shared<siconos::algebra::BlockVector>();
  ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::p1)] =
      std::make_shared<siconos::algebra::BlockVector>();
  ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::p2)] =
      std::make_shared<siconos::algebra::BlockVector>();
  bool has_two_ds = &ds1 != &ds2;
  auto* neds1 = dynamic_cast<NewtonEulerDS*>(&ds1);

  ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::q0)]->insertPtr(neds1->q());
  ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::velocity)]->insertPtr(
      neds1->twist());
  //  ds_vars[NewtonEulerR::ds_var::deltaq]->insertPtr(deltaq());
  ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::dotq)]->insertPtr(neds1->dotq());
  if (neds1->p(0))
    ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::p0)]->insertPtr(neds1->p(0));
  if (neds1->p(1))
    ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::p1)]->insertPtr(neds1->p(1));
  if (neds1->p(2))
    ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::p2)]->insertPtr(neds1->p(2));

  if (has_two_ds) {
    auto* neds2 = dynamic_cast<NewtonEulerDS*>(&ds2);

    ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::q0)]->insertPtr(neds2->q());
    ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::velocity)]->insertPtr(
        neds2->twist());
    //  ds_vars[NewtonEulerR::ds_var::deltaq]->insertPtr(deltaq());
    ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::dotq)]->insertPtr(
        neds2->dotq());
    if (neds2->p(0))
      ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::p0)]->insertPtr(neds2->p(0));
    if (neds2->p(1))
      ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::p1)]->insertPtr(neds2->p(1));
    if (neds2->p(2))
      ds_vars[siconos::tools::enum_to_index(NewtonEulerR::ds_var::p2)]->insertPtr(neds2->p(2));
  }
}
