/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*! \file CouplerJointRR.cpp
 */

#include "CouplerJointR.hpp"

#include <Interaction.hpp>
#include <optional>

#include "BlockVector.hpp"
#include "NewtonEulerDS.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #include <boost/math/quaternion.hpp>
// #include <cfloat>
// #include <iostream>

// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #include "siconos_debug.h"

/** Initialize a coupler joint. For use with EqualityConditionNSL to
 * bind two degrees of freedom into one by a ratio and offset. */
// Warning : shared_ptr --> can't be wrapped with pb11
siconos::joints::CouplerJointR::CouplerJointR(
    std::shared_ptr<NewtonEulerJointR> joint1, siconos::algebra::Index dof1,
    std::shared_ptr<NewtonEulerJointR> joint2, siconos::algebra::Index dof2, double ratio,
    siconos::algebra::Index ref1_index, std::shared_ptr<siconos::algebra::SiconosVector> ref1,
    siconos::algebra::Index ref2_index, std::shared_ptr<siconos::algebra::SiconosVector> ref2)
    : NewtonEulerJointR{},
      joint1_(joint1),
      joint2_(joint2),
      ref1_{ref1},
      ref2_{ref2},
      dof1_(dof1),
      dof2_(dof2),
      ref1_index_(ref1_index),
      ref2_index_(ref2_index),
      ratio_(ratio),
      offset_(0.0) {
  if (ref1_) assert(ref1_->size() == 7);
  if (ref2_) assert(ref2_->size() == 7);
  assert(dof1_ < joint1_->numberOfDoF());
  assert(dof2_ < joint2_->numberOfDoF());
}

/* A constructor taking a DS exists just because it's hard to pass
 * ds->q() through Python without it automatically converting to numpy
 * array and back, which messes up the shared_ptr reference.
 *
 * NOTE that using q() as the reference is not quite right, in fact it
 * should be using the reference ds's temporary work vector in order
 * to perform correctly in the Newton loop. (TODO) */
siconos::joints::CouplerJointR::CouplerJointR(
    std::shared_ptr<NewtonEulerJointR> joint1, siconos::algebra::Index dof1,
    std::shared_ptr<NewtonEulerJointR> joint2, siconos::algebra::Index dof2, double ratio,
    std::shared_ptr<siconos::modeling::NewtonEulerDS> refds1,
    siconos::algebra::Index ref1_index,
    std::shared_ptr<siconos::modeling::NewtonEulerDS> refds2,
    siconos::algebra::Index ref2_index)
    : NewtonEulerJointR(),
      joint1_(joint1),
      joint2_(joint2),
      dof1_(dof1),
      dof2_(dof2),
      ref1_index_(ref1_index),
      ref2_index_(ref2_index),
      ratio_(ratio) {
  if (refds1) ref1_ = refds1->q();
  if (refds2) ref2_ = refds2->q();
  assert(dof1_ < joint1_->numberOfDoF());
  assert(dof2_ < joint2_->numberOfDoF());
}

void siconos::joints::CouplerJointR::build_vectors(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2) {
  if (ref1_)  // v1_=[ref1, q1] or [q1, ref1]
  {
    v1_[ref1_index_].emplace(ref1_->data());
    v1_[1 - ref1_index_].emplace(q1.data());

    // v1_[ref1_index_] = Eigen::Map<const siconos::algebra::SiconosVector7>(ref1_->data());
    // //        Eigen::Ref<const siconos::algebra::SiconosVector7>(*ref1_);  // view on ref1
    // v1_[1 - ref1_index_] = Eigen::Map<const siconos::algebra::SiconosVector7>(q1.data());
  } else {  // no ref_1
    v1_[0].emplace(q1.data());
    //    v1_[0] = Eigen::Map<const siconos::algebra::SiconosVector7>(q1.data());

    if (q2) {  // v1_= [q1, q2]
               //      v1_[1] = Eigen::Map<const siconos::algebra::SiconosVector7>(q2->data());
      v1_[1].emplace(q2->data());

    } else {  // v1_ = [q1, nullopt]
      v1_[1] = std::nullopt;
    }
  }

  if (q2) {
    if (ref2_) {
      v2_[ref2_index_].emplace(ref2_->data());
      v2_[1 - ref2_index_].emplace(q2->data());

      // v2_[ref2_index_] = Eigen::Map<const siconos::algebra::SiconosVector7>(ref2_->data());
      // v2_[1 - ref2_index_] = Eigen::Map<const siconos::algebra::SiconosVector7>(q2->data());

    } else {
      v2_[0].emplace(q1.data());
      v2_[1].emplace(q2->data());

      // v2_[0] = Eigen::Map<const siconos::algebra::SiconosVector7>(q1.data());

      // v2_[1] = Eigen::Map<const siconos::algebra::SiconosVector7>(q2->data());
    }
  } else {
    for (int i = 0; i < 2; ++i) {
      if (v1_[i])
        v2_[i].emplace(v1_[i]->data());
      else
        v2_[i] = std::nullopt;
    }

    //    v2_ = v1_;
  }
}

// void siconos::joints::CouplerJointR::resolveVectors(
//     const siconos::algebra::SiconosVector7* q1, const siconos::algebra::SiconosVector7* q2,
//     const siconos::algebra::SiconosVector7*& v1_1,
//     const siconos::algebra::SiconosVector7*& v1_2,
//     const siconos::algebra::SiconosVector7*& v2_1,
//     const siconos::algebra::SiconosVector7*& v2_2) const {
//   // Case 1: no references, q2 present → coupling q1/q2
//   if (!ref1_ && !ref2_ && q2) {
//     v1_1 = q1;
//     v1_2 = q2;
//     v2_1 = q1;
//     v2_2 = q2;
//   }
//   // Case 2: no q2, maybe ref1
//   else if (!q2) {
//     if (ref1_) {
//       v1_1 = (ref1_index_ == 0) ? ref1_.get() : q1;
//       v1_2 = (ref1_index_ == 0) ? q1 : ref1_.get();
//       v2_1 = v1_1;
//       v2_2 = v1_2;
//     } else {
//       v1_1 = q1;
//       v1_2 = nullptr;
//       v2_1 = q1;
//       v2_2 = nullptr;
//     }
//   }
//   // Case 3: q2 present, potentially ref1 and ref2
//   else {
//     // v1
//     if (ref1_) {
//       v1_1 = (ref1_index_ == 0) ? ref1_.get() : q1;
//       v1_2 = (ref1_index_ == 0) ? q1 : ref1_.get();
//     } else {
//       v1_1 = q1;
//       v1_2 = q2;
//     }

//     // v2
//     if (ref2_) {
//       v2_1 = (ref2_index_ == 0) ? ref2_.get() : q2;
//       v2_2 = (ref2_index_ == 0) ? q2 : ref2_.get();
//     } else {
//       v2_1 = q1;
//       v2_2 = q2;
//     }
//   }
// }

void siconos::joints::CouplerJointR::setBasePositions(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2) {
  offset_ = 0;
  siconos::algebra::SiconosVector y(1);
  y << 0.;
  computeh(q1, q2, y);
  offset_ = -y(0);
}

void siconos::joints::CouplerJointR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  siconos::algebra::SiconosVector y1(y), y2(y);
  build_vectors(q1, q2);
  // Compute hDoF for both joints

  joint1_->computehDoF(*v1_[0], v1_[1], y1, dof1_);
  joint2_->computehDoF(*v2_[0], v2_[1], y2, dof2_);
  // Constraint is the linear relation between them
  y(0) = y2(0) - y1(0) * ratio_ + offset_;
}

void siconos::joints::CouplerJointR::computeH_NE_(double time,
                                                  siconos::modeling::Interaction& inter,
                                                  const siconos::algebra::BlockVector& q0) {
  siconos::algebra::SiconosDenseMatrix jachq1{1, q0.size()};
  siconos::algebra::SiconosDenseMatrix jachq2{1, q0.size()};
  jachq1.setZero();
  jachq2.setZero();
  // Get jacobians for the implicated degrees of freedom

  if (q0.numberOfBlocks() == 2)
    build_vectors(*q0.vector(0), *q0.vector(1));
  else
    build_vectors(*q0.vector(0), std::nullopt);

  joint1_->computeJachqDoF(inter, *v1_[0], v1_[1], jachq1, dof1_);
  joint2_->computeJachqDoF(inter, *v2_[0], v2_[1], jachq2, dof2_);
  // Constraint is the linear relation between them
  (*H_NE_view_) = jachq2 - ratio_ * jachq1;
}

void siconos::joints::CouplerJointR::computehDoF(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y, siconos::algebra::Index axis) {
  // The DoF of the constraint is the same as the constraint itself
  assert(axis == 0);
  computeh(q1, q2, y);
}

void siconos::joints::CouplerJointR::computeJachqDoF(
    siconos::modeling::Interaction& inter,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> jachq, siconos::algebra::Index axis) {
  // The Jacobian of the DoF of the constraint is the same as the
  // Jacobian of the constraint itself. (Same as computeJacobianhOver_q(), but
  // don't store result in member object.)
  assert(axis == 0);

  build_vectors(q1, q2);
  auto cols = q1.size() + (q2 ? q2->size() : 0);
  siconos::algebra::SiconosDenseMatrix jachq1{1, cols};
  siconos::algebra::SiconosDenseMatrix jachq2{1, cols};
  jachq1.setZero();
  jachq2.setZero();
  joint1_->computeJachqDoF(inter, *v1_[0], v1_[1], jachq1, dof1_);
  joint2_->computeJachqDoF(inter, *v2_[0], v2_[1], jachq2, dof2_);
  jachq = jachq2 - ratio_ * jachq1;
}
