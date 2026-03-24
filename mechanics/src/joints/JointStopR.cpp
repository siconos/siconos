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
/*! \file JointStopR.cpp

*/

#include "JointStopR.hpp"

#include <vector>

#include "BlockVector.hpp"
#include "NewtonEulerJointR.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #include "siconos_debug.h"

/** Initialize a joint stop for a common case: a single axis with a
 * single stop, either positive or negative. For use with
 * NewtonImpactNSL. */
siconos::joints::JointStopR::JointStopR(std::shared_ptr<NewtonEulerJointR> joint, double pos,
                                        bool dir, siconos::algebra::Index axis)
    : NewtonEulerR{}, joint_(joint) {
  axis_.push_back(axis);
  pos_.resize(1);
  pos_(0) = pos;
  dir_.resize(1);
  dir_(0) = dir ? -1 : 1;
  axisMin_ = axis;
  axisMax_ = axis;
  assert(siconos::algebra::to_index(axisMax_ - axisMin_ + 1) <= joint_->numberOfDoF());
}

/** Initialize a multidimensional joint stop, e.g. the cone stop on
 * a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
siconos::joints::JointStopR::JointStopR(
    std::shared_ptr<NewtonEulerJointR> joint,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& pos,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& dir,
    const std::vector<siconos::algebra::Index>& axes)
    : NewtonEulerR{}, joint_{joint}, axis_{axes}, pos_{pos}, dir_{dir} {
  axisMin_ = 100;
  axisMax_ = 0;

  for (size_t i = 0; i < axis_.size(); i++) {
    if (axis_[i] > axisMax_) axisMax_ = axis_[i];
    if (axis_[i] < axisMin_) axisMin_ = axis_[i];
  }
  assert(siconos::algebra::to_index(axisMax_ - axisMin_ + 1) <= joint_->numberOfDoF());
}

void siconos::joints::JointStopR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  // Common cases optimisation
  bool case_onestop = y.size() == 1;
  bool case_posneg = y.size() == 2 && axis_[0] == axis_[1];
  if (case_onestop || case_posneg) {
    if (q2)
      joint_->computehDoF(q1, q2, y, axis_[0]);
    else
      joint_->computehDoF(q1, std::nullopt, y, axis_[0]);
    y(0) = (y(0) - pos_(0)) * dir_(0);
    if (case_posneg) y(1) = (y(0) - pos_(1)) * dir_(1);
    return;
  }

  // Get h for each relevant axis
  siconos::algebra::SiconosVector tmp_y(axisMax_ - axisMin_ + 1);
  tmp_y.setZero();
  if (q2)
    joint_->computehDoF(q1, q2, tmp_y, axisMin_);
  else
    joint_->computehDoF(q1, std::nullopt, tmp_y, axisMin_);

  // Copy and scale each stop for its axis/position/direction
  auto size_y = y.size();
  y = tmp_y.head(size_y) - pos_.head(size_y).cwiseProduct(dir_.head(size_y));
}

void siconos::joints::JointStopR::computeH_NE_(double time,
                                               siconos::modeling::Interaction& inter,
                                               const siconos::algebra::BlockVector& q0) {
  auto n = siconos::algebra::to_index(axisMax_ - axisMin_ + 1);

  if (jacobian_buffer_.cols() != q0.size() || jacobian_buffer_.rows() != n) {
    jacobian_buffer_.resize(n, q0.size());
  }
  jacobian_buffer_.setZero();

  // Compute the jacobian for the required range of axes
  if (q0.numberOfBlocks() > 1) {
    joint_->computeJachqDoF(inter, *q0.vector(0), *q0.vector(1), jacobian_buffer_, axisMin_);
  } else {
    joint_->computeJachqDoF(inter, *q0.vector(0), std::nullopt, jacobian_buffer_, axisMin_);
  }

  // Copy indicated axes into the stop jacobian, possibly flipped for negative stops
  for (siconos::algebra::Index i = 0; i < H_NE_view_->rows(); i++)
    for (siconos::algebra::Index j = 0; j < H_NE_view_->cols(); j++)
      H_NE_view_->setValue(
          i, j,
          (jacobian_buffer_)(axis_[siconos::algebra::to_unsigned<size_t>(i)] - axisMin_, j) *
              dir_(i));
}

std::size_t siconos::joints::JointStopR::numberOfConstraints() const { return axis_.size(); }

auto siconos::joints::JointStopR::axis(size_t _index) { return axis_.at(_index); }

double siconos::joints::JointStopR::position(siconos::algebra::Index _index) {
  return pos_(_index);
}

double siconos::joints::JointStopR::direction(siconos::algebra::Index _index) {
  return dir_(_index);
}

auto siconos::joints::JointStopR::numberOfAxes() { return axis_.size(); }
