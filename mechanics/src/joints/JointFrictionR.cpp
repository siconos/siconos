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
/*! \file JointFrictionR.cpp
 */

#include "JointFrictionR.hpp"

#include "BlockVector.hpp"
#include "NewtonEulerJointR.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

/** Initialize a joint friction for a common case: a single axis with a
 * single friction, either positive or negative. For use with
 * NewtonImpactNSL. */
siconos::joints::JointFrictionR::JointFrictionR(std::shared_ptr<NewtonEulerJointR> joint,
                                                siconos::algebra::Index axis)
    : NewtonEulerR{},
      _joint(joint),
      _axis(std::make_shared<std::vector<siconos::algebra::Index>>()) {
  _axis->push_back(axis);
  _axisMin = axis;
  _axisMax = axis;
  assert(siconos::algebra::to_index(_axisMax - _axisMin + 1) <= _joint->numberOfDoF());
}

/** Initialize a multidimensional joint friction, e.g. the cone friction on
 * a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
siconos::joints::JointFrictionR::JointFrictionR(
    std::shared_ptr<NewtonEulerJointR> joint,
    std::shared_ptr<std::vector<siconos::algebra::Index>> axes)
    : NewtonEulerR{}, _joint(joint), _axis(axes) {
  if (axes) {
    _axisMin = 100;
    _axisMax = 0;
    for (size_t i = 0; i < _axis->size(); i++) {
      if ((*_axis)[i] > _axisMax) _axisMax = (*_axis)[i];
      if ((*_axis)[i] < _axisMin) _axisMin = (*_axis)[i];
    }
  } else {
    _axisMin = _axisMax = 0;
    _axis = std::make_shared<std::vector<siconos::algebra::Index>>();
    _axis->push_back(0);
  }

  assert(siconos::algebra::to_index(_axisMax - _axisMin + 1) <= _joint->numberOfDoF());
}

void siconos::joints::JointFrictionR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  // Velocity-level constraint, no position-level h
  y.setZero();
}

void siconos::joints::JointFrictionR::computeH_NE_(double time,
                                                   siconos::modeling::Interaction& inter,
                                                   const siconos::algebra::BlockVector& q0) {
  auto n = siconos::algebra::to_index(_axisMax - _axisMin + 1);
  assert(n == 1);  // For now, multi-axis support TODO

  if (!jacobianhOver_q_Tmp ||
      !(jacobianhOver_q_Tmp->cols() == q0.size() && jacobianhOver_q_Tmp->rows() == n)) {
    jacobianhOver_q_Tmp = std::make_shared<siconos::algebra::SiconosMatrix>(n, q0.size());
  }

  // Compute the jacobian for the required range of axes
  if (q0.numberOfBlocks() > 1) {
    _joint->computeJachqDoF(inter, *q0.vector(0), *q0.vector(1), *jacobianhOver_q_Tmp,
                            _axisMin);
  } else {
    _joint->computeJachqDoF(inter, *q0.vector(0), std::nullopt, *jacobianhOver_q_Tmp,
                            _axisMin);
  }

  // Copy indicated axes into the friction jacobian, negative and positive sides
  // NOTE trying ==1 using Relay, maybe don't need LCP formulation
  assert(H_NE_view_->rows() == 1);
  for (siconos::algebra::Index j = 0; j < H_NE_view_->cols(); j++) {
    (*H_NE_view_)(0, j) = -(*jacobianhOver_q_Tmp)((*_axis)[0] - _axisMin, j);
  }
}

size_t siconos::joints::JointFrictionR::numberOfConstraints() const { return _axis->size(); }
siconos::algebra::Index siconos::joints::JointFrictionR::axis(size_t _index) {
  return _axis->at(_index);
}

size_t siconos::joints::JointFrictionR::numberOfAxes() { return _axis->size(); }
