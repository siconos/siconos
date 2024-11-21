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
#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

/** Initialize a joint friction for a common case: a single axis with a
 * single friction, either positive or negative. For use with
 * NewtonImpactNSL. */
siconos::joints::JointFrictionR::JointFrictionR(std::shared_ptr<NewtonEulerJointR> joint,
                                                unsigned int axis)
    : NewtonEulerR{}, _joint(joint), _axis(std::make_shared<std::vector<unsigned int>>()) {
  _axis->push_back(axis);
  _axisMin = axis;
  _axisMax = axis;
  assert((_axisMax - _axisMin + 1) <= _joint->numberOfDoF());
}

/** Initialize a multidimensional joint friction, e.g. the cone friction on
 * a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
siconos::joints::JointFrictionR::JointFrictionR(
    std::shared_ptr<NewtonEulerJointR> joint, std::shared_ptr<std::vector<unsigned int>> axes)
    : NewtonEulerR{}, _joint(joint), _axis(axes) {
  if (axes) {
    _axisMin = 100;
    _axisMax = 0;
    for (unsigned int i = 0; i < _axis->size(); i++) {
      if ((*_axis)[i] > _axisMax) _axisMax = (*_axis)[i];
      if ((*_axis)[i] < _axisMin) _axisMin = (*_axis)[i];
    }
  } else {
    _axisMin = _axisMax = 0;
    _axis = std::make_shared<std::vector<unsigned int>>();
    _axis->push_back(0);
  }

  assert((_axisMax - _axisMin + 1) <= _joint->numberOfDoF());
}

void siconos::joints::JointFrictionR::computeh(double time,
                                               const siconos::algebra::BlockVector& q0,
                                               siconos::algebra::SiconosVector& y) {
  // Velocity-level constraint, no position-level h
  y.setZero();
}

void siconos::joints::JointFrictionR::computeJacobianhOver_q(
    double time, siconos::modeling::Interaction& inter,
    std::shared_ptr<siconos::algebra::BlockVector> q0) {
  unsigned int n = _axisMax - _axisMin + 1;
  assert(n == 1);  // For now, multi-axis support TODO

  if (!jacobianhOver_q_Tmp || !(jacobianhOver_q_Tmp->size(1) == q0->size() && jacobianhOver_q_Tmp->size(0) == n)) {
    jacobianhOver_q_Tmp = std::make_shared<siconos::algebra::SiconosMatrix>(n, q0->size());
  }

  // Compute the jacobian for the required range of axes
  _joint->computeJachqDoF(time, inter, q0, *jacobianhOver_q_Tmp, _axisMin);

  // Copy indicated axes into the friction jacobian, negative and positive sides
  // NOTE trying ==1 using Relay, maybe don't need LCP formulation
  assert(jacobianhOver_q_->size(0) == 1);
  for (unsigned int i = 0; i < 1; i++)
    for (unsigned int j = 0; j < jacobianhOver_q_->size(1); j++) {
      jacobianhOver_q_->setValue(i, j,
                       jacobianhOver_q_Tmp->getValue((*_axis)[i] - _axisMin, j) * (i == 1 ? 1 : -1));
    }
}

unsigned int siconos::joints::JointFrictionR::numberOfConstraints() { return _axis->size(); }
unsigned int siconos::joints::JointFrictionR::axis(unsigned int _index) {
  return _axis->at(_index);
}

unsigned int siconos::joints::JointFrictionR::numberOfAxes() { return _axis->size(); }
