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
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

/** Initialize a joint stop for a common case: a single axis with a
 * single stop, either positive or negative. For use with
 * NewtonImpactNSL. */
siconos::joints::JointStopR::JointStopR(std::shared_ptr<NewtonEulerJointR> joint, double pos,
                                        bool dir, unsigned int axis)
    : NewtonEulerR{},
      _joint(joint),
      _axis(std::make_shared<std::vector<unsigned int>>()),
      _pos(std::make_shared<siconos::algebra::SiconosVector>(1)),
      _dir(std::make_shared<siconos::algebra::SiconosVector>(1)) {
  _axis->push_back(axis);
  _pos->setValue(0, pos);
  _dir->setValue(0, dir ? -1 : 1);
  _axisMin = axis;
  _axisMax = axis;
  assert((_axisMax - _axisMin + 1) <= _joint->numberOfDoF());
}

/** Initialize a multidimensional joint stop, e.g. the cone stop on
 * a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
siconos::joints::JointStopR::JointStopR(std::shared_ptr<NewtonEulerJointR> joint,
                                        std::shared_ptr<siconos::algebra::SiconosVector> pos,
                                        std::shared_ptr<siconos::algebra::SiconosVector> dir,
                                        std::shared_ptr<std::vector<unsigned int>> axes)
    : NewtonEulerR{}, _joint(joint), _axis(axes), _pos(pos), _dir(dir) {
  _axisMin = 100;
  _axisMax = 0;
  for (unsigned int i = 0; i < _axis->size(); i++) {
    if ((*_axis)[i] > _axisMax) _axisMax = (*_axis)[i];
    if ((*_axis)[i] < _axisMin) _axisMin = (*_axis)[i];
  }
  assert((_axisMax - _axisMin + 1) <= _joint->numberOfDoF());
}

#if 0  // Disabled, see JointStopR.hpp.  Use multiple JointStopR instead.
/** Initialize a joint stop for a common case: a single axis with a
 * double stop, one positive and one negative. */
siconos::joints::JointStopR::JointStopR(std::shared_ptr<NewtonEulerJointR> joint, double pos, double neg,
                       unsigned int axis)
  : NewtonEulerR()
  , _joint(joint)
  , _axis(std::make_shared< std::vector<unsigned int> >())
  , _pos(std::make_shared<siconos::algebra::SiconosVector>(2))
  , _dir(std::make_shared<siconos::algebra::SiconosVector>(2))
{
  _axis->push_back(axis);
  _axis->push_back(axis);
  _pos->setValue(0, pos);
  _pos->setValue(1, neg);
  _dir->setValue(0, 1);
  _dir->setValue(1, -1);
  _axisMin = axis;
  _axisMax = axis;
  assert((_axisMax - _axisMin + 1) <= _joint->numberOfDoF());
}
#endif

void siconos::joints::JointStopR::computeh(double time,
                                           const siconos::algebra::BlockVector& q0,
                                           siconos::algebra::SiconosVector& y) {
  // Common cases optimisation
  bool case_onestop = y.size() == 1;
  bool case_posneg = y.size() == 2 && (*_axis)[0] == (*_axis)[1];
  if (case_onestop || case_posneg) {
    _joint->computehDoF(time, q0, y, (*_axis)[0]);

    y.setValue(0, (y.getValue(0) - _pos->getValue(0)) * _dir->getValue(0));
    if (case_posneg) y.setValue(1, (y.getValue(0) - _pos->getValue(1)) * _dir->getValue(1));
    return;
  }

  // Get h for each relevant axis
  siconos::algebra::SiconosVector tmp_y(_axisMax - _axisMin + 1);
  _joint->computehDoF(time, q0, tmp_y, _axisMin);

  // Copy and scale each stop for its axis/position/direction
  for (unsigned int i = 0; i < y.size(); i++) {
    y.setValue(i, (tmp_y.getValue((*_axis)[i]) - _pos->getValue(i)) * _dir->getValue(i));
  }
}

void siconos::joints::JointStopR::computeJachq(
    double time, siconos::modeling::Interaction& inter,
    std::shared_ptr<siconos::algebra::BlockVector> q0) {
  unsigned int n = _axisMax - _axisMin + 1;

  if (!_jachqTmp || !(_jachqTmp->size(1) == q0->size() && _jachqTmp->size(0) == n)) {
    _jachqTmp = std::make_shared<siconos::algebra::SimpleMatrix>(n, q0->size());
  }

  // Compute the jacobian for the required range of axes
  _joint->computeJachqDoF(time, inter, q0, *_jachqTmp, _axisMin);

  // Copy indicated axes into the stop jacobian, possibly flipped for negative stops
  for (unsigned int i = 0; i < _jachq->size(0); i++)
    for (unsigned int j = 0; j < _jachq->size(1); j++)
      _jachq->setValue(i, j,
                       _jachqTmp->getValue((*_axis)[i] - _axisMin, j) * _dir->getValue(i));
}

unsigned int siconos::joints::JointStopR::numberOfConstraints() { return _axis->size(); }

unsigned int siconos::joints::JointStopR::axis(unsigned int _index) {
  return _axis->at(_index);
}

double siconos::joints::JointStopR::position(unsigned int _index) {
  return _pos->getValue(_index);
}

double siconos::joints::JointStopR::direction(unsigned int _index) {
  return _dir->getValue(_index);
}

unsigned int siconos::joints::JointStopR::numberOfAxes() { return _axis->size(); }
