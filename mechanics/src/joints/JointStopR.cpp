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
  (*_pos)(0) = pos;
  (*_dir)(0) = dir ? -1 : 1;
  _axisMin = axis;
  _axisMax = axis;
  assert((_axisMax - _axisMin + 1) <= _joint->numberOfDoF());
}

/** Initialize a multidimensional joint stop, e.g. the cone stop on
 * a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
siconos::joints::JointStopR::JointStopR(
    std::shared_ptr<NewtonEulerJointR> joint,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& pos,
    const Eigen::Ref<const siconos::algebra::SiconosVector>& dir,
    std::shared_ptr<std::vector<unsigned int>> axes)
    : NewtonEulerR{}, _joint(joint), _axis(axes) {
  _axisMin = 100;
  _axisMax = 0;

  _pos = std::make_shared<siconos::algebra::SiconosVector>(pos);
  _dir = std::make_shared<siconos::algebra::SiconosVector>(dir);

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
  (*_pos)(0) = pos;
  (*_pos)(1) = neg;
  (*_dir)(0) = 1;
  (*_dir)(1) = -1;
  _axisMin = axis;
  _axisMax = axis;
  assert((_axisMax - _axisMin + 1) <= _joint->numberOfDoF());
}
#endif

void siconos::joints::JointStopR::computeh(const siconos::algebra::BlockVector& q,
                                           Eigen::Ref<siconos::algebra::SiconosVector> y) {
  // Common cases optimisation
  bool case_onestop = y.size() == 1;
  bool case_posneg = y.size() == 2 && (*_axis)[0] == (*_axis)[1];
  if (case_onestop || case_posneg) {
    if (q.numberOfBlocks() == 2)
      _joint->computehDoF(*q.vector(0), *q.vector(1), y, (*_axis)[0]);
    else
      _joint->computehDoF(*q.vector(0), std::nullopt, y, (*_axis)[0]);
    y(0) = (y(0) - (*_pos)(0)) * (*_dir)(0);
    if (case_posneg) y(1) = (y(0) - (*_pos)(1)) * (*_dir)(1);
    return;
  }

  // Get h for each relevant axis
  siconos::algebra::SiconosVector tmp_y(_axisMax - _axisMin + 1);

  if (q.numberOfBlocks() == 2)
    _joint->computehDoF(*q.vector(0), *q.vector(1), tmp_y, _axisMin);
  else
    _joint->computehDoF(*q.vector(0), std::nullopt, tmp_y, _axisMin);

  // Copy and scale each stop for its axis/position/direction
  for (unsigned int i = 0; i < y.size(); i++) {
    y(i) = (tmp_y((*_axis)[i]) - (*_pos)(i)) * (*_dir)(i);
  }
}

void siconos::joints::JointStopR::computeH_NE_(
    double time, siconos::modeling::Interaction& inter,
    const siconos::algebra::BlockVector& q0) {
  unsigned int n = _axisMax - _axisMin + 1;

  if (!jacobianhOver_q_Tmp ||
      !(jacobianhOver_q_Tmp->cols() == q0.size() && jacobianhOver_q_Tmp->rows() == n)) {
    jacobianhOver_q_Tmp = std::make_shared<siconos::algebra::SiconosMatrix>(n, q0.size());
  }

  // Compute the jacobian for the required range of axes
  _joint->computeJachqDoF(inter, q0, *jacobianhOver_q_Tmp, _axisMin);

  // Copy indicated axes into the stop jacobian, possibly flipped for negative stops
  for (unsigned int i = 0; i < H_NE_view_->rows(); i++)
    for (unsigned int j = 0; j < H_NE_view_->cols(); j++)
      H_NE_view_->setValue(
          i, j, (*jacobianhOver_q_Tmp)((*_axis)[i] - _axisMin, j) * (*_dir)(i));
}

unsigned int siconos::joints::JointStopR::numberOfConstraints() const { return _axis->size(); }

auto siconos::joints::JointStopR::axis(unsigned int _index) { return _axis->at(_index); }

double siconos::joints::JointStopR::position(unsigned int _index) {
  return (*_pos)(_index);
}

double siconos::joints::JointStopR::direction(unsigned int _index) {
  return (*_dir)(_index);
}

auto siconos::joints::JointStopR::numberOfAxes() { return _axis->size(); }
