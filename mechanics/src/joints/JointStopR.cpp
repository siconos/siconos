/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include <NewtonEulerDS.hpp>
#include <Interaction.hpp>
#include <boost/math/quaternion.hpp>
#include <BlockVector.hpp>
#include <cfloat>
#include <iostream>

// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

/** Initialize a joint stop for a common case: a single axis with a
 * single stop, either positive or negative. For use with
 * NewtonImpactNSL. */
JointStopR::JointStopR(SP::NewtonEulerJointR joint, double pos, bool dir,
                       unsigned int axis)
  : NewtonEulerR()
  , _joint(joint)
  , _axis(std11::make_shared< std::vector<unsigned int> >())
  , _pos(std11::make_shared<SiconosVector>(1))
  , _dir(std11::make_shared<SiconosVector>(1))
{
  _axis->push_back(axis);
  _pos->setValue(0, pos);
  _dir->setValue(0, dir ? -1 : 1);
  _axisMin = axis;
  _axisMax = axis;
  assert( (_axisMax - _axisMin + 1) <= _joint->numberOfDoF() );
}

/** Initialize a multidimensional joint stop, e.g. the cone stop on
 * a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
JointStopR::JointStopR(SP::NewtonEulerJointR joint, SP::SiconosVector pos,
                       SP::SiconosVector dir, SP::UnsignedIntVector axes)
  : NewtonEulerR()
  , _joint(joint)
  , _axis(axes)
  , _pos(pos)
  , _dir(dir)
{
  _axisMin = 100;
  _axisMax = 0;
  for (unsigned int i=0; i < _axis->size(); i++)
  {
    if ((*_axis)[i] > _axisMax) _axisMax = (*_axis)[i];
    if ((*_axis)[i] < _axisMin) _axisMin = (*_axis)[i];
  }
  assert( (_axisMax - _axisMin + 1) <= _joint->numberOfDoF() );
}

#if 0  // Disabled, see JointStopR.hpp.  Use multiple JointStopR instead.
/** Initialize a joint stop for a common case: a single axis with a
 * double stop, one positive and one negative. */
JointStopR::JointStopR(SP::NewtonEulerJointR joint, double pos, double neg,
                       unsigned int axis)
  : NewtonEulerR()
  , _joint(joint)
  , _axis(std11::make_shared< std::vector<unsigned int> >())
  , _pos(std11::make_shared<SiconosVector>(2))
  , _dir(std11::make_shared<SiconosVector>(2))
{
  _axis->push_back(axis);
  _axis->push_back(axis);
  _pos->setValue(0, pos);
  _pos->setValue(1, neg);
  _dir->setValue(0, 1);
  _dir->setValue(1, -1);
  _axisMin = axis;
  _axisMax = axis;
  assert( (_axisMax - _axisMin + 1) <= _joint->numberOfDoF() );
}
#endif

void JointStopR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  // Common cases optimisation
  bool case_onestop = y.size()==1;
  bool case_posneg = y.size()==2 && (*_axis)[0] == (*_axis)[1];
  if (case_onestop || case_posneg)
  {
    _joint->computehDoF(time, q0, y, (*_axis)[0]);

    y.setValue(0, (y.getValue(0) - _pos->getValue(0)) * _dir->getValue(0));
    if (case_posneg)
      y.setValue(1, (y.getValue(0) - _pos->getValue(1)) * _dir->getValue(1));
    return;
  }

  // Get h for each relevant axis
  SiconosVector tmp_y(_axisMax - _axisMin + 1);
  _joint->computehDoF(time, q0, tmp_y, _axisMin);

  // Copy and scale each stop for its axis/position/direction
  for (unsigned int i=0; i < y.size(); i++) {
    y.setValue(i, (tmp_y.getValue((*_axis)[i])
                   - _pos->getValue(i))*_dir->getValue(i));
  }
}

void JointStopR::computeJachq(double time, Interaction& inter, SP::BlockVector q0)
{
  unsigned int n = _axisMax - _axisMin + 1;

  if (!_jachqTmp || !(_jachqTmp->size(1) == q0->size() &&
                      _jachqTmp->size(0) == n))
  {
    _jachqTmp = std11::make_shared<SimpleMatrix>(n, q0->size());
  }

  // Compute the jacobian for the required range of axes
  _joint->computeJachqDoF(time, inter, q0, *_jachqTmp, _axisMin);

  // Copy indicated axes into the stop jacobian, possibly flipped for negative stops
  for (unsigned int i=0; i<_jachq->size(0); i++)
    for (unsigned int j=0; j<_jachq->size(1); j++)
      _jachq->setValue(i,j,
                       _jachqTmp->getValue((*_axis)[i]-_axisMin,j) * _dir->getValue(i));
}

unsigned int JointStopR::numberOfConstraints()
{
  return _axis->size();
}
