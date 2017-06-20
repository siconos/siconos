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
/*! \file KneeJointR.cpp

*/

#include "JointFrictionR.hpp"
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

/** Initialize a joint friction for a common case: a single axis with a
 * single friction, either positive or negative. For use with
 * NewtonImpactNSL. */
JointFrictionR::JointFrictionR(SP::NewtonEulerJointR joint, unsigned int axis)
  : NewtonEulerR()
  , _joint(joint)
  , _axis(std11::make_shared< std::vector<unsigned int> >())
{
  _axis->push_back(axis);
  _axisMin = axis;
  _axisMax = axis;
  assert( (_axisMax - _axisMin + 1) <= _joint->numberOfDoF() );
}

/** Initialize a multidimensional joint friction, e.g. the cone friction on
 * a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
JointFrictionR::JointFrictionR(SP::NewtonEulerJointR joint, SP::UnsignedIntVector axes)
  : NewtonEulerR()
  , _joint(joint)
  , _axis(axes)
{
  if (axes)
  {
    _axisMin = 100;
    _axisMax = 0;
    for (unsigned int i=0; i < _axis->size(); i++)
    {
      if ((*_axis)[i] > _axisMax) _axisMax = (*_axis)[i];
      if ((*_axis)[i] < _axisMin) _axisMin = (*_axis)[i];
    }
  }
  else
  {
    _axisMin = _axisMax = 0;
    _axis = std11::make_shared< std::vector<unsigned int> >();
    _axis->push_back(0);
  }

  assert( (_axisMax - _axisMin + 1) <= _joint->numberOfDoF() );
}

void JointFrictionR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  // Velocity-level constraint, no position-level h
  y.zero();
}

void JointFrictionR::computeJachq(double time, Interaction& inter, SP::BlockVector q0)
{
  unsigned int n = _axisMax - _axisMin + 1;
  assert(n==1); // For now, multi-axis support TODO

  if (!_jachqTmp || !(_jachqTmp->size(1) == q0->size() &&
                      _jachqTmp->size(0) == n))
  {
    _jachqTmp = std11::make_shared<SimpleMatrix>(n, q0->size());
  }

  // Compute the jacobian for the required range of axes
  _joint->computeJachqDoF(time, inter, q0, *_jachqTmp, _axisMin);

  // Copy indicated axes into the friction jacobian, negative and positive sides
  // NOTE trying ==1 using Relay, maybe don't need LCP formulation
  assert(_jachq->size(0)==1);
  for (unsigned int i=0; i<1; i++)
    for (unsigned int j=0; j<_jachq->size(1); j++) {
      _jachq->setValue(i,j,_jachqTmp->getValue((*_axis)[i]-_axisMin,j) * (i==1?1:-1));
    }
}

unsigned int JointFrictionR::numberOfConstraints()
{
  return _axis->size();
}
