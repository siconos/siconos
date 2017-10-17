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
/*! \file CouplerJointRR.cpp
*/

#include "CouplerJointR.hpp"
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

CouplerJointR::CouplerJointR() : NewtonEulerJointR() {}

/** Initialize a coupler joint. For use with EqualityConditionNSL to
 * bind two degrees of freedom into one by a ratio and offset. */
CouplerJointR::CouplerJointR(SP::NewtonEulerJointR joint, unsigned int dof1,
                             unsigned int dof2, double ratio)
  : NewtonEulerJointR()
  , _joint(joint)
  , _dof1(dof1)
  , _dof2(dof2)
  , _ratio(ratio)
  , _offset(0.0)
{
  assert( _dof1 < _joint->numberOfDoF() );
  assert( _dof2 < _joint->numberOfDoF() );
}

#include "CylindricalJointR.hpp"
void CouplerJointR::setBasePositions(SP::SiconosVector q1, SP::SiconosVector q2)
{
  // Get current positions of the implicated degrees of freedom
  SiconosVector y1(1), y2(1);
  if (q2)
  {
    BlockVector q0( q1, q2 );
    // TODO: time=0?
    _joint->computehDoF(0, q0, y1, _dof1);
    _joint->computehDoF(0, q0, y2, _dof2);
  }
  else
  {
    BlockVector q0(1); q0.setVectorPtr(0, q1);
    // TODO: time=0?
    _joint->computehDoF(0, q0, y1, _dof1);
    _joint->computehDoF(0, q0, y2, _dof2);
  }

  // Compute initial offset between the DoFs
  _offset = y1(0) * _ratio - y2(0);
}

void CouplerJointR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  // Get current positions of the implicated degrees of freedom
  SiconosVector y1(y), y2(y);
  _joint->computehDoF(time, q0, y1, _dof1);
  _joint->computehDoF(time, q0, y2, _dof2);

  // Constraint is the linear relation between them
  y(0) = y2(0) - y1(0)*_ratio + _offset;
}

void CouplerJointR::computeJachq(double time, Interaction& inter, SP::BlockVector q0)
{
  SP::SimpleMatrix jachq1 = std11::make_shared<SimpleMatrix>(1, q0->size());
  SP::SimpleMatrix jachq2 = std11::make_shared<SimpleMatrix>(1, q0->size());

  // Get jacobians for the implicated degrees of freedom
  // Compute the jacobian for the required range of axes
  _joint->computeJachqDoF(time, inter, q0, *jachq1, _dof1);
  _joint->computeJachqDoF(time, inter, q0, *jachq2, _dof2);

  // Constraint is the linear relation between them
  for (unsigned int i=0; i<1; i++)
    for (unsigned int j=0; j<_jachq->size(1); j++)
      _jachq->setValue(i,j, jachq2->getValue(i,j) - jachq1->getValue(i,j)*_ratio);
}

void CouplerJointR::_normalDoF(SiconosVector& ans, const BlockVector& q0, int axis,
                               bool absoluteRef)
{
  // We define the normal of this constraint as simply the cross
  // product of the normals of the two coupled DoFs.
  SiconosVector n1(3), n2(3);
  _joint->normalDoF(n1, q0, _dof1, absoluteRef);
  _joint->normalDoF(n2, q0, _dof2, absoluteRef);
  cross_product(n1, n2, ans);
}

void CouplerJointR::computehDoF(double time, BlockVector& q0, SiconosVector& y,
                                unsigned int axis)
{
  // The DoF of the constraint is the same as the constraint itself
  assert(axis==0);
  this->computeh(time, q0, y);
}

void CouplerJointR::computeJachqDoF(double time, Interaction& inter,
                                    SP::BlockVector q0, SimpleMatrix& jachq,
                                    unsigned int axis)
{
  // The Jacobian of the DoF of the constraint is the same as the
  // Jacobian of the constraint itself. (Same as computeJachq(), but
  // don't store result in member object.)
  assert(axis==0);

  SP::SimpleMatrix jachq1 = std11::make_shared<SimpleMatrix>(1, q0->size());
  SP::SimpleMatrix jachq2 = std11::make_shared<SimpleMatrix>(1, q0->size());

  // Get jacobians for the implicated degrees of freedom
  // Compute the jacobian for the required range of axes
  _joint->computeJachqDoF(time, inter, q0, *jachq1, _dof1);
  _joint->computeJachqDoF(time, inter, q0, *jachq2, _dof2);

  // Constraint is the linear relation between them
  for (unsigned int i=0; i<1; i++)
    for (unsigned int j=0; j<_jachq->size(1); j++)
      jachq.setValue(i,j, jachq2->getValue(i,j) - jachq1->getValue(i,j)*_ratio);
}
