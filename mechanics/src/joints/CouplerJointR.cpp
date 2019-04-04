/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
CouplerJointR::CouplerJointR(SP::NewtonEulerJointR joint1, unsigned int dof1,
                             SP::NewtonEulerJointR joint2, unsigned int dof2,
                             double ratio,
                             SP::SiconosVector ref1, unsigned int ref1_index,
                             SP::SiconosVector ref2, unsigned int ref2_index)
  : NewtonEulerJointR()
  , _joint1(joint1)
  , _joint2(joint2)
  , _ref1(ref1)
  , _ref2(ref2)
  , _dof1(dof1)
  , _dof2(dof2)
  , _ref1_index(ref1_index)
  , _ref2_index(ref2_index)
  , _ratio(ratio)
  , _offset(0.0)
{
  assert( _dof1 < _joint1->numberOfDoF() );
  assert( _dof2 < _joint2->numberOfDoF() );
}

/* A constructor taking a DS exists just because it's hard to pass
 * ds->q() through Python without it automatically converting to numpy
 * array and back, which messes up the shared_ptr reference.
 *
 * NOTE that using q() as the reference is not quite right, in fact it
 * should be using the reference ds's temporary work vector in order
 * to perform correctly in the Newton loop. (TODO) */
CouplerJointR::CouplerJointR(SP::NewtonEulerJointR joint1, unsigned int dof1,
                             SP::NewtonEulerJointR joint2, unsigned int dof2,
                             double ratio,
                             SP::NewtonEulerDS refds1, unsigned int ref1_index,
                             SP::NewtonEulerDS refds2, unsigned int ref2_index)
  : NewtonEulerJointR()
  , _joint1(joint1)
  , _joint2(joint2)
  , _ref1(refds1->q())
  , _ref2(refds2->q())
  , _dof1(dof1)
  , _dof2(dof2)
  , _ref1_index(ref1_index)
  , _ref2_index(ref2_index)
  , _ratio(ratio)
  , _offset(0.0)
{
  assert( _dof1 < _joint1->numberOfDoF() );
  assert( _dof2 < _joint2->numberOfDoF() );
}

void CouplerJointR::setReferences(SP::NewtonEulerJointR joint1, unsigned int dof1,
                                  SP::NewtonEulerJointR joint2, unsigned int dof2,
                                  SP::SiconosVector ref1, unsigned int ref1_index,
                                  SP::SiconosVector ref2, unsigned int ref2_index)
{
  _joint1 = joint1;
  _joint2 = joint2;
  _dof1 = dof1;
  _dof2 = dof2;
  _ref1 = ref1;
  _ref2 = ref2;
  _ref1_index = ref1_index;
  _ref2_index = ref2_index;
}

void CouplerJointR::setReferences(SP::NewtonEulerJointR joint1, unsigned int dof1,
                                  SP::NewtonEulerJointR joint2, unsigned int dof2,
                                  SP::NewtonEulerDS refds1, unsigned int ref1_index,
                                  SP::NewtonEulerDS refds2, unsigned int ref2_index)
{
  _joint1 = joint1;
  _joint2 = joint2;
  _dof1 = dof1;
  _dof2 = dof2;
  _ref1 = refds1->q();
  _ref2 = refds2->q();
  _ref1_index = ref1_index;
  _ref2_index = ref2_index;
}

void CouplerJointR::setRatio(double ratio)
{
  _ratio = ratio;
}

void CouplerJointR::makeBlockVectors(SP::SiconosVector q1, SP::SiconosVector q2,
                                     BlockVector& q01, BlockVector& q02)
{
  /* Decide what to do for each combination of presence of references
   * and/or q2. Only the existence of q1 is absolutely required.
   * Logic: if no references, then the coupling is between DoFs of q1
   * and q2; in most cases here, joint1 and joint2 will be the same,
   * but they don't have to be.  If there are references, then the
   * coupling is between the measurements of ref1-q1 and ref2-q2, with
   * respect to joint1 and joint2.  The refX_index variables allow to
   * specify the position of the references with respect to the joints
   * (i.e. whether joint1 is defined as q1-ref1 or ref1-q1.)  This has
   * an effect on the sign of the measurement, and therefore must be
   * specified by the user.  In mechanics_run.py, the user only needs
   * to specify the two reference joints: we assume the reference body
   * is that which is not ds1 or ds2, and reference indexes are
   * calculated automatically.  Here we do not have that luxury since
   * we only have SP::SiconosVectors. */
  VectorOfVectors vects1(2), vects2(2);
  if (!_ref1 && !_ref2 && q2)
  {
    vects1[0] = q1;
    vects1[1] = q2;
    vects2 = vects1;
  }
  else if (!q2)
  {
    if (_ref1) {
      vects1[_ref1_index] = _ref1;
      vects1[1-_ref1_index] = q1;
      vects2 = vects1;
    }
    else {
      vects1.resize(1);
      vects1[0] = q1;
      vects2 = vects1;
    }
  }
  else if (q2)
  {
    if (_ref1)
    {
      vects1[_ref1_index] = _ref1;
      vects1[1-_ref1_index] = q1;
    }
    else
    {
      vects1[0] = q1;
      vects1[1] = q2;
    }

    if (_ref2)
    {
      vects2[_ref2_index] = _ref2;
      vects2[1-_ref2_index] = q2;
    }
    else
    {
      vects2[0] = q1;
      vects2[1] = q2;
    }
  }

  q01.setAllVect(vects1);
  q02.setAllVect(vects2);
}

void CouplerJointR::setBasePositions(SP::SiconosVector q1, SP::SiconosVector q2)
{
  // TODO: time=0?

  // Get current positions of the implicated degrees of freedom
  SiconosVector y1(1), y2(1);
  BlockVector q01, q02;
  makeBlockVectors(q1, q2, q01, q02);
  _joint1->computehDoF(0, q01, y1, _dof1);
  _joint2->computehDoF(0, q02, y2, _dof2);

  // Compute initial offset between the DoFs
  _offset = y1(0) * _ratio - y2(0);
}

void CouplerJointR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  SiconosVector y1(y), y2(y);

  // Get current positions of the implicated degrees of freedom
  BlockVector q01, q02;
  makeBlockVectors(q0.vector(0),
                   q0.getAllVect().size()>1 ? q0.vector(1) : SP::SiconosVector(),
                   q01, q02);
  _joint1->computehDoF(time, q01, y1, _dof1);
  _joint2->computehDoF(time, q02, y2, _dof2);

  // Constraint is the linear relation between them
  y(0) = y2(0) - y1(0)*_ratio + _offset;
}

void CouplerJointR::computeJachq(double time, Interaction& inter, SP::BlockVector q0)
{
  SP::SimpleMatrix jachq1 = std11::make_shared<SimpleMatrix>(1, q0->size());
  SP::SimpleMatrix jachq2 = std11::make_shared<SimpleMatrix>(1, q0->size());

  // Get jacobians for the implicated degrees of freedom
  // Compute the jacobian for the required range of axes
  SP::BlockVector q01 = std11::make_shared<BlockVector>();
  SP::BlockVector q02 = std11::make_shared<BlockVector>();
  makeBlockVectors(q0->vector(0),
                   q0->getAllVect().size()>1 ? q0->vector(1) : SP::SiconosVector(),
                   *q01, *q02);
  _joint1->computeJachqDoF(time, inter, q01, *jachq1, _dof1);
  _joint2->computeJachqDoF(time, inter, q02, *jachq2, _dof2);

  // Constraint is the linear relation between them
  for (unsigned int i=0; i<1; i++)
    for (unsigned int j=0; j<_jachq->size(1); j++)
      _jachq->setValue(i,j, jachq2->getValue(i,j) - jachq1->getValue(i,j)*_ratio);
}

void CouplerJointR::_normalDoF(SiconosVector& ans, const BlockVector& q0, int axis,
                               bool absoluteRef)
{
  assert("CouplerJointR::_normalDoF is not defined.");

  // We could define the normal of this constraint as simply the cross
  // product of the normals of the two coupled DoFs, as calculated
  // below, however this is not really sensible as the constraint is
  // defined in terms of a line through a manifold, and would need to
  // be projected down to world coordinates somehow here.  For
  // example, what is the "direction of the constraint" in world
  // coordinates for a coupling between rotational and lineawr motion?

  /*
  SiconosVector n1(3), n2(3);
  BlockVector q01, q02;
  makeBlockVectors(q0.getAllVect()[0],
                   q0.getAllVect().size()>1 ? q0.getAllVect()[1] : SP::SiconosVector(),
                   q01, q02);
  _joint1->normalDoF(n1, q01, _dof1, absoluteRef);
  _joint2->normalDoF(n2, q02, _dof2, absoluteRef);
  cross_product(n1, n2, ans);
  */
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
  SP::BlockVector q01 = std11::make_shared<BlockVector>();
  SP::BlockVector q02 = std11::make_shared<BlockVector>();
  makeBlockVectors(q0->vector(0),
                   q0->getAllVect().size()>1 ? q0->vector(1) : SP::SiconosVector(),
                   *q01, *q02);
  _joint1->computeJachqDoF(time, inter, q01, *jachq1, _dof1);
  _joint2->computeJachqDoF(time, inter, q02, *jachq2, _dof2);

  // Constraint is the linear relation between them
  for (unsigned int i=0; i<1; i++)
    for (unsigned int j=0; j<_jachq->size(1); j++)
      jachq.setValue(i,j, jachq2->getValue(i,j) - jachq1->getValue(i,j)*_ratio);
}
