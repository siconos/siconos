/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
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
/*! \file JointFrictionR.hpp

*/
#ifndef JointFrictionRELATION_H
#define JointFrictionRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerR.hpp>
#include <NewtonEulerJointR.hpp>
#include <Tools.hpp>

/** \class JointFrictionR
 *  \brief This class implements a friction on a DoF for any NewtonEulerJointR.
 */
class JointFrictionR : public NewtonEulerR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(JointFrictionR);
  JointFrictionR() : NewtonEulerR() {}

  SP::NewtonEulerJointR _joint;

  SP::UnsignedIntVector _axis;

  unsigned int _axisMin, _axisMax;
  SP::SimpleMatrix _jachqTmp;

public:

  /** Initialize a joint friction for a common case: a single axis with a
   * single friction, either positive or negative. For use with
   * NewtonImpactNSL. */
  JointFrictionR(SP::NewtonEulerJointR joint, unsigned int axis);

  /** Initialize a multidimensional joint friction, e.g. the cone friction on
   * a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
  JointFrictionR(SP::NewtonEulerJointR joint,
                 SP::UnsignedIntVector axes=SP::UnsignedIntVector());

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  virtual void computeJachq(double time, Interaction& inter, SP::BlockVector q0);

  virtual unsigned int numberOfConstraints();

  /* Return the joint axis number assigned to a friction axis. */
  unsigned int axis(unsigned int _index) { return _axis->at(_index); }

  /* Return the joint assigned to this friction relation. */
  SP::NewtonEulerJointR joint() { return _joint; }

  /* Return the number of joint axes indexed by this relation. */
  unsigned int numberOfAxes() { return _axis->size(); }

  /** destructor
   */
  virtual ~JointFrictionR() {};
};
#endif  //JointFrictionRELATION_H
