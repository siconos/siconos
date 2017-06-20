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
/*! \file JointStopR.hpp

*/
#ifndef JointStopRELATION_H
#define JointStopRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerR.hpp>
#include <NewtonEulerJointR.hpp>
#include <Tools.hpp>

/** \class JointStopR
 *  \brief This class implements a stop on a DoF for any NewtonEulerJointR.
 */
class JointStopR : public NewtonEulerR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(JointStopR);
  JointStopR() : NewtonEulerR() {}

  SP::NewtonEulerJointR _joint;

  SP::UnsignedIntVector _axis;
  SP::SiconosVector _pos;
  SP::SiconosVector _dir;

  unsigned int _axisMin, _axisMax;
  SP::SimpleMatrix _jachqTmp;

public:

  /** Initialize a joint stop for a common case: a single axis with a
   * single stop, either positive or negative. For use with
   * NewtonImpactNSL. */
  JointStopR(SP::NewtonEulerJointR joint, double pos, bool dir,
             unsigned int axis=0);

  /** Initialize a multidimensional joint stop, e.g. the cone stop on
   * a ball joint. For use with NewtonImpactFrictionNSL size 2 or 3. */
  JointStopR(SP::NewtonEulerJointR joint, SP::SiconosVector pos,
             SP::SiconosVector dir, SP::UnsignedIntVector axes);

  /* The following constructor is disabled for now.  In fact
   * JointStopR is designed to support multiple stops, but it does not
   * work in practice since only the first value of y is examined to
   * determine whether to put the Interaction into indexSet1,
   * c.f. MoreauJeanOSI::addInteractionInIndexSet(). Use an individual
   * instance of JointStopR for each stop instead.*/
#if 0
  /** Initialize a joint stop for a common case: a single axis with a
   * double stop, one positive and one negative. */
  JointStopR(SP::NewtonEulerJointR joint, double pos, double neg,
             unsigned int axis=0);
#endif

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  virtual void computeJachq(double time, Interaction& inter, SP::BlockVector q0);

  virtual unsigned int numberOfConstraints();

  /** destructor
   */
  virtual ~JointStopR() {};
};
#endif  //JointStopRELATION_H
