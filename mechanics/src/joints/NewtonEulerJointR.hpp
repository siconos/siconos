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
/*! \file NewtonEulerJointR.hpp

*/
#ifndef NewtonEulerJointRELATION_H
#define NewtonEulerJointRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerR.hpp>

/** \class NewtonEulerJointR
 *  \brief This class implements an abstract Joint relation (articulation) between one or two Newton/Euler dynamical systems.
 */
class NewtonEulerJointR : public NewtonEulerR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(NewtonEulerJointR);
  NewtonEulerJointR(): NewtonEulerR(), _allowSelfCollide(false) {};

  /** A flag determining whether this joint should block
   * "self-collision", i.e., if true, bodies connected by this joint
   * will not enter into unilateral contact. */
  bool _allowSelfCollide;

public:

  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(double time, BlockVector& q0, SiconosVector& y,
                           unsigned int axis=0) {}

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(double time, Interaction& inter,
                               SP::BlockVector q0, SimpleMatrix& jachq,
                               unsigned int axis=0) {}

  /** Compute the vector of linear and angular velocities of the free axes */
  virtual void computeVelDoF(double time, BlockVector& q0, SiconosVector& v) {}

  /** Project a vector (assumed to be in q1 frame) onto the given
   * 0-indexed free axis. Useful for calculating velocities in the
   * axis, or for calculating axis-aligned forces applied to connected
   * bodies. */
  virtual void projectOntoAxis(SP::SiconosVector v, SP::SiconosVector ans, int axis=0)
    {}

  /** Return the value of the _allowSelfCollide flag. */
  bool allowSelfCollide() { return _allowSelfCollide; }

  /** Set the value of the _allowSelfCollide flag. */
  void setAllowSelfCollide(bool x) { _allowSelfCollide = x; }

  /** Get the number of constraints defined in the joint
      \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() = 0;

  /** Return the number of degrees of freedom of this joint.
      \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() = 0;

  typedef enum {
    DOF_TYPE_INVALID=0,
    DOF_TYPE_LINEAR=1,
    DOF_TYPE_ANGULAR=2,
  } DoF_Type;

  /** Return the type of a degree of freedom of this joint.
      \return the type of the degree of freedom (DoF)
   */
  virtual DoF_Type typeOfDoF(unsigned int axis) = 0;

  /** destructor
   */
  virtual ~NewtonEulerJointR() {};
};
#endif  //NewtonEulerJointRELATION_H
