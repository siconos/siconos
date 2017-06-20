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
  NewtonEulerJointR(): NewtonEulerR()
                     , _allowSelfCollide(false)
                     , _absoluteRef(true) {};

  /** A flag determining whether this joint should block
   * "self-collision", i.e., if true, bodies connected by this joint
   * will not enter into unilateral contact. */
  bool _allowSelfCollide;

  /** Points used to defined the joint constraint. */
  VectorOfVectors _points;

  /** Axes used to defined the joint constraint. */
  VectorOfVectors _axes;

  /** Defines whether points and axes are specified in absolute or
   * relative frame. */
  bool _absoluteRef;

  /** Private version of normalDoF for subclasses to override. */
  virtual void _normalDoF(const BlockVector& q0, SiconosVector& ans, int axis,
                          bool absoluteRef=true) {}

public:

  /** Set a point for this joint. The role of each point is specific
   * to the joint subclass. Won't take effect until
   * setInitialConditions is called.
   *
   * \param index The index of the points.
   * \param point A SiconosVector of size 3.
   */
  void setPoint(unsigned int index, SP::SiconosVector point)
    { _points[index] = point; }

  /** Set an axis for this joint. The role of each axis is specific to
   * the joint subclass. Won't take effect until setInitialConditions
   * is called.
   *
   * \param index The index of the points.
   * \param axis A SiconosVector of size 3.
   */
  void setAxis(unsigned int index, SP::SiconosVector axis)
    { _axes[index] = axis; }

  /** Set whether points and axes should be interpreted in absolute or
   * relative frame. Won't take effect until setInitialConditions is
   * called.
   *
   * \param absoluteRef true for absolute frame, false for relative frame.
   */
  void setAbsolute(bool absoluteRef)
    { _absoluteRef = absoluteRef; }

  /** Initialize the joint constants based on the provided initial positions. */
  virtual void setInitialConditions(SP::SiconosVector q1,
                                    SP::SiconosVector q2=SP::SiconosVector()) = 0;

  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(double time, BlockVector& q0, SiconosVector& y,
                           unsigned int axis=0) {}

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(double time, Interaction& inter,
                               SP::BlockVector q0, SimpleMatrix& jachq,
                               unsigned int axis=0) {}

  /** Project a vector onto the given 0-indexed free axis. Useful for
   * calculating velocities in the axis, or for calculating
   * axis-aligned forces applied to connected bodies.  If axis is of
   * angular type (see typeOfDoF), then the projection is onto the
   * axis of rotation.
   *
   * \param v The vector to project
   * \param q0 The state q of one or more NewtonEulerDS
   * \param ans The vector to receive the projection.
   * \param absoluteRef If true, v and ans are in the inertial frame,
   *                    otherwise the q1 frame is assumed.
   */
  void projectVectorDoF(const SiconosVector& v, const BlockVector& q0,
                        SiconosVector& ans, int axis,
                        bool absoluteRef=true);

  SP::SiconosVector projectVectorDoF(const SiconosVector& v,
                                     const BlockVector& q0, int axis,
                                     bool absoluteRef=true);

  /** Retrieve a normal in the direction of a 0-indexed free
   * axis. Useful for calculating velocities in the axis, or for
   * calculating axis-aligned forces applied to connected bodies.  If
   * axis is of angular type (see typeOfDoF), then the returned normal
   * is the axis of rotation.
   *
   * \param q0 The state q of one or more NewtonEulerDS
   * \param ans The vector to receive the projection.
   * \param absoluteRef If true, ans is in the inertial frame,
   *                    otherwise the q1 frame is assumed.
   */
  void normalDoF(const BlockVector& q0, SiconosVector& ans, int axis,
                 bool absoluteRef=true);

  SP::SiconosVector normalDoF(const BlockVector& q0, int axis,
                              bool absoluteRef=true);

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
