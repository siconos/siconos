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

/*! \file CylindricalJointR.hpp
*/

#ifndef CylindricalJointRELATION_H
#define CylindricalJointRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerJointR.hpp>

/** \class CylindricalJointR
 *  \brief This class implements a cylindrical joint between one or
 *  two Newton/Euler Dynamical system.  It is similar to a
 *  PrismaticJointR but allows for rotation around the axis.
 *
 * From a given axis, we construct two unit othorgonal vectors to the
 *  axis V1 and V2 such that (axis,V1,V2) is an orthogonal frame
 */
class CylindricalJointR : public NewtonEulerJointR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(CylindricalJointR);
  CylindricalJointR(): NewtonEulerJointR() {};

public:
  /** Axis of the cylindrical point in the inertial frame of reference
   */
  SP::SiconosVector _axis0;

  /** _V1 is an unit vector that is orthogonal to the cylindrical axis
   * _axis0.  It forms with _V2 and _axis0 a base such that
   * (_axis0,_V1,_v2) is an orthogonal frame
   */
  SP::SiconosVector _V1;

  /** _V2 is an unit vector that is orthogonal to the cylindrical axis
   * _axis0.  It forms with _V2 and _axis0 a base such that
   * (_axis0,_V1,_v2) is an orthogonal frame
   */
  SP::SiconosVector _V2;

  double _cq2q101;
  double _cq2q102;
  double _cq2q103;
  double _cq2q104;

  /** P is the point defining the location of the line created by
   * _axis0.  It is stored in the q1 frame, i.e. the vector from
   * initial G1 to P, called _G1P0. */
  SP::SiconosVector _G1P0;

  /** _G2P0 is the vector from initial G1 to P */
  SP::SiconosVector _G2P0;

  /** Cumulative number of twists around the joint relative to initial
   * angular difference. */
  int _twistCount;    // TODO: Should be in a graph work vector?
  double _previousAngle; // Needed to track _twistCount, TODO: work vector?
  double _initialAngle;

  /** constructor from two dynamical systems and an axis
   *  \param d1 first  DynamicalSystem link by the  joint
   *  \param d2 second  DynamicalSystem link by the joint
   *  \param A SiconosVector of size 3 that defines the cylindrical axis
   *           in the body frame of d1
   */
  CylindricalJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2,
                    SP::SiconosVector P, SP::SiconosVector A,
                    bool absoluteRef = false);

  /** constructor from one dynamical systems and an axis
   *  \param d1 the  DynamicalSystem link by the  joint
   *  \param P SiconosVector of size 3 that defines the point around
   *           which rotation is allowed
   *  \param A SiconosVector of size 3 that defines the cylindrical axis
   *           in the inertial frame of reference
   *  \param absoluteRef if true, P is in the absolute frame,
   *                     otherwise P is in d1 frame
   */
  CylindricalJointR(SP::NewtonEulerDS d1,
                    SP::SiconosVector P, SP::SiconosVector A,
                    bool absoluteRef = false);

  void computeFromInitialPosition(SP::SiconosVector q2,
                                  SP::SiconosVector q1=SP::SiconosVector());

  void computeV1V2FromAxis();

  /** destructor
   */
  virtual ~CylindricalJointR() {};

  virtual void computeJachq(double time, Interaction& inter, SP::BlockVector q0 );

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(double time, BlockVector& q0, SiconosVector& y,
                           unsigned int axis);

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(double time, Interaction& inter,
                               SP::BlockVector q0, SimpleMatrix& jachq,
                               unsigned int axis);

  void Jd1d2(
    double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
    double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  void Jd1(
    double X1, double Y1, double Z1, double q10, double q11, double q12, double q13);

  /** Get the number of constraints defined in the joint
      \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() { return 4; }

  /** Return the number of degrees of freedom of this joint.
      \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() { return 2; }

  /** Return the type of a degree of freedom of this joint.
      \return the type of the degree of freedom (DoF)
  */
  virtual DoF_Type typeOfDoF(unsigned int axis) {
    if (axis==0) return DOF_TYPE_LINEAR;
    else if (axis==1) return DOF_TYPE_ANGULAR;
    else return DOF_TYPE_INVALID;
  }
};
#endif  //CylindricalJointRELATION_H
