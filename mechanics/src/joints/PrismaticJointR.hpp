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
/*! \file PrismaticJointR.hpp

*/
#ifndef PrismaticJointRELATION_H
#define PrismaticJointRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerJointR.hpp>

/** \class PrismaticJointR
 *  \brief This class implements a prismatic joint between one or two Newton/Euler Dynamical system
 *
 * From a given axis, we construct two unit othorgonal vectors to the axis V1 and V2 such that
 *  (axis,V1,V2) is an orthogonal frame
 *
 *
 *
 */
class PrismaticJointR : public NewtonEulerJointR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(PrismaticJointR);
  PrismaticJointR(): NewtonEulerJointR() {};

  /** Axis of the prismatic point in the q1 inertial frame of reference
   */
  SP::SiconosVector _axis0;

  /** _V1 is an unit vector that is orthogonal to the prismatic axis _axis0.
   * It forms with _V2 and _axis0 a base such that (_axis0,_V1,_v2) is an orthogonal
   * frame
   */
  SP::SiconosVector _V1;

  /** _V2 is an unit vector that is orthogonal to the prismatic axis _axis0.
   * It forms with _V2 and _axis0 a base such that (_axis0,_V1,_v2) is an orthogonal
   * frame
   */
  SP::SiconosVector _V2;

  /** Convenient storage of the components of _V1 and _V2
   */
  double _V1x;
  double _V1y;
  double _V1z;
  double _V2x;
  double _V2y;
  double _V2z;

  /**
   */
  double _G10G20d1x;
  double _G10G20d1y;
  double _G10G20d1z;

  double _cq2q101;
  double _cq2q102;
  double _cq2q103;
  double _cq2q104;

  /** Return the normal of the linear DoF axis.  \param axis must be 0 */
  virtual void _normalDoF(const BlockVector& q0, SiconosVector& ans, int axis,
                          bool absoluteRef=true);

public:

  /** Constructor based on one or two dynamical systems and an axis.
   *  \param d1 first DynamicalSystem linked by the joint.
   *  \param d2 second DynamicalSystem linked by the joint, or NULL
   *            for absolute frame.
   *  \param A SiconosVector of size 3 that defines the prismatic axis.
   *  \param absoluteRef if true, A is in the absolute frame,
   *                     otherwise A is in d1 frame.
   */
  PrismaticJointR(SP::SiconosVector axis, bool absoluteRef,
                  SP::NewtonEulerDS d1 = SP::NewtonEulerDS(),
                  SP::NewtonEulerDS d2 = SP::NewtonEulerDS());

  /** Initialize the joint constants based on the provided initial positions. */
  virtual void setInitialConditions(SP::SiconosVector q1,
                                    SP::SiconosVector q2 = SP::SiconosVector());

  void displayInitialPosition();

  void computeV1V2FromAxis();

  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(double time, BlockVector& q0, SiconosVector& y,
                           unsigned int axis);

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(double time, Interaction& inter,
                               SP::BlockVector q0, SimpleMatrix& jachq,
                               unsigned int axis);

  /** destructor
   */
  virtual ~PrismaticJointR() {};

  virtual void computeJachq(double time, Interaction& inter, SP::BlockVector q0 );

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  virtual void computeDotJachq(double time, BlockVector& workQ, BlockVector& workZ, BlockVector& workQdot);

  /* The options were    : operatorarrow */
  double H1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  /* The options were    : operatorarrow */
  double H2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  /* The options were    : operatorarrow */
  double H3(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  /* The options were    : operatorarrow */
  double H5(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  /* The options were    : operatorarrow */
  double H4(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
             double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13);

  void DotJd1d2(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11, double qdot12, double qdot13,
                double Xdot2, double Ydot2, double Zdot2, double qdot20, double qdot21, double qdot22, double qdot23);

  void DotJd2(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11, double qdot12, double qdot13,
              double X2, double Y2, double Z2, double qdot20, double qdot21, double qdot22, double qdot23);


  /** Get the number of constraints defined in the joint
      \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() { return 5; }

  /** Return the number of degrees of freedom of this joint.
      \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() { return 2; }

  /** Return the type of a degree of freedom of this joint.
      \return the type of the degree of freedom (DoF)
  */
  virtual DoF_Type typeOfDoF(unsigned int axis) {
    if (axis==0) return DOF_TYPE_LINEAR;
    else return DOF_TYPE_INVALID;
  };
};
#endif  //PrismaticJointRELATION_H
