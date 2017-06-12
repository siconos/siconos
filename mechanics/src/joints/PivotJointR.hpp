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
/*! \file PivotJointR.hpp
*/
#ifndef PivotJointRELATION_H
#define PivotJointRELATION_H

#include "KneeJointR.hpp"

/** \class PivotJointR
 * \brief This class implements a pivots joint between one or two Newton/Euler Dynamical system. - Inherits from KneeJointR
 *
 */
class PivotJointR : public KneeJointR
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(PivotJointR);
  PivotJointR() : KneeJointR() {};

  /*Axis coordonates*/
  SP::SiconosVector _A;
  double _A1x, _A1y, _A1z;
  double _A2x, _A2y, _A2z;

  /*Initial conditions*/
  double _cq2q101, _cq2q102, _cq2q103, _cq2q104;
  double _initial_AscalA, _initial_AscalA1, _initial_AscalA2;

  /** Cumulative number of twists around the joint relative to initial
   * angular difference. */
  int _twistCount;    // TODO: Should be in a graph work vector?
  int _previousAngle; // Needed to track _twistCount, TODO: work vector?

  void buildA1A2();

  /** Compute the vector of linear and angular positions of the free axes */
  virtual void computehDoF(double time, BlockVector& q0, SiconosVector& y,
                           unsigned int axis);

  /** Compute the jacobian of linear and angular DoF with respect to some q */
  virtual void computeJachqDoF(double time, Interaction& inter,
                               SP::BlockVector q0, SimpleMatrix& jachq,
                               unsigned int axis);

  virtual void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  virtual void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13);

  void rot2to1(double q10, double q11, double q12, double q13,
               double q20, double q21, double q22, double q23,
               double *q2to1w, double *q2to1x,
               double *q2to1y, double *q2to1z);

  double AscalA1(double q2to1x, double q2to1y, double q2to1z);
  double AscalA2(double q2to1x, double q2to1y, double q2to1z);
  double AscalA(double q2to1x, double q2to1y, double q2to1z);

  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& work);

public:
  /* constructor,
     \param a SP::NewtonEulerDS d1, a dynamical system containing the intial position
     \param a SP::NewtonEulerDS d2, a dynamical system containing the intial position
     \param a SP::SiconosVector P, see KneeJointR documentation.
     \param a SP::SiconosVector A, Axis of the pivot in the frame of d1.
  */
  PivotJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SiconosVector P, SP::SiconosVector A);

  /* constructor,
     \param a SP::NewtonEulerDS d1, a dynamical system containing the intial position
     \param a SP::SiconosVector P, see KneeJointR documentation.
     \param a SP::SiconosVector A, axis in the frame of the object.
     \param a bool, used only by the KneeJointR constructor see KneeJointR documentation.
  */
  PivotJointR(SP::NewtonEulerDS d1, SP::SiconosVector P, SP::SiconosVector A, bool absoluteRef = true);

  /** destructor
   */
  virtual ~PivotJointR() {};

  SP::SiconosVector A() { return _A; }

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  /** Get the number of constraints defined in the joint
      \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() { return 5; }

  /** Return the number of degrees of freedom of this joint.
      \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() { return 1; }

  /** Return the type of a degree of freedom of this joint.
      \return the type of the degree of freedom (DoF)
  */
  virtual DoF_Type typeOfDoF(unsigned int axis) {
    if (axis==0) return DOF_TYPE_ANGULAR;
    else return DOF_TYPE_INVALID;
  };
};
#endif // PivotJointRELATION_H
