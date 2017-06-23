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
/** \file FixedJointR.hpp
 */
#ifndef FixedJointRELATION_H
#define FixedJointRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerJointR.hpp>

/** \class FixedJointR
 *  \brief This class implements a fixed joint between one or two Newton/Euler Dynamical system
 *
 */
class FixedJointR : public NewtonEulerJointR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(FixedJointR);

  /*Initial conditions*/
  double _G10G20d1x, _G10G20d1y, _G10G20d1z;
  double _cq2q101, _cq2q102, _cq2q103, _cq2q104;

public:
  /** Empty constructor. The relation may be initialized later by
   * setInitialConditions. */
  FixedJointR() : NewtonEulerJointR() {};

  /* constructor,
     \param a SP::NewtonEulerDS d1, a dynamical system containing the initial position
     \param a SP::NewtonEulerDS d2, a dynamical system containing the initial position
  */
  FixedJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2 = SP::NewtonEulerDS());

  /** destructor
   */
  virtual ~FixedJointR() {};

  /** Initialize the joint constants based on the provided initial positions. */
  virtual void setInitialConditions(SP::SiconosVector q1,
                                    SP::SiconosVector q2 = SP::SiconosVector());

  /** Get the number of constraints defined in the joint
      \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() { return 6; }

  virtual void computeJachq(double time, Interaction& inter, SP::BlockVector q0);

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  virtual unsigned int numberOfDoF() { return 0; }

  virtual DoF_Type typeOfDoF(unsigned int axis) { return DOF_TYPE_INVALID; }

protected:

  virtual void Jd1d2(double X1, double Y1, double Z1,
                     double q10, double q11, double q12, double q13,
                     double X2, double Y2, double Z2,
                     double q20, double q21, double q22, double q23);

  virtual void Jd1(double X1, double Y1, double Z1,
                   double q10, double q11, double q12, double q13);
};

#endif // FixedJointRELATION_H
