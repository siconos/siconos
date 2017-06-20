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
/** \file KneeJointR.hpp
 */
#ifndef KneeJointRELATION_H
#define KneeJointRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerJointR.hpp>

/** \class KneeJointR
 *  \brief This class implements a knee joint between one or two Newton/Euler Dynamical system
 *
 */
class KneeJointR : public NewtonEulerJointR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(KneeJointR);
  KneeJointR(): NewtonEulerJointR() {};

  /** Coordinate of the knee point in the body frame of the first dynamical system _d1
   */
  SP::SiconosVector _P0;

  /**Absolute coodinates of the vector  G1P0 when d1 is located in q=(0,0,0,1,0,0,0)
   * i.e. P0 in the body frame of d1.
   * These values are computed when the constructor is called.
   */
  double _G1P0x;
  double _G1P0y;
  double _G1P0z;

  /** Absolute coodinates of the vector G2P0 when d2 is located in q=(0,0,0,1,0,0,0)
   *  i.e. P0 in the body frame of d2.
   * These values are computed when the constructor is called.
   */
  double _G2P0x;
  double _G2P0y;
  double _G2P0z;
  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

public:
  /* constructor,
     \param a SP::NewtonEulerDS d1, a dynamical system containing the intial position
     \param a SP::NewtonEulerDS d2, a dynamical system containing the intial position
     \param a SP::SiconosVector P, P contains the coordinates of the Knee point, in the frame of d1 where the origin is G1.
                                  ie P contains the coordinates of the Knee point, in the object frame G1.
  */
  KneeJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SiconosVector P);

  /* constructor,
     \param a SP::NewtonEulerDS d1, a dynamical system containing the intial position
     \param a SP::SiconosVector P, P contains the coordinates of the Knee point
     \param bool indicating whether P is in the absolute frame (=true, default)
            default) or the frame of the DS (=false)
  */
  KneeJointR(SP::NewtonEulerDS d1, SP::SiconosVector P, bool absoluteRef = true);

  /** destructor
   */
  void checkInitPos(SP::SiconosVector q1, SP::SiconosVector q2);
  virtual ~KneeJointR() {};

  /** Get the number of constraints defined in the joint
      \return the number of constraints
   */
  virtual unsigned int numberOfConstraints() { return 3; }

  /** Get the number of degrees of freedom defined in the joint
      \return the number of degrees of freedom (DoF)
   */
  virtual unsigned int numberOfDoF() { return 3; }

  /** Return the type of a degree of freedom of this joint.
      \return the type of the degree of freedom (DoF)
  */
  virtual DoF_Type typeOfDoF(unsigned int axis) {
    if (axis<3) return DOF_TYPE_ANGULAR;
    else return DOF_TYPE_INVALID;
  };

  virtual void computeJachq(double time, Interaction& inter, SP::BlockVector q0);

  
  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  virtual void computeDotJachq(double time, BlockVector& workQ, BlockVector& workZ, BlockVector& workQdot);

  virtual void computeDotJachq(double time, SP::SiconosVector qdot1, SP::SiconosVector qdot2=SP::SiconosVector());

  SP::SiconosVector P() { return _P0; }

protected:

  virtual void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
                     double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  virtual void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13);

  /* \warning, the following function should also depend on q */
  virtual void DotJd1d2(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11, double qdot12, double qdot13,
                        double Xdot2, double Ydot2, double Zdot2, double qdot20, double qdot21, double qdot22, double qdot23);
  virtual void DotJd1(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11, double qdot12, double qdot13);


//public:
  double Hx(double X1, double Y1, double Z1, double  q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  double Hy(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  double Hz(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
            double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
};
#endif // KneeJointRELATION_H
