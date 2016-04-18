/* Siconos-Kernel  Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY siconos-team@lists.gforge.inria.fr
 */
/** \file KneeJointR.hpp
 */
#ifndef KneeJointRELATION_H
#define KneeJointRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerR.hpp>

/** \class KneeJointR
 *  \brief This class implements a knee joint between one or two Newton/Euler Dynamical system
 *
 */
class KneeJointR : public NewtonEulerR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(KneeJointR);
  KneeJointR(): NewtonEulerR() {};

  /** Coordinate of the knee point in the body frame of the first dynamical system _d1
   */
  SP::SiconosVector _P0;

  /** Pointers on the first concerned dynamical system*/
  SP::NewtonEulerDS _d1;
  /** Pointers on the second concerned dynamical system*/
  SP::NewtonEulerDS _d2;

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
     \param a SP::SiconosVector P, P contains the coordinates of the Knee point, in the absolute frame.
  */
  KneeJointR(SP::NewtonEulerDS d1, SP::SiconosVector P0, bool absolutRef = true);
  /** destructor
   */
  void checkInitPos(SP::SiconosVector q1, SP::SiconosVector q2);
  virtual ~KneeJointR() {};

  /** Get the number of constraints defined in the joint
      \return the number of constraints
   */
  static unsigned int numberOfConstraints() { return 3; }

  virtual void computeJachq(double time, Interaction& inter, VectorOfBlockVectors& DSlink);

  virtual void computeJachq(double time, Interaction& inter, SP::SiconosVector q1, SP::SiconosVector q2=SP::SiconosVector());

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  virtual void computeDotJachq(double time, SiconosVector& workQ, SiconosVector& workZ, SiconosVector& workQdot);

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
