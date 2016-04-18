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
/*! \file PrismaticJointR.hpp

*/
#ifndef PrismaticJointRELATION_H
#define PrismaticJointRELATION_H

#include <MechanicsFwd.hpp>
#include <SiconosFwd.hpp>
#include <NewtonEulerR.hpp>

/** \class PrismaticJointR
 *  \brief This class implements a prismatic joint between one or two Newton/Euler Dynamical system
 *
 * From a given axis, we construct two unit othorgonal vectors to the axis V1 and V2 such that
 *  (axis,V1,V2) is an orthogonal frame
 *
 *
 *
 */
class PrismaticJointR : public NewtonEulerR
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(PrismaticJointR);
  PrismaticJointR(): NewtonEulerR() {};

public:
  /** Axis of the prismatic point in the inertial frame of reference
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

  double _q1cq202;
  double _q1cq203;
  double _q1cq204;

  /** constructor from two dynamical systems and an axis
   *  \param d1 first  DynamicalSystem link by the  joint
   *  \param d2 second  DynamicalSystem link by the joint
   * \param axis SiconosVector of size 3 that defines the prismatic axis
   *  in the body  frame of d1 ?
   */
  PrismaticJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SiconosVector axis);


  /** constructor from one dynamical systems and an axis
   *  \param d1 the  DynamicalSystem link by the  joint
   * \param axis SiconosVector of size 3 that defines the prismatic axis
   *  in the inertial frame of reference
   */
  PrismaticJointR(SP::NewtonEulerDS d1, SP::SiconosVector axis);


  void computeFromInitialPosition(SP::SiconosVector q2, SP::SiconosVector q1=SP::SiconosVector());

  void displayInitialPosition();

  void computeV1V2FromAxis();

  /** destructor
   */
  virtual ~PrismaticJointR() {};

  virtual void computeJachq(double time, Interaction& inter, VectorOfBlockVectors& DSlink);
  
  virtual void computeJachq(double time, Interaction& inter,
                            SP::SiconosVector q2, SP::SiconosVector q1=SP::SiconosVector() );

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

  void Jd2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13,
           double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  void DotJd1d2(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11, double qdot12, double qdot13,
                double Xdot2, double Ydot2, double Zdot2, double qdot20, double qdot21, double qdot22, double qdot23);

  void DotJd2(double Xdot1, double Ydot1, double Zdot1, double qdot10, double qdot11, double qdot12, double qdot13,
              double X2, double Y2, double Z2, double qdot20, double qdot21, double qdot22, double qdot23);


  /** Get the number of constraints defined in the joint
      \return the number of constraints
   */
  static unsigned int numberOfConstraints() { return 5; }

};
#endif  //PrismaticJointRELATION_H
