/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
/*! \file NewtonEulerR.h

*/
#ifndef PivotJointRELATION_H
#define PivotJointRELATION_H

#include "SiconosKernel.hpp"
#include "KneeJoint.hpp"


class PivotJointR : public KneeJointR
{
public:
  static int _sNbEqualities;
protected:

  /*Axis coordonates*/
  double _Ax, _Ay, _Az;
  double _A1x, _A1y, _A1z;
  double _A2x, _A2y, _A2z;
  void buildA1A2();
public:
  /* constructor,
     \param a SP::NewtonEulerDS d1, a dynamical system containing the intial position
     \param a SP::NewtonEulerDS d2, a dynamical system containing the intial position
     \param a SP::SimpleVector P, P contains the coordinates of the Pivot point, in the frame of d1 where the origine is G1.
                                  ie P contains the coordinates of the Pivot point, in the object frame.
     \param a SP::SimpleVector A, Axis of the pivot in the frame of d1.
  */
  PivotJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SimpleVector P, SP::SimpleVector A);
  PivotJointR(SP::NewtonEulerDS d1, SP::SimpleVector P0, SP::SimpleVector A);
  /** destructor
   */
  virtual ~PivotJointR() {};


  virtual void computeh(double t);
protected:
  virtual void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  virtual void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13);
  double AscalA1(double q10, double q11, double q12, double q13, double q20, double q21, double q22, double q23);
  double AscalA2(double q10, double q11, double q12, double q13, double q20, double q21, double q22, double q23);

};
TYPEDEF_SPTR(PivotJointR);
#endif // BEAMRELATION_H
