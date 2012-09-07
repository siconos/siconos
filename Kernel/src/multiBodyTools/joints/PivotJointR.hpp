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
/*! \file NewtonEulerR.hpp

*/
#ifndef PivotJointRELATION_H
#define PivotJointRELATION_H

#include "SiconosKernel.hpp"
#include "KneeJointR.hpp"


class PivotJointR : public KneeJointR
{
public:
  static int _sNbEqualities;
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(PivotJointR);


  /*Axis coordonates*/
  double _Ax, _Ay, _Az;
  double _A1x, _A1y, _A1z;
  double _A2x, _A2y, _A2z;
  void buildA1A2();

  virtual void initComponents(Interaction& inter);
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
     \param a SP::SiconosVector P0, see KneeJointR documentation.
     \param a SP::SiconosVector A, axis in the frame of the object.
     \param a bool, used only by the KneeJointR constructor see KneeJointR documentation.
  */
  PivotJointR(SP::NewtonEulerDS d1, SP::SiconosVector P0, SP::SiconosVector A, bool absolutRef = true);
  /** destructor
   */
  virtual ~PivotJointR() {};

  virtual void computeh(const double time, Interaction& inter);
protected:

  virtual void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  virtual void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13);
  double AscalA1(double q10, double q11, double q12, double q13, double q20, double q21, double q22, double q23);
  double AscalA2(double q10, double q11, double q12, double q13, double q20, double q21, double q22, double q23);

};
TYPEDEF_SPTR(PivotJointR)
#endif // BEAMRELATION_H
