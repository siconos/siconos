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
#ifndef KneeJointRELATION_H
#define KneeJointRELATION_H

#include "SiconosKernel.hpp"

class KneeJointR : public NewtonEulerR
{
public:
  static int _sNbEqualities;
protected:
  /*coodinate of the Knee point in the frame of d1*/
  SP::SimpleVector _P0;

  SP::NewtonEulerDS _d1;
  SP::NewtonEulerDS _d2;
  /*absolut coodinate G1P0 when d1 is located in q=(0,0,0,1,0,0,0). ie P0 in the frame of the object d1.*/
  double _G1P0x;
  double _G1P0y;
  double _G1P0z;
  /*absolut coodinate G2P0 when d2 is located in q=(0,0,0,1,0,0,0). ie P0 in the frame of the object d2.*/
  double _G2P0x;
  double _G2P0y;
  double _G2P0z;

public:
  /* constructor,
     \param a SP::NewtonEulerDS d1, a dynamical system containing the intial position
     \param a SP::NewtonEulerDS d2, a dynamical system containing the intial position
     \param a SP::SimpleVector P, P contains the coordinates of the Knee point, in the frame of d1 where the origine is G1.
                                  ie P contains the coordinates of the Knee point, in the object frame G1.
  */
  KneeJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SimpleVector P);
  /* constructor,
     \param a SP::NewtonEulerDS d1, a dynamical system containing the intial position
     \param a SP::SimpleVector P, P contains the coordinates of the Knee point, in the absolute frame.
  */
  KneeJointR(SP::NewtonEulerDS d1, SP::SimpleVector P0, bool absolutRef = true);
  /** destructor
   */
  void checkInitPos();
  virtual ~KneeJointR() {};

  virtual void computeJachq(double t);
  virtual void computeh(double t);
protected:
  virtual void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  virtual void Jd1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13);
public:
  double Hx(double X1, double Y1, double Z1, double  q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  double Hy(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  double Hz(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
};
TYPEDEF_SPTR(KneeJointR);
#endif // BEAMRELATION_H
