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
#ifndef PrismaticJointRELATION_H
#define PrismaticJointRELATION_H

#include "SiconosKernel.hpp"

class PrismaticJointR : public NewtonEulerR
{
public:
  static int _sNbEqualities;



public:
  SP::SiconosVector _axe0;

  SP::SiconosVector _V1;
  SP::SiconosVector _V2;
  double _G10G20d1x;
  double _G10G20d1y;
  double _G10G20d1z;
  double _V1x;
  double _V1y;
  double _V1z;
  double _V2x;
  double _V2y;
  double _V2z;
  double _q1cq202;
  double _q1cq203;
  double _q1cq204;

  SP::NewtonEulerDS _d1;
  SP::NewtonEulerDS _d2;

  /*axe is the axis of the prismatic joint, in the frame of the first DS, d1.*/
  PrismaticJointR(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2, SP::SiconosVector axe);
  /*axe is the axis of the prismatic joint, in the absolute frame.*/
  PrismaticJointR(SP::NewtonEulerDS d2, SP::SiconosVector axe);

  void computeFromInitialPosition();
  void displayInitialPosition();

  void computeV1V2FromAxe();

  /** destructor
   */
  virtual ~PrismaticJointR() {};



  virtual void computeJachq(const double time, Interaction& inter);


  virtual void computeh(const double time, Interaction& inter);


  /* The options were    : operatorarrow */
  double H1(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  /* The options were    : operatorarrow */
  double H2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  /* The options were    : operatorarrow */
  double H3(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);
  /* The options were    : operatorarrow */
  double H5(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  /* The options were    : operatorarrow */
  double H4(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  void Jd1d2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);

  void Jd2(double X1, double Y1, double Z1, double q10, double q11, double q12, double q13, double X2, double Y2, double Z2, double q20, double q21, double q22, double q23);






};


TYPEDEF_SPTR(PrismaticJointR)
#endif
