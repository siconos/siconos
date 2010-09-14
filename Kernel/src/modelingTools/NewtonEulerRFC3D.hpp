/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
/*! \file NewtonEulerR.h

*/
#ifndef NEWTONEULERRELATIONFC3D_H
#define NEWTONEULERRELATIONFC3D_H

#include "NewtonEulerR.hpp"


class NewtonEulerRFC3D : public NewtonEulerR
{

protected:
  /*Point of contact*/
  SP::SimpleVector _Pc;
  /*Normal to the contact */
  SP::SimpleVector _Nc;
public:
  NewtonEulerRFC3D(): _Pc(new SimpleVector(3)), _Nc(new SimpleVector(3)), NewtonEulerR() {}

  /** destructor
   */
  virtual ~NewtonEulerRFC3D() {};

  /*default implementation consists in multiplying jachq and T*/
  virtual void computeJachqT();

};
TYPEDEF_SPTR(NewtonEulerRFC3D);
#endif // NEWTONEULERRELATIONFC3D_H
