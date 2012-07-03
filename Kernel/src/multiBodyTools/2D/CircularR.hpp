/* Siconos-Example version 3.0.0, Copyright INRIA 2005-2011.
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
 * Foundation, Inc., 51 Franklin St, Fifth FLOOR, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 *
 */

/*! \file CircleCircleR.hpp
  \brief Two disks relation - Inherits from LagrangianScleronomousR
*/

#ifndef CircularR_h
#define CircularR_h

#include "Interaction.hpp"
#include "LagrangianScleronomousR.hpp"

class CircularR : public LagrangianScleronomousR
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(CircularR);

  double _r1, _r2;

  CircularR() : LagrangianScleronomousR() {};

public:

  /** Constructor

  \param disk1 radius
  \param disk2 radius
  */
  CircularR(double r1, double r2): _r1(r1), _r2(r2) {};

  double getRadius1()
  {
    return _r1;
  };

  double getRadius2()
  {
    return _r2;
  };

  virtual double distance(double, double, double, double, double, double)
  {
    assert(0);
    return(0);
  };

  ~CircularR() {};

};

TYPEDEF_SPTR(CircularR);
#endif /* CircularR_h */
