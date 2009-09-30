/* Siconos-Example version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 *
 */

/*! \file SphereSphereR.h
  \brief Two spheres relation - Inherits from LagrangianScleronomousR
*/

#ifndef SphereSphereR_h
#define SphereSphereR_h

#include "LagrangianScleronomousR.h"

class SphereSphereR : public LagrangianScleronomousR, public boost::enable_shared_from_this<SphereSphereR>
{
private:
  double r1, r2, r1pr2;

  SphereSphereR();

public:

  /** Constructor

  \param disk1 radius
  \param disk2 radius
  */
  SphereSphereR(double, double);

  double distance(double, double, double, double, double, double, double, double);

  void computeH(double);

  void computeJacQH(double);

  /** visitors hook
   */
  ACCEPT_VISITORS();

};

TYPEDEF_SPTR(SphereSphereR);

#endif /* SphereSphereR_h */
