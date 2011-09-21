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

/*! \file SphereNEDSSphereNEDSR.hpp
  \brief Two spheres relation with the Newton Euler formalism and quaternions.
*/

#ifndef SphereNEDSSphereNEDSR_h
#define SphereNEDSSphereNEDSR_h

#include "NewtonEulerR.hpp"
#include "NewtonEulerFrom3DLocalFrameR.hpp"
class SphereNEDSSphereNEDSR : public NewtonEulerFrom3DLocalFrameR,
  public boost::enable_shared_from_this<SphereNEDSSphereNEDSR>
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SphereNEDSSphereNEDSR);

  double r1, r2, r1pr2;

  SphereNEDSSphereNEDSR();

public:

  /** Constructor

  \param disk1 radius
  \param disk2 radius
  */
  SphereNEDSSphereNEDSR(double, double);

  double distance(double, double, double, double, double, double, double, double);

  void computeh(double);

  //void computeJachq(double);

  /** visitors hook
   */
  ACCEPT_VISITORS();

};

TYPEDEF_SPTR(SphereNEDSSphereNEDSR);

#endif /* SphereNEDSSphereNEDSR_h */
