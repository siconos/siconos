/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

/** \file CircleCircleR.hpp
 *
 */

#ifndef CircleCircleR_h
#define CircleCircleR_h

#include "CircularR.hpp"

/** \class CircleCircleR
 *  \brief Two disks relation - Inherits from LagrangianScleronomousR
 */
class CircleCircleR : public CircularR, public std11::enable_shared_from_this<CircleCircleR>
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(CircleCircleR);

  CircleCircleR() {};

public:

  /** Constructor
  \param rdisk1 radius
  \param rdisk2 radius
  */
  CircleCircleR(double rdisk1, double rdisk2);

  /** compute distance between 2 disks
      \param x1 x position of first disk
      \param y1 y position of first disk
      \param r1 radius of first disk
      \param x2 x position of second disk
      \param y2 y position of second disk
      \param r2 radius of second disk
      \return distance 
  */
  double distance(double x1, double y1, double r1,
                  double x2, double y2, double r2);

  using LagrangianScleronomousR::computeh;

  /** computeh implementation, see LagrangianScleronomousR
      \param time
      \param inter
  */
  void computeh(SiconosVector& q, SiconosVector& z, SiconosVector& y);

  /** computeh implementation, see LagrangianScleronomousR
      \param time
      \param inter
  */
  void computeJachq(SiconosVector& q, SiconosVector& z);

  /** visitors hook
   */
  ACCEPT_VISITORS();

  ~CircleCircleR() {};

};
#endif /* CircleCircleR_h */
