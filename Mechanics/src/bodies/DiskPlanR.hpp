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


/** \file DiskPlanR.hpp
 */

#ifndef DiskPlanR_h
#define DiskPlanR_h

#include "MechanicsFwd.hpp"
#include "LagrangianScleronomousR.hpp"

/** \class DiskPlanR
 *  \brief disk - plan relation - Inherits from LagrangianScleronomousR
 */
class DiskPlanR : public LagrangianScleronomousR, public std11::enable_shared_from_this<DiskPlanR>
{
private:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(DiskPlanR);

  double r, A, B, C, sqrA2pB2,
    AC, B2, A2, AB, BC, xCenter, yCenter, width, halfWidth, x1, x2, y1, y2;
  bool finite;

  void init(double, double, double, double, double, double, double);

  DiskPlanR() : LagrangianScleronomousR() {};

public:

  /** Infinite Plan

  \param r disk radius
  \param A component of line equation Ax + By + C = 0
  \param B component of line equation Ax + By + C = 0
  \param C component of line equation Ax + By + C = 0
  */
  DiskPlanR(double r, double A, double B, double C);

  /** Finite or infinite Plan (segment)

    \param disk radius
    \param A
    \param B
    \param C
    \param xCenter
    \param yCenter
    \param width
    */
  DiskPlanR(double disk, double A, double B, double C,
            double xCenter, double yCenter, double width);

  /** Finite Plan
  */
  DiskPlanR(double, double, double, double, double);

  /* distance between disk and plan */
  double distance(double x, double y, double r) const;

  double getRadius() const
  {
    return r;
  };

  double getA() const
  {
    return A;
  };

  double getB() const
  {
    return B;
  };

  double getC() const
  {
    return C;
  };

  double gethypotAB() const
  {
    return sqrA2pB2;
  };

  double getXCenter() const
  {
    return xCenter;
  };

  double getYCenter() const
  {
    return yCenter;
  };

  double getWidth() const
  {
    return width;
  };

  using LagrangianScleronomousR::computeh;
  void computeh(SiconosVector& q, SiconosVector& z, SiconosVector& y);

  void computeJachq(SiconosVector& q, SiconosVector& z);

  bool equal(double, double, double, double) const;

  bool equal(double, double, double, double, double, double, double) const;

  bool equal(const DiskPlanR&) const;

  bool isFinite() const
  {
    return finite;
  };

  /** visitor hooks
   */
  ACCEPT_VISITORS();

  ~DiskPlanR() {};

};
#endif /* DiskPlanR */

