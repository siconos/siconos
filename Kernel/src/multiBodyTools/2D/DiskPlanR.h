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


/*! \file DiskPlanR.h
  \brief disk - plan relation - Inherits from LagrangianScleronomousR
*/

#ifndef DiskPlanR_h
#define DiskPlanR_h

#include "Interaction.h"
#include "LagrangianScleronomousR.h"

class DiskPlanR : public LagrangianScleronomousR, public boost::enable_shared_from_this<DiskPlanR>
{
private:
  double r, A, B, C, sqrA2pB2,
         AC, B2, A2, AB, BC, xCenter, yCenter, width, halfWidth, x1, x2, y1, y2;
  bool finite;

  void init(double, double, double, double, double, double, double);

  DiskPlanR();

public:

  /** Infinite Plan

  \param disk radius
  \param A
  \param B
  \param C
  */
  DiskPlanR(double, double, double, double);

  /** Finite or infinite Plan (segment)

    \param disk radius
    \param A
    \param B
    \param C
    \param xCenter
    \param yCenter
    \param width
    */
  DiskPlanR(double, double, double, double, double, double, double);

  /** Finite Plan
  */
  DiskPlanR(double, double, double, double, double);

  /* distance between disk and plan */
  double distance(double x, double y, double r);

  double getRadius()
  {
    return r;
  };

  double getA()
  {
    return A;
  };

  double getB()
  {
    return B;
  };

  double getC()
  {
    return C;
  };

  double getHypotAB()
  {
    return sqrA2pB2;
  };

  double getXCenter()
  {
    return xCenter;
  };

  double getYCenter()
  {
    return yCenter;
  };

  double getWidth()
  {
    return width;
  };

  void computeH(double);

  void computeJacH(double, unsigned int);

  bool equal(double, double, double, double);

  bool equal(double, double, double, double, double, double, double);

  bool equal(DiskPlanR&);

  bool isFinite()
  {
    return finite;
  };

  /** visitor hooks
   */
  virtual void accept(SiconosVisitor& tourist)
  {
    tourist.visit(*this);
  }
  virtual void accept(SP::SiconosVisitor tourist)
  {
    tourist->visit(shared_from_this());
  }

};

TYPEDEF_SPTR(DiskPlanR);

#endif /* DiskPlanR */

