/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
class DiskPlanR : public LagrangianScleronomousR, public std::enable_shared_from_this<DiskPlanR>
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

