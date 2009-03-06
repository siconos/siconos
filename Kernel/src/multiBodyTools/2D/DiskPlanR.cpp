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

#include <math.h>
#include "DiskPlanR.h"
#include "Disk.h"

void DiskPlanR::init(double r, double A, double B, double C,
                     double xC, double yC, double w)
{
  sqrA2pB2 = hypot(A, B);

  if (width == 0 || width == std::numeric_limits<double>::infinity())
    finite = false;
  else
    finite = true;
  if (fabs(A * xC + B * yC + C) > std::numeric_limits<double>::epsilon())
  {
    if (A == 0)
      // By+C=0
    {
      yCenter = -C / B;
    }
    else if (B == 0)
      // Ax+C=0
    {
      xCenter = -C / A;
    }
    else
      // Ax+By+C=0
    {
      if (xCenter != 0)
        yCenter = - (A * xCenter + C) / B;
      else
        xCenter = - (B * yCenter + C) / A;
    }
  }

  AC = A * C;
  B2 = B * B;
  A2 = A * A;
  AB = A * B;
  BC = B * C;
  x1 = xCenter + halfWidth * fabs(B) / sqrA2pB2;
  y1 =  yCenter - copysign(1, B) * halfWidth * A / sqrA2pB2;
  x2 = xCenter - halfWidth * fabs(B) / sqrA2pB2;
  y2 = yCenter + copysign(1, B) * halfWidth * A / sqrA2pB2;

  halfWidth = 0.5 * width;

}

DiskPlanR::DiskPlanR(double r, double A, double B, double C) :
  LagrangianScleronomousR(),
  r(r), A(A), B(B), C(C)
{
  init(r, A, B, C, 0, 0, std::numeric_limits<double>::infinity());
}

DiskPlanR::DiskPlanR(double r, double A, double B, double C,
                     double xC, double yC, double w) :
  LagrangianScleronomousR(),
  r(r), A(A), B(B), C(C),
  xCenter(xC), yCenter(yC), width(w)
{
  init(r, A, B, C, xCenter, yCenter, width);
}

DiskPlanR::DiskPlanR(double r, double xa, double ya, double xb, double yb)
{

  double d = xa * yb - xb * ya;
  init(r, (yb - ya) / d, (xb - xa) / d, 1, (xa + xb) / 2, (ya + yb) / 2, hypot(xa - xb, ya - yb));
}

double DiskPlanR::distance(double x, double y, double rad)
{
  if (finite)
  {
    double x0 = - (AC - B2 * x + AB * y) / (A2 + B2);
    double y0 = - (BC - A2 * y + AB * x) / (A2 + B2);

    if (hypot(xCenter - x0, yCenter - y0) >= halfWidth)
    {
      // ... and no jacH...
      double r = fmin(hypot(x - x1, y - y1), hypot(x - x2, y - y2));
      return r - rad;
    }
  }

  return (fabs(A * x + B * y + C) / sqrA2pB2 - rad);

}

void DiskPlanR::computeH(double)
{
  SP::SiconosVector y = getInteractionPtr()->getYPtr(0);

  double q_0 = boost::static_pointer_cast<Disk>
               (*getInteractionPtr()->dynamicalSystemsBegin())->getQ(0);
  double q_1 = boost::static_pointer_cast<Disk>
               (*(getInteractionPtr()->dynamicalSystemsBegin()))->getQ(1);

  y->setValue(0, distance(q_0, q_1, r));

}

void DiskPlanR::computeJacH(double, unsigned int)
{

  double *g = &(*(JacH[0]))(0, 0);

  double x0 = boost::static_pointer_cast<Disk>
              (*getInteractionPtr()->dynamicalSystemsBegin())->getQ(0);
  double y0 = boost::static_pointer_cast<Disk>
              (*(getInteractionPtr()->dynamicalSystemsBegin()))->getQ(1);
  double D1 = A * x0 + B * y0 + C;
  double D2 = sqrA2pB2 * fabs(D1);

  g[0] = A * D1 / D2;
  g[1] = -B * D1 / D2;
  g[2] = B * D1 / D2;
  g[3] = A * D1 / D2;
  g[4] = 0.;
  g[5] = -r;
}

bool DiskPlanR::equal(double pA, double pB, double pC, double pr)
{
  return (A == pA && B == pB && C == pC && r == pr);
}

bool DiskPlanR::equal(double pA, double pB, double pC, double pr,
                      double pXc, double pYc, double pw)
{
  if (finite)
    return (A == pA && B == pB && C == pC && r == pr &&
            pXc == xCenter && pYc == yCenter && pw == width);
  else
    return equal(pA, pB, pC, pr);
}

bool DiskPlanR::equal(DiskPlanR& odpr)
{
  if (finite)
    return (equal(odpr.getA(), odpr.getB(), odpr.getC(),
                  odpr.getRadius()));
  else
    return (equal(odpr.getA(), odpr.getB(), odpr.getC(),
                  odpr.getRadius(), odpr.getXCenter(),
                  odpr.getYCenter(), odpr.getWidth()));
}
