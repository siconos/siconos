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
#include "SpherePlanR.h"

SpherePlanR::SpherePlanR(double r, double A, double B, double C, double D)
  : LagrangianScleronomousR(), r(r), A(A), B(B), C(C), D(D)
{
  sqrA2pB2pC2 = sqrt(A * A + B * B + C * C);

}

double SpherePlanR::distance(double x, double y, double z, double rad)
{
  return (fabs(A * x + B * y + C * z + D) / sqrA2pB2pC2 - rad);
}


void SpherePlanR::computeH(double)
{

  // Warning: temporary method to have contiguous values in memory,
  // copy of block to simple.
  *workX = *data[q0];
  double *q = &(*workX)(0);

  SP::SiconosVector y = getInteractionPtr()->getYPtr(0);

  y->setValue(0, distance(q[0], q[1], q[2], r));
  y->setValue(1, 0.);

};

void SpherePlanR::computeJacH(double, unsigned int)
{

  *workX = *data[q0];
  double *q = &(*workX)(0);
  double *g = &(*(JacH[0]))(0, 0);

  double D1 = A * q[0] + B * q[1] + C * q[2] + D;
  double D2 = sqrA2pB2pC2 * fabs(D1);

  double D1oD2 = D1 / D2;
  double AD1oD2 = A * D1oD2;
  double BD1oD2 = B * D1oD2;
  double CD1oD2 = C * D1oD2;

  g[0] = AD1oD2;
  g[3] = BD1oD2;
  g[6] = CD1oD2;
  g[11] = r;
}

