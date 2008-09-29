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

#include "DiskPlanR.h"

DiskPlanR::DiskPlanR(double r, double A, double B, double C) :
  LagrangianScleronomousR(),
  r(r), A(A), B(B), C(C)
{
  sqrA2pB2 = sqrt(A * A + B * B);
}

void DiskPlanR::computeH(double)
{

  SP::SiconosVector y = interaction->getYPtr(0);

  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *workX = *data["q0"];
  double *q = &(*workX)(0);

  y->setValue(0, fabs(A * q[0] + B * q[1] + C) / sqrA2pB2 - r);

}

void DiskPlanR::computeG(double, unsigned int)
{
  *workX = *data["q0"];
  double *q = &(*workX)(0);
  double *g = &(*(G[0]))(0, 0);

  double x0 = q[0];
  double y0 = q[1];
  double D1 = A * x0 + B * y0 + C;
  double D2 = sqrA2pB2 * fabs(D1);

  g[0] = A * D1 / D2;
  g[1] = -B * D1 / D2;
  g[2] = B * D1 / D2;
  g[3] = A * D1 / D2;
  g[4] = 0.;
  g[5] = -r;
}
