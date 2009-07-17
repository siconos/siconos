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
#include "SphereSphereR.h"

SphereSphereR::SphereSphereR(double r, double rr) : LagrangianScleronomousR()
{
  r1 = r;
  r2 = rr;
  r1pr2 = r1 + r2;
}

double SphereSphereR::distance(double x1, double y1, double z1, double r1, double x2, double y2, double z2, double r2)
{
  double dx = x1 - x2;
  double dy = y1 - y2;
  double dz = z1 - z2;

  return (sqrt(dx * dx + dy * dy + dz * dz) - r1pr2);
}


void SphereSphereR::computeH(double)
{

  // Warning: temporary method to have contiguous values in memory,
  // copy of block to simple.
  *workX = *data[q0];
  double *q = &(*workX)(0);

  SP::SiconosVector y = getInteractionPtr()->getYPtr(0);

  y->setValue(0, distance(q[0], q[1], q[2], r1, q[6], q[7], q[8], r2));
  y->setValue(1, 0.);

};

void SphereSphereR::computeJacH(double, unsigned int)
{

  *workX = *data[q0];
  double *q = &(*workX)(0);
  double *g = &(*(JacH[0]))(0, 0);

  double dx = q[6] - q[0];
  double dy = q[7] - q[1];
  double dz = q[8] - q[2];

  double d = sqrt(dx * dx + dy * dy + dz * dz);

  double dmr1pr2 = d - r1pr2;

  double dxsd = dx / d;
  double dysd = dy / d;
  double dzsd = dz / d;

  g[0] = -dxsd;
  g[3] = -dysd;
  g[6] = -dzsd;
  g[18] = dxsd;
  g[21] = dysd;
  g[24] = dzsd;
  g[11] = r1;
  g[13] = -r2;
  g[29] = -r1;
  g[31] = r2;
}

