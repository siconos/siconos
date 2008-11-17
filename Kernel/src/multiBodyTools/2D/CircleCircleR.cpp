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
#include "CircleCircleR.h"

CircleCircleR::CircleCircleR(double r, double rr) : CircularR()
{
  r1 = r;
  r2 = rr;
  ar1mr2 = fabs(r1 - r2);
}

double CircleCircleR::distance(double x1, double y1, double r1, double x2, double y2, double r2)
{

  return (fabs(r1 - r2) - hypot(x1 - x2, y1 - y2));

}

void CircleCircleR::computeH(double)
{
  // Warning: temporary method to have contiguous values in memory,
  // copy of block to simple.
  *workX = *data[q0];
  double *q = &(*workX)(0);
  SP::SiconosVector y = getInteractionPtr()->getYPtr(0);

  y->setValue(0, distance(q[0], q[1], r1, q[3], q[4], r2));
  y->setValue(1, 0.);

};

void CircleCircleR::computeJacH(double, unsigned int)
{

  *workX = *data[q0];
  double *q = &(*workX)(0);
  double *g = &(*(JacH[0]))(0, 0);

  double dx = q[3] - q[0];
  double dy = q[4] - q[1];

  double d = hypot(dx, dy);

  double dxsd = dx / d;
  double dysd = dy / d;

  g[0] = dxsd;       // dh[0]/dq[0]
  g[1] = dysd;       // dh[1]/dq[0]
  g[2] = dysd;       // dh[0]/dq[1]
  g[3] = -dxsd;        // dh[1]/dq[1]
  g[4] = 0.;          // dh[0]/dq[2]
  g[5] = -r2;          // dh[1]/dq[2]
  g[6] = -dxsd;        // ...
  g[7] = -dysd;
  g[8] = -dysd;
  g[9] = dxsd;
  g[10] = 0.;
  g[11] = r1;

}

