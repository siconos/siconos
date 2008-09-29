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

#include "DiskDiskR.h"

DiskDiskR::DiskDiskR(double r, double rr) : LagrangianScleronomousR()
{
  r1 = r;
  r2 = rr;
  r1pr2 = r1 + r2;
}

void DiskDiskR::computeH(double)
{

  SP::SiconosVector y = interaction->getYPtr(0);

  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *workX = *data["q0"];
  double *q = &(*workX)(0);

  double dx = q[3] - q[0];
  double dy = q[4] - q[1];

  double d =  sqrt(dx * dx + dy * dy);

  y->setValue(0, d - r1pr2);
  y->setValue(1, 0.);

};

void DiskDiskR::computeG(double, unsigned int)
{

  *workX = *data["q0"];
  double *q = &(*workX)(0);
  double *g = &(*(G[0]))(0, 0);

  double dx = q[3] - q[0];
  double dy = q[4] - q[1];

  double d = sqrt(dx * dx + dy * dy);

  double dxsd = dx / d;
  double dysd = dy / d;

  g[0] = -dxsd;       // dh[0]/dq[0]
  g[1] = -dysd;       // dh[1]/dq[0]
  g[2] = -dysd;       // dh[0]/dq[1]
  g[3] = dxsd;        // dh[1]/dq[1]
  g[4] = 0.;          // dh[0]/dq[2]
  g[5] = r2;          // dh[1]/dq[2]
  g[6] = dxsd;        // ...
  g[7] = dysd;
  g[8] = dysd;
  g[9] = -dxsd;
  g[10] = 0.;
  g[11] = r1;

}

