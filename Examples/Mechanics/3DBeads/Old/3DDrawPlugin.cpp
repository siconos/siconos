/* Siconos-sample version 2.0.0, Copyright INRIA 2005-2008.
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

const double R = 0.1; // beads radius
const double m = 1; // beads mass
const double g = 9.81; // gravity

extern "C" void gravity(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* param)
{
  for (unsigned int i = 0; i < sizeOfq; i++)
    fExt[i] = 0.0;

  fExt[2] = -m * g;
}


extern "C" void h0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double d = sqrt((q[6] - q[0]) * (q[6] - q[0]) + (q[7] - q[1]) * (q[7] - q[1]) + (q[8] - q[2]) * (q[8] - q[2]));
  y[0] =  d - 2 * R;
}

extern "C" void G0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{

  double d = sqrt((q[6] - q[0]) * (q[6] - q[0]) + (q[7] - q[1]) * (q[7] - q[1]) + (q[8] - q[2]) * (q[8] - q[2]));

  G[0]  = -(q[6] - q[0]) / d;
  G[3]  = -(q[7] - q[1]) / d;
  G[6]  = -(q[8] - q[2]) / d;
  G[18] = (q[6] - q[0]) / d;
  G[21] = (q[7] - q[1]) / d;
  G[24] = (q[8] - q[2]) / d;


  G[11] =  R;
  G[13] = -R;

  G[29] = -R;
  G[31] =  R;
}
