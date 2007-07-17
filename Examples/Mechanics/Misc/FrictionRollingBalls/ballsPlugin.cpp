/* Siconos-sample version 2.1.1, Copyright INRIA 2005-2006.
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
const double omega = 1.4; // force pulsation
const double g = 10;
const double m = 1;
const double k = 1;
const double theta = 0;//0.78539816; //(pi/4)
const double F = 1; // force amplitude;

extern "C" void FExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* z)
{
  fExt[0] = F * sin(omega * time) * cos(theta);
  fExt[1] = F * sin(omega * time) * sin(theta);
  fExt[2] = -m * g;
}

extern "C" void FInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *fInt, unsigned int sizeZ, double * z)
{
  fInt[0] = k * sqrt((q[0] * q[0] + q[1] * q[1])) * cos(theta);
  fInt[1] = k * sqrt((q[0] * q[0] + q[1] * q[1])) * sin(theta);
  fInt[2] = 0;
}

extern "C" void NNL(unsigned int sizeOfq, const double *q, const double *velocity, double *NNL, unsigned int sizeZ, double *z)
{
  NNL[0] = 0.0;
  NNL[1] = 0.0;
  NNL[2] = 0.0;
}


extern "C" void Mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfq * sizeOfq; i++)
    mass[i] = 0.0;

  mass[0] = m;
  mass[4] = m;
  mass[8] = m;
}

extern "C" void jacobianQFInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  jacob[0] = k * cos(theta) * q[0] / (sqrt((q[0] * q[0] + q[1] * q[1])));
  jacob[1] = k * cos(theta) * q[1] / (sqrt((q[0] * q[0] + q[1] * q[1])));
  jacob[2] = k * sin(theta) * q[0] / (sqrt((q[0] * q[0] + q[1] * q[1])));
  jacob[3] = k * sin(theta) * q[1] / (sqrt((q[0] * q[0] + q[1] * q[1])));
}

extern "C" void jacobianVelocityFInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfq * sizeOfq; i++)
    jacob[i] = 0.0;
}

extern "C" void jacobianQNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfq * sizeOfq; i++)
    jacob[i] = 0.0;
}

extern "C" void jacobianVelocityNNL(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfq * sizeOfq; i++)
    jacob[i] = 0.0;
}
