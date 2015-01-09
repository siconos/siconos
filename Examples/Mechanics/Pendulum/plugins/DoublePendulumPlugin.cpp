/* Siconos-sample , Copyright INRIA 2005-2011.
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

#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <stdio.h>
#include <math.h>

double gravity = 10.0;
double m1 = 1.0;
double m2 = 1.0 ;
double l1 = 1.0 ;
double l2 = 1.0 ;

SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  mass[0] = (m1 + m2) * l1;
  mass[1] = m2 * l1 * cos(q[0] - q[1]);
  mass[2] = m2 * l2 * cos(q[0] - q[1]);
  mass[3] = m2 * l2;
}

SICONOS_EXPORT void FGyr(unsigned int sizeOfq, const double *q, const double *velocity, double *FGyr, unsigned int sizeZ, double* z)
{
  FGyr[0] =  m2 * l2 * velocity[1] * velocity[1] * sin(q[0] - q[1]) ;
  FGyr[1] = -m2 * l1 * velocity[0] * velocity[0] * sin(q[0] - q[1]);
}

SICONOS_EXPORT void jacobianFGyrq(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  jacob[0] =  m2 * l2 * velocity[1] * velocity[0] * cos(q[0] - q[1]);
  jacob[1] =  -m2 * l1 * velocity[0] * velocity[0] * cos(q[0] - q[1]);
  jacob[2] =  -m2 * l2 * velocity[1] * velocity[1] * cos(q[0] - q[1]);
  jacob[3] =  m2 * l1 * velocity[0] * velocity[0] * cos(q[0] - q[1]);
}

SICONOS_EXPORT void jacobianVFGyr(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  jacob[0] =  0.0;
  jacob[1] = -2.0 * m2 * l1 * velocity[0] * sin(q[0] - q[1]);
  jacob[2] =  2.0 * m2 * l2 * velocity[1] * sin(q[0] - q[1]) ;
  jacob[3] =  0.0;
}

SICONOS_EXPORT void FInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *fInt, unsigned int sizeZ, double* z)
{
  fInt[0] = sin(q[0]) * gravity * (m1 + m2);
  fInt[1] = sin(q[1]) * gravity * m2;
}
SICONOS_EXPORT void jacobianFIntq(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  jacob[0] = cos(q[0]) * gravity * (m1 + m2);
  jacob[1] = 0.0;
  jacob[2] = 0.0;
  jacob[3] = cos(q[1]) * gravity * (m2);;
}

SICONOS_EXPORT void jacobianVFInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  jacob[0] = 0.0;
  jacob[1] = 0.0;
  jacob[2] = 0.0;
  jacob[3] = 0.0;
}

SICONOS_EXPORT void h0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* z)
{
  y[0] = l1 * sin(q[0]);
}

SICONOS_EXPORT void G0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* z)
{
  G[0] = l1 * cos(q[0]);
  G[1] = 0.0;
}

SICONOS_EXPORT void h1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* z)
{
  y[0] = l1 * sin(q[0]) + l2 * sin(q[1]);
}

SICONOS_EXPORT void G1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* z)
{
  G[0] = l1 * cos(q[0]);
  G[1] = l2 * cos(q[1]);
}

