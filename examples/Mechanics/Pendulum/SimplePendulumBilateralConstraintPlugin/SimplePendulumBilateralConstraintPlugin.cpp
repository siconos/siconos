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
using namespace std;

#undef restrict
#define restrict __restrict

extern "C" double L;

SICONOS_EXPORT void h0(unsigned int sizeOfq,  double* restrict  q, unsigned int sizeOfY, double* restrict  y, unsigned int sizeZ, double* restrict  z)
{
  y[0] = pow(L, 2) - (pow(q[0], 2) + pow(q[1], 2));
}

SICONOS_EXPORT void G0(unsigned int sizeOfq,  double* restrict  q, unsigned int sizeOfY, double* restrict  G, unsigned int sizeZ, double* restrict  z)
{
  G[0] = -2.0 * q[0];
  G[1] = -2.0 * q[1];
}

SICONOS_EXPORT void G0dot(unsigned int sizeOfq,  double* restrict  q, unsigned int sizeOfqdot,  double* restrict  qdot, double* restrict  S, unsigned int sizeOfZ, double* restrict  z)
{
  S[0] = -2.0 * qdot[0];
  S[1] = -2.0 * qdot[1];
}
