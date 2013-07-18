/* Siconos-sample version 2.0.0, Copyright INRIA 2005-2011.
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

#if defined(_MSC_VER)
#define _USE_MATH_DEFINES
#endif
#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <math.h>
#include <iostream>
#include <stdio.h>

//const double PI = 3.14159;


// ===== Dynamical System =====

// function to compute u
SICONOS_EXPORT void computeU(double time, unsigned int sizeU, double *U, unsigned int sizeZ, double* z)
{
  double f = 55000.0;
  if (time == 0.0)
    U[0] = 0;
  else
  {
    U[0] = z[0] * sin(2.0 * M_PI * f * time) / fabs(sin(2.0 * M_PI * f * time));
  }
}

