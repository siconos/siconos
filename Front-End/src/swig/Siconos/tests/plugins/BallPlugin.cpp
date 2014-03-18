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

const double R = 0.1; // ball radius
const double m = 1; // ball mass
const double g = 9.8; // gravity
const double h = -0.1; // altitude of the plan

extern "C" double FextFunction(double time)
{
  double res = -0.0;
  return res;
}


SICONOS_EXPORT void ballFExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfq; i++)
    fExt[i] = 0.0;
  fExt[0] = -m * g + FextFunction(time);
}

