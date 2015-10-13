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

// for M_PI
#if defined(_MSC_VER)
#define _USE_MATH_DEFINES
#endif
#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <math.h>
#include <stdio.h>

SICONOS_EXPORT void prescribedvelocity(double time, unsigned int sizeofprescribedvelocity, double *pv)
{
  /* the plugin implements v(t) = C + A cos(omega *t) */

  double C = 0.0 ;
  double omega = M_PI / 2.0;
  double A = 1.0;

  pv[0] =  C +  A * cos(omega * time);
  //printf("prescribed velocity = %e\n", pv[0]);
}
