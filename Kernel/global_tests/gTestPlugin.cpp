/* Siconos-sample version 1.1.3, Copyright INRIA 2005-2006.
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

// BOUNCING BALL

const double m = 1; // ball mass
const double g = 9.8; // gravity
extern "C" void ballFExt(const unsigned int *sizeOfq, const double *time, double *fExt, double* param)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector);
   * output parameter : fExt (pointer to Fext vector)
   */

  //  printf("Call of the function 'ballFExt' of the basic plugin.\nYou have to implement this function.\n");

  int i;
  int n = *sizeOfq;
  for (i = 0; i < n; i++)
  {
    fExt[i] = 0.0;
  }

  fExt[0] = -m * g;
}

// BALLBOWL

const double R = 0.5; // ball radius


extern "C" void groundFExt(const unsigned int *sizeOfq, const double *time, double *fExt, double* param)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector);
   * output patarmeter : fExt (pointer to Fext vector)
   */

  //  printf("Call of the function 'ballFExt' of the basic plugin.\nYou have to implement this function.\n");

  int i;
  int n = *sizeOfq;

  for (i = 0; i < n; i++)
    fExt[i] = 0.0;
}

extern "C" void h0(const unsigned int* sizeOfq, const double* q, const unsigned int* sizeOfY, double* y, double* param)
{
  double R0 = 0.0;
  y[0] = q[0] + sqrt(R * R - q[1] * q[1]) - R0;
}

extern "C" void G0(const unsigned int* sizeOfq, const double* q, const unsigned int* sizeOfY, double* G, double* param)
{
  G[0] = 1.0;
  G[1] = -q[1] / (sqrt(R * R - q[1] * q[1]));
}


