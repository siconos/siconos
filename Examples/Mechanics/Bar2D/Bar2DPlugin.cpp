/* Siconos-sample version 2.1.0, Copyright INRIA 2005-2006.
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

const double R = 1.0; // bar2D radius
const double m = 5.5135; // bar2D mass
const double g = 9.8; // gravity
const double h = -0.1; // altitude of the plan


extern "C" double FextFunction(double time)
{
  double res = -0.0;
  return res;
}


extern "C" void bar2DFExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* z)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector);
   * output parameter : fExt (pointer to Fext vector)
   */

  //  printf("Call of the function 'bar2DFExt' of the basic plugin.\nYou have to implement this function.\n");

  int i;
  int n = sizeOfq;

  for (i = 0; i < n; i++)
    fExt[i] = 0.0;

  fExt[0] = -m * g + FextFunction(time);
}

extern "C" void groundFExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* z)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector);
   * output patarmeter : fExt (pointer to Fext vector)
   */

  //  printf("Call of the function 'bar2DFExt' of the basic plugin.\nYou have to implement this function.\n");

  int i;
  int n = sizeOfq;

  for (i = 0; i < n; i++)
    fExt[i] = 0.0;
}

