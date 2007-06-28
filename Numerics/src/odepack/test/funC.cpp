/* Siconos-Numerics version 2.1.0, Copyright INRIA 2005-2006.
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
#include "blaslapack.h"

#define F77NAME(x) x##_


extern "C" void F77NAME(f1)(integer *sizeOfX, doublereal *time, doublereal *x, doublereal *xdot)
{
  /* input parameter : sizeOfX (size of the vector X); time ; x (pointer to X vector);
   * output parameter : xdot (pointer to Xdot vector)
   */

  //        printf("Call of the function 'f1' of the file funC.\n");

  if (sizeOfX[0] == 2)
  {
    xdot[0] = x[1];
    xdot[1] = 3.0 * (1.0 - x[0] * x[0]) * x[1] - x[0];
  }
  else
  {
    printf("OscillatorPlugin:vectorField --- Bad size of x. %i -- %i\n", *sizeOfX, sizeOfX);
  }

}

extern "C"void F77NAME(jac1)(integer *sizeOfX, doublereal *time, doublereal *x, integer* ml, integer *mu,  doublereal *jacob, integer *nrowpd)
{
  /* input parameter : sizeOfX (size of the vector X); time; x (pointer to x vector);
   * output parameter : jacob (pointer to JacobianX matrix)
   */

  //      printf("Call of the function 'jac1' of the  the file funC.c .\n");

  if (*sizeOfX == 2)
  {
    jacob[0] = 0.0;
    jacob[1] = -6.0 * x[0] * x[1] - 1.0;
    jacob[2] = 1.0;
    jacob[3] =  3.0 * (1.0 - x[0] * x[0]);
  }
  else
  {
    printf("OscillatorPlugin:computeJacobianX --- Bad size of x. \n");
  }
}
