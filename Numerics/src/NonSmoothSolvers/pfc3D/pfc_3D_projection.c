
/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include <time.h>
#include "pfc_3D_Solvers.h"

void pfc_3D_projection(int n , double *C , double *b , double *zz , double *ww , double coef , pfc3D_fPtr* Compute_G, pfc3D_fPtr* Compute_JacG, double* param1, double *param2, double *param3, int *iparam_local , double *dparam_local)
{

  double mrn, num, coef2;

  coef2 = coef * coef;

  if (b[0] > 0.)
  {
    zz[0] = 0.;
    zz[1] = 0.;
    zz[2] = 0.;
  }
  else
  {
    zz[0] = -b[0] / C[0 * n + 0];
    zz[1] = -b[1] / C[1 * n + 1];
    zz[2] = -b[2] / C[2 * n + 2];

    mrn = zz[1] * zz[1] + zz[2] * zz[2];

    if (mrn > coef2 * zz[0]*zz[0])
    {
      num = coef * zz[0] / sqrt(mrn);
      zz[1] = zz[1] * num;
      zz[2] = zz[2] * num;
    }
  }
}

