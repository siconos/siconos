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

void pfc_3D_projc(int nc , double* mu , double *z , double *p , int *status)
{

  int i;

  for (i = 0 ; i < nc ; ++i)
  {

    /* No contact case */

    if (z[2 * i] < 0.0)
    {

      z[2 * i  ]  = 0.0;
      z[2 * i + 1]  = 0.0;
      status[i] = 0;
    }
    else
    {
      /* contact case */
      if (p[2 * i + 1] == 0.0)
      {
        /* sliding contact */
        if (z[2 * i + 1] > 0.0)
        {
          z[2 * i + 1] = mu[i] * z[2 * i];
          status[i] = 2;
        }
        else
        {
          z[2 * i + 1] = -mu[i] * z[2 * i];
          status[i] = 3;
        }
      }
      else
      {
        /* slide forward */
        if (z[2 * i + 1] < -mu[i]*z[2 * i])
        {
          z[2 * i + 1]  = -mu[i] * z[2 * i];
          status[i] = 3;
        }
        /* slide backward */
        else if (z[2 * i + 1] > mu[i]*z[2 * i])
        {
          z[2 * i + 1] = mu[i] * z[2 * i];
          status[i] = 2;
        }
        /* sticking contact */
        else
        {
          status[i] = 1;
        }
      }
    }
  }
}

