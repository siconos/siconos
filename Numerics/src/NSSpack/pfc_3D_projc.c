/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
/*!\file pfc_2D_projc.c
 *
 * \fn  pfc_2D_projc( int n , double mu , double *z , double *p , int *status )
 *
 * pfc_2D_projc is a specific projection operator related to CPG (conjugated projected gradient) algorithm
 *              for primal contact problem with friction.\n
 *
 * Ref: Renouf, M. and Alart, P. "" Comp. Method Appl. Mech. Engrg. (2004).
 *
 * \param n       Unchanged parameter which represents the half dimension of the system.
 * \param mu      Unchanged parameter which represents the friction coefficient
 * \param z       Modified parameter which retruns the corrected iterate.
 * \param p       Unchanged parameter which contains the components of the descent direction.
 * \param status  Unchanged parameter which contains the vector status
 *
 * \author Mathieu Renouf.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void pfc_3D_projc(int nc , double mu , double *z , double *p , int *status)
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
          z[2 * i + 1] = mu * z[2 * i];
          status[i] = 2;
        }
        else
        {
          z[2 * i + 1] = -mu * z[2 * i];
          status[i] = 3;
        }
      }
      else
      {
        /* slide forward */
        if (z[2 * i + 1] < -mu * z[2 * i])
        {
          z[2 * i + 1]  = -mu * z[2 * i];
          status[i] = 3;
        }
        /* slide backward */
        else if (z[2 * i + 1] > mu * z[2 * i])
        {
          z[2 * i + 1] = mu * z[2 * i];
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

