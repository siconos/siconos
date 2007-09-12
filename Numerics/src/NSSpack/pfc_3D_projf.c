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
/*!\file pfc_2D_projf.c
 *
 * \fn  pfc_2D_projf( int n , double mu , double *ww , double *zz , double *rr , double *pp , int *status )
 *
 * pfc_2D_projf is a specific projection operator related to CPG (conjugated projected gradient) algorithm
 *              for primal contact problem with friction.\n
 *
 * Ref: Renouf, M. and Alart, P. "" Comp. Method Appl. Mech. Engrg. (2004).
 *
 * \param n       Unchanged parameter which represents the half dimension of the system.
 * \param ww      Modified parameter which returns the projected residue.
 * \param zz      Modified parameter which retruns the projected descent direction.
 * \param rr      Unchanged parameter which contains the components of the residue vector.
 * \param pp      Unchanged parameter which contains the components of the descent direction.
 * \param status  Unchanged parameter which contains the vector status
 *
 * \author Mathieu Renouf.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "NSSpack.h"

void pfc_3D_projf(int nc , double *ww , double *zz , double *rr , double *pp , int *status)
{

  int i;

  for (i = 0; i < nc; ++i)
  {

    /* no contact case */
    if (status[i] == 0)
    {
      if (rr[2 * i] < 0)
      {
        ww[2 * i]   = 0.0;
        zz[2 * i]   = 0.0;
        ww[2 * i + 1] = 0.0;
        zz[2 * i + 1] = 0.0;
      }
      else
      {
        ww[2 * i]   = rr[2 * i];
        ww[2 * i + 1] = 0.0;
        if (pp[2 * i] < 0)
        {
          zz[2 * i]   = 0.0;
          zz[2 * i + 1] = 0.0;
        }
        else
        {
          zz[2 * i]   = pp[2 * i];
          zz[2 * i + 1] = 0.0;
        }
      }
    }
    /* backward contact case */
    /* zt = mu*zn            */
    /* wt < 0                */
    else if (status[i] == 2)
    {
      ww[2 * i]   = rr[2 * i];
      zz[2 * i]   = pp[2 * i];
      ww[2 * i + 1] = fmin(0.0 , rr[2 * i + 1]);
      if (ww[2 * i + 1] == 0.0) zz[2 * i + 1] = 0.0;
      else zz[2 * i + 1] = fmin(0.0 , pp[2 * i + 1]);
    }
    /* forward contact case : */
    /* zt = -mu*zn            */
    /* wt > 0                 */
    else if (status[i] == 3)
    {
      ww[2 * i]   = rr[2 * i];
      zz[2 * i]   = pp[2 * i];
      ww[2 * i + 1] = fmax(0.0 , rr[2 * i + 1]);
      if (ww[2 * i + 1] == 0.0) zz[2 * i + 1] = 0.0;
      else zz[2 * i + 1] = fmax(0.0 , pp[2 * i + 1]);
    }
    /* sticking contact case */
    /* zt in [-mu*zn,mu*zn]  */
    /* wt = 0                */
    else
    {
      ww[2 * i]   = rr[2 * i];
      zz[2 * i]   = pp[2 * i];
      ww[2 * i + 1] = rr[2 * i + 1];
      zz[2 * i + 1] = pp[2 * i + 1];
    }
  }

}
