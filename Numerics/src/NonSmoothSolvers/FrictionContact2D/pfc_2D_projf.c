
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
/*!\file pfc_2D_projf.c
 *
 * pfc_2D_projf is a specific projection operator related to CPG (conjugated projected gradient) algorithm
 *              for primal contact problem with friction.\n
 * \author Shéhérazade Nineb.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void pfc_2D_projf(int etat[], int *nn, double y[], double fric[], double projf1[])
{

  int      i, nc, n = *nn;

  double   mina, maxa, bb;




  nc = n / 2;
  bb = 0.0;

  for (i = 0; i < nc; i++)
  {

    if (etat[i] == 0)                      /*  No contact status       */
    {
      if (y[2 * i] <=  0.0)
      {

        projf1[2 * i]   = 0.0;
        projf1[2 * i + 1] = 0.0;
      }
      else
      {
        projf1[2 * i]   = y[2 * i];
        projf1[2 * i + 1] = y[2 * i + 1];
      }
    }
    else if (etat[i] == 3)                  /*   Etat de contact glissant+ */
    {
      projf1[2 * i] = y[2 * i];

      if (y[2 * i + 1] > bb)
      {
        mina = bb;
      }
      else
      {
        mina = y[2 * i + 1];
      }

      projf1[2 * i + 1] = mina;
    }
    else if (etat[i] == 1)                  /*   Etat de contact glissant-  */
    {
      projf1[2 * i] = y[2 * i];

      if (y[2 * i + 1] < bb)
      {
        maxa = bb;
      }
      else
      {
        maxa = y[2 * i + 1];
      }

      projf1[2 * i + 1] = maxa;
    }
    else                                /*   Etat de contact adhérent   */
    {
      projf1[2 * i]   = y[2 * i];
      projf1[2 * i + 1] = y[2 * i + 1];
    }
  }

}




















