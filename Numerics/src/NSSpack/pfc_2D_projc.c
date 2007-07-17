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
 * \fn pfc_2D_projc ( double xi[], int *nn, int statusi[], double pi[], double fric[], double *projc1, int *projc2 )

 *
 * pfc_2D_projc is a specific projection operator related to CPG (conjugated projected gradient) algorithm for primal contact problem with friction.\n
 *
 *
 * \param xi        On enter, the intermediate iterate which goes to be projected (projc1).
 * \param nn        On enter, the dimension of the system.
 * \param statusi   On enter, a vector which contains the initial status.
 * \param pi        On enter, a vector which contains the components of the descent direction.
 * \param fric      On enter, a vector which contains the friction coefficient.
 * \param projc1    On return, the corrected iterate.
 * \param projc2    On return, the new status.
 *
 * \author Sheherazade Nineb.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void pfc_2D_projc(double xi[], int *nn, int statusi[], double pi[], double fric[], double *projc1, int *projc2)

{

  int     i, nc, n = *nn, stat;

  double  mu1;




  nc  = n / 2;



  for (i = 0 ; i < nc ; i++)
  {
    mu1  = fric[i];
    stat = statusi[i];

    if (xi[2 * i] <= 0.0)                       /* No contact status  */
    {

      projc1[2 * i]   = 0.0;
      projc1[2 * i + 1] = 0.0;
      projc2[i]     = 0;

    }
    else
    {
      projc1[2 * i] = xi[2 * i];

      if (xi[2 * i + 1] <= -mu1 * xi[2 * i])   /*  Slide backward     */
      {

        projc1[2 * i + 1] = -mu1 * xi[2 * i] ;
        projc2[i]     = 1;

      }
      else if (xi[2 * i + 1] >= mu1 * xi[2 * i]) /*  Slide forward      */
      {
        projc1[2 * i + 1] = mu1 * xi[2 * i];
        projc2[i]     = 3;
      }
      else
      {
        if (pi[2 * i + 1] == 0.0)
        {
          if (stat == 1)                     /*  Slide backward     */
          {
            projc1[2 * i + 1] = -mu1 * xi[2 * i];
            projc2[i]     = 1;
          }
          else if (stat == 3)                  /*  Slide forward        */
          {
            projc1[2 * i + 1] = mu1 * xi[2 * i];
            projc2[i]     = 3;
          }
          else
            /*   Stick contact        */
          {
            projc1[2 * i + 1] = xi[2 * i + 1];
            projc2[i]     = 2;
          }
        }
        else
          /*   Stick contact      */
        {
          projc1[2 * i + 1]   = xi[2 * i + 1];
          projc2[i]       = 2;
        }
      }

    }

  }

}

