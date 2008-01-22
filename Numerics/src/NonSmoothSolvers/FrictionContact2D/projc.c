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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void projc(double xi[], int *nn, int statusi[], double pi[], double fric[], double *projc1, int *projc2)
{

  int i, nc, n = *nn, stat;
  double  mu1, eps;




  nc = n / 2;
  eps = 1.e-10;


  for (i = 0; i < nc; i++)
  {
    mu1 = fric[i];
    stat = statusi[i];
    if (xi[2 * i] <= 0.0) /*/ !status de non contact*/
    {
      projc1[2 * i] = 0.0;
      projc1[2 * i + 1] = 0.0;
      projc2[i] = 0;
    }
    else
    {
      projc1[2 * i] = xi[2 * i];
      /*/       printf("else\n");  */
      if (xi[2 * i + 1] <= -mu1 * xi[2 * i]) /*/!slide backward*/
      {
        projc1[2 * i + 1] = -mu1 * xi[2 * i] ;
        projc2[i] = 1;
      }
      else if (xi[2 * i + 1] >= mu1 * xi[2 * i]) /*/!slide forward*/
      {
        projc1[2 * i + 1] = mu1 * xi[2 * i];
        projc2[i] = 3;
      }
      else
      {
        if (pi[2 * i + 1] == 0.0)
        {
          if (stat == 1) /*/ !slide backward*/
          {
            projc1[2 * i + 1] = -mu1 * xi[2 * i];
            projc2[i] = 1;
          }
          else if (stat == 3) /*/!slide forward*/
          {
            projc1[2 * i + 1] = mu1 * xi[2 * i];
            projc2[i] = 3;
          }
          else
            /*/     !stick contact*/
          {
            projc1[2 * i + 1] = xi[2 * i + 1];
            projc2[i] = 2;
          }
        }
        else
          /*/      !stick contact*/
        {
          projc1[2 * i + 1] = xi[2 * i + 1];
          projc2[i] = 2;
        }
      }

    }


  }

}

