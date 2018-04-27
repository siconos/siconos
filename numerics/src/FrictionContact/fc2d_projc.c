/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
/*!\file fc2d_projc.c
 *
 * fc2d_projc is a specific projection operator related to CPG (conjugated projected gradient) algorithm for global contact problem with friction.\n
 *
 * \author Sheherazade Nineb.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fc2d_Solvers.h"


void fc2d_projc(double xi[], int *nn, int statusi[], double pi[], double fric[], double *projc1, int *projc2)
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

