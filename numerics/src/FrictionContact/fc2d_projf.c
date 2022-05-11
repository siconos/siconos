
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "fc2d_Solvers.h"  // for fc2d_projf

void fc2d_projf(int etat[], int *nn, double y[], double fric[], double projf1[])
{

  int      i, nc, n = *nn;

  double   mina, maxa, bb;




  nc = n / 2;
  bb = 0.0;

  for(i = 0; i < nc; i++)
  {

    if(etat[i] == 0)                       /*  No contact status       */
    {
      if(y[2 * i] <=  0.0)
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
    else if(etat[i] == 3)                   /*   Etat de contact glissant+ */
    {
      projf1[2 * i] = y[2 * i];

      if(y[2 * i + 1] > bb)
      {
        mina = bb;
      }
      else
      {
        mina = y[2 * i + 1];
      }

      projf1[2 * i + 1] = mina;
    }
    else if(etat[i] == 1)                   /*   Etat de contact glissant-  */
    {
      projf1[2 * i] = y[2 * i];

      if(y[2 * i + 1] < bb)
      {
        maxa = bb;
      }
      else
      {
        maxa = y[2 * i + 1];
      }

      projf1[2 * i + 1] = maxa;
    }
    else                                /*  Etat de contact adhérent   */
    {
      projf1[2 * i]   = y[2 * i];
      projf1[2 * i + 1] = y[2 * i + 1];
    }
  }

}




















