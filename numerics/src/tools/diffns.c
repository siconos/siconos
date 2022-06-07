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
#include "NSSTools.h"  // for diffns

void diffns(int *na, int *a, int *nb, int * b, int *nc, int *c)
{

  int pta, ptb, ptc;
  int aa, i;

  pta = 0;
  ptb = 0;
  ptc = 0;

  if(*nb == 0)
  {

    for(i = 0 ; i < *na ; i++)
      c[i] = a[i];
    *nc  = *na;

  }

  else
  {

    for(i = 0 ; i < *na ; i++)
      c[i] = -1;

    while((pta < *na) && (ptb < *nb))
    {

      aa  = a[pta];

      if(b[ptb] > aa)
      {

        c[ptc] = aa ;
        ptc    = ptc + 1 ;
        pta = pta + 1;
      }
      else if(b[ptb] == aa)
      {

        pta = pta + 1;

      }
      else
      {

        while((b[ptb] < aa) && (ptb < *nb))
        {


          ptb = ptb + 1;

          if(ptb >= *nb)
          {

            c[ptc] = aa;
            ptc    = ptc + 1;

            break;

          }
        }

      }



    }



    for(i = pta + 1; i < *na ; i++)
    {


      c[ptc] = a[i];
      ptc = ptc + 1;
    }

    *nc = ptc;

  }

}
