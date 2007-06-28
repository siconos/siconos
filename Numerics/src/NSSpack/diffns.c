/* Siconos-Numerics version 2.1.0, Copyright INRIA 2005-2006.
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

/*
Input na, a, nb, b
Output nc, c
a and b: interger vectors in increasing order
c : vector of integers of a that are not in b.
\author Nineb Sheherazade & Dureisseix David.

 */

void diffns(int *na, int *a, int *nb, int * b, int *nc, int *c)
{

  int pta, ptb, ptc;
  int aa, i;





  pta = 0;
  ptb = 0;
  ptc = 0;


  if (*nb == 0)
  {

    for (i = 0 ; i < *na ; i++)
      c[i] = a[i];
    *nc  = *na;


  }
  else
  {

    for (i = 0 ; i < *na ; i++)
      c[i] = -1;


    while ((pta < *na) && (ptb < *nb))
    {

      aa  = a[pta];

      if (b[ptb] > aa)
      {

        c[ptc] = aa ;
        ptc    = ptc + 1 ;
        pta = pta + 1;
      }
      else if (b[ptb] == aa)
      {

        pta = pta + 1;

      }
      else
      {

        while ((b[ptb] < aa) && (ptb < *nb))
        {


          ptb = ptb + 1;

          if (ptb >= *nb)
          {

            c[ptc] = aa;
            ptc    = ptc + 1;

            break;

          }
        }

      }



    }



    for (i = pta + 1; i < *na ; i++)
    {


      c[ptc] = a[i];
      ptc = ptc + 1;
    }

    *nc = ptc;

  }

}
