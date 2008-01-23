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
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include <time.h>


int pr_driver(double *vec, double *q, int *nn, method *pt, double *z, double *w)
{


  int info = -1, it_end;

  char prkey1[10] = "NLGS", prkey2[10] = "Latin";

  double res;

  clock_t t1, t2;


  t1 = clock();


  if (strcmp(pt->pr.name , prkey1) == 0)
  {
    pr_nlgs(vec, q, nn, pt->pr.a, pt->pr.b, & pt->pr.itermax, & pt->pr.tol, &pt->pr.chat, z, w, &it_end, &res, &info);

    pt->pr.err = res;
    pt->pr.iter = it_end;

  }
  else if (strcmp(pt->pr.name, prkey2) == 0)
  {
    pr_latin(vec, q, nn, &pt->pr.k_latin, pt->pr.a, pt->pr.b, &pt->pr.itermax, &pt->pr.tol, &pt->pr.chat, z, w, &it_end, &res, &info);

    pt->pr.err = res;
    pt->pr.iter = it_end;
  }

  else printf("Warning : Unknown solving method : %s\n", pt->pr.name);

  t2 = clock();


  printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

  return info;
}
