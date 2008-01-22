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
#include <time.h>
#ifndef MEXFLAG
#include "NSSpack.h"
#endif

int pfc_2D_solver(int n, double *vec, double *q, method *pt , double *z , double *w, double* mu)
{

  char pfckey1[10] = "NLGS", pfckey2[10] = "CPG", pfckey3[10] = "Latin";

  int i, info;

  int     iparam[5];
  double  dparam[4];

  clock_t t1, t2;

  for (i = 0 ; i < 5 ; ++i) iparam[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparam[i] = 0.0;

  info    = -1;

  t1 = clock();

  if (strcmp(pt->pfc_2D.name , pfckey1) == 0)
  {

    iparam[0] = pt->pfc_2D.itermax;
    iparam[1] = pt->pfc_2D.chat;
    dparam[0] = pt->pfc_2D.tol;

    pfc_2D_nlgs(n , vec , q , z , w , mu, &info , iparam , dparam);

    pt->pfc_2D.iter = iparam[2];
    pt->pfc_2D.err  = dparam[1];

  }
  else if (strcmp(pt->pfc_2D.name , pfckey2) == 0)
  {

    iparam[0] = pt->pfc_2D.itermax;
    iparam[1] = pt->pfc_2D.chat;
    dparam[0] = pt->pfc_2D.tol;

    pfc_2D_cpg(n , vec , q , z , w , mu, &info , iparam , dparam);

    pt->pfc_2D.iter = iparam[2];
    pt->pfc_2D.err  = dparam[1];

  }
  else if (strcmp(pt->pfc_2D.name , pfckey3) == 0)
  {

    iparam[0] = pt->pfc_2D.itermax;
    iparam[1] = pt->pfc_2D.chat;
    dparam[0] = pt->pfc_2D.tol;
    dparam[1] = pt->pfc_2D.k_latin;

    pfc_2D_latin(n , vec , q , z , w , mu, &info , iparam , dparam);

    pt->pfc_2D.iter = iparam[2];
    pt->pfc_2D.err  = dparam[2];

  }
  else printf("Warning : Unknown solving method : %s\n", pt->pfc_2D.name);

  t2 = clock();

  return info;

}
