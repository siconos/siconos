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
#include "LA.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NonSmoothDrivers.h"

int filter_result_pfc(int nc, int Gsize, double * w, double* M, double *z , double *q, pfc3D_fPtr* G, double* mu, int* iparam, double* dparam, double** work)
{

  /*
     work[0] = qLocal
     work[1] = zBlock
     work[2] = temp copy of z
     work[3] = res
     work[4] = wLocal
     work[5] = currentBlock
  */
  //  double * res = (double*)malloc(Gsize*nc*sizeof(*res));
  //  double * wLocal = (double*)malloc(Gsize*sizeof(*wLocal));
  //  double * qLocal = (double*)malloc(Gsize*sizeof(*qLocal));
  /*   if ( res == NULL  || wLocal == NULL) */
  /*     { */
  /*       fprintf(stderr, "filter_result_pfc, memory allocation for MBlock failed.\n"); */
  /*       exit(EXIT_FAILURE); */
  /*     }  */

  double * p1, *p2, *p3; /*Useless parameters of G function ...*/
  double an, at, alpha, beta, det;
  double tol = dparam[1];
  //  double *zBlock   = (double*)malloc(Gsize*sizeof(*zBlock));
  int n = 3 * nc;
  //  double * zTmp = (double*)malloc(n*sizeof(*zTmp));

  int blockPos = 0;
  int i, j;

  //  double * currentBlock =(double*)malloc(Gsize*Gsize*sizeof(*currentBlock ));
  int blockNum = 0;
  int it, is;
  DCOPY(n , z , 1 , work[2] , 1);
  /* Loop through the contact points */
  for (i = 0 ; i < nc ; ++i)
  {
    blockPos = 3 * i;
    blockNum = i * nc;
    it = blockPos + 1;
    is = it + 1;
    work[5][0] = M[(blockPos) * n + blockPos];
    work[5][1] = M[(blockPos) * n + it];
    work[5][2] = M[(blockPos) * n + is];
    work[5][3] = M[(it) * n + blockPos];
    work[5][4] = M[(it) * n + it];
    work[5][5] = M[(it) * n + is];
    work[5][6] = M[(is) * n + blockPos];
    work[5][7] = M[(is) * n + it];
    work[5][8] = M[(is) * n + is];

    work[2][blockPos] = 0.0;
    work[2][it] = 0.0;
    work[2][is] = 0.0;
    work[0][0] = q[blockPos] + DDOT(n , &M[blockPos] , n , work[2] , 1);
    work[0][1] = q[it] + DDOT(n , &M[it] , n , work[2], 1);
    work[0][2] = q[is] + DDOT(n , &M[is] , n , work[2] , 1);
    work[2][blockPos] = z[blockPos];
    work[2][it] = z[it];
    work[2][is] = z[is];
    for (j = 0 ; j < Gsize ; ++j)
      work[1][j] = z[Gsize * i + j];

    an = 1. / work[5][0];
    alpha = work[5][Gsize + 1] + work[5][2 * Gsize + 2];
    det = work[5][Gsize + 1] * work[5][2 * Gsize + 2] - work[5][2 * Gsize + 1] + work[5][Gsize + 2];
    beta = alpha * alpha - 4 * det;
    if (beta > 0.)
      beta = sqrt(beta);
    else
      beta = 0.;

    at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));

    for (j = 0; j < Gsize; j++)
      work[4][j] = w[j + blockPos];
    (*G)(Gsize, &work[3][blockPos], work[1], work[5], work[4],  work[0], p1, p2, p3, an, at, mu[i]);
  }

  double error = DNRM2(nc * Gsize, work[3], 1);

  //  free (res);
  // free (wLocal);
  //  free (zBlock);
  //  free (zTmp);
  //  free (currentBlock);
  if (iparam[1] > 0)
    printf("PFC filter, error = %f\n" , error);

  if (error > tol)
    return 1;
  else return 0;
}
