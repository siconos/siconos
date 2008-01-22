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
#include "NSSpack.h"

int filter_result_pfc_block2(int nc, int Gsize, double * w, SparseBlockStructuredMatrix *M, double *z , double *qLocal , pfc3D_fPtr* G, double* mu, int* iparam, double* dparam, double *work)
{

  // work[0... Gsize*nc-1] is used to save result form G function call.
  if (work == NULL)
  {
    fprintf(stderr, "filter_result_pfc_block, no memory allocated for work vector.\n");
    exit(EXIT_FAILURE);
  }
  double * p1, *p2, *p3; /*Useless parameters of G function ...*/
  double an, at, alpha, beta, det;
  double tol = dparam[1];

  int blockPos = 0;
  int i, j;
  int diagBlockPos = 0; /* index of the current diagonal block of M */
  double * currentBlock;
  int blockNum = 0;
  /* Loop through the contact points */
  int pos;
  for (i = 0 ; i < nc ; ++i)
  {
    blockPos = 3 * i;
    blockNum = i * nc + i + 1;
    currentBlock = M->block[diagBlockPos];

    /* Loop through the columns(blocks) of M to compute qLocal */
    pos = 3 * (i + 1);
    for (j = i + 1; j < nc ; ++j)
    {
      DGEMV(LA_NOTRANS, 3, 3, 1.0, M->block[blockNum], Gsize, &z[pos], 1 , 1.0, &qLocal[blockPos], 1);
      pos += 3;
      blockNum = blockNum + 1;
    }

    an = 1. / currentBlock[0];
    alpha = currentBlock[Gsize + 1] + currentBlock[2 * Gsize + 2];
    det = currentBlock[Gsize + 1] * currentBlock[2 * Gsize + 2] - currentBlock[2 * Gsize + 1] + currentBlock[Gsize + 2];
    beta = alpha * alpha - 4 * det;
    if (beta > 0.)
      beta = sqrt(beta);
    else
      beta = 0.;

    at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));

    (*G)(Gsize, &work[blockPos], &z[Gsize * i], currentBlock, &w[blockPos],  &qLocal[blockPos], p1, p2, p3, an, at, mu[i]);
    diagBlockPos += 1 + nc;
  }

  double error = DNRM2(nc * Gsize, work , 1);
  if (iparam[1] > 0)
    printf("PFC filter, error = %f\n" , error);

  if (error > tol) return 1;
  else return 0;
}
