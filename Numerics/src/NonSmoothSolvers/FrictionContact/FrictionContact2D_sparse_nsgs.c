/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include "LA.h"
#include "FrictionContact2D_compute_error.h"


#define SGN(x) ((x) < 0 ? -1 : (x) > 0 ? 1 : 0)

/* in lcp_nsgs_SBM.c */
void buildLocalProblem(int rowNumber, const SparseBlockStructuredMatrix* const blmat,
                       LinearComplementarityProblem* local_problem,
                       double* q, double* z);

int frictionContact2DLocalSolve(double *W, double *q, double mu, double *P, double *U)
{
  double D, muPn;

  /* | Wnn Wnt |
     | Wtn Wtt | */

#define Wnn W[0]
#define Wtn W[1]
#define Wnt W[2]
#define Wtt W[3]

  if (q[0] > 0)
  {
    P[0] = 0;
    P[1] = 0;
  }
  else
  {
    /* solve WP + q = 0  */
    D = Wnn * Wtt - Wnt * Wtn;
    if (D < DBL_EPSILON) return(1);

    P[0] = - (Wtt * q[0] - Wnt * q[1]) / D;
    P[1] = - (-Wtn * q[0] + Wnn * q[1]) / D;

    muPn = mu * P[0];

    if (fabs(P[1]) > muPn)
      /* outside cone */
    {

      if (P[1] + muPn < 0)
      {

        P[0] = - q[0] / (Wnn - mu * Wnt);
        P[1] = - mu * P[0];
      }
      else
      {

        P[0] = - q[0] / (Wnn + mu * Wnt);
        P[1] = mu * P[0];

      }
    }
  }


#undef Wnn
#undef Wnt
#undef Wtn
#undef Wtt

  return(0);
}


void frictionContact2D_sparse_nsgs(FrictionContactProblem* problem, double *z, double *w,
                                   int *info, SolverOptions* options)
{
  /* Notes:
     - we suppose that the trivial solution case has been checked before,
     and that all inputs differs from NULL since this function is
     supposed to be called from lcp_driver_global().

     - Input matrix M of the problem is supposed to be sparse-block
     with no null row (ie no rows with all blocks equal to null)
  */

  assert(problem->M->matrix1);

  /*
    The options for the global "block" solver are defined in
    options[0].
   */

  /* Global Solver parameters*/
  int itermax = options[0].iparam[0];
  double tolerance = options[0].dparam[0];
  int erritermax = options[0].iparam[7];

  /* warning needed */
  if (erritermax > itermax || erritermax < 1) erritermax = 0;

  /* Matrix M/vector q of the LCP */
  SparseBlockStructuredMatrix* blmat = problem->M->matrix1;
  double * q = problem->q;

  assert(blmat->nbblocks >= 1);

  /* Local problem initialization */

  LinearComplementarityProblem * local_problem = (LinearComplementarityProblem *)
      malloc(sizeof(*local_problem));
  local_problem->M = (NumericsMatrix *)malloc(sizeof(*local_problem->M));
  local_problem->M->storageType = 0; // dense storage
  local_problem->M->matrix0 = NULL;
  local_problem->M->matrix1 = NULL;

  /* Memory allocation for q. Size of q = blsizemax, size of the
     largest square-block in blmat */

  int blsizemax = blmat->blocksize0[0];
  int k;
  for (unsigned int i = 1 ; i < blmat->blocknumber0 ; i++)
  {
    k = blmat->blocksize0[i] - blmat->blocksize0[i - 1];
    if (k > blsizemax) blsizemax = k;
  }
  local_problem->q = (double*)malloc(blsizemax * sizeof(double));

  /* Current row (of blocks) number */
  unsigned int rowNumber;

  /*****  Gauss-Seidel iterations *****/
  int iter = 0; /* Current iteration number */
  int erriter = 0; /* Current error iteration number */
  double error = INFINITY; /* Current error */
  int hasNotConverged = 1;

  /* Output from local solver */
  options[0].iparam[2] = 0;
  options[0].dparam[2] = 0.0;

  int pos = 0;

  /* Output from local solver */
  int infoLocal = -1;

  while ((iter < itermax) && hasNotConverged)
  {
    ++iter;
    ++erriter;

    /* Loop over the rows of blocks in blmat */
    for (pos = 0, rowNumber = 0; rowNumber < blmat->blocknumber0; ++rowNumber, ++pos, ++pos)
    {
      /* Local problem formalization */
      buildLocalProblem(rowNumber, blmat, local_problem, q, z);

      /* Solve local problem */
      infoLocal = frictionContact2DLocalSolve(local_problem->M->matrix0,
                                              local_problem->q,
                                              problem->mu[rowNumber],
                                              &z[pos], &w[pos]);


      if (infoLocal)
      {
        free(local_problem->q);
        free(local_problem->M);
        free(local_problem);
        /* Number of GS iterations */
        options[0].iparam[1] = iter;
        fprintf(stderr, "FrictionContact2D_nsgs error: local LCP solver failed at global iteration %d.\n for block-row number %d. Output info equal to %d.\n", iter, rowNumber, infoLocal);
        printf("FrictionContact2D_nsgs warning: local LCP solver failed at global iteration %d.\n for block-row number %d. Output info equal to %d.\n", iter, rowNumber, infoLocal);
        return;
      }

    }

    if (erriter >= erritermax)
    {
      erriter = 0;
      FrictionContact2D_compute_error(problem, z, w, tolerance, &error);
      hasNotConverged = error > tolerance  ;
    }
  }

  if (verbose > 0)
    printf("Siconos Numerics : problem size=%d, nb iterations=%d, error=%g\n",
           blmat->blocknumber0,
           iter,
           error);

  *info = hasNotConverged;

  /* Number of GS iterations */
  options[0].iparam[1] = iter;

  /* Resulting error */
  options[0].dparam[1] = error;

  free(local_problem->q);
  free(local_problem->M);
  free(local_problem);
}


int frictionContact2D_sparse_nsgs_setDefaultSolverOptions(SolverOptions *options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGS Solver\n");
  }

  /*strcpy(options->solverName,"NSGS");*/
  options->solverId = SICONOS_FRICTION_2D_NSGS;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 8;
  options->dSize = 8;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 8; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;
  options ->internalSolvers = NULL;
  return 0;
}
