/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include "LA.h"
#include <assert.h>

void buildLocalProblem(int rowNumber, const SparseBlockStructuredMatrix* const blmat, LinearComplementarityProblem* local_problem, double* q, double* z)
{

  assert(blmat->blocksize0[rowNumber] > 0);

  /* Position in vector blmat->block of the required diagonal block */
  int diagPos = getDiagonalBlockPos(blmat, rowNumber);
  /* Gets diagonal block = MLocal  */
  local_problem->M->matrix0 = blmat->block[diagPos];
  local_problem->size = blmat->blocksize0[rowNumber];
  int pos = 0;
  if (rowNumber != 0)
  {
    local_problem->size -= blmat->blocksize0[rowNumber - 1];
    pos =  blmat->blocksize0[rowNumber - 1];
  }

  local_problem->M->size0 = local_problem->size;
  local_problem->M->size1 = local_problem->size;

  /* Computes qLocal */
  /* qLocal = q[rowNumber] + sum (rowM.z),
     sum over all non-diagonal blocks, rowM being the current
     row of blocks of M
  */
  DCOPY(local_problem->size, &q[pos], 1, local_problem->q, 1);
  rowProdNoDiagSBM(blmat->blocksize0[blmat->blocknumber0 - 1], local_problem->size, rowNumber, blmat, z, local_problem->q, 0);

}

void lcp_nsgs_SBM(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  /* Notes:

     - we suppose that the trivial solution case has been checked
     before, and that all inputs differs from NULL since this function
     is supposed to be called from lcp_driver_global().

     - Input matrix M of the problem is supposed to be sparse-block
       with no null row (ie no rows with all blocks equal to null)
  */
  if (problem->M->matrix1 == NULL)
  {
    fprintf(stderr, "lcp_NSGS_SBM error: wrong storage type for input matrix M of the LCP.\n");
    exit(EXIT_FAILURE);
  }

  /*
    The options for the global "block" solver are defined in options[0].\n
    options[i], for 0<i<numberOfSolvers-1 correspond to local solvers.
   */

  /* Global Solver parameters*/
  int itermax = options[0].iparam[0];
  double tolerance = options[0].dparam[0];

  /* Matrix M/vector q of the LCP */
  SparseBlockStructuredMatrix* blmat = problem->M->matrix1;
  double * q = problem->q;

  /* Number of non-null blocks in blmat */
  int nbOfNonNullBlocks = blmat->nbblocks;
  if (nbOfNonNullBlocks < 1)
  {
    fprintf(stderr, "Numerics::lcp_NSGS_SBM error: empty M matrix (all blocks = NULL).\n");
    exit(EXIT_FAILURE);
  }

  /* Local problem initialization */

  LinearComplementarityProblem * local_problem = malloc(sizeof(*local_problem));
  local_problem->M = malloc(sizeof(*local_problem->M));
  local_problem->M->storageType = 0; // dense storage
  local_problem->M->matrix0 = NULL;
  local_problem->M->matrix1 = NULL;

  /* Memory allocation for q. Size of q = blsizemax, size of the largest square-block in blmat */
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
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  /* Output from local solver */
  options[0].iparam[2] = 0;
  options[0].dparam[2] = 0.0;

  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("lcp_nsgs_SBM", "The NSGS_SBM method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }



  /*Number of the local solver */
  int localSolverNum = options->numberOfInternalSolvers ;
  SolverOptions * internalSolvers = options->internalSolvers ;

  int pos = 0;
  /* Output from local solver */
  int infoLocal = -1;

  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    /* Loop over the rows of blocks in blmat */
    localSolverNum = 0;
    pos = 0;
    /*       DCOPY(problem->size,w,1,wBackup,1); */
    for (rowNumber = 0; rowNumber < blmat->blocknumber0; ++rowNumber)
    {
      /* Local problem formalization */
      buildLocalProblem(rowNumber, blmat, local_problem, q, z);
      /* Solve local problem */
      infoLocal = lcp_driver_DenseMatrix(local_problem, &z[pos], &w[pos], &internalSolvers[localSolverNum]);
      pos += local_problem->size;
      /* sum of local number of iterations (output from local_driver)*/
      if (options[localSolverNum].iparam != NULL)
        options[0].iparam[2] += internalSolvers[localSolverNum].iparam[1];
      /* sum of local errors (output from local_driver)*/
      options[0].dparam[2] += internalSolvers[localSolverNum].dparam[1];

      if (infoLocal > 0)
      {
        //free(local_problem->q);
        //free(local_problem->M);
        //free(local_problem);
        /* Number of GS iterations */
        options[0].iparam[1] = iter;
        fprintf(stderr, "lcp_NSGS_SBM error: Warning local LCP solver  at global iteration %d.\n for block-row number %d. Output info equal to %d.\n", iter, rowNumber, infoLocal);
        //exit(EXIT_FAILURE);

        break;

      }

      while (localSolverNum < options->numberOfInternalSolvers - 1)
        localSolverNum++;
    }

    /*       DCOPY(problem->size , problem->q , 1 , w , 1); */
    /*       prod(problem->size,problem->size, 1.0, problem->M,z,1.0,w); */
    /*       DAXPY(problem->size, -1.0, w,1,wBackup, 1); */
    /*       num = DNRM2(problem->size,wBackup,1); */
    /*       error = num*den; */
    /* Criterium convergence */
    hasNotConverged = lcp_compute_error(problem, z, w, tolerance, &error);
    /*       if(error<tolerance) hasNotConverged = 0; */
  }
  *info = hasNotConverged;
  /* Number of GS iterations */
  options[0].iparam[1] = iter;
  /* Resulting error */
  options[0].dparam[1] = error;

  free(local_problem->q);
  free(local_problem->M);
  free(local_problem);
  /*   free(wBackup); */
}



int linearComplementarity_nsgs_SBM_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGS Solver\n");
  }

  options->solverId = SICONOS_LCP_NSGS_SBM;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-6;
  options->internalSolvers = (SolverOptions*)malloc(options->numberOfInternalSolvers * sizeof(SolverOptions));

  linearComplementarity_pgs_setDefaultSolverOptions(options->internalSolvers);

  return 0;
}
