/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>
/* #ifndef MEXFLAG */
/* #include "NonSmoothDrivers.h" */
/* #endif */
#include "SparseBlockMatrix.h"
#include "LinearComplementarityProblem.h"
#include "MixedLinearComplementarityProblem.h"
#include "NumericsMatrix.h"
#include "SolverOptions.h"
#include "SiconosBlas.h"
#include "LCP_Solvers.h"
#include "MLCP_Solvers.h"
#include "SiconosCompat.h"
#include <assert.h>
#include "numerics_verbose.h"

static void mlcp_pgs_sbm_buildLocalProblem(int rowNumber, const SparseBlockStructuredMatrix* const blmat, LinearComplementarityProblem* local_problem, double* q, double* z);

void mlcp_pgs_sbm_buildLocalProblem(int rowNumber, const SparseBlockStructuredMatrix* const blmat, LinearComplementarityProblem* local_problem, double* q, double* z)
{

  assert(blmat->blocksize0[rowNumber] > 0);

  /* Position in vector blmat->block of the required diagonal block */
  int diagPos = SBM_get_position_diagonal_block(blmat, rowNumber);
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
  cblas_dcopy(local_problem->size, &q[pos], 1, local_problem->q, 1);
  SBM_row_prod_no_diag(blmat->blocksize0[blmat->blocknumber0 - 1], local_problem->size, rowNumber, blmat, z, local_problem->q, 0);

}

void mlcp_pgs_SBM(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
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
    fprintf(stderr, "mlcp_NSGS_SBM error: wrong storage type for input matrix M of the LCP.\n");
    exit(EXIT_FAILURE);
  }

  /*
    The options for the global "block" solver are defined in options[0].\n
    options[i], for 0<i<numberOfSolvers-1 correspond to local solvers.
   */

  /* Global Solver parameters*/
  int itermax = options[0].iparam[0];
  double tolerance = options[0].dparam[0];

  /* Matrix M/vector q of the MLCP */
  SparseBlockStructuredMatrix* blmat = problem->M->matrix1;
  double * q = problem->q;

  /* Number of non-null blocks in blmat */
  int nbOfNonNullBlocks = blmat->nbblocks;
  if (nbOfNonNullBlocks < 1)
  {
    fprintf(stderr, "Numerics::mlcp_NSGS_SBM error: empty M matrix (all blocks = NULL).\n");
    exit(EXIT_FAILURE);
  }

  /* Local problem initialization */
  /* We choose to set a LCP problem which can be solved by a Linear System solver if the
     block corresponds to equalities */
  LinearComplementarityProblem * local_problem = (LinearComplementarityProblem *)malloc(sizeof(*local_problem));
  local_problem->M = (NumericsMatrix *)malloc(sizeof(*local_problem->M));
  local_problem->M->storageType = 0; // dense storage
  local_problem->M->matrix0 = NULL;
  local_problem->M->matrix1 = NULL;
  local_problem->M->matrix2 = NULL;
  local_problem->M->internalData = NULL;

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
    numerics_error("mlcp_nsgs_SBM", "The MLCP_PGS_SBM method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
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
    /*       cblas_dcopy(problem->size,w,1,wBackup,1); */
    for (rowNumber = 0; rowNumber < blmat->blocknumber0; ++rowNumber)
    {
      /* Local problem formalization */
      mlcp_pgs_sbm_buildLocalProblem(rowNumber, blmat, local_problem, q, z);
      /* Solve local problem */
      if (problem->blocksIsComp[rowNumber])
      {
      infoLocal = lcp_driver_DenseMatrix(local_problem, &z[pos], &w[pos], &internalSolvers[localSolverNum]);
      }
      else /* Solve a linear system*/
      {
        if (local_problem->size ==1)
        {
          double d =  local_problem->M->matrix0[0];
          if (fabs(d) < DBL_EPSILON)
          {
            printf("Numerics::mlcp_pgs_sbm, error: vanishing diagonal term \n");
            printf(" The problem cannot be solved with this method \n");
            exit(EXIT_FAILURE);
          }
          w[pos] = 0.0;
          z[pos] = local_problem->q[0]/d;
        }
        else
        {
          printf("Numerics::mlcp_pgs_sbm, error: nontrivial diagonal term no yet implemented\n");
          exit(EXIT_FAILURE);
        }
      }
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
        fprintf(stderr, "MCLP_PGS_SBM error: Warning local LCP solver  at global iteration %d.\n for block-row number %d. Output info equal to %d.\n", iter, rowNumber, infoLocal);
        //exit(EXIT_FAILURE);

        break;

      }

      while (localSolverNum < options->numberOfInternalSolvers - 1)
        localSolverNum++;
    }

    /*       cblas_dcopy(problem->size , problem->q , 1 , w , 1); */
    /*       prod(problem->size,problem->size, 1.0, problem->M,z,1.0,w); */
    /*       cblas_daxpy(problem->size, -1.0, w,1,wBackup, 1); */
    /*       num = cblas_dnrm2(problem->size,wBackup,1); */
    /*       error = num*den; */
    /* Criterium convergence */
    hasNotConverged = mlcp_compute_error(problem, z, w, tolerance, &error);
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



int mixedLinearComplementarity_pgs_SBM_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the MLCP_PGS_SBM Solver\n");
  }

  options->solverId = SICONOS_MLCP_PGS_SBM;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
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
