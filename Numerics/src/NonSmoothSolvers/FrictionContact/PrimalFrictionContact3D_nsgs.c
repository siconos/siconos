/* Siconos-Numerics, Copyright INRIA 2005-2012.
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
#include "FrictionContact3D_projection.h"
//#include "PrimalFrictionContact3D_projection.h"
#include "PrimalFrictionContact3D_Solvers.h"
#include "PrimalFrictionContact3D_compute_error.h"
#include "projectionOnCone.h"
#include "LA.h"
#include "SparseBlockMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

extern int *Primal_ipiv;
extern int  Primal_MisInverse;
extern int  Primal_MisLU;

void PrimalfrictionContact3D_projection_free(PrimalFrictionContactProblem* problem)
{
  assert(problem->M);
  if (problem->M->storageType == 0)
  {
    free(Primal_ipiv);
  }
}


void initializePrimalLocalSolver(int n, SolverPrimalPtr* solve, FreeSolverPrimalPtr* freeSolver, ComputeErrorPrimalPtr* computeError, const NumericsMatrix* const M, const double* const q, const double* const mu, int* iparam)
{
  /** Connect to local solver */
  /* Projection */
  if (iparam[4] == 0)
  {
    /*       *solve = &frictionContact3D_projectionOnCone_solve; */
    *freeSolver = &PrimalfrictionContact3D_projection_free;
    *computeError = (ComputeErrorPrimalPtr)&PrimalFrictionContact3D_compute_error;
    /*       frictionContact3D_projection_initialize(n,M,q,mu); */
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown local solver set by iparam[4]\n");
    exit(EXIT_FAILURE);
  }
}


void primalFrictionContact3D_nsgs(PrimalFrictionContactProblem* problem, double *reaction, double *velocity, double *globalVelocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;
  NumericsMatrix* M = problem->M;
  NumericsMatrix* H = problem->H;
  double* q = problem->q;
  double* b = problem->b;
  double* mu = problem->mu;

  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];

  /* Check for trivial case */
  *info = checkTrivialCasePrimal(n, q, velocity, reaction, globalVelocity, options);

  if (*info == 0)
    return;

  SolverPrimalPtr local_solver = NULL;
  FreeSolverPrimalPtr freeSolver = NULL;
  ComputeErrorPrimalPtr computeError = NULL;

  /* Connect local solver */
  initializePrimalLocalSolver(n, &local_solver, &freeSolver, &computeError, M, q, mu, iparam);

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  int contact; /* Number of the current row of blocks in M */
  SparseBlockStructuredMatrix *Htrans = (SparseBlockStructuredMatrix*)malloc(sizeof(SparseBlockStructuredMatrix));

  if (H->storageType != M->storageType)
  {
    //     if(verbose==1)
    fprintf(stderr, "Numerics, PrimalFrictionContact3D_nsgs. H->storageType != M->storageType :This case is not taken into account.\n");
    exit(EXIT_FAILURE);
  }
  else if (M->storageType == 1)
  {
    inverseDiagSBM(M->matrix1);
    Primal_MisInverse = 1;
    transposeSBM(H->matrix1, Htrans);
  }
  else if (M->storageType == 0)
  {
    /*  Assume that M is not already LU */
    int infoDGETRF = -1;
    Primal_ipiv = (int *)malloc(n * sizeof(int));
    assert(!Primal_MisLU);
    DGETRF(n, n, M->matrix0, n, Primal_ipiv, &infoDGETRF);
    Primal_MisLU = 1;
    assert(!infoDGETRF);
  }

  dparam[0] = dparam[2]; // set the tolerance for the local solver
  double* qtmp = (double*)malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) qtmp[i] = 0.0;



  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    /* Solve the first part with the current reaction */

    /* qtmp <--q */
    DCOPY(n, q, 1, qtmp, 1);

    double alpha = 1.0;
    double beta = 1.0;
    /*qtmp = H reaction +qtmp */
    prodNumericsMatrix(m, n, alpha, H, reaction , beta, qtmp);

    if (M->storageType == 1)
    {
      beta = 0.0;
      assert(Primal_MisInverse);
      /*  globalVelocity = M^-1 qtmp */
      prodNumericsMatrix(n, n, alpha, M, qtmp , beta, globalVelocity);
    }
    else if (M->storageType == 0)
    {
      int infoDGETRS = -1;
      DCOPY(n, qtmp, 1, globalVelocity, 1);
      assert(Primal_MisLU);
#ifdef USE_MKL
      DGETRS(CLA_NOTRANS, n, 1,  M->matrix0, n, Primal_ipiv, globalVelocity , n, infoDGETRS);
#else
      DGETRS(LA_NOTRANS, n, 1,  M->matrix0, n, Primal_ipiv, globalVelocity , n, &infoDGETRS);
#endif
      assert(!infoDGETRS);
    }
    /* Compute current local velocity */
    /*      velocity <--b */
    DCOPY(m, b, 1, velocity, 1);

    if (H->storageType == 1)
    {
      /* velocity <-- H^T globalVelocity + velocity*/
      beta = 1.0;
      prodSBM(n, m, alpha, Htrans, globalVelocity , beta, velocity);
    }
    else if (H->storageType == 0)
    {
      DGEMV(LA_TRANS, n, m, 1.0, H->matrix0 , n, globalVelocity , 1, 1.0, velocity, 1);
    }

    /* Loop through the contact points */

    for (contact = 0 ; contact < nc ; ++contact)
    {
      /*    (*local_solver)(contact,n,reaction,iparam,dparam); */
      int pos = contact * 3;
      double normUT = sqrt(velocity[pos + 1] * velocity[pos + 1] + velocity[pos + 2] * velocity[pos + 2]);
      double an = 1.0;
      reaction[pos] -= an * (velocity[pos] + mu[contact] * normUT);
      reaction[pos + 1] -= an * velocity[pos + 1];
      reaction[pos + 2] -= an * velocity[pos + 2];
      projectionOnCone(&reaction[pos], mu[contact]);
    }
    /*       int k; */
    /*       printf("\n"); */
    /*       for (k = 0 ; k < m; k++) printf("velocity[%i] = %12.8e \t \t reaction[%i] = %12.8e \n ", k, velocity[k], k , reaction[k]); */
    /*       for (k = 0 ; k < n; k++) printf("globalVelocity[%i] = %12.8e \t \n ", k, globalVelocity[k]); */
    /*       printf("\n"); */



    /* **** Criterium convergence **** */
    (*computeError)(problem, reaction , velocity, globalVelocity, tolerance, &error);

    if (verbose > 0)
      printf("----------------------------------- FC3D - NSGS - Iteration %i Error = %14.7e\n", iter, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }

  if (H->storageType == 1)
  {
    freeSBM(Htrans);
  }
  free(Htrans);
  free(qtmp);
  /*   free(Primal_ipiv); */
  dparam[0] = tolerance;
  dparam[1] = error;


  /***** Free memory *****/
  (*freeSolver)(problem);
}

