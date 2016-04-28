/* Siconos-Numerics, Copyright INRIA 2005-2016
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
#include "fc3d_projection.h"
//#include "gfc3d_projection.h"
#include "gfc3d_Solvers.h"
#include "gfc3d_compute_error.h"
#include "projectionOnCone.h"
#include "SiconosLapack.h"
#include "SparseBlockMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "sanitizer.h"

static void Globalfc3d_projection_free(GlobalFrictionContactProblem* problem)
{
  assert(problem->M);
}

static void initializeGlobalLocalSolver(int n, SolverGlobalPtr* solve, FreeSolverGlobalPtr* freeSolver, ComputeErrorGlobalPtr* computeError, const NumericsMatrix* const M, const double* const q, const double* const mu, int* iparam)
{
  /** Connect to local solver */
  /* Projection */
  if (iparam[4] == 0)
  {
    /*       *solve = &fc3d_projectionOnCone_solve; */
    *freeSolver = &Globalfc3d_projection_free;
    *computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error;
    /*       fc3d_projection_initialize(n,M,q,mu); */
  }
  else
  {
    fprintf(stderr, "Numerics, fc3d_nsgs failed. Unknown local solver set by iparam[4]\n");
    exit(EXIT_FAILURE);
  }
}


void gfc3d_nsgs(GlobalFrictionContactProblem* restrict problem, double* restrict reaction, double* restrict velocity, double* restrict globalVelocity, int* restrict info, SolverOptions* restrict options)
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
  unsigned int erritermax = options->iparam[7];
  if (erritermax == 0)
  {
    fprintf(stderr, "gfc3d_nsgs :: erritermax is 0, something is wrong\n");
  }
  /* Tolerance */
  double tolerance = dparam[0];

  /* Check for trivial case */
  *info = checkTrivialCaseGlobal(n, q, velocity, reaction, globalVelocity, options);

  if (*info == 0)
    return;

  gfc3d_init_workspace(problem);

  NumericsMatrix* factorized_M = problem->workspace->factorized_M;
  double* qtmp = problem->workspace->globalVelocity;

  SolverGlobalPtr local_solver = NULL;
  FreeSolverGlobalPtr freeSolver = NULL;
  ComputeErrorGlobalPtr computeError = NULL;

  /* Connect local solver */
  initializeGlobalLocalSolver(n, &local_solver, &freeSolver, &computeError, M, q, mu, iparam);

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  int contact; /* Number of the current row of blocks in M */

  if (H->storageType != M->storageType)
  {
    //     if(verbose==1)
    fprintf(stderr, "Numerics, gfc3d_nsgs. H->storageType != M->storageType :This case is not taken into account.\n");
    exit(EXIT_FAILURE);
  }

  dparam[0] = dparam[2]; // set the tolerance for the local solver

  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    /* Solve the first part with the current reaction */

    /* qtmp <--q */
    cblas_dcopy_msan(n, q, 1, qtmp, 1);

    /*qtmp = H reaction +qtmp */
    NM_gemv(1., H, reaction, 1., qtmp);

    cblas_dcopy(n, qtmp, 1, globalVelocity, 1);

    CHECK_RETURN(!NM_gesv_expert(factorized_M, globalVelocity, true));

    /* Compute current local velocity */
    /*      velocity <--b */
    cblas_dcopy(m, b, 1, velocity, 1);

    /* velocity <-- H^T globalVelocity + velocity*/
    NM_tgemv(1., H, globalVelocity, 1., velocity);

    /* Loop through the contact points */

    for (contact = 0 ; contact < nc ; ++contact)
    {
      /*    (*local_solver)(contact,n,reaction,iparam,dparam); */
      int pos = contact * 3;
      double normUT = sqrt(velocity[pos + 1] * velocity[pos + 1] + velocity[pos + 2] * velocity[pos + 2]);
      double an = 1.0;
      reaction[pos] -= an * (velocity[pos] + mu[contact] * normUT);
//      reaction[pos] -= an * mu[contact] * normUT;
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
    /* this is very expensive to check, you better do it only once in a while  */
    if (!(iter % erritermax))
    {
      (*computeError)(problem, reaction , velocity, globalVelocity, tolerance, &error);

      if (verbose > 0)
        printf("----------------------------------- FC3D - NSGS - Iteration %i Residual = %14.7e; Tol = %g\n", iter, error, tolerance);
    }

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }

  /*  One last error computation in case where are at the very end */
  if (iter == itermax)
  {
    (*computeError)(problem, reaction , velocity, globalVelocity, tolerance, &error);
  }

  dparam[0] = tolerance;
  dparam[1] = error;


  /***** Free memory *****/
  (*freeSolver)(problem);
}

