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

#include "FrictionContact3D_projection.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContact3D_compute_error.h"
#include "LA.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "pinv.h"

void initializeLocalSolver_nsgs_velocity(int n, SolverPtr* solve, FreeSolverPtr* freeSolver, ComputeErrorPtr* computeError, const NumericsMatrix* const M, const double* const q, const double* const mu, Solver_Options* localsolver_options)
{
  /** Connect to local solver */
  /* Projection */
  if (strcmp(localsolver_options->solverName, "ProjectionOnCone_velocity") == 0)
  {
    *solve = &frictionContact3D_projectionOnCone_velocity_solve;
    *freeSolver = &frictionContact3D_projection_free;
    *computeError = &FrictionContact3D_compute_error_velocity;
    frictionContact3D_projection_initialize(n, M, q, mu);
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs_velocity failed. Unknown internal solver : %s.\n", localsolver_options->solverName);
    exit(EXIT_FAILURE);
  }
}

void frictionContact3D_nsgs_velocity(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, Solver_Options* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* Number of contacts */
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  NumericsMatrix* M = problem->M;
  double* mu = problem->mu;
  /* Dimension of the problem */
  int n = 3 * nc;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];

  /* Check for trivial case */
  /*   *info = checkTrivialCase(n, q,velocity, reaction, iparam, dparam); */

  if (*info == 0)
    return;

  SolverPtr local_solver = NULL;
  FreeSolverPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;





  if (M->storageType == 0)
  {
    /*  /\* Inversion of the matrix M *\/ */
    /*   int* ipiv = (int *)malloc(n*sizeof(*ipiv));  */
    /*   int infoDGETRF=0; */
    /*   DGETRF(n,n,M->matrix0,n, ipiv,infoDGETRF ); */
    /*   assert(!infoDGETRF); */
    /*   int infoDGETRI; */
    /*   DGETRI(n,M->matrix0,n, ipiv,infoDGETRI ); */
    /*   assert(!infoDGETRI); */
    /*   double* qtmp = (double*)malloc(n*sizeof(double)); */
    /*   DCOPY(n,  q, 1, qtmp, 1); */
    /*   DGEMV(LA_NOTRANS, n, n, -1.0, M->matrix0 , n, qtmp,1,0.0,q, 1); */
    /*   free(ipiv); */
    /*   free(qtmp); */
    double tolpinv = 1e-07;
    pinv(M->matrix0, n, n, tolpinv);
  }
  else
  {
    fprintf(stderr, "Numerics, frictionContact3D_nsgs_velocity. Not yet implemented for storageType %i\n", M->storageType);
    exit(EXIT_FAILURE);
  }
  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("frictionContact3D_nsgs_velocity", "The NSGS method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }
  assert(&options[1]);

  Solver_Options * localsolver_options = &options[1];
  if (verbose > 0)
  {
    printf("Local solver data :");
    printSolverOptions(localsolver_options);
  }


  /* Connect local solver */
  initializeLocalSolver_nsgs_velocity(n, &local_solver, &freeSolver, &computeError, M, q, mu, localsolver_options);

  /*****  NSGS_VELOCITY Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */

  dparam[0] = dparam[2]; // set the tolerance for the local solver
  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    /* Loop through the contact points */
    //DCOPY( n , q , incx , velocity , incy );
    for (contact = 0 ; contact < nc ; ++contact)
    {
      (*local_solver)(contact, n, velocity, localsolver_options);
      for (int ncc = 0; ncc < 3; ncc ++)
      {
        printf("velocity[%i]=%14.7e\t", ncc, velocity[contact * 3 + ncc]);
      }
      printf("\n");
    }



    /* **** Criterium convergence **** */
    (*computeError)(problem, reaction , velocity, tolerance, &error);

    if (verbose > 0)
      printf("----------------------------------- FC3D - NSGS_VELOCITY - Iteration %i Error = %14.7e\n", iter, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }
  printf("----------------------------------- FC3D - NSGS_VELOCITY - # Iteration %i Final Error = %14.7e\n", iter, error);
  dparam[0] = tolerance;
  dparam[1] = error;


  /***** Free memory *****/
  (*freeSolver)();
}

