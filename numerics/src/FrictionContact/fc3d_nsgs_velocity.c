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

#include "fc3d_projection.h"
#include "fc3d_Solvers.h"
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"
#include "fc3d_compute_error.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "pinv.h"
#include "Friction_cst.h"

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void initializeLocalSolver_nsgs_velocity(SolverPtr* solve, FreeSolverPtr* freeSolver, ComputeErrorPtr* computeError, FrictionContactProblem* problem, FrictionContactProblem* localproblem, SolverOptions* localsolver_options)
{




  /** Connect to local solver */
  /* Projection */
  if (localsolver_options->solverId == SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity)
  {
    *solve = &fc3d_projectionOnCone_velocity_solve;
    *freeSolver = (FreeSolverPtr)&fc3d_projection_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error_velocity;
    fc3d_projection_initialize(problem, localproblem);
  }
  else
  {
    fprintf(stderr, "Numerics, fc3d_nsgs_velocity failed. Unknown internal solver : %s.\n", idToName(localsolver_options->solverId));
    exit(EXIT_FAILURE);
  }
}

void fc3d_nsgs_velocity(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* Number of contacts */
  int nc = problem->numberOfContacts;
  NumericsMatrix* M = problem->M;
  /* Dimension of the problem */
  int n = 3 * nc;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];

  /* Check for trivial case */
  /*   *info = checkTrivialCase(n, q,velocity, reaction, options); */

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
    /*   cblas_dcopy(n,  q, 1, qtmp, 1); */
    /*   cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, -1.0, M->matrix0 , n, qtmp,1,0.0,q, 1); */
    /*   free(ipiv); */
    /*   free(qtmp); */
    double tolpinv = 1e-07;
    pinv(M->matrix0, n, n, tolpinv);
  }
  else
  {
    fprintf(stderr, "Numerics, fc3d_nsgs_velocity. Not yet implemented for storageType %i\n", M->storageType);
    exit(EXIT_FAILURE);
  }
  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("fc3d_nsgs_velocity", "The NSGS method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }


  SolverOptions * localsolver_options = options->internalSolvers;

  /* Connect local solver */


  FrictionContactProblem* localproblem = 0;
  initializeLocalSolver_nsgs_velocity(&local_solver, &freeSolver, &computeError, problem, localproblem, localsolver_options);

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
    //cblas_dcopy( n , q , incx , velocity , incy );
    for (contact = 0 ; contact < nc ; ++contact)
    {

      (*local_solver)(localproblem, velocity, localsolver_options);
      for (int ncc = 0; ncc < 3; ncc ++)
      {
        printf("velocity[%i]=%14.7e\t", ncc, velocity[contact * 3 + ncc]);
      }
      printf("\n");
    }



    /* **** Criterium convergence **** */
    (*computeError)(problem, reaction , velocity, tolerance, options, &error);

    if (verbose > 0)
      printf("----------------------------------- FC3D - NSGS_VELOCITY - Iteration %i Residual = %14.7e\n", iter, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }
  printf("----------------------------------- FC3D - NSGS_VELOCITY - # Iteration %i Final Residual = %14.7e\n", iter, error);
  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;

  /***** Free memory *****/
  (*freeSolver)();
}
int fc3d_nsgs_velocity_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGSV Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_NSGSV;

  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 8;
  options->dSize = 8;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  null_SolverOptions(options);
  for (i = 0; i < 8; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_onecontact_nonsmooth_Newtow_setDefaultSolverOptions(options->internalSolvers);

  return 0;
}
