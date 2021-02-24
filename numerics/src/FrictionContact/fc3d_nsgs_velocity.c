/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#include <assert.h>                  // for assert
#include <stdio.h>                   // for printf, fprintf, NULL, stderr
#include <stdlib.h>                  // for exit, EXIT_FAILURE
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "Friction_cst.h"            // for SICONOS_FRICTION_3D_ONECONTACT_NSN
#include "NumericsFwd.h"             // for SolverOptions, FrictionContactPr...
#include "NumericsMatrix.h"          // for NumericsMatrix
#include "SiconosBlas.h"             // for cblas_dnrm2
#include "SolverOptions.h"           // for SolverOptions, SICONOS_DPARAM_TOL
#include "fc3d_Solvers.h"            // for ComputeErrorPtr, FreeSolverPtr
#include "fc3d_compute_error.h"      // for fc3d_compute_error_velocity
#include "fc3d_projection.h"         // for fc3d_projection_initialize, fc3d...
#include "numerics_verbose.h"        // for numerics_error, verbose
#include "pinv.h"                    // for pinv

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void fc3d_nsgs_initialize_local_solver_velocity(SolverPtr* solve, FreeSolverPtr* freeSolver, ComputeErrorPtr* computeError, FrictionContactProblem* problem, FrictionContactProblem* localproblem, SolverOptions* localsolver_options)
{




  /** Connect to local solver */
  /* Projection */
  if(localsolver_options->solverId == SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity)
  {
    *solve = &fc3d_projectionOnCone_velocity_solve;
    *freeSolver = (FreeSolverPtr)&fc3d_projection_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error_velocity;
    fc3d_projection_initialize(problem, localproblem);
  }
  else
  {
    fprintf(stderr, "Numerics, fc3d_nsgs_velocity failed. Unknown internal solver : %s.\n", solver_options_id_to_name(localsolver_options->solverId));
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
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  double norm_q = cblas_dnrm2(nc*3, problem->q, 1);
  /* Check for trivial case */
  /*   *info = fc3d_checkTrivialCase(n, q,velocity, reaction, options); */

  if(*info == 0)
    return;

  SolverPtr local_solver = NULL;
  FreeSolverPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;





  if(M->storageType == 0)
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
  if(options->numberOfInternalSolvers < 1)
  {
    numerics_error("fc3d_nsgs_velocity", "The NSGS method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }


  SolverOptions * localsolver_options = options->internalSolvers[0];

  /* Connect local solver */


  FrictionContactProblem* localproblem = 0;
  fc3d_nsgs_initialize_local_solver_velocity(&local_solver, &freeSolver, &computeError, problem, localproblem, localsolver_options);

  /*****  NSGS_VELOCITY Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */



  dparam[SICONOS_DPARAM_TOL] = dparam[2]; // set the tolerance for the local solver
  while((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    /* Loop through the contact points */
    //cblas_dcopy( n , q , incx , velocity , incy );
    for(contact = 0 ; contact < nc ; ++contact)
    {

      (*local_solver)(localproblem, velocity, localsolver_options);
      for(int ncc = 0; ncc < 3; ncc ++)
      {
        printf("velocity[%i]=%14.7e\t", ncc, velocity[contact * 3 + ncc]);
      }
      printf("\n");
    }



    /* **** Criterium convergence **** */
    (*computeError)(problem, reaction, velocity, tolerance, options, norm_q,  &error);

    if(verbose > 0)
      printf("--------------- FC3D - NSGS_VELOCITY - Iteration %i Residual = %14.7e\n", iter, error);

    if(error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }
  printf("--------------- FC3D - NSGS_VELOCITY - # Iteration %i Final Residual = %14.7e\n", iter, error);
  dparam[SICONOS_DPARAM_TOL] = tolerance;
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[7] = iter;

  /***** Free memory *****/
  (*freeSolver)();
}

void fc3d_nsgs_velocity_set_default(SolverOptions* options)
{
  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] = solver_options_create(SICONOS_FRICTION_3D_ONECONTACT_NSN);
}
