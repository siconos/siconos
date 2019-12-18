/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <assert.h>                        // for assert
#include <math.h>                          // for sqrt
#include <stdio.h>                         // for fprintf, NULL, stderr
#include <stdlib.h>                        // for exit, EXIT_FAILURE
#include "Friction_cst.h"                  // for SICONOS_FRICTION_3D_IPARAM...
#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_dnrm2
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NumericsFwd.h"                   // for GlobalFrictionContactProblem
#include "NumericsMatrix.h"                // for NumericsMatrix, NM_gemv
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_DPA...
#include "debug.h"                         // for DEBUG_EXPR, DEBUG_PRINTF
#include "gfc3d_Solvers.h"                 // for ComputeErrorGlobalPtr, gfc...
#include "gfc3d_compute_error.h"           // for gfc3d_compute_error
#include "numerics_verbose.h"              // for numerics_printf_verbose
#include "projectionOnCone.h"              // for projectionOnCone
#include "sanitizer.h"                     // for cblas_dcopy_msan

static void gfc3d_nsgs_local_solver_projection_free(GlobalFrictionContactProblem* problem)
{
  assert(problem->M);
}

static void gfc3d_nsgs_initialize_local_solver(int n, SolverGlobalPtr* solve, FreeSolverGlobalPtr* freeSolver, ComputeErrorGlobalPtr* computeError, const NumericsMatrix* const M, const double* const q, const double* const mu, int* iparam)
{
  /** Connect to local solver */
  /* Projection */
  if(iparam[4] == 0)
  {
    /*       *solve = &fc3d_projectionOnCone_solve; */
    *freeSolver = &gfc3d_nsgs_local_solver_projection_free;
    *computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error;
    /*       fc3d_projection_initialize(n,M,q,mu); */
  }
  else
  {
    fprintf(stderr, "Numerics, gfc3d_nsgs failed. Unknown local solver set by iparam[4]\n");
    exit(EXIT_FAILURE);
  }
}


void gfc3d_nsgs(GlobalFrictionContactProblem* restrict problem, double* restrict reaction,
                double* restrict velocity, double* restrict globalVelocity,
                int* restrict info, SolverOptions* restrict options)
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

  assert((int)H->size1 == problem->numberOfContacts * problem->dimension);
  assert((int)M->size0 == M->size1);
  assert((int)M->size0 == H->size0); /* size(velocity) ==
                                      * Htrans*globalVelocity */



  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];

  /* Check for trivial case */
  *info = gfc3d_checkTrivialCaseGlobal(n, q, velocity, reaction, globalVelocity, options);

  if(*info == 0)
    return;

  SolverGlobalPtr local_solver = NULL;
  FreeSolverGlobalPtr freeSolver = NULL;
  ComputeErrorGlobalPtr computeError = NULL;

  /* Connect local solver */
  gfc3d_nsgs_initialize_local_solver(n, &local_solver, &freeSolver, &computeError, M, q, mu, iparam);

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  int contact; /* Number of the current row of blocks in M */

  if(H->storageType != M->storageType)
  {
    //     if(verbose==1)
    fprintf(stderr, "Numerics, gfc3d_nsgs. H->storageType != M->storageType :This case is not taken into account.\n");
    exit(EXIT_FAILURE);
  }

  double norm_q = cblas_dnrm2(n, problem->q, 1);
  double norm_b = cblas_dnrm2(m, problem->b, 1);
  /* verbose=1; */
  while((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    /* Solve the first part with the current reaction */
    DEBUG_PRINTF("--------- iter = %i\n", iter);
    /* globalVelocity <--q */
    cblas_dcopy_msan(n, q, 1, globalVelocity, 1);
    /* globalVelocity = H reaction + globalVelocity */
    if(nc > 0) NM_gemv(1., H, reaction, 1., globalVelocity);

    CHECK_RETURN(!NM_gesv_expert(problem->M, globalVelocity, NM_KEEP_FACTORS));

    DEBUG_EXPR(NM_vector_display(reaction,m));
    DEBUG_EXPR(NM_vector_display(globalVelocity,n));

    if(nc > 0)
    {
      /* Compute current local velocity */
      /*      velocity <--b */
      cblas_dcopy(m, b, 1, velocity, 1);

      /* velocity <-- H^T globalVelocity + velocity*/
      NM_tgemv(1., H, globalVelocity, 1., velocity);
      DEBUG_EXPR(NM_vector_display(velocity,m););

      /* Loop through the contact points */
      for(contact = 0 ; contact < nc ; ++contact)
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
      DEBUG_EXPR(NM_vector_display(reaction,m););
    }


    /* **** Criterium convergence **** */
    /* this is very expensive to check, you better do it only once in a while  */
    if(options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY]>0)
    {
      if(!(iter % options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY]))
      {
        /* computeGlobalVelocity(problem, reaction, globalVelocity); */
        (computeError)(problem, reaction, velocity, globalVelocity, tolerance, options, norm_q, norm_b, &error);
      }
    }
    else
    {
      (computeError)(problem, reaction, velocity, globalVelocity, tolerance, options, norm_q, norm_b, &error);

      numerics_printf_verbose(1,"----- GFC3D - NSGS - Iteration %i Residual = %14.7e; Tol = %g", iter, error, tolerance);
    }
    if(error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }

  /*  One last error computation in case where are at the very end */
  if(iter == itermax)
  {
    (*computeError)(problem, reaction, velocity, globalVelocity, tolerance, options, norm_q, norm_b, &error);
  }

  dparam[SICONOS_DPARAM_TOL] = tolerance;
  dparam[SICONOS_DPARAM_RESIDU] = error;


  /***** Free memory *****/
  (*freeSolver)(problem);
}
