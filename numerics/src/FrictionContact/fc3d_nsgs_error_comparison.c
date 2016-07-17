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
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"
#include "fc3d_Path.h"
#include "fc3d_NCPGlockerFixedPoint.h"
#include "fc3d_projection.h"
#include "fc3d_unitary_enumerative.h"
#include "fc3d_compute_error.h"
#include "NCP_Solvers.h"
#include "SiconosBlas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <alloca.h>
#include "op3x3.h"
#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void snPrintf(int level, SolverOptions* opts, const char *fmt, ...);

void fc3d_nsgs_error_comparison(FrictionContactProblem* problem, double *reaction,
                               double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  unsigned int nc = problem->numberOfContacts;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];

  if (*info == 0)
    return;

  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("fc3d_nsgs_redblack_openmp", "The NSGS method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >= 1");
  }
  assert(options->internalSolvers);

  SolverOptions * localsolver_options = options->internalSolvers;

  SolverPtr local_solver = NULL;
  UpdatePtr update_localproblem = NULL;
  FreeSolverNSGSPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;

  if (verbose > 0) printf("----------------------------------- number of contacts %i\n", nc );


  FrictionContactProblem *localproblem = malloc(sizeof(FrictionContactProblem));
  localproblem->numberOfContacts = 1;
  localproblem->dimension = 3;
  localproblem->q = (double*)malloc(3 * sizeof(double));
  localproblem->mu = (double*)malloc(sizeof(double));

  if (problem->M->storageType == NM_DENSE || problem->M->storageType == NM_SPARSE)
  {
    localproblem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3,
                                                   malloc(9 * sizeof(double)));
  }
  else /* NM_SPARSE_BLOCK */
  {
    localproblem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3, NULL);
  }


  fc3d_nsgs_initialize_local_solver(&local_solver, &update_localproblem,
                                    (FreeSolverNSGSPtr *)&freeSolver, &computeError,
                                    problem, localproblem,
                                    options, localsolver_options);


  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  double error_delta_reaction=0.0;
  double error_delta_velocity=0.0;
  double error_delta_velocity_accurate =0.0;
  unsigned int *scontacts = NULL;
  double normq = cblas_dnrm2(nc*3 , problem->q , 1);



  double * velocity_accurate = (double*) malloc(3*nc*sizeof(double));
  double * velocity_accurate_k = (double*) malloc(3*nc*sizeof(double));

  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    error_delta_reaction=0.0;
    error_delta_velocity=0.0;
    error_delta_velocity_accurate =0.0;

    /* Loop through the contact points */
    //cblas_dcopy( n , q , incx , velocity , incy );

    for (unsigned int kk=0; kk < 3*nc; kk++ ) velocity_accurate_k[kk]=velocity_accurate[kk];

    for ( unsigned int contact = 0 ; contact < nc ; contact+=1)
    {


        if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
        (*update_localproblem)(contact, problem, localproblem,
                               reaction, localsolver_options);

        localsolver_options->iparam[4] = contact;

        /* version without localreaction */
        double localreaction[3], localvelocity[3];
        /* double  worktmp[3]; */
        {
          localreaction[0] = reaction[3 * contact+0];
          localreaction[1] = reaction[3 * contact+1];
          localreaction[2] = reaction[3 * contact+2];

          localvelocity[0] = localproblem->q[0];
          localvelocity[1] = localproblem->q[1];
          localvelocity[2] = localproblem->q[2];

        };

        (*local_solver)(localproblem, localreaction,
                        localsolver_options);
        {

         error_delta_reaction += pow(reaction[3 * contact] - localreaction[0], 2) +
            pow(reaction[3 * contact + 1] - localreaction[1], 2) +
            pow(reaction[3 * contact + 2] - localreaction[2], 2);



          mvp3x3(localproblem->M->matrix0,localreaction,localvelocity);

          error_delta_velocity += pow(velocity[3 * contact] - localvelocity[0], 2) +
            pow(velocity[3 * contact + 1] - localvelocity[1], 2) +
            pow(velocity[3 * contact + 2] - localvelocity[2], 2);

          /* double localerror =0.0; */
          /* fc3d_unitary_compute_and_add_error(localreaction, localvelocity, */
          /*                                    localproblem->mu[0], */
          /*                                    &localerror, worktmp); */
          /* error_delta_velocity += localerror; */


          reaction[3 * contact+0] = localreaction[0];
          reaction[3 * contact+1] = localreaction[1];
          reaction[3 * contact+2] = localreaction[2];
          velocity[3 * contact+0] = localvelocity[0];
          velocity[3 * contact+1] = localvelocity[1];
          velocity[3 * contact+2] = localvelocity[2];
        }

        /* version without localreaction */
        /* (*local_solver)(localproblems[tid], &(reaction[3 * contact]),
           localsolvoptions[tid]);*/


        //printf("### tid = %i\n",tid);
      }



    /* /\**** Criterium convergence ****\/ */
    if (iparam[14] > 0)
    {
      if (iter % iparam[14] == 0)
      {
        (*computeError)(problem, reaction , velocity_accurate, tolerance, options, normq, &error);
      }
    }
    else
    {
      (*computeError)(problem, reaction , velocity_accurate, tolerance, options, normq,  &error);
    }
    for ( unsigned int contact = 0 ; contact < nc ; contact+=1)
    {
      error_delta_velocity_accurate += pow(velocity_accurate[3 * contact] - velocity_accurate_k[3*contact], 2) +
        pow(velocity_accurate[3 * contact + 1] - velocity_accurate_k[3*contact+1], 2) +
        pow(velocity_accurate[3 * contact + 2] - velocity_accurate_k[3*contact+2], 2);
    }

    error_delta_velocity = sqrt(error_delta_velocity);

    error_delta_velocity_accurate = sqrt(error_delta_velocity_accurate);

    error_delta_reaction = sqrt(error_delta_reaction);



    //error = error_delta_reaction;

    if (error < tolerance)
    {
      hasNotConverged = 0;
      if (verbose > 0)
      {
        printf("----------------------------------- FC3D - NSGS - Iteration %i Residual = %14.7e < %7.3e\n", iter, error, options->dparam[0]);

      }
    }
    else
    {
      if (verbose > 0)
      {
        printf("----------------------------------- FC3D - NSGS - Iteration %i Residual = %14.7e > %7.3e\n", iter, error, options->dparam[0]);
        printf(" error  = %e\n", error);

        printf("norm v = %e\n",cblas_dnrm2(nc*3 , velocity , 1));
        printf("norm v_accurate = %e\n",cblas_dnrm2(nc*3 , velocity_accurate , 1));
        printf("norm r = %e\n",cblas_dnrm2(nc*3 , reaction , 1));
        printf("norm q = %e\n", normq);
        printf("error_delta_velocity  = %e\t", error_delta_velocity);
        printf("rel error_delta_velocity  = %e\n", error_delta_velocity/cblas_dnrm2(nc*3 , velocity , 1));
        printf("error_delta_velocity_accurate  = %e\t", error_delta_velocity_accurate);
        printf("rel error_delta_velocity_accurate  = %e\n", error_delta_velocity_accurate/cblas_dnrm2(nc*3 , velocity_accurate , 1));
        printf("error_delta_reaction  = %e\t", error_delta_reaction);

        printf("rel error_delta_reaction  = %e\n", error_delta_reaction/cblas_dnrm2(nc*3 , reaction , 1));
      }
      /* test of consistency */
      double c  = 10.0;
      if   ( (error/c >=  error_delta_reaction/cblas_dnrm2(nc*3 , reaction , 1)) ||
             (error_delta_reaction/cblas_dnrm2(nc*3 , reaction , 1) >= c *error))
      {
        printf("%e %e %e   \n",error/c, error_delta_reaction/cblas_dnrm2(nc*3 , reaction , 1),c *error     );
        printf(" WARNING: rel error_delta_reaction is really different from error  \n");
      }
      
      if   ( (error_delta_velocity/c >=  error_delta_velocity_accurate/cblas_dnrm2(nc*3 , velocity_accurate , 1)) ||
             (error_delta_velocity/cblas_dnrm2(nc*3 , reaction , 1) >= c *error_delta_velocity))
      {
        printf("%e %e %e   \n",error_delta_velocity/c,
               error_delta_velocity_accurate/cblas_dnrm2(nc*3 , velocity_accurate , 1),
               c * error_delta_velocity    );
        printf(" WARNING: rel error_delta_velocity_accurate is really different from error_delta_velocity  \n");
      }


    }

    *info = hasNotConverged;

    if (options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, 3 * nc,
                                               reaction, velocity,
                                               error, NULL);
    }
  }

  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;

  /***** Free memory *****/
  (*freeSolver)(problem,localproblem,localsolver_options);
  if (problem->M->storageType == NM_DENSE && localproblem->M->matrix0)
  {
    free(localproblem->M->matrix0);
  }
  localproblem->M->matrix0 = NULL;
  freeFrictionContactProblem(localproblem);
  solver_options_delete(localsolver_options);
  free(localsolver_options);


  if (scontacts) /* shuffle */
  {
    free(scontacts);
  }

}
