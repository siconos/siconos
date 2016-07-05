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

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

#include <SiconosConfig.h>
#if defined(WITH_OPENMP) && defined(_OPENMP)
#define USE_OPENMP 1
#include <omp.h>
#endif

void initializeLocalSolver_nsgs(SolverPtr* solve, UpdatePtr* update, FreeSolverNSGSPtr* freeSolver, ComputeErrorPtr* computeError,
                                FrictionContactProblem* problem, FrictionContactProblem* localproblem,
                                SolverOptions * options, SolverOptions * localsolver_options);
void fc3d_nsgs_redblack_openmp(FrictionContactProblem* problem, double *reaction,
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

  /* Connect local solver and local problem*/
#if defined(USE_OPENMP) && defined(_OPENMP)
  const unsigned int max_threads = omp_get_max_threads();
  FrictionContactProblem **localproblems = alloca(max_threads*sizeof(void*));
  SolverOptions **localsolvoptions = alloca(max_threads*sizeof(void*));
#else
  const unsigned int max_threads = 1;
  FrictionContactProblem *localproblems[1];
  SolverOptions *localsolvoptions[1];
#endif

  for (unsigned int i=0; i < max_threads; i++)
  {
    FrictionContactProblem *localproblem = malloc(sizeof(FrictionContactProblem));
    localproblems[i] = localproblem;
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

    localsolvoptions[i] = malloc(sizeof(SolverOptions));
    null_SolverOptions(localsolvoptions[i]);
    localsolvoptions[i]->dparam = NULL;
    localsolvoptions[i]->iparam = NULL;
    copy_SolverOptions(localsolver_options,localsolvoptions[i]);

    initializeLocalSolver_nsgs(&local_solver, &update_localproblem,
                               (FreeSolverNSGSPtr *)&freeSolver, &computeError,
                               problem, localproblems[i],
                               options, localsolvoptions[i]);
  }

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  unsigned int *scontacts = NULL;

  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    /* Loop through the contact points */
    //cblas_dcopy( n , q , incx , velocity , incy );
    for (int c=0; c<2; c++)
    {
    #if defined(_OPENMP) && defined(USE_OPENMP)
    #pragma omp parallel for private(contact)
    #endif
      for ( unsigned int contact = c ; contact < nc ; contact+=2)
      {
      #if defined(_OPENMP) && defined(USE_OPENMP)
        unsigned int tid = omp_get_thread_num();
      #else
        unsigned int tid = 0;
      #endif

        if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
        (*update_localproblem)(contact, problem, localproblems[tid],
                               reaction, localsolvoptions[tid]);
        localsolvoptions[tid]->iparam[4] = contact;
        double localreaction[3];
        #pragma omp critical(localreaction)
        {
          localreaction[0] = reaction[3 * contact+0];
          localreaction[1] = reaction[3 * contact+1];
          localreaction[2] = reaction[3 * contact+2];
        };
        (*local_solver)(localproblems[tid], localreaction,
                        localsolvoptions[tid]);
        #pragma omp critical(localreaction)
        {
          reaction[3 * contact+0] = localreaction[0];
          reaction[3 * contact+1] = localreaction[1];
          reaction[3 * contact+2] = localreaction[2];
        }
        //printf("### tid = %i\n",tid);
      }
    }

    /* **** Criterium convergence **** */
    (*computeError)(problem, reaction , velocity, tolerance, options, &error);

    if (error < tolerance)
    {
      hasNotConverged = 0;
      if (verbose > 0)
        printf("----------------------------------- FC3D - NSGS - Iteration %i Residual = %14.7e < %7.3e\n", iter, error, options->dparam[0]);
    }
    else
    {
      if (verbose > 0)
        printf("----------------------------------- FC3D - NSGS - Iteration %i Residual = %14.7e > %7.3e\n", iter, error, options->dparam[0]);
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
  for (unsigned int i=0; i < max_threads; i++)
  {
    (*freeSolver)(problem,localproblems[i],localsolvoptions[i]);
    if (problem->M->storageType == NM_DENSE && localproblems[i]->M->matrix0)
    {
      free(localproblems[i]->M->matrix0);
    }
    localproblems[i]->M->matrix0 = NULL;
    freeFrictionContactProblem(localproblems[i]);
    deleteSolverOptions(localsolvoptions[i]);
    free(localsolvoptions[i]);
  }

  if (scontacts) /* shuffle */
  {
    free(scontacts);
  }

}
