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
#include "rolling_fc3d_projection.h"
#include "rolling_fc3d_local_problem_tools.h"
#include "rolling_fc3d_compute_error.h"
#include "SiconosBlas.h"
#include "NumericsArrays.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <float.h>
#include <string.h>
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"
#include "numerics_verbose.h"


//#define FCLIB_OUTPUT

#ifdef FCLIB_OUTPUT
static int fccounter = -1;
#include "fclib_interface.h"
#endif



#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void rolling_fc3d_nsgs_update(int contact, RollingFrictionContactProblem* problem, RollingFrictionContactProblem* localproblem, double * reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  rolling_fc3d_local_problem_fill_M(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
         excluding the block corresponding to the current contact. ****/
  rolling_fc3d_local_problem_compute_q(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];


}

void rolling_fc3d_nsgs_initialize_local_solver(RollingSolverPtr* solve, RollingUpdatePtr* update,
                                               RollingFreeSolverNSGSPtr* freeSolver,
                                               RollingComputeErrorPtr* computeError,
                                               RollingFrictionContactProblem* problem,
                                               RollingFrictionContactProblem* localproblem,
                                               SolverOptions * options)
{
  SolverOptions * localsolver_options = options->internalSolvers;
  /** Connect to local solver */
  switch (localsolver_options->solverId)
  {
  case SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration:
  {
    *solve = &rolling_fc3d_projectionOnConeWithLocalIteration_solve;
    *update = &rolling_fc3d_projection_update;
    *freeSolver = (RollingFreeSolverNSGSPtr)&rolling_fc3d_projectionOnConeWithLocalIteration_free;
    *computeError = (RollingComputeErrorPtr)&rolling_fc3d_compute_error;
    rolling_fc3d_projectionOnConeWithLocalIteration_initialize(problem, localproblem,localsolver_options );
    break;
  }
  case SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone:
  {
    *solve = &rolling_fc3d_projectionOnCone_solve;
    *update = &rolling_fc3d_projection_update;
    *freeSolver = (RollingFreeSolverNSGSPtr)&rolling_fc3d_projection_free;
    *computeError = (RollingComputeErrorPtr)&rolling_fc3d_compute_error;
    rolling_fc3d_projection_initialize(problem, localproblem);
    break;
  }
  default:
  {
    numerics_error("rolling_fc3d_nsgs_initialize_local_solver", "Numerics, rolling_fc3d_nsgs failed. Unknown internal solver : %s.\n", solver_options_id_to_name(localsolver_options->solverId));
  }
  }
}


static
unsigned int* allocShuffledContacts(RollingFrictionContactProblem *problem,
                                    SolverOptions *options)
{
  unsigned int *scontacts = 0;
  unsigned int nc = problem->numberOfContacts;
  if (options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE||
      options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)
  {
    if (options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED] >0)
    {
      srand((unsigned int)options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED]);
    }
    else
      srand(1);
    scontacts = (unsigned int *) malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i < nc ; ++i)
    {
      scontacts[i] = i;
    }
    uint_shuffle(scontacts, nc);
  }
  return scontacts;
}

static
int solveLocalReaction(RollingUpdatePtr update_localproblem, RollingSolverPtr local_solver,
                       unsigned int contact, RollingFrictionContactProblem *problem,
                       RollingFrictionContactProblem *localproblem, double *reaction,
                       SolverOptions *localsolver_options, double localreaction[5])
{
  (*update_localproblem)(contact, problem, localproblem,
                         reaction, localsolver_options);

  localsolver_options->iparam[SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_CONTACTNUMBER] = contact;

  localreaction[0] = reaction[contact*5 + 0];
  localreaction[1] = reaction[contact*5 + 1];
  localreaction[2] = reaction[contact*5 + 2];
  localreaction[3] = reaction[contact*5 + 3];
  localreaction[4] = reaction[contact*5 + 4];

  return (*local_solver)(localproblem, localreaction, localsolver_options);
}

static
void performRelaxation(double localreaction[5], double *oldreaction, double omega)
{
  localreaction[0] = omega*localreaction[0]+(1.0-omega)*oldreaction[0];
  localreaction[1] = omega*localreaction[1]+(1.0-omega)*oldreaction[1];
  localreaction[2] = omega*localreaction[2]+(1.0-omega)*oldreaction[2];
  localreaction[3] = omega*localreaction[3]+(1.0-omega)*oldreaction[3];
  localreaction[4] = omega*localreaction[4]+(1.0-omega)*oldreaction[4];
}

static
void accumulateLightErrorSum(double *light_error_sum, double localreaction[5],
                             double *oldreaction)
{
  *light_error_sum += ( pow(oldreaction[0] - localreaction[0], 2) +
                        pow(oldreaction[1] - localreaction[1], 2) +
                        pow(oldreaction[2] - localreaction[2], 2) +
                        pow(oldreaction[3] - localreaction[3], 2) +
                        pow(oldreaction[4] - localreaction[4], 2) );
}

static
void acceptLocalReactionFiltered(RollingFrictionContactProblem *localproblem,
                                 SolverOptions *localsolver_options,
                                 unsigned int contact, unsigned int iter,
                                 double *reaction, double localreaction[5])
{
  if (isnan(localsolver_options->dparam[SICONOS_DPARAM_RESIDU])
      || isinf(localsolver_options->dparam[SICONOS_DPARAM_RESIDU])
      || localsolver_options->dparam[SICONOS_DPARAM_RESIDU] > 1.0)
  {
    DEBUG_EXPR(rollingFrictionContact_display(localproblem));
    DEBUG_PRINTF("Discard local reaction for contact %i at iteration %i "
                 "with local_error = %e\n",
                 contact, iter, localsolver_options->dparam[SICONOS_DPARAM_RESIDU]);
    
    numerics_printf("Discard local reaction for contact %i at iteration %i "
                    "with local_error = %e",
                    contact, iter, localsolver_options->dparam[SICONOS_DPARAM_RESIDU]);
  }
  else
    memcpy(&reaction[contact*5], localreaction, sizeof(double)*5);
}

static
void acceptLocalReactionUnconditionally(unsigned int contact,
                                        double *reaction, double localreaction[5])
{
  memcpy(&reaction[contact*5], localreaction, sizeof(double)*5);
}

static
double calculateLightError(double light_error_sum, unsigned int nc, double *reaction)
{
  DEBUG_BEGIN("calculateLightError(...)\n");
  double error = sqrt(light_error_sum);
  double norm_r = cblas_dnrm2(nc*5 , reaction , 1);
  if (fabs(norm_r) > DBL_EPSILON)
    error /= norm_r;
  DEBUG_PRINTF("error = %f\n", error);
  DEBUG_END("calculateLightError(...)\n");
  return error;
}

static
double calculateFullErrorAdaptiveInterval(RollingFrictionContactProblem *problem,
                                          RollingComputeErrorPtr computeError,
                                          SolverOptions *options, int iter,
                                          double *reaction, double *velocity,
                                          double tolerance, double norm_q)
{
  double error=1e+24;
  if (options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] >0)
  {
    if (iter % options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] == 0) {
      (*computeError)(problem, reaction , velocity, tolerance, options, norm_q,  &error);
      if (error > tolerance
          && options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE)
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] *= 2;
    }
    numerics_printf("--------------- RFC3D - NSGS - Iteration %i "
                    "options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = %i, options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = % i",
                    iter, options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY], options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION]);
  }
  else
    (*computeError)(problem, reaction , velocity, tolerance, options, norm_q,  &error);

  return error;
}



static
double calculateFullErrorFinal(RollingFrictionContactProblem *problem, SolverOptions *options,
                               RollingComputeErrorPtr computeError,
                               double *reaction, double *velocity, double tolerance,
                               double norm_q)
{
  double absolute_error;
  (*computeError)(problem, reaction , velocity, tolerance,
                  options, norm_q, &absolute_error);



  if (verbose > 0)
  {
    if (absolute_error > options->dparam[SICONOS_DPARAM_TOL])
    {
      numerics_printf("------- RFC3D - NSGS - Warning absolute "
                      "Residual = %14.7e is larger than required precision = %14.7e",
                      absolute_error, options->dparam[SICONOS_DPARAM_TOL]);
    }
    else
    {
      numerics_printf("------- RFC3D - NSGS - absolute "
                      "Residual = %14.7e is smaller than required precision = %14.7e",
                      absolute_error, options->dparam[SICONOS_DPARAM_TOL]);
    }
  }
  return absolute_error;
}


static
int determine_convergence(double error, double tolerance, int iter,
                          SolverOptions *options)
{
  int hasNotConverged = 1;
  if (error < tolerance)
  {
    hasNotConverged = 0;
    numerics_printf("--------------- RFC3D - NSGS - Iteration %i "
                    "Residual = %14.7e < %7.3e\n", iter, error, tolerance);
  }
  else
  {
    numerics_printf("--------------- RFC3D - NSGS - Iteration %i "
                    "Residual = %14.7e > %7.3e\n", iter, error, tolerance);
  }
  return hasNotConverged;
}

static
int determine_convergence_with_full_final(RollingFrictionContactProblem *problem, SolverOptions *options,
                                          RollingComputeErrorPtr computeError,
                                          double *reaction, double *velocity,
                                          double *tolerance, double norm_q, double error,
                                          int iter)
{
  int hasNotConverged = 1;
  if (error < *tolerance)
  {
    hasNotConverged = 0;
    numerics_printf("--------------- RFC3D - NSGS - Iteration %i "
                    "Residual = %14.7e < %7.3e", iter, error, *tolerance);

    double absolute_error = calculateFullErrorFinal(problem, options,
                                                    computeError,
                                                    reaction, velocity,
                                                    options->dparam[SICONOS_DPARAM_TOL],
                                                    norm_q);
    if (absolute_error > options->dparam[SICONOS_DPARAM_TOL])
    {
      *tolerance = error/absolute_error*options->dparam[SICONOS_DPARAM_TOL];
      assert(*tolerance > 0.0 && "tolerance has to be positive");
      /* if (*tolerance < DBL_EPSILON) */
      /* { */
      /*   numerics_warning("determine_convergence_with_full_final", "We try to set a very smal tolerance"); */
      /*   *tolerance = DBL_EPSILON; */
      /* } */
      numerics_printf("------- RFC3D - NSGS - We modify the required incremental precision to reach accuracy to %e", *tolerance);
      hasNotConverged = 1;
    }
    else
    {
      numerics_printf("------- RFC3D - NSGS - The incremental precision is sufficient to reach accuracy to %e", *tolerance);
    }




  }
  else
  {
    numerics_printf("--------------- RFC3D - NSGS - Iteration %i "
                    "Residual = %14.7e > %7.3e", iter, error, *tolerance);
  }
  return hasNotConverged;
}


static
void statsIterationCallback(RollingFrictionContactProblem *problem,
                            SolverOptions *options,
                            double *reaction, double *velocity, double error)
{
  if (options->callback)
  {
    options->callback->collectStatsIteration(options->callback->env,
                                             problem->numberOfContacts * 5,
                                             reaction, velocity,
                                             error, NULL);
  }
}




void rolling_fc3d_nsgs(RollingFrictionContactProblem* problem, double *reaction,
               double *velocity, int* info, SolverOptions* options)
{

  /* problem->mu_r[0]=0.1; */
  /* problem->mu[0]=1.0; */

  
  /* verbose=1; */
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* Number of contacts */
  unsigned int nc = problem->numberOfContacts;

  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];

  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  double norm_q = cblas_dnrm2(nc*5 , problem->q , 1);
  double omega = dparam[SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE];

  SolverOptions * localsolver_options = options->internalSolvers;
  RollingSolverPtr local_solver = NULL;
  RollingUpdatePtr update_localproblem = NULL;
  RollingFreeSolverNSGSPtr freeSolver = NULL;
  RollingComputeErrorPtr computeError = NULL;

  RollingFrictionContactProblem* localproblem;
  double localreaction[5];
  
  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  unsigned int contact; /* Number of the current row of blocks in M */
  unsigned int *scontacts = NULL;

  if (*info == 0)
    return;

  if (options->numberOfInternalSolvers < 1)
  {
    numerics_error("rolling_fc3d_nsgs",
                   "The NSGS method needs options for the internal solvers, "
                   "options[0].numberOfInternalSolvers should be >= 1");
  }
  assert(options->internalSolvers);

  /*****  Initialize various solver options *****/
  localproblem = rolling_fc3d_local_problem_allocate(problem);

  rolling_fc3d_nsgs_initialize_local_solver(&local_solver, &update_localproblem,
                                            (RollingFreeSolverNSGSPtr *)&freeSolver, &computeError,
                                            problem, localproblem, options);

  scontacts = allocShuffledContacts(problem, options);

  /*****  Check solver options *****/
  if (! (iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE
         || iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE
         || iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP))
  {
    numerics_error(
      "rolling_fc3d_nsgs", "iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] must be equal to "
      "SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE (0), "
      "SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE (1) or "
      "SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP (2)");
    return;
  }

  if (! (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL
         || iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL
         || iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT
         || iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE))
  {
    numerics_error(
      "rolling_fc3d_nsgs", "iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] must be equal to "
      "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL (0), "
      "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (1), "
      "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT (2) or "
      "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE (3)");
    return;
  }

  /*****  NSGS Iterations *****/

  /* A special case for the most common options (should correspond
   * with mechanics_run.py **/
  if (iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE
      && iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE
      && iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] == SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE
      && iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT)
  {
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      double light_error_sum = 0.0;

      rolling_fc3d_set_internalsolver_tolerance(problem, options, &localsolver_options[0], error);

      for (unsigned int i = 0 ; i < nc ; ++i)
      {
        contact = i;


        solveLocalReaction(update_localproblem, local_solver, contact,
                           problem, localproblem, reaction, localsolver_options,
                           localreaction);

        accumulateLightErrorSum(&light_error_sum, localreaction, &reaction[contact*5]);

        /* #if 0 */
        acceptLocalReactionFiltered(localproblem, localsolver_options,
                                    contact, iter, reaction, localreaction);
      }

      error = calculateLightError(light_error_sum, nc, reaction);

      hasNotConverged = determine_convergence(error, tolerance, iter, options);

      statsIterationCallback(problem, options, reaction, velocity, error);
    }
  }

  /* All other cases, we put all the ifs inline.. otherwise, too many
   * variations to have dedicated loops, but add more if there are
   * common cases to avoid checking booleans on every iteration. **/
  else
  {
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      double light_error_sum = 0.0;
      rolling_fc3d_set_internalsolver_tolerance(problem, options, &localsolver_options[0], error);

      for (unsigned int i = 0 ; i < nc ; ++i)
      {
        if (iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE
            || iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)
        {
          if (iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)
            uint_shuffle(scontacts, nc);
          contact = scontacts[i];
        }
        else
          contact = i;


        solveLocalReaction(update_localproblem, local_solver, contact,
                           problem, localproblem, reaction, localsolver_options,
                           localreaction);

        if (iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] == SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE)
          performRelaxation(localreaction, &reaction[contact*3], omega);

        accumulateLightErrorSum(&light_error_sum, localreaction, &reaction[contact*5]);

        /* int test =100; */
        /* if (contact == test) */
        /* { */
        /*   printf("reaction[%i] = %16.8e\t",3*contact-1,reaction[3*contact]); */
        /*   printf("localreaction[%i] = %16.8e\n",2,localreaction[0]); */
        /* } */


        if (iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] == SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE)
          acceptLocalReactionFiltered(localproblem, localsolver_options,
                                      contact, iter, reaction, localreaction);
        else
          acceptLocalReactionUnconditionally(contact, reaction, localreaction);



      }

      if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT)
      {
        error = calculateLightError(light_error_sum, nc, reaction);
        hasNotConverged = determine_convergence(error, tolerance, iter, options);
      }
      else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
      {
        error = calculateLightError(light_error_sum, nc, reaction);
        hasNotConverged = determine_convergence_with_full_final(problem,  options, computeError,
                                                                reaction, velocity,
                                                                &tolerance, norm_q, error,
                                                                iter);

      }
      else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL)
      {
        error = calculateFullErrorAdaptiveInterval(problem, computeError, options,
                                                   iter, reaction, velocity,
                                                   tolerance, norm_q);
        hasNotConverged = determine_convergence(error, tolerance, iter, options);
      }

      statsIterationCallback(problem, options, reaction, velocity, error);
    }
  }


  /* Full criterium */
  if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
  {
    error = calculateFullErrorFinal(problem, options, computeError, reaction, velocity,
                                    tolerance, norm_q);

    hasNotConverged = determine_convergence(error,  dparam[SICONOS_DPARAM_TOL] , iter, options);


  }

  *info = hasNotConverged;



  /** return parameter values */
  /* dparam[SICONOS_DPARAM_TOL] = tolerance; */
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

  /** Free memory **/
  (*freeSolver)(problem,localproblem,localsolver_options);
  rolling_fc3d_local_problem_free(localproblem, problem);
  if (scontacts) free(scontacts);
}

int rolling_fc3d_nsgs_setDefaultSolverOptions(SolverOptions* options)
{
  numerics_printf_verbose(1,"rolling_fc3d_nsgs_setDefaultSolverOptions");

  /*  strcpy(options->solverName,"NSGS");*/
  options->solverId = SICONOS_ROLLING_FRICTION_3D_NSGS;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE;
  /* options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE; */
  /* options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT; */

  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] = 10.0;


  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));

  rolling_fc3d_projectionOnConeWithLocalIteration_setDefaultSolverOptions(options->internalSolvers);

  return 0;
}
