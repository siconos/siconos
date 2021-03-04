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
#include <assert.h>                            // for assert
#include <float.h>                             // for DBL_EPSILON
#include <math.h>                              // for pow, fabs, sqrt, isinf
#include <stdlib.h>                            // for calloc, malloc, srand
#include <string.h>                            // for NULL, memcpy
#include "SiconosBlas.h"                             // for cblas_dnrm2
#include "Friction_cst.h"                      // for SICONOS_FRICTION_3D_IP...
#include "NumericsArrays.h"                    // for uint_shuffle
#include "NumericsFwd.h"                       // for SolverOptions, Rolling...
#include "RollingFrictionContactProblem.h"     // for RollingFrictionContact...
#include "SolverOptions.h"                     // for SolverOptions, SICONOS...
#include "siconos_debug.h"                             // for DEBUG_PRINTF, DEBUG_BEGIN
#include "numerics_verbose.h"                  // for numerics_printf, numer...
#include "rolling_fc_Solvers.h"              // for RollingComputeErrorPtr
#include "rolling_fc3d_compute_error.h"        // for rolling_fc3d_compute_e...
#include "rolling_fc3d_local_problem_tools.h"  // for rolling_fc3d_local_pro...
#include "rolling_fc3d_projection.h"           // for rolling_fc3d_projectio...

//#define FCLIB_OUTPUT

#ifdef FCLIB_OUTPUT
static int fccounter = -1;
#include "fclib_interface.h"
#endif



#pragma GCC diagnostic ignored "-Wmissing-prototypes"

/* static void rolling_fc3d_nsgs_update(int contact, RollingFrictionContactProblem* problem, RollingFrictionContactProblem* localproblem, double * reaction, SolverOptions* options) */
/* { */
/*   /\* Build a local problem for a specific contact */
/*      reaction corresponds to the global vector (size n) of the global problem. */
/*   *\/ */
/*   /\* Call the update function which depends on the storage for MGlobal/MBGlobal *\/ */
/*   /\* Build a local problem for a specific contact */
/*      reaction corresponds to the global vector (size n) of the global problem. */
/*   *\/ */

/*   /\* The part of MGlobal which corresponds to the current block is copied into MLocal *\/ */
/*   rolling_fc3d_local_problem_fill_M(problem, localproblem, contact); */

/*   /\****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock, */
/*          excluding the block corresponding to the current contact. ****\/ */
/*   rolling_fc3d_local_problem_compute_q(problem, localproblem, reaction, contact); */

/*   /\* Friction coefficient for current block*\/ */
/*   localproblem->mu[0] = problem->mu[contact]; */
/*   /\* Rolling Friction coefficient for current block*\/ */
/*   localproblem->mu_r[0] = problem->mu_r[contact]; */


/* } */

void rolling_fc3d_nsgs_initialize_local_solver(RollingSolverPtr* solve, RollingUpdatePtr* update,
    RollingFreeSolverNSGSPtr* freeSolver,
    RollingComputeErrorPtr* computeError,
    RollingFrictionContactProblem* problem,
    RollingFrictionContactProblem* localproblem,
    SolverOptions * options)
{
  SolverOptions * localsolver_options = options->internalSolvers[0];
  /** Connect to local solver */
  switch(localsolver_options->solverId)
  {
  case SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration:
  {
    *solve = &rolling_fc3d_projectionOnConeWithLocalIteration_solve;
    *update = &rolling_fc3d_projection_update;
    *freeSolver = (RollingFreeSolverNSGSPtr)&rolling_fc3d_projectionOnConeWithLocalIteration_free;
    *computeError = (RollingComputeErrorPtr)&rolling_fc3d_compute_error;
    rolling_fc3d_projectionOnConeWithLocalIteration_initialize(problem, localproblem,localsolver_options);
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

// Never used ...
/* static */
/* unsigned int* allocShuffledContacts(RollingFrictionContactProblem *problem, */
/*                                     SolverOptions *options) */
/* { */
/*   unsigned int *scontacts = 0; */
/*   unsigned int nc = problem->numberOfContacts; */
/*   if(options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE|| */
/*       options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP) */
/*   { */
/*     if(options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED] >0) */
/*     { */
/*       srand((unsigned int)options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED]); */
/*     } */
/*     else */
/*       srand(1); */
/*     scontacts = (unsigned int *) malloc(nc * sizeof(unsigned int)); */
/*     for(unsigned int i = 0; i < nc ; ++i) */
/*     { */
/*       scontacts[i] = i; */
/*     } */
/*     uint_shuffle(scontacts, nc); */
/*   } */
/*   return scontacts; */
/* } */

static
unsigned int* allocfreezingContacts(RollingFrictionContactProblem *problem,
                                    SolverOptions *options)
{
  unsigned int *fcontacts = 0;
  unsigned int nc = problem->numberOfContacts;
  if(options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0)
  {
    fcontacts = (unsigned int *) malloc(nc * sizeof(unsigned int));
    for(unsigned int i = 0; i < nc ; ++i)
    {
      fcontacts[i] = 0;
    }
  }
  return fcontacts;
}

static
int solveLocalReaction(RollingUpdatePtr update_localproblem, RollingSolverPtr local_solver,
                       unsigned int contact, RollingFrictionContactProblem *problem,
                       RollingFrictionContactProblem *localproblem, double *reaction,
                       SolverOptions *localsolver_options, double localreaction[5])
{
  (*update_localproblem)(contact, problem, localproblem,
                         reaction, localsolver_options);

  localsolver_options->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] = contact;

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
double light_error_squared( double localreaction[5],
                            double *oldreaction)
{
  return (pow(oldreaction[0] - localreaction[0], 2) +
          pow(oldreaction[1] - localreaction[1], 2) +
          pow(oldreaction[2] - localreaction[2], 2) +
          pow(oldreaction[3] - localreaction[3], 2) +
          pow(oldreaction[4] - localreaction[4], 2));
}

static
double squared_norm(double localreaction[5])
{
  return (pow(localreaction[0], 2) +
          pow(localreaction[1], 2) +
          pow(localreaction[2], 2) +
          pow(localreaction[3], 2) +
          pow(localreaction[4], 2));
}

static
void acceptLocalReactionFiltered(RollingFrictionContactProblem *localproblem,
                                 SolverOptions *localsolver_options,
                                 unsigned int contact, unsigned int iter,
                                 double *reaction, double localreaction[5])
{
  if(isnan(localsolver_options->dparam[SICONOS_DPARAM_RESIDU])
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
double calculateLightError(double light_error_sum, unsigned int nc, double *reaction, double * norm_r)
{
  DEBUG_BEGIN("calculateLightError(...)\n");
  double error = sqrt(light_error_sum);
  *norm_r = cblas_dnrm2(nc*5, reaction, 1);
  if(fabs(*norm_r) > DBL_EPSILON)
    error /= (*norm_r);
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
  if(options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] >0)
  {
    if(iter % options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] == 0)
    {
      (*computeError)(problem, reaction, velocity, tolerance, options, norm_q,  &error);
      if(error > tolerance
          && options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE)
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] *= 2;
    }
    numerics_printf("--------------- RFC3D - NSGS - Iteration %i "
                    "options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = %i, options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = % i",
                    iter, options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY], options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION]);
  }
  else
    (*computeError)(problem, reaction, velocity, tolerance, options, norm_q,  &error);

  return error;
}



static
double calculateFullErrorFinal(RollingFrictionContactProblem *problem, SolverOptions *options,
                               RollingComputeErrorPtr computeError,
                               double *reaction, double *velocity, double tolerance,
                               double norm_q)
{
  double absolute_error;
  (*computeError)(problem, reaction, velocity, tolerance,
                  options, norm_q, &absolute_error);



  if(verbose > 0)
  {
    if(absolute_error > options->dparam[SICONOS_DPARAM_TOL])
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
  if(error < tolerance)
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
  if(error < *tolerance)
  {
    hasNotConverged = 0;
    numerics_printf("--------------- RFC3D - NSGS - Iteration %i "
                    "Residual = %14.7e < %7.3e", iter, error, *tolerance);

    double absolute_error = calculateFullErrorFinal(problem, options,
                            computeError,
                            reaction, velocity,
                            options->dparam[SICONOS_DPARAM_TOL],
                            norm_q);
    if(absolute_error > options->dparam[SICONOS_DPARAM_TOL])
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
  if(options->callback)
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
  double norm_q = cblas_dnrm2(nc*5, problem->q, 1);
  double omega = dparam[SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE];
 
  double  norm_r[] = {1e24};
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
  unsigned int *freeze_contacts = NULL;

  if(*info == 0)
    return;

  if(options->numberOfInternalSolvers < 1)
  {
    numerics_error("rolling_fc3d_nsgs",
                   "The NSGS method needs options for the internal solvers, "
                   "options[0].numberOfInternalSolvers should be >= 1");
  }
  SolverOptions * localsolver_options = options->internalSolvers[0];

  /*****  Initialize various solver options *****/
  localproblem = rolling_fc3d_local_problem_allocate(problem);

  rolling_fc3d_nsgs_initialize_local_solver(&local_solver, &update_localproblem,
      (RollingFreeSolverNSGSPtr *)&freeSolver,
      &computeError,
      problem, localproblem, options);

  freeze_contacts = allocfreezingContacts(problem, options);

  /*****  Check solver options *****/
  if(!(iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE
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

  if(!(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL
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
  if(iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE
     && iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] == 0
     && iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE
     && iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] == SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE
     && iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT)
  {
    while((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      double light_error_sum = 0.0;

      rolling_fc3d_set_internalsolver_tolerance(problem, options, localsolver_options, error);

      for(unsigned int i = 0 ; i < nc ; ++i)
      {
        contact = i;


        solveLocalReaction(update_localproblem, local_solver, contact,
                           problem, localproblem, reaction, localsolver_options,
                           localreaction);

        light_error_sum += light_error_squared(localreaction, &reaction[contact*5]);

        /* #if 0 */
        acceptLocalReactionFiltered(localproblem, localsolver_options,
                                    contact, iter, reaction, localreaction);
      }

      error = calculateLightError(light_error_sum, nc, reaction, norm_r);

      hasNotConverged = determine_convergence(error, tolerance, iter, options);

      statsIterationCallback(problem, options, reaction, velocity, error);
    }
  }

  /* All other cases, we put all the ifs inline.. otherwise, too many
   * variations to have dedicated loops, but add more if there are
   * common cases to avoid checking booleans on every iteration. **/
  else
  {
    while((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      double light_error_sum = 0.0;
      double light_error_2 = 0.0;

      rolling_fc3d_set_internalsolver_tolerance(problem, options, localsolver_options, error);

      for(unsigned int i = 0 ; i < nc ; ++i)
      {
        if(iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE
            || iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)
        {
          if(iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)
            uint_shuffle(scontacts, nc);
          contact = scontacts[i];
        }
        else
          contact = i;


        if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0)
        {
          if (freeze_contacts[contact] >0)
          {
            /* we skip freeze contacts */
            freeze_contacts[contact] -=  1;
            continue;
          }
        }

        solveLocalReaction(update_localproblem, local_solver, contact,
                           problem, localproblem, reaction, localsolver_options,
                           localreaction);

        if(iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] == SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE)
          performRelaxation(localreaction, &reaction[contact*5], omega);

        light_error_2= light_error_squared(localreaction, &reaction[contact*5]);
        light_error_sum += light_error_2;


        /* int test =100; */
        /* if (contact == test) */
        /* { */
        /*   printf("reaction[%i] = %16.8e\t",3*contact-1,reaction[3*contact]); */
        /*   printf("localreaction[%i] = %16.8e\n",2,localreaction[0]); */
        /* } */

        if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0)
        {
          if ((light_error_2*squared_norm(localreaction) <= tolerance*tolerance/(nc*nc*10)
               || squared_norm(localreaction) <=  (*norm_r* *norm_r/(nc*nc*1000)))
              && iter >=10)
          {
            /* we  freeze the contact for n iterations*/
            //printf("first criteria : light_error_2*squared_norm(localreaction) <= tolerance*tolerance/(nc*nc*10) ==> %e <= %e\n", light_error_2*squared_norm(localreaction), tolerance*tolerance/(nc*nc*10));
            //printf("second criteria :  squared_norm(localreaction) <=  (*norm_r* *norm_r/(nc*nc))/1000. ==> %e <= %e\n",  squared_norm(localreaction) ,  (*norm_r* *norm_r/(nc*nc))/1000.);
            //printf("Contact % i is freezed for %i steps\n", contact,  iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]);
            freeze_contacts[contact] = iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] ;
          }
        }

        if(iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] == SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE)
          acceptLocalReactionFiltered(localproblem, localsolver_options,
                                      contact, iter, reaction, localreaction);
        else
          acceptLocalReactionUnconditionally(contact, reaction, localreaction);

      }


      if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT)
      {
        error = calculateLightError(light_error_sum, nc, reaction, norm_r);
        hasNotConverged = determine_convergence(error, tolerance, iter, options);
      }
      else if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
      {
        error = calculateLightError(light_error_sum, nc, reaction, norm_r);
        hasNotConverged = determine_convergence_with_full_final(problem,  options, computeError,
                          reaction, velocity,
                          &tolerance, norm_q, error,
                          iter);

      }
      else if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL)
      {
        error = calculateFullErrorAdaptiveInterval(problem, computeError, options,
                iter, reaction, velocity,
                tolerance, norm_q);
        hasNotConverged = determine_convergence(error, tolerance, iter, options);
      }

      statsIterationCallback(problem, options, reaction, velocity, error);

    }
    /* if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0) */
    /* { */
    /*   int frozen_contact=0; */
    /*   for(unsigned int i = 0 ; i < nc ; ++i) */
    /*   { */
    /*     if (freeze_contacts[contact] >0) */
    /*     { */
    /*       frozen_contact++; */
    /*     } */
    /*   } */
    /*   printf("number of frozen contacts : %i\n", frozen_contact ); */
    /* } */
  }

  /* Full criterium */
  if(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
  {
    error = calculateFullErrorFinal(problem, options, computeError, reaction, velocity,
                                    tolerance, norm_q);

    hasNotConverged = determine_convergence(error,  dparam[SICONOS_DPARAM_TOL], iter, options);


  }

  *info = hasNotConverged;



  /** return parameter values */
  /* dparam[SICONOS_DPARAM_TOL] = tolerance; */
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

  /** Free memory **/
  (*freeSolver)(problem,localproblem,localsolver_options);
  rolling_fc3d_local_problem_free(localproblem, problem);
  if(scontacts) free(scontacts);
}

void rfc3d_nsgs_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE;
  /* options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE; */
  /* options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT; */
  options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] =  SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE;
  options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED] = 0;
  options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 0;
  options->iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] = SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_FALSE;
  options->iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] = SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 0;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] = 10.0;

  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] = solver_options_create(SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
}
