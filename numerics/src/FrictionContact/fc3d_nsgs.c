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
#include <float.h>

#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include "debug.h"

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void fake_compute_error_nsgs(FrictionContactProblem* problem, double *reaction, double *velocity, double tolerance, SolverOptions  *options,  double* error)
{
  int n = 3 * problem->numberOfContacts;
  *error = 0.;
  int i, m;
  m = 5 * n / 3;
  double err = INFINITY;
  for (i = 0 ; i < m ; ++i)
  {
    *error += Compute_NCP_error1(i, err);
  }
}
void fc3d_nsgs_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem, double * reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  fc3d_nsgs_fillMLocal(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/
  fc3d_nsgs_computeqLocal(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];


}
void fc3d_nsgs_initialize_local_solver(SolverPtr* solve, UpdatePtr* update, FreeSolverNSGSPtr* freeSolver, ComputeErrorPtr* computeError,
                                FrictionContactProblem* problem, FrictionContactProblem* localproblem,
                                SolverOptions * options, SolverOptions * localsolver_options)
{

  /** Connect to local solver */
  switch (localsolver_options->solverId)
  {
    /* Projection */
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization:
  {
    *solve = &fc3d_projectionWithDiagonalization_solve;
    *update = &fc3d_projectionWithDiagonalization_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projection_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_projection_initialize(problem, localproblem);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone:
  {
    *solve = &fc3d_projectionOnCone_solve;
    *update = &fc3d_projection_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projection_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_projection_initialize(problem, localproblem);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration:
  {
    *solve = &fc3d_projectionOnConeWithLocalIteration_solve;
    *update = &fc3d_projection_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projectionOnConeWithLocalIteration_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_projectionOnConeWithLocalIteration_initialize(problem, localproblem,localsolver_options );
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization:
  {
    *solve = &fc3d_projectionOnCone_solve;
    *update = &fc3d_projection_update_with_regularization;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projection_with_regularization_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_projection_initialize_with_regularization(problem, localproblem);
    break;
  }
  /* Newton solver (Alart-Curnier) */
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_AC:
  {
    *solve = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
    *update = &fc3d_onecontact_nonsmooth_Newton_AC_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_onecontact_nonsmooth_Newton_solvers_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP:
  {
    *solve = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
    *update = &fc3d_onecontact_nonsmooth_Newton_AC_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_onecontact_nonsmooth_Newton_solvers_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP_P:
  {
    *solve = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
    *update = &fc3d_onecontact_nonsmooth_Newton_AC_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_onecontact_nonsmooth_Newton_solvers_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
    break;
  }  /* Newton solver (Glocker-Fischer-Burmeister)*/
  case SICONOS_FRICTION_3D_NCPGlockerFBNewton:
  {
    *solve = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
    *update = &NCPGlocker_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_onecontact_nonsmooth_Newton_solvers_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    // *computeError = &fake_compute_error;
    fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
    break;
  }
  /* Path solver (Glocker Formulation) */
  case SICONOS_FRICTION_3D_NCPGlockerFBPATH:
  {
    *solve = &fc3d_Path_solve;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_Path_free;
    *update = &NCPGlocker_update;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    // *computeError = &fake_compute_error;
    fc3d_Path_initialize(problem, localproblem, localsolver_options);
    break;
  }

  /* Fixed Point solver (Glocker Formulation) */
  case SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint:
  {
    *solve = &fc3d_FixedP_solve;
    *update = &NCPGlocker_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_FixedP_free;
    /* *computeError = &fake_compute_error_nsgs; */
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_FixedP_initialize(problem, localproblem, localsolver_options);
    break;
  }
  /*iparam[4] > 10 are reserved for Tresca resolution */
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder:
  {
    *solve = &fc3d_projectionOnCylinder_solve;
    *update = &fc3d_projectionOnCylinder_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projectionOnCylinder_free;
    *computeError = (ComputeErrorPtr)&fc3d_Tresca_compute_error;
    fc3d_projectionOnCylinder_initialize(problem, localproblem, options );
    break;
  }
    /*iparam[4] > 10 are reserved for Tresca resolution */
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration:
  {
    *solve = &fc3d_projectionOnCylinderWithLocalIteration_solve;
    *update = &fc3d_projectionOnCylinder_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_projectionOnCylinderWithLocalIteration_free;
    *computeError = (ComputeErrorPtr)&fc3d_Tresca_compute_error;
    fc3d_projectionOnCylinderWithLocalIteration_initialize(problem, localproblem, options, localsolver_options );
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC:
  {
    *solve = &fc3d_unitary_enumerative_solve;
    *update = &fc3d_nsgs_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_unitary_enumerative_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_unitary_enumerative_initialize(localproblem);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU:
  {
    *solve = &fc3d_unitary_enumerative_solve;
    *update = &fc3d_nsgs_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_unitary_enumerative_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_unitary_enumerative_initialize(localproblem);
    break;
  }
  default:
  {
    fprintf(stderr, "Numerics, fc3d_nsgs failed. Unknown internal solver : %s.\n", solver_options_id_to_name(localsolver_options->solverId));
    exit(EXIT_FAILURE);
  }
  }
}
void fc3d_nsgs_computeqLocal(FrictionContactProblem * problem, FrictionContactProblem * localproblem, double *reaction, int contact)
{

  double *qLocal = localproblem->q;
  int n = 3 * problem->numberOfContacts;


  int in = 3 * contact, it = in + 1, is = it + 1;

  /* qLocal computation*/



  qLocal[0] = problem->q[in];
  qLocal[1] =  problem->q[it];
  qLocal[2] =  problem->q[is];

  if (problem->M->storageType == 0)
  {
    double * MM = problem->M->matrix0;
    int incx = n, incy = 1;
    /* reaction current block set to zero, to exclude current contact block */
    double rin = reaction[in] ;
    double rit = reaction[it] ;
    double ris = reaction[is] ;
    reaction[in] = 0.0;
    reaction[it] = 0.0;
    reaction[is] = 0.0;
    qLocal[0] += cblas_ddot(n , &MM[in] , incx , reaction , incy);
    qLocal[1] += cblas_ddot(n , &MM[it] , incx , reaction , incy);
    qLocal[2] += cblas_ddot(n , &MM[is] , incx , reaction , incy);
    reaction[in] = rin;
    reaction[it] = rit;
    reaction[is] = ris;
  }
  else if (problem->M->storageType == 1)
  {
    /* qLocal += rowMB * reaction
     * with rowMB the row of blocks of MGlobal which corresponds
     * to the current contact
     */
    SBM_row_prod_no_diag_3x3(n, 3, contact, problem->M->matrix1, reaction, qLocal);
  }



}

void fc3d_nsgs_fillMLocal(FrictionContactProblem * problem, FrictionContactProblem * localproblem, int contact)
{

  NumericsMatrix * MGlobal = problem->M;

  int n = 3 * problem->numberOfContacts;


  int storageType = MGlobal->storageType;
  if (storageType == 0)
    // Dense storage
  {
    int in = 3 * contact, it = in + 1, is = it + 1;
    int inc = n * in;
    double * MM = MGlobal->matrix0;
    double * MLocal =  localproblem->M->matrix0;

    /* The part of MM which corresponds to the current block is copied into MLocal */
    MLocal[0] = MM[inc + in];
    MLocal[1] = MM[inc + it];
    MLocal[2] = MM[inc + is];
    inc += n;
    MLocal[3] = MM[inc + in];
    MLocal[4] = MM[inc + it];
    MLocal[5] = MM[inc + is];
    inc += n;
    MLocal[6] = MM[inc + in];
    MLocal[7] = MM[inc + it];
    MLocal[8] = MM[inc + is];
  }
  else if (storageType == 1)
  {
    int diagPos = getDiagonalBlockPos(MGlobal->matrix1, contact);
    localproblem->M->matrix0 = MGlobal->matrix1->block[diagPos];

  }
  else
    numericsError("fc3d_projection -", "unknown storage type for matrix M");
}


/* swap two indices */
void uint_swap (unsigned int *a, unsigned int *b)
{
    unsigned int temp = *a;
    *a = *b;
    *b = temp;
}

/* shuffle an unsigned array */
void uint_shuffle (unsigned int *a, unsigned int n) {

  for (unsigned int i = 0; i < n - 1; i++)
  {
    uint_swap  (&a[i], &a[i + rand()%(n - i)]);
  }
}



void fc3d_nsgs(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
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
  double normq = cblas_dnrm2(nc*3 , problem->q , 1);

  if (*info == 0)
    return;

  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("fc3d_nsgs", "The NSGS method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >= 1");
  }
  assert(options->internalSolvers);

  SolverOptions * localsolver_options = options->internalSolvers;


  SolverPtr local_solver = NULL;
  UpdatePtr update_localproblem = NULL;
  FreeSolverNSGSPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;

  /* Connect local solver and local problem*/
  FrictionContactProblem* localproblem = (FrictionContactProblem*)malloc(sizeof(FrictionContactProblem));
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
                             problem , localproblem,
                             options, localsolver_options);

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  unsigned int contact; /* Number of the current row of blocks in M */

  unsigned int *scontacts = NULL;
  if (iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE||
      iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP) /* shuffle */
  {
    if (iparam[6] >0)
    {
      srand((unsigned int)iparam[6]);
    }
    else
      srand(1);
    scontacts = (unsigned int *) malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i<nc ; ++i)
    {
      scontacts[i] = i;
    }
    uint_shuffle(scontacts, nc);
    /* printf("scontacts :\t"); */
    /* for (unsigned int i = 0; i<nc ; ++i) printf("%i\t", scontacts[i]); */
    /* printf("\n"); */
  }
  /*  dparam[0]= dparam[2]; // set the tolerance for the local solver */
  if (iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL ||
      iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT)
  {
    double local_reaction[3];
    while ((iter < itermax) && (hasNotConverged > 0))
    {

      ++iter;
      /* Loop through the contact points */
      //cblas_dcopy( n , q , incx , velocity , incy );
      error = 0.0;

      for (unsigned int i= 0 ; i < nc ; ++i)
      {
        if (iparam[5])
        {
          contact = scontacts[i];
        }
        else
        {
          contact = i;
        }
        local_reaction[0] = reaction[3 * contact];
        local_reaction[1] = reaction[3 * contact + 1];
        local_reaction[2] = reaction[3 * contact + 2];
        if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
        (*update_localproblem)(contact, problem, localproblem, reaction, localsolver_options);

        localsolver_options->iparam[4] = contact;
        (*local_solver)(localproblem, local_reaction , localsolver_options);

        error += pow(reaction[3 * contact] - local_reaction[0], 2) +
                 pow(reaction[3 * contact + 1] - local_reaction[1], 2) +
                 pow(reaction[3 * contact + 2] - local_reaction[2], 2);

        if ((localsolver_options->dparam[1] > 1.0) && iparam[14] == SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE)
        {
          DEBUG_EXPR(
            frictionContact_display(localproblem);
            printf("Discard local reaction for contact %i at iteration %i with local_error = %e\n", contact, iter, localsolver_options->dparam[1]);
            local_reaction[0] = reaction[3 * contact];
            local_reaction[1] = reaction[3 * contact + 1];
            local_reaction[2] = reaction[3 * contact + 2];
            fc3d_projectionOnConeWithLocalIteration_initialize(problem, localproblem, localsolver_options );
            fc3d_projectionOnConeWithLocalIteration_solve (localproblem, local_reaction , localsolver_options);
            printf("Try local fc3d_projectionOnConeWithLocalIteration_solve with local_error = %e\n", localsolver_options->dparam[1]);
            (*local_solver)(localproblem, local_reaction , localsolver_options);
            printf("Try local another newton solve with local_error = %e\n", localsolver_options->dparam[1]);
            if (localsolver_options->dparam[1] <= localsolver_options->dparam[0])
            {
              DEBUG_PRINTF("Finally keep the local solution = %e\n", localsolver_options->dparam[1]);
              reaction[3 * contact]     = local_reaction[0];
              reaction[3 * contact + 1] = local_reaction[1];
              reaction[3 * contact + 2] = local_reaction[2];
              getchar();
            }
            );

        }
        else
        {
          reaction[3 * contact]     = local_reaction[0];
          reaction[3 * contact + 1] = local_reaction[1];
          reaction[3 * contact + 2] = local_reaction[2];
        }
      }


      /* **** Criterium convergence **** */
      error = sqrt(error);
      double norm_r = cblas_dnrm2(nc*3 , reaction , 1);
      if (fabs(norm_r) > DBL_EPSILON)
        error /= norm_r;


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
    }

    if (iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) /* Full criterium */
    {
      double absolute_error;
      (*computeError)(problem, reaction , velocity, tolerance, options, normq, &absolute_error);
      if (verbose > 0)
      {
        if (absolute_error > error)
        {
          printf("--------------------------- FC3D - NSGS - Warning absolute Residual = %14.7e is larger than incremental error = %14.7e\n", absolute_error, error);
        }
      }
    }
  }
  else if (iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL
           || iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE)
  {
    if (iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE) /* shuffle */
    {
      if (iparam[4] == SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE) /* relaxation */
      {
        double reactionold[3];

        double omega = dparam[8];
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );
          for (unsigned int i= 0 ; i < nc ; ++i)
          {

            contact = scontacts[i];

            reactionold[0] = reaction[3 * contact];
            reactionold[1] = reaction[3 * contact + 1];
            reactionold[2] = reaction[3 * contact + 2];

            if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
            (*update_localproblem)(contact, problem, localproblem, reaction, localsolver_options);
            localsolver_options->iparam[4] = contact;
            (*local_solver)(localproblem, &(reaction[3 * contact]), localsolver_options);

            reaction[3 * contact] = omega*reaction[3 * contact]+(1.0-omega)*reactionold[0];
            reaction[3 * contact+1] = omega*reaction[3 * contact+1]+(1.0-omega)*reactionold[1];
            reaction[3 * contact+2] = omega*reaction[3 * contact+2]+(1.0-omega)*reactionold[2];
          }

          /* **** Criterium convergence **** */
          if (iparam[8] >0)
          {
            if (iter % iparam[8] ==0 )
              (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
          }
          else
            (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);

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
        } /* end while loop */
      }
      else if (iparam[4] == SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE)
      {
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );

          for (unsigned int i= 0 ; i < nc ; ++i)
          {
            contact = scontacts[i];

            if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
            (*update_localproblem)(contact, problem, localproblem, reaction, localsolver_options);
            localsolver_options->iparam[4] = contact;
            (*local_solver)(localproblem, &(reaction[3 * contact]), localsolver_options);

          }

          /* **** Criterium convergence **** */
          if (iparam[8] >0)
          {
            if (iter % iparam[8] ==0 )
              (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
          }
          else
            (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);


          if (verbose > 0)
            printf("----------------------------------- FC3D - NSGS - Iteration %i Residual = %14.7e <= %7.3e\n", iter, error, options->dparam[0]);
          if (error < tolerance) hasNotConverged = 0;
          *info = hasNotConverged;

          if (options->callback)
          {
            options->callback->collectStatsIteration(options->callback->env, 3 * nc,
                                                     reaction, velocity,
                                                     error, NULL);
          }
        }
      }
      else
      {
        numericsError("fc3d_nsgs", "iparam[4] must be equal to SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE (0) or SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE (1)");
      }
    }
    else if (iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP) /* shuffle in each loop */
    {
      if (iparam[4] == SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE)
      {
        double reactionold[3];

        double omega = dparam[8];
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          uint_shuffle(scontacts, nc);

          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );

          for (unsigned int i= 0 ; i < nc ; ++i)
          {

            contact = scontacts[i];

            reactionold[0] = reaction[3 * contact];
            reactionold[1] = reaction[3 * contact + 1];
            reactionold[2] = reaction[3 * contact + 2];

            if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
            (*update_localproblem)(contact, problem, localproblem, reaction, localsolver_options);

            localsolver_options->iparam[4] = contact; /* We write in dWork always with respect to the initial index i*/

            (*local_solver)(localproblem, &(reaction[3 * contact]), localsolver_options);

            reaction[3 * contact] = omega*reaction[3 * contact]+(1.0-omega)*reactionold[0];
            reaction[3 * contact+1] = omega*reaction[3 * contact+1]+(1.0-omega)*reactionold[1];
            reaction[3 * contact+2] = omega*reaction[3 * contact+2]+(1.0-omega)*reactionold[2];
          }

          /* **** Criterium convergence **** */

          if (iparam[8] >0)
          {
            if (iter % iparam[8] ==0 )
              (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
          }
          else
            (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);

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

      }
      else if (iparam[4] == SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE)
      {
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          uint_shuffle(scontacts, nc);
          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );
          for (unsigned int i= 0 ; i < nc ; ++i)
          {
            contact = scontacts[i];
            if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
            (*update_localproblem)(contact, problem, localproblem, reaction, localsolver_options);
            localsolver_options->iparam[4] = contact; /* We write in dWork always with respect to the initial index i*/
            (*local_solver)(localproblem, &(reaction[3 * contact]), localsolver_options);
          }
          /* **** Criterium convergence **** */
          if (iparam[8] >0)
          {
            if (iter % iparam[8] ==0 )
              (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
          }
          else
            (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);


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
      }
      else
      {
        numericsError("fc3d_nsgs", "iparam[4] must be equal to SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE (0) or SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE (1)");
      }
    }
    else if (iparam[5] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE) /* no shufle */
    {

      if(iparam[4] == SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE)
      {
        double reactionold[3];
        double omega = dparam[8];
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );
          for (contact= 0 ; contact < nc ; ++contact)
          {
            reactionold[0] = reaction[3 * contact];
            reactionold[1] = reaction[3 * contact + 1];
            reactionold[2] = reaction[3 * contact + 2];


            if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
            (*update_localproblem)(contact, problem, localproblem, reaction, localsolver_options);
            localsolver_options->iparam[4] = contact;
            (*local_solver)(localproblem, &(reaction[3 * contact]), localsolver_options);
            reaction[3 * contact] = omega*reaction[3 * contact]+(1.0-omega)*reactionold[0];
            reaction[3 * contact+1] = omega*reaction[3 * contact+1]+(1.0-omega)*reactionold[1];
            reaction[3 * contact+2] = omega*reaction[3 * contact+2]+(1.0-omega)*reactionold[2];
          }

          /* **** Criterium convergence **** */

          if (iparam[8] >0)
          {
            if (iter % iparam[8] ==0 )
              (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
          }
          else
            (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);


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

      }
      else if (iparam[4] == SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE)
      {
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );

          for (contact= 0 ; contact < nc ; ++contact)
          {


            if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);

            (*update_localproblem)(contact, problem, localproblem, reaction, localsolver_options);
            localsolver_options->iparam[4] = contact;
            (*local_solver)(localproblem, &(reaction[3 * contact]), localsolver_options);
          }

          /* **** Criterium convergence **** */
          if (iparam[8] >0)
          {
            if (iter % iparam[8] == 0) {
              (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
              if (error > tolerance && iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE)
                iparam[8] *= 2;
            }
            if (verbose > 0)
              printf("----------------------------------- FC3D - NSGS - Iteration %i iparam[8] = %i, iparam[1] = % i \n", iter, iparam[8], iparam[1]);
          }
          else
            (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);

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
        if (iparam[1] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE && (iter == iparam[8] || iparam[8] > itermax))
        {
          iparam[8] = iparam[8]==1 ? 1 : iparam[8]/2;
          if (iter == itermax)
            (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
        }
        if (verbose > 0)
              printf("----------------------------------- FC3D - NSGS - Iteration %i iparam[8] = %i, iparam[1] = % i \n", iter, iparam[8], iparam[1]);

      }
      else
      {
        numericsError("fc3d_nsgs", "iparam[4] must be equal to SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE (0) or SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE (1)");
      }
    }
    else
    {
      numericsError("fc3d_nsgs", "iparam[5] must be equal to SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE (0), SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE (1) or  SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP (2)");
    }
  }
  else
  {
    numericsError("fc3d_nsgs", "iparam[1] must be equal to SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL (0), SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (1), SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT (2) or SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE (3)");
  }

  /** return parameter values */
  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;

  /** Free memory **/
  (*freeSolver)(problem,localproblem,localsolver_options);
  if (problem->M->storageType == NM_DENSE && localproblem->M->matrix0)
  {
    free(localproblem->M->matrix0);
  }
  localproblem->M->matrix0 = NULL;
  freeFrictionContactProblem(localproblem);

  if (scontacts) /* shuffle */
  {
    free(scontacts);
  }

}

int fc3d_nsgs_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGS Solver\n");
  }

  /*  strcpy(options->solverName,"NSGS");*/
  options->solverId = SICONOS_FRICTION_3D_NSGS;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  for (i = 0; i < 10; i++)
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
