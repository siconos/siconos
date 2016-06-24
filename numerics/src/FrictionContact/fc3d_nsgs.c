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
void initializeLocalSolver_nsgs(SolverPtr* solve, UpdatePtr* update, FreeSolverNSGSPtr* freeSolver, ComputeErrorPtr* computeError,
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
  /* Newton solver (Glocker-Fischer-Burmeister)*/
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
    fprintf(stderr, "Numerics, fc3d_nsgs failed. Unknown internal solver : %s.\n", idToName(localsolver_options->solverId));
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
    rowProdNoDiagSBM(n, 3, contact, problem->M->matrix1, reaction, qLocal, 0);
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
  unsigned int contact; /* Number of the current row of blocks in M */

  unsigned int *scontacts = NULL;
  if (iparam[5]) /* shuffle */
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
  if (iparam[1] == 1 || iparam[1] == 2)
  {
    while ((iter < itermax) && (hasNotConverged > 0))
    {

      ++iter;
      /* Loop through the contact points */
      //cblas_dcopy( n , q , incx , velocity , incy );
      error = 0.0;

      #if defined(_OPENMP) && defined(USE_OPENMP)
      #pragma omp parallel for private(contact) reduction(+:error)
      #endif
      for (unsigned int i= 0 ; i < nc ; ++i)
      {
        double reactionold[3];
        if (iparam[5])
        {
          contact = scontacts[i];
        }
        else
        {
          contact = i;
        }

        #if defined(_OPENMP) && defined(USE_OPENMP)
        unsigned int tid = omp_get_thread_num();
        #else
        unsigned int tid = 0;
        #endif

        reactionold[0] = reaction[3 * contact];
        reactionold[1] = reaction[3 * contact + 1];
        reactionold[2] = reaction[3 * contact + 2];
        if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
        (*update_localproblem)(contact, problem, localproblems[tid],
                               reaction, localsolvoptions[tid]);

        localsolvoptions[tid]->iparam[4] = contact;
        (*local_solver)(localproblems[tid], &(reaction[3 * contact]),
                        localsolvoptions[tid]);

        error += pow(reaction[3 * contact] - reactionold[0], 2) +
                 pow(reaction[3 * contact + 1] - reactionold[1], 2) +
                 pow(reaction[3 * contact + 2] - reactionold[2], 2);

      }


      /* **** Criterium convergence **** */
      error = sqrt(error);

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

    if (iparam[1] == 1) /* Full criterium */
    {
      double absolute_error;
      (*computeError)(problem, reaction , velocity, tolerance, options, &absolute_error);
      if (verbose > 0)
      {
        if (absolute_error > error)
        {
          printf("----------------------------------- FC3D - NSGS - Warning absolute Residual = %14.7e is larger than incremental error = %14.7e\n", absolute_error, error);
        }
      }
    }
  }
  else
  {
    if (iparam[5] == 1) /* shuffle */
    {
      int withRelaxation=iparam[4];
      if (withRelaxation)
      {
        double omega = dparam[8];
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );
          #if defined(_OPENMP) && defined(USE_OPENMP)
          #pragma omp parallel for private(contact)
          #endif
          for (unsigned int i= 0 ; i < nc ; ++i)
          {
            double reactionold[3];

            contact = scontacts[i];

            #if defined(_OPENMP) && defined(USE_OPENMP)
            unsigned int tid = omp_get_thread_num();
            #else
            unsigned int tid = 0;
            #endif

            reactionold[0] = reaction[3 * contact];
            reactionold[1] = reaction[3 * contact + 1];
            reactionold[2] = reaction[3 * contact + 2];

            if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
            (*update_localproblem)(contact, problem, localproblems[tid],
                                   reaction, localsolvoptions[tid]);
            localsolvoptions[tid]->iparam[4] = contact;
            (*local_solver)(localproblems[tid], &(reaction[3 * contact]),
                            localsolvoptions[tid]);

            reaction[3 * contact] = omega*reaction[3 * contact]+(1.0-omega)*reactionold[0];
            reaction[3 * contact+1] = omega*reaction[3 * contact+1]+(1.0-omega)*reactionold[1];
            reaction[3 * contact+2] = omega*reaction[3 * contact+2]+(1.0-omega)*reactionold[2];
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
        } /* end while loop */
      }
      else /* without relaxation but shuffle */
      {
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );

          #if defined(_OPENMP) && defined(USE_OPENMP)
          #pragma omp parallel for private(contact)
          #endif
          for (unsigned int i= 0 ; i < nc ; ++i)
          {
            contact = scontacts[i];
            #if defined(_OPENMP) && defined(USE_OPENMP)
            unsigned int tid = omp_get_thread_num();
            #else
            unsigned int tid = 0;
            #endif

            if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
            (*update_localproblem)(contact, problem, localproblems[tid],
                                   reaction, localsolvoptions[tid]);
            localsolvoptions[tid]->iparam[4] = contact;
            (*local_solver)(localproblems[tid], &(reaction[3 * contact]),
                            localsolvoptions[tid]);
          }

          /* **** Criterium convergence **** */
          (*computeError)(problem, reaction , velocity, tolerance, options, &error);

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
    }
    else if (iparam[5] == 2) /* shuffle in each loop */
    {
      int withRelaxation=iparam[4];
      if (withRelaxation)
      {
        double omega = dparam[8];
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          uint_shuffle(scontacts, nc);

          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );

          #if defined(_OPENMP) && defined(USE_OPENMP)
          #pragma omp parallel for private(contact)
          #endif
          for (unsigned int i= 0 ; i < nc ; ++i)
          {
            double reactionold[3];

            contact = scontacts[i];

            #if defined(_OPENMP) && defined(USE_OPENMP)
            unsigned int tid = omp_get_thread_num();
            #else
            unsigned int tid = 0;
            #endif

            reactionold[0] = reaction[3 * contact];
            reactionold[1] = reaction[3 * contact + 1];
            reactionold[2] = reaction[3 * contact + 2];

            if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
            (*update_localproblem)(contact, problem, localproblems[tid],
                                   reaction, localsolvoptions[tid]);

            localsolvoptions[tid]->iparam[4] = contact; /* We write in dWork always with respect to the initial index i*/

            (*local_solver)(localproblems[tid], &(reaction[3 * contact]),
                            localsolvoptions[tid]);

            reaction[3 * contact] = omega*reaction[3 * contact]+(1.0-omega)*reactionold[0];
            reaction[3 * contact+1] = omega*reaction[3 * contact+1]+(1.0-omega)*reactionold[1];
            reaction[3 * contact+2] = omega*reaction[3 * contact+2]+(1.0-omega)*reactionold[2];
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

      }
      else
      {
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          uint_shuffle(scontacts, nc);
          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );
          #if defined(_OPENMP) && defined(USE_OPENMP)
          #pragma omp parallel for private(contact)
          #endif
          for (unsigned int i= 0 ; i < nc ; ++i)
          {
            #if defined(_OPENMP) && defined(USE_OPENMP)
            unsigned int tid = omp_get_thread_num();
            #else
            unsigned int tid = 0;
            #endif
            contact = scontacts[i];
            if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
            (*update_localproblem)(contact, problem, localproblems[tid],
                                   reaction, localsolvoptions[tid]);
            localsolvoptions[tid]->iparam[4] = contact; /* We write in dWork always with respect to the initial index i*/
            (*local_solver)(localproblems[tid], &(reaction[3 * contact]),
                            localsolvoptions[tid]);
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
      }
    }
    else
    {
      int withRelaxation=iparam[4];
      if(withRelaxation)
      {
        double omega = dparam[8];
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );
          #if defined(_OPENMP) && defined(USE_OPENMP)
          #pragma omp parallel for private(contact)
          #endif
          for (contact= 0 ; contact < nc ; ++contact)
          {
            double reactionold[3];
            if (withRelaxation)
            {
              reactionold[0] = reaction[3 * contact];
              reactionold[1] = reaction[3 * contact + 1];
              reactionold[2] = reaction[3 * contact + 2];
            }

            #if defined(_OPENMP) && defined(USE_OPENMP)
            unsigned int tid = omp_get_thread_num();
            #else
            unsigned int tid = 0;
            #endif

            if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
            (*update_localproblem)(contact, problem, localproblems[tid],
                                   reaction, localsolvoptions[tid]);
            localsolvoptions[tid]->iparam[4] = contact;
            (*local_solver)(localproblems[tid], &(reaction[3 * contact]),
                            localsolvoptions[tid]);
            if(withRelaxation)
            {
              reaction[3 * contact] = omega*reaction[3 * contact]+(1.0-omega)*reactionold[0];
              reaction[3 * contact+1] = omega*reaction[3 * contact+1]+(1.0-omega)*reactionold[1];
              reaction[3 * contact+2] = omega*reaction[3 * contact+2]+(1.0-omega)*reactionold[2];
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

      }
      else
      {
        while ((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          /* Loop through the contact points */
          //cblas_dcopy( n , q , incx , velocity , incy );

          #if defined(_OPENMP) && defined(USE_OPENMP)
          #pragma omp parallel for private(contact)
          #endif
          for (contact= 0 ; contact < nc ; ++contact)
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
            (*local_solver)(localproblems[tid], &(reaction[3 * contact]),
                            localsolvoptions[tid]);
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
      }
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
    free(localsolvoptions[i]);
  }

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
  options->iSize = 10;
  options->dSize = 10;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  null_SolverOptions(options);
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
