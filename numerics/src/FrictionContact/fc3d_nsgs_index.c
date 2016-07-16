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

#pragma GCC diagnostic ignored "-Wmissing-prototypes"


#include "fc3d_nsgs_openmp.h"
void fc3d_nsgs_index_computeqLocal(FrictionContactProblem * problem,
                                   double *reaction, int contact,
                                   unsigned int * index, unsigned int index_size,
                                   double * qLocal)
{
  int n = 3 * problem->numberOfContacts;
  int in = 3 * contact, it = in + 1, is = it + 1;
  /* qLocal computation*/
  qLocal[0] = problem->q[in];
  qLocal[1] =  problem->q[it];
  qLocal[2] =  problem->q[is];

  if (problem->M->storageType == 0)
  {
    fprintf(stderr, "Numerics, fc3d_nsgs failed.  : %s.\n");
    exit(EXIT_FAILURE);
  }
  else if (problem->M->storageType == 1)
  {
    /* qLocal += rowMB * reaction
     * with rowMB the row of blocks of MGlobal which corresponds
     * to the current contact
     */
    rowProdNoDiagSBM3x3_index_block(n, 3, contact, problem->M->matrix1, reaction, qLocal, index, index_size);
    /* rowProdNoDiagSBM3x3(n, 3, contact, problem->M->matrix1, reaction, qLocal); */
  }

  /* for(unsigned int i =0; i < 3; i++) printf("qLocal[%i]= %e\t", i, qLocal[i]); */
  /* printf("\n"); */

  
}






void fc3d_index_nsgs_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem,
                            double * reaction, SolverOptions* options,
                            unsigned int * index, unsigned int index_size)
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
  fc3d_nsgs_index_computeqLocal(problem, reaction, contact, index, index_size, localproblem->q);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];


}
void fc3d_nsgs_index_initialize_local_solver(SolverPtr* solve, Update_indexPtr* update,
                                             FreeSolverNSGSPtr* freeSolver, ComputeErrorPtr* computeError,
                                             FrictionContactProblem* problem,
                                             FrictionContactProblem* localproblem,
                                             SolverOptions * options, SolverOptions * localsolver_options)
{

  /** Connect to local solver */
  switch (localsolver_options->solverId)
  {
  /* Newton solver (Alart-Curnier) */
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_AC:
  {
    *solve = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
    *update = &fc3d_index_nsgs_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_onecontact_nonsmooth_Newton_solvers_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
    break;
  }
   /* Newton solver (Alart-Curnier) */
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP:
  {
    *solve = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
    *update = &fc3d_index_nsgs_update;
    *freeSolver = (FreeSolverNSGSPtr)&fc3d_onecontact_nonsmooth_Newton_solvers_free;
    *computeError = (ComputeErrorPtr)&fc3d_compute_error;
    fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
    break;
  }
 
  default:
  {
    fprintf(stderr, "Numerics, fc3d_nsgs failed. Unknown internal solver : %s.\n", idToName(localsolver_options->solverId));
    exit(EXIT_FAILURE);
  }
  }
}


void fc3d_nsgs_index(FrictionContactProblem* problem,
                     double *reaction, double *velocity,
                     int* info, SolverOptions* options,
                     unsigned int* index, unsigned int index_size)
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


  
  /* printf("start fc3d_nsgs_index\n \n" ); */
  /* for (int ii=0; ii < index_size; ii++) */
  /* { */
  /*   printf("index[%i] = %i\n", ii,index[ii] ); */
  /* } */
  
  if (*info == 0)
    return;

  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("fc3d_nsgs", "The NSGS method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >= 1");
  }
  assert(options->internalSolvers);

  SolverOptions * localsolver_options = options->internalSolvers;
  
  SolverPtr local_solver = NULL;
  Update_indexPtr update_problem = NULL;
  FreeSolverNSGSPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;

  /* Connect local solver and local problem*/
  FrictionContactProblem* localproblem = (FrictionContactProblem*)malloc(sizeof(FrictionContactProblem));
  localproblem->numberOfContacts = 1;
  localproblem->dimension = 3;
  localproblem->q = (double*)malloc(3 * sizeof(double));
  localproblem->mu = (double*)malloc(sizeof(double));
  /* printf("start of initialization\n" ); */

  //display(problem->M);
  if (problem->M->storageType == NM_DENSE || problem->M->storageType == NM_SPARSE)
  {
    localproblem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3,
                                           malloc(9 * sizeof(double)));
  }
  else /* NM_SPARSE_BLOCK */
  {
    localproblem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3, NULL);
  }

  fc3d_nsgs_index_initialize_local_solver(&local_solver, &update_problem,
                             (FreeSolverNSGSPtr *)&freeSolver, &computeError,
                             problem , localproblem,
                             options, localsolver_options);
  /* printf("end of initialization\n" ); */
  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  unsigned int contact; /* Number of the current row of blocks in M */

  unsigned int *scontacts = NULL;
  /*  dparam[0]= dparam[2]; // set the tolerance for the local solver */
  if (iparam[1] == SICONOS_FRICTION_3D_NSGS_LIGHT_ERROR_EVALUATION_WITH_FULL_FINAL ||
      iparam[1] == SICONOS_FRICTION_3D_NSGS_LIGHT_ERROR_EVALUATION)
  {
    double reactionold[3];
    while ((iter < itermax) && (hasNotConverged > 0))
    {

      ++iter;
      /* Loop through the contact points */
      //cblas_dcopy( n , q , incx , velocity , incy );
      error = 0.0;

      for (unsigned int i= 0 ; i < index_size ; ++i)
      {
        contact = index[i];

        reactionold[0] = reaction[3 * contact];
        reactionold[1] = reaction[3 * contact + 1];
        reactionold[2] = reaction[3 * contact + 2];
        if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);

        (*update_problem)(contact, problem, localproblem,
                                   reaction, localsolver_options,
                                   index, index_size);
        /* frictionContact_display(localproblem); */
        localsolver_options->iparam[4] = contact;
        (*local_solver)(localproblem, &(reaction[3 * contact]) , localsolver_options);

        error += pow(reaction[3 * contact] - reactionold[0], 2) +
                 pow(reaction[3 * contact + 1] - reactionold[1], 2) +
                 pow(reaction[3 * contact + 2] - reactionold[2], 2);

      }

      dparam[2] = error;
      /* **** Criterium convergence **** */
      error = sqrt(error);
      double norm_r = cblas_dnrm2(nc*3 , reaction , 1);
      if (fabs(norm_r) > DBL_EPSILON)
        error /= norm_r;

      if (error < tolerance)
      {
        hasNotConverged = 0;
        if (verbose > 0)
          printf("----------------------------------- FC3D - NSGS INDEX - Iteration %i Residual = %14.7e < %7.3e\n", iter, error, options->dparam[0]);
      }
      else
      {
        if (verbose > 0)
          printf("----------------------------------- FC3D - NSGS INDEX - Iteration %i Residual = %14.7e > %7.3e\n", iter, error, options->dparam[0]);
      }
      *info = hasNotConverged;
    }
    if (iparam[1] == SICONOS_FRICTION_3D_NSGS_LIGHT_ERROR_EVALUATION_WITH_FULL_FINAL) /* Full criterium */
    {
      double absolute_error;
      (*computeError)(problem, reaction , velocity, tolerance, options, normq, &absolute_error);
      if (verbose > 0)
      {
        if (absolute_error > error)
        {
          printf("----------------------------------- FC3D - NSGS INDEX - WARNING absolute Residual = %14.7e is larger than incremental error = %14.7e\n", absolute_error, error);
        }
      }
    }
  }
  else
  {
    numericsError("fc3d_nsgs", "iparam[1] must be equal to SICONOS_FRICTION_3D_NSGS_LIGHT_ERROR_EVALUATION_WITH_FULL_FINAL (1) or SICONOS_FRICTION_3D_NSGS_LIGHT_ERROR_EVALUATION  (2)");
    exit(1);
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

