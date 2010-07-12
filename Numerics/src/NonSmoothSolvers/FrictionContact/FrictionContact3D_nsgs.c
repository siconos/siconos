/* Siconos-Numerics, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "FrictionContact3D_Newton.h"
#include "FrictionContact3D_Path.h"
#include "FrictionContact3D_NCPGlockerFixedPoint.h"
#include "FrictionContact3D_projection.h"

#include "FrictionContact3D_compute_error.h"
#include "NCP_Solvers.h"
#include "LA.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

void fake_compute_error_nsgs(FrictionContactProblem* problem, double *reaction, double *velocity, double tolerance, SolverOptions  *options,  double* error)
{
  int n = 3 * problem->numberOfContacts;
  *error = 0.;
  int i, m;
  m = 5 * n / 3;
  double err;
  for (i = 0 ; i < m ; ++i)
  {
    *error += Compute_NCP_error1(i, err);
  }
}

void initializeLocalSolver_nsgs(SolverPtr* solve, UpdatePtr* update, FreeSolverPtr* freeSolver, ComputeErrorPtr* computeError, FrictionContactProblem* problem, FrictionContactProblem* localproblem, SolverOptions * localsolver_options)
{


  /** Connect to local solver */
  switch (localsolver_options->solverId)
  {
    /* Projection */
  case SICONOS_FRICTION_3D_ProjectionOnConeWithDiagonalization:
  {
    *solve = &frictionContact3D_projectionWithDiagonalization_solve;
    *update = &frictionContact3D_projectionWithDiagonalization_update;
    *freeSolver = &frictionContact3D_projection_free;
    *computeError = &FrictionContact3D_compute_error;
    frictionContact3D_projection_initialize(problem, localproblem);
    break;
  }
  case SICONOS_FRICTION_3D_ProjectionOnCone:
  {
    *solve = &frictionContact3D_projectionOnCone_solve;
    *update = &frictionContact3D_projection_update;
    *freeSolver = &frictionContact3D_projection_free;
    *computeError = &FrictionContact3D_compute_error;
    frictionContact3D_projection_initialize(problem, localproblem);
    break;
  }
  case SICONOS_FRICTION_3D_ProjectionOnConeWithLocalIteration:
  {
    *solve = &frictionContact3D_projectionOnConeWithLocalIteration_solve;
    *update = &frictionContact3D_projection_update;
    *freeSolver = &frictionContact3D_projection_free;
    *computeError = &FrictionContact3D_compute_error;
    frictionContact3D_projection_initialize(problem, localproblem);
    break;
  }
  case SICONOS_FRICTION_3D_projectionOnConeWithRegularization:
  {
    *solve = &frictionContact3D_projectionOnCone_solve;
    *update = &frictionContact3D_projection_update_with_regularization;
    *freeSolver = &frictionContact3D_projection_with_regularization_free;
    *computeError = &FrictionContact3D_compute_error;
    frictionContact3D_projection_initialize_with_regularization(problem, localproblem);
    break;
  }
  /* Newton solver (Alart-Curnier) */
  case SICONOS_FRICTION_3D_AlartCurnierNewton:
  {
    *solve = &frictionContact3D_Newton_solve;
    *update = &frictionContact3D_AC_update;
    *freeSolver = &frictionContact3D_Newton_free;
    *computeError = &FrictionContact3D_compute_error;
    frictionContact3D_Newton_initialize(problem, localproblem, localsolver_options);
    break;
  }
  /* Newton solver (Glocker-Fischer-Burmeister)*/
  case SICONOS_FRICTION_3D_NCPGlockerFBNewton:
  {
    *solve = &frictionContact3D_Newton_solve;
    *update = &NCPGlocker_update;
    *freeSolver = &frictionContact3D_Newton_free;
    *computeError = &FrictionContact3D_compute_error;
    // *computeError = &fake_compute_error;
    frictionContact3D_Newton_initialize(problem, localproblem, localsolver_options);
    break;
  }
  /* Path solver (Glocker Formulation) */
  case SICONOS_FRICTION_3D_NCPGlockerFBPATH:
  {
    *solve = &frictionContact3D_Path_solve;
    *freeSolver = &frictionContact3D_Path_free;
    *update = &NCPGlocker_update;
    *computeError = &FrictionContact3D_compute_error;
    // *computeError = &fake_compute_error;
    frictionContact3D_Path_initialize(problem, localproblem, localsolver_options);
    break;
  }

  /* Fixed Point solver (Glocker Formulation) */
  case SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint:
  {
    *solve = &frictionContact3D_FixedP_solve;
    *update = &NCPGlocker_update;
    *freeSolver = &frictionContact3D_FixedP_free;
    *computeError = &fake_compute_error_nsgs;
    frictionContact3D_FixedP_initialize(problem, localproblem, localsolver_options);
    break;
  }
  /*iparam[4] > 10 are reserved for Tresca resolution */
  case SICONOS_FRICTION_3D_projectionOnCylinder:
  {
    *solve = &frictionContact3D_projectionOnCylinder_solve;
    *update = &frictionContact3D_projectionOnCylinder_update;
    *freeSolver = &frictionContact3D_projection_free;
    *computeError = &FrictionContact3D_Tresca_compute_error;
    frictionContact3D_projection_initialize(problem, localproblem);
    break;
  }

  default:
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown internal solver : %s.\n", idToName(localsolver_options->solverId));
    exit(EXIT_FAILURE);
  }
  }
}
void frictionContact3D_nsgs_computeqLocal(FrictionContactProblem * problem, FrictionContactProblem * localproblem, double *reaction, int contact)
{

  double *qLocal = localproblem->q;
  int n = 3 * problem->numberOfContacts;


  int in = 3 * contact, it = in + 1, is = it + 1;
  /* reaction current block set to zero, to exclude current contact block */
  double rin = reaction[in] ;
  double rit = reaction[it] ;
  double ris = reaction[is] ;
  reaction[in] = 0.0;
  reaction[it] = 0.0;
  reaction[is] = 0.0;
  /* qLocal computation*/



  qLocal[0] = problem->q[in];
  qLocal[1] =  problem->q[it];
  qLocal[2] =  problem->q[is];

  if (problem->M->storageType == 0)
  {
    double * MM = problem->M->matrix0;
    int incx = n, incy = 1;
    qLocal[0] += DDOT(n , &MM[in] , incx , reaction , incy);
    qLocal[1] += DDOT(n , &MM[it] , incx , reaction , incy);
    qLocal[2] += DDOT(n , &MM[is] , incx , reaction , incy);
  }
  else if (problem->M->storageType == 1)
  {
    /* qLocal += rowMB * reaction
    with rowMB the row of blocks of MGlobal which corresponds to the current contact
    */
    rowProdNoDiagSBM(n, 3, contact, problem->M->matrix1, reaction, qLocal, 0);
  }
  reaction[in] = rin;
  reaction[it] = rit;
  reaction[is] = ris;


}

void frictionContact3D_nsgs_fillMLocal(FrictionContactProblem * problem, FrictionContactProblem * localproblem, int contact)
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
    numericsError("FrictionContact3D_projection -", "unknown storage type for matrix M");
}




void frictionContact3D_nsgs(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];

  if (*info == 0)
    return;

  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("frictionContact3D_nsgs", "The NSGS method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }
  assert(&options[1]);

  SolverOptions * localsolver_options = options->internalSolvers;


  SolverPtr local_solver = NULL;
  UpdatePtr update_localproblem = NULL;
  FreeSolverPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;

  /* Connect local solver and local problem*/
  FrictionContactProblem* localproblem = (FrictionContactProblem*)malloc(sizeof(FrictionContactProblem));
  localproblem->numberOfContacts = 1;
  localproblem->dimension = 3;
  localproblem->isComplete = 1;
  localproblem->q = (double*)malloc(3 * sizeof(double));
  localproblem->mu = (double*)malloc(sizeof(double));

  localproblem->M = (NumericsMatrix*)malloc(sizeof(NumericsMatrix));
  localproblem->M->storageType = 0;
  localproblem->M->size0 = 3;
  localproblem->M->size1 = 3;

  if (problem->M->storageType == 0)
  {
    localproblem->M->matrix0 = (double*)malloc(9 * sizeof(double));
    localproblem->M->matrix1 = NULL;
  }
  else
  {
    localproblem->M->matrix0 = NULL;
    localproblem->M->matrix1 = NULL;
  }

  initializeLocalSolver_nsgs(&local_solver, &update_localproblem, &freeSolver, &computeError, problem , localproblem, localsolver_options);

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */




  /*  dparam[0]= dparam[2]; // set the tolerance for the local solver */

  if (iparam[1] == 1 || iparam[1] == 2)
  {

    while ((iter < itermax) && (hasNotConverged > 0))
    {
      double reactionold[3];
      ++iter;
      /* Loop through the contact points */
      //DCOPY( n , q , incx , velocity , incy );
      error = 0.0;

      for (contact = 0 ; contact < nc ; ++contact)
      {

        reactionold[0] = reaction[3 * contact];
        reactionold[1] = reaction[3 * contact + 1];
        reactionold[2] = reaction[3 * contact + 2];

        (*update_localproblem)(contact, problem, localproblem, reaction, localsolver_options);

        (*local_solver)(localproblem, &(reaction[3 * contact]) , localsolver_options);

        error += pow(reaction[3 * contact] - reactionold[0], 2) +
                 pow(reaction[3 * contact + 1] - reactionold[1], 2) +
                 pow(reaction[3 * contact + 2] - reactionold[2], 2);

      }

      /* **** Criterium convergence **** */
      error = sqrt(error);
      if (verbose > 0)
        printf("----------------------------------- FC3D - NSGS - Iteration %i Error = %14.7e\n", iter, error);
      if (error < tolerance) hasNotConverged = 0;
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
          printf("----------------------------------- FC3D - NSGS - Warning absolute Error = %14.7e is larger than incremental error = %14.7e\n", absolute_error, error);
        }
      }
    }
  }
  else
  {
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      /* Loop through the contact points */
      //DCOPY( n , q , incx , velocity , incy );
      for (contact = 0 ; contact < nc ; ++contact)
      {

        (*update_localproblem)(contact, problem, localproblem, reaction, localsolver_options);

        (*local_solver)(localproblem, &(reaction[3 * contact]), localsolver_options);
      }

      /* **** Criterium convergence **** */
      (*computeError)(problem, reaction , velocity, tolerance, options, &error);

      if (verbose > 0)
        printf("----------------------------------- FC3D - NSGS - Iteration %i Error = %14.7e\n", iter, error);

      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }
  dparam[0] = tolerance;
  dparam[1] = error;


  /***** Free memory *****/
  (*freeSolver)(localproblem);
  localproblem->M->matrix0 = NULL;
  freeFrictionContact_problem(localproblem);

}

int frictionContact3D_nsgs_setDefaultSolverOptions(SolverOptions* options)
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
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  frictionContact3D_AlartCurnierNewton_setDefaultSolverOptions(options->internalSolvers);

  return 0;
}
