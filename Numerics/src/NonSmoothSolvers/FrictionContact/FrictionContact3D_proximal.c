/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#include "FrictionContact3D_Solvers.h"
#include "FrictionContact3D_compute_error.h"
#include "NCP_Solvers.h"
#include "SiconosBlas.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"


void frictionContact3D_proximal(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
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
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];


  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("frictionContact3D_proximal", "The PROX method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }
  SolverOptions *internalsolver_options = options->internalSolvers;


  /*****  PROXIMAL Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  double rho = dparam[3];
  if (dparam[3] < 1e-12)
  {
    dparam[3] = 1.0;
  }

  double sigma = 5.0;
  double * reactionold = (double *)malloc(n * sizeof(double));
  cblas_dcopy(n , reaction , 1 , reactionold , 1);


  internalSolverPtr internalsolver =0;
  int iter_iparam =7;
  options->iparam[6]= 0;


  if (internalsolver_options->solverId == SICONOS_FRICTION_3D_NSGS)
  {
    frictionContact3D_nsgs_setDefaultSolverOptions(options->internalSolvers);
    internalsolver = &frictionContact3D_nsgs;
  }
  else if (internalsolver_options->solverId == SICONOS_FRICTION_3D_DeSaxceFixedPoint)
  {
    frictionContact3D_DeSaxceFixedPoint_setDefaultSolverOptions(options->internalSolvers);
    internalsolver = &frictionContact3D_DeSaxceFixedPoint;
  }
  else if (internalsolver_options->solverId == SICONOS_FRICTION_3D_EG)
  {
    frictionContact3D_ExtraGradient_setDefaultSolverOptions(options->internalSolvers);
    internalsolver = &frictionContact3D_ExtraGradient;
  }
  else if (internalsolver_options->solverId == SICONOS_FRICTION_3D_LOCALAC)
  {
    frictionContact3D_AlartCurnier_setDefaultSolverOptions(options->internalSolvers);
    internalsolver = &frictionContact3D_sparseLocalAlartCurnier;
    iter_iparam =1;
    options->internalSolvers->iparam[3]=1000000;
  }
  else
  {
    frictionContact3D_nsgs_setDefaultSolverOptions(options->internalSolvers);
    internalsolver = &frictionContact3D_nsgs;
  }
  options->internalSolvers->iparam[0] = options->iparam[3]+1; // Default Maximum iteration

  FrictionContact3D_compute_error(problem, reaction , velocity, tolerance, options, &error);


  internalsolver_options->dparam[0] = options->dparam[0];
  internalsolver_options->dparam[0] = error;

  rho = sigma*error;


  DEBUG_PRINTF("options->iparam[2] = %i\n",options->iparam[2]);
  DEBUG_PRINTF("options->iparam[3] = %i\n",options->iparam[3]);



  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    cblas_dcopy(n , reaction , 1 , reactionold , 1);
    //Add proximal regularization on q
    cblas_daxpy(n, -rho, reactionold, 1, problem->q , 1) ;

    //Add proximal regularization on M
    if (M->storageType == 0)
    {
      for (int i = 0 ; i < n ; i++) M->matrix0[i + i * n] += rho;
    }
    else if (M->storageType == 1)
    {
      for (int ic = 0 ; ic < nc ; ic++)
      {
        int diagPos = getDiagonalBlockPos(M->matrix1, ic);
        for (int i = 0 ; i < 3 ; i++) M->matrix1->block[diagPos][i + 3 * i] += rho ;
      }
    }
    // internal solver for the regularized problem
    /*       frictionContact3D_nsgs(problem, reaction , velocity , info , internalsolver_options); */

    /* internalsolver_options->dparam[0] = max(error/10.0, options->dparam[0]); */
    //internalsolver_options->dparam[0] = options->dparam[0];

    internalsolver_options->dparam[0] = rho*error;
    DEBUG_PRINTF("internal solver tolerance = %21.8e \n",internalsolver_options->dparam[0]);
    (*internalsolver)(problem, reaction , velocity , info , internalsolver_options);



    /* **** Criterium convergence **** */
    //substract proximal regularization on q
    cblas_daxpy(n, rho, reactionold, 1, problem->q, 1) ;
    //substract proximal regularization on M
    if (M->storageType == 0)
    {
      for (int i = 0 ; i < n ; i++) M->matrix0[i + i * n] -= rho;
    }
    else if (M->storageType == 1)
    {
      for (int ic = 0 ; ic < nc ; ic++)
      {
        int diagPos = getDiagonalBlockPos(M->matrix1, ic);
        for (int i = 0 ; i < 3 ; i++) M->matrix1->block[diagPos][i + 3 * i] -= rho ;
      }
    }

    FrictionContact3D_compute_error(problem, reaction , velocity, tolerance, options, &error);
    //update the rho with respect to the number of internal iterations.

    int iter_internalsolver = internalsolver_options->iparam[iter_iparam];
    options->iparam[6] +=iter_internalsolver;
    DEBUG_PRINTF("iter_internalsolver = %i\n",iter_internalsolver);
    DEBUG_PRINTF("info = %i\n",*info);
    DEBUG_PRINTF("options->iparam[1] = %i\n",options->iparam[1]);
    DEBUG_PRINTF("options->iparam[2] = %i\n",options->iparam[2]);
    DEBUG_PRINTF("options->iparam[3] = %i\n",options->iparam[3]);

    /* if (iter_internalsolver < options->iparam[2])// || (*info == 0))  */
    /* { */
    /*   rho = rho/sigma; */
    /*   DEBUG_PRINTF("We decrease rho = %8.4e\n",rho); */
    /* } */
    /* else if (iter_internalsolver > options->iparam[3]) */
    /* { */
    /*   rho = sigma *rho; */
    /*   DEBUG_PRINTF("We increase rho = %8.4e\n",rho); */
    /* } */

    rho = sigma*error;


    DEBUG_PRINTF("rho = %8.4e\n",rho);

    if (options->callback)
    {
      options->callback->endIteration(options->callback->env, nc * 3,
                                      reaction, velocity,
                                      error);
    }

    if (verbose > 0)
      printf("------------------------ FC3D - PROXIMAL - Iteration %i Error = %14.7e with rho = %12.8e\n\n", iter, error, rho);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }
  if (verbose > 0)
  {
    printf("------------------------ FC3D - PROXIMAL - # Iteration %i Final Error = %14.7e  \n", iter, error);
    printf("------------------------ FC3D - PROXIMAL - # Iteration of internal solver %i \n", options->iparam[6]);
  }

  iparam[7] = iter;
  dparam[0] = tolerance;
  dparam[1] = error;

  free(reactionold);



}


int frictionContact3D_proximal_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the PROX Solver\n");
  }

  /*  strcpy(options->solverName,"PROX");*/
  options->solverId = SICONOS_FRICTION_3D_PROX;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 8;
  options->dSize = 8;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;   options->callback = NULL; options->numericsOptions = NULL;
  for (i = 0; i < 8; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->iparam[2] = 5;   // Default Mimimun iteration of the internal solver for decreasing rho
  options->iparam[3] = 15;  // Default Maximum iteration of the internal solver for increasing rho

  options->dparam[0] = 1e-4;
  options->dparam[3] = 1.e4; // default value for proximal parameter;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  options->internalSolvers->solverId = SICONOS_FRICTION_3D_LOCALAC;
  return 0;
}
