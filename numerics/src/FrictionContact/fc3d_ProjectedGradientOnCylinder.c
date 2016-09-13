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

#include "projectionOnCylinder.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Friction_cst.h"
#include "misc.h"

//#define VERBOSE_DEBUG

void fc3d_ProjectedGradientOnCylinder(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  NumericsMatrix* M = problem->M;
  /* Dimension of the problem */
  int n = 3 * nc;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];




  /*****  Projected Gradient iterations *****/
  int j, iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */
  int nLocal = 3;
  dparam[0] = dparam[2]; // set the tolerance for the local solver
  double * velocitytmp = (double *)malloc(n * sizeof(double));

  double rho = 0.0;
  int isVariable = 0;
  double rhoinit, rhomin;
  if (dparam[3] > 0.0)
  {
    rho = dparam[3];
  }
  else
  {
    /* Variable step in fixed*/
    isVariable = 1;
    printf("Variable step (line search) in Projected Gradient iterations\n");
    rhoinit = dparam[3];
    rhomin = dparam[4];
  }

  double * reactionold;
  double * direction;
  if (isVariable)
  {
    reactionold = (double *)malloc(n * sizeof(double));
    direction = (double *)malloc(n * sizeof(double));
  }
  double alpha = 1.0;
  double beta = 1.0;

  /*   double minusrho  = -1.0*rho; */

  if (!isVariable)
  {
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      cblas_dcopy(n , q , 1 , velocitytmp, 1);
      prodNumericsMatrix(n, n, alpha, M, reaction, beta, velocitytmp);
      // projection for each contact
      cblas_daxpy(n, -1.0, velocitytmp, 1, reaction , 1);
      for (contact = 0 ; contact < nc ; ++contact)
        projectionOnCylinder(&reaction[ contact * nLocal],
                             options->dWork[contact]);

#ifdef VERBOSE_DEBUG

      printf("reaction before LS\n");
      for (contact = 0 ; contact < nc ; ++contact)
      {
        for (j = 0; j < 3; j++)
          printf("reaction[%i] = %le\t", 3 * contact + j, reaction[3 * contact + j]);
        printf("\n");
      }
      printf("velocitytmp before LS\n");
      for (contact = 0 ; contact < nc ; ++contact)
      {
        for (j = 0; j < 3; j++)
          printf("velocitytmp[%i] = %le\t", 3 * contact + j, velocitytmp[3 * contact + j]);
        printf("\n");
      }
#endif
      /* **** Criterium convergence **** */
      fc3d_Tresca_compute_error(problem, reaction , velocity, tolerance, options, &error);

      if (options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, nc * 3, 
                                        reaction, velocity, 
                                        error, NULL);
      }

      if (verbose > 0)
        printf("----------------------------------- FC3D - Projected Gradient On Cylinder (PGoC) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);

      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }
  else
  {
    rho =  rhoinit;


    cblas_dcopy(n , q , 1 , velocitytmp, 1);
    prodNumericsMatrix(n, n, 1.0, M, reaction, 1.0, velocitytmp);

    cblas_daxpy(n, rho, velocitytmp, 1, reaction, 1);

    for (contact = 0 ; contact < nc ; ++contact)
      projectionOnCylinder(&reaction[contact * nLocal],
                           options->dWork[contact]);
    cblas_dcopy(n , q , 1 , velocitytmp, 1);
    prodNumericsMatrix(n, n, 1.0, M, reaction, 1.0, velocitytmp);

    double oldcriterion = cblas_ddot(n, reaction, 1, velocitytmp, 1);
#ifdef VERBOSE_DEBUG
    printf("oldcriterion =%le \n", oldcriterion);
#endif


    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      // store the old reaction
      cblas_dcopy(n , reaction , 1 , reactionold , 1);
      // compute the direction
      cblas_dcopy(n , q , 1 , velocitytmp, 1);
      prodNumericsMatrix(n, n, 1.0, M, reaction, 1.0, velocitytmp);
      cblas_dcopy(n, velocitytmp, 1, direction, 1);

      // start line search
      j = 0;

      if (rho <= 100 * rhoinit) rho = 10.0 * rho;

      double newcriterion = 1e24;
      do
      {





        cblas_dcopy(n , reactionold , 1 , reaction , 1);
        cblas_daxpy(n, rho, direction, 1, reaction , 1) ;
#ifdef VERBOSE_DEBUG
        printf("LS iteration %i step 0 \n", j);
        printf("rho = %le \n", rho);
        for (contact = 0 ; contact < nc ; ++contact)
        {
          for (int k = 0; k < 3; k++)
            printf("reaction[%i] = %le\t",
                   3 * contact + k, reaction[3 * contact + k]);
          printf("\n");
        }
#endif
        for (contact = 0 ; contact < nc ; ++contact)
          projectionOnCylinder(&reaction[contact * nLocal],
                               options->dWork[contact]);
        /*          printf("options->dWork[%i] = %le\n",contact, options->dWork[contact]  );} */
#ifdef VERBOSE_DEBUG
        printf("LS iteration %i step 1 after projection\n", j);
        for (contact = 0 ; contact < nc ; ++contact)
        {
          for (int k = 0; k < 3; k++)
            printf("reaction[%i] = %le\t",
                   3 * contact + k, reaction[3 * contact + k]);
          printf("\n");
        }
#endif
        cblas_dcopy(n , q , 1 , velocitytmp, 1);
        prodNumericsMatrix(n, n, 1.0, M, reaction, 1.0, velocitytmp);

#ifdef VERBOSE_DEBUG
        printf("LS iteration %i step 3 \n", j);
        for (contact = 0 ; contact < nc ; ++contact)
        {
          for (int k = 0; k < 3; k++)
            printf("velocitytmp[%i] = %le\t", 3 * contact + k, velocitytmp[3 * contact + k]);
          printf("\n");
        }
#endif


        newcriterion = cblas_ddot(n, reaction, 1, velocitytmp, 1);

#ifdef VERBOSE_DEBUG
        printf("LS iteration %i newcriterion =%le\n", j, newcriterion);
#endif
        if (rho > rhomin)
        {
          rho = rhomin;
          break;
        }

        rho = 0.5 * rho;
      }
      while (newcriterion > oldcriterion &&
             ++j <= options->iparam[2]);
      oldcriterion = newcriterion;

      /* **** Criterium convergence **** */
      fc3d_Tresca_compute_error(problem, reaction , velocity, tolerance, options, &error);

      if (verbose > 0)
        printf("----------------------------------- FC3D - Projected Gradient On Cylinder (PGoC) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);

      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }





  printf("----------------------------------- FC3D - Projected Gradient On Cylinder (PGoC)- #Iteration %i Final Residual = %14.7e\n", iter, error);
  dparam[0] = tolerance;
  dparam[1] = error;
  free(velocitytmp);
  if (isVariable)
  {
    free(reactionold);
    free(direction);
  }
}


int fc3d_ProjectedGradientOnCylinder_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the PGoC Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_PGoC;

  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 10;
  options->dSize = 10;
  options->iparam = (int *)realloc(options->iparam, options->iSize*sizeof(int));
  options->dparam = (double *)realloc(options->dparam, options->dSize*sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  memset(options->iparam, 0, options->iSize*sizeof(int));
  memset(options->dparam, 0, options->dSize*sizeof(int));
  options->iparam[0] = 20000;
  options->dparam[0] = 1e-6;
  options->dparam[3] = 1e-3;

  if (options->internalSolvers != NULL)
  {
    solver_options_delete(options->internalSolvers);
    free(options->internalSolvers);
  }

  options->internalSolvers = NULL;


  return 0;
}
