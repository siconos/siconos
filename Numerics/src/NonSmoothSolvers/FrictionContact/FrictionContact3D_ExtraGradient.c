/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
#include "projectionOnCone.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContact3D_compute_error.h"
#include "LA.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void frictionContact3D_ExtraGradient(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  NumericsMatrix* M = problem->M;
  double* mu = problem->mu;
  /* Dimension of the problem */
  int n = 3 * nc;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];





  /*****  Fixed point iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */
  int nLocal = 3;
  dparam[0] = dparam[2]; // set the tolerance for the local solver
  double * velocitytmp = (double *)malloc(n * sizeof(double));
  double * reactiontmp = (double *)malloc(n * sizeof(double));

  double rho = 0.0;

  if (dparam[3] > 0.0)
  {
    rho = dparam[3];
  }
  else
  {
    printf("Extragradient method. rho <=0  is not well defined\n");

  }

  double alpha = 1.0;
  double beta = 1.0;

  /*   double minusrho  = -1.0*rho; */
  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;

    DCOPY(n , q , 1 , velocitytmp, 1);
    DCOPY(n , reaction , 1 , reactiontmp, 1);

    prodNumericsMatrix(n, n, alpha, M, reactiontmp, beta, velocitytmp);
    // projection for each contact
    for (contact = 0 ; contact < nc ; ++contact)
    {
      int pos = contact * nLocal;
      double  normUT = sqrt(velocitytmp[pos + 1] * velocitytmp[pos + 1] + velocitytmp[pos + 2] * velocitytmp[pos + 2]);
      reactiontmp[pos] -= rho * (velocitytmp[pos] + mu[contact] * normUT);
      reactiontmp[pos + 1] -= rho * velocitytmp[pos + 1];
      reactiontmp[pos + 2] -= rho * velocitytmp[pos + 2];
      projectionOnCone(&reactiontmp[pos], mu[contact]);
    }
    DCOPY(n , q , 1 , velocitytmp, 1);
    prodNumericsMatrix(n, n, alpha, M, reactiontmp, beta, velocitytmp);
    // projection for each contact
    for (contact = 0 ; contact < nc ; ++contact)
    {
      int pos = contact * nLocal;
      double  normUT = sqrt(velocitytmp[pos + 1] * velocitytmp[pos + 1] + velocitytmp[pos + 2] * velocitytmp[pos + 2]);
      reaction[pos] -= rho * (velocitytmp[pos] + mu[contact] * normUT);
      reaction[pos + 1] -= rho * velocitytmp[pos + 1];
      reaction[pos + 2] -= rho * velocitytmp[pos + 2];
      projectionOnCone(&reaction[pos], mu[contact]);
    }

    /* **** Criterium convergence **** */
    FrictionContact3D_compute_error(problem, reaction , velocity, tolerance, options, &error);

    if (verbose > 0)
      printf("----------------------------------- FC3D - Extra Gradient (EG) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }
  printf("----------------------------------- FC3D -  Extra Gradient (EG) - #Iteration %i Final Error = %14.7e\n", iter, error);
  dparam[0] = tolerance;
  dparam[1] = error;
  free(velocitytmp);
  free(reactiontmp);

}


int frictionContact3D_ExtraGradient_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the ExtraGradient Solver\n");
  }

  /*strcpy(options->solverName,"DSFP");*/
  options->solverId = SICONOS_FRICTION_3D_EG;
  options->numberOfInternalSolvers = 0;
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
  options->iparam[0] = 20000;
  options->dparam[0] = 1e-3;
  options->dparam[3] = 1e-3;

  options->internalSolvers = NULL;

  return 0;
}
