/* Siconos-Numerics version 2.0.1, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include "FrictionContact3D_projection.h"
//#include "PrimalFrictionContact3D_projection.h"
#include "PrimalFrictionContact3D_Solvers.h"
#include "PrimalFrictionContact3D_compute_error.h"
#include "LA.h"
#include <stdio.h>
#include <stdlib.h>


void initializePrimalLocalSolver(int n, SolverPrimalPtr* solve, FreeSolverPrimalPtr* freeSolver, ComputeErrorPrimalPtr* computeError, const NumericsMatrix* const M, const double* const q, const double* const mu, int* iparam)
{
  /** Connect to local solver */
  /* Projection */
  if (iparam[4] == 0)
  {
    *solve = &frictionContact3D_projectionOnCone_solve;
    *freeSolver = &frictionContact3D_projection_free;
    *computeError = &PrimalFrictionContact3D_compute_error;
    frictionContact3D_projection_initialize(n, M, q, mu);
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown solver.\n");
    exit(EXIT_FAILURE);
  }
}

void primalFrictionContact3D_nsgs(PrimalFrictionContact_Problem* problem, double *reaction, double *velocity, double *globalVelocity, int* info, Solver_Options* options)
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

  /* Check for trivial case */
  *info = checkTrivialCasePrimal(n, q, velocity, reaction, globalVelocity, iparam, dparam);

  if (*info == 0)
    return;

  SolverPrimalPtr local_solver = NULL;
  FreeSolverPrimalPtr freeSolver = NULL;
  ComputeErrorPrimalPtr computeError = NULL;

  /* Connect local solver */
  initializePrimalLocalSolver(n, &local_solver, &freeSolver, &computeError, M, q, mu, iparam);

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  int contact; /* Number of the current row of blocks in M */

  dparam[0] = dparam[2]; // set the tolerance for the local solver
  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;



    /* Loop through the contact points */
    //DCOPY( n , q , incx , velocity , incy );
    for (contact = 0 ; contact < nc ; ++contact)
      (*local_solver)(contact, n, reaction, iparam, dparam);

    /* **** Criterium convergence **** */
    (*computeError)(problem, reaction , velocity, globalVelocity, tolerance, &error);

    if (verbose > 0)
      printf("----------------------------------- FC3D - NSGS - Iteration %i Error = %14.7e\n", iter, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }
  dparam[0] = tolerance;
  dparam[1] = error;


  /***** Free memory *****/
  (*freeSolver)();
}

