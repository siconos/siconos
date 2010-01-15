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
#include "FrictionContact3D_FixedP.h"
#include "FrictionContact3D_projection.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContact3D_compute_error.h"
#include "NCP_Solvers.h"
#include "LA.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void fake_compute_error_nsgs(FrictionContact_Problem* problem, double *reaction, double *velocity, double tolerance, double* error)
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

void initializeLocalSolver_nsgs(int n, SolverPtr* solve, FreeSolverPtr* freeSolver, ComputeErrorPtr* computeError, const NumericsMatrix* const M, const double* const q, const double* const mu, int* iparam)
{
  /** Connect to local solver */
  /* Projection */
  if (iparam[4] == 8)
  {
    *solve = &frictionContact3D_projectionWithDiagonalization_solve;
    *freeSolver = &frictionContact3D_projection_free;
    *computeError = &FrictionContact3D_compute_error;
    frictionContact3D_projection_initialize(n, M, q, mu);
  }
  else  if (iparam[4] == 3)
  {
    *solve = &frictionContact3D_projectionWithDiagonalization_solve;
    *freeSolver = &frictionContact3D_projection_free;
    *computeError = &FrictionContact3D_compute_error;
    frictionContact3D_projection_initialize(n, M, q, mu);
  }
  else if (iparam[4] == 0)
  {
    *solve = &frictionContact3D_projectionOnCone_solve;
    *freeSolver = &frictionContact3D_projection_free;
    *computeError = &FrictionContact3D_compute_error;
    frictionContact3D_projection_initialize(n, M, q, mu);
  }
  else if (iparam[4] == 2)
  {
    *solve = &frictionContact3D_projectionOnConeWithLocalIteration_solve;
    *freeSolver = &frictionContact3D_projection_free;
    *computeError = &FrictionContact3D_compute_error;
    frictionContact3D_projection_initialize(n, M, q, mu);
  }
  /* Newton solver (Alart-Curnier) */
  else if (iparam[4] == 1)
  {
    *solve = &frictionContact3D_Newton_solve;
    *freeSolver = &frictionContact3D_Newton_free;
    *computeError = &FrictionContact3D_compute_error;
    frictionContact3D_Newton_initialize(n, M, q, mu, iparam);
  }
  /* Newton solver (Glocker-Fischer-Burmeister)*/
  else if (iparam[4] == 6)
  {
    *solve = &frictionContact3D_Newton_solve;
    *freeSolver = &frictionContact3D_Newton_free;
    *computeError = &FrictionContact3D_compute_error;
    // *computeError = &fake_compute_error;
    frictionContact3D_Newton_initialize(n, M, q, mu, iparam);
  }
  /* Path solver (Glocker Formulation) */
  else if (iparam[4] == 4)
  {
    *solve = &frictionContact3D_Path_solve;
    *freeSolver = &frictionContact3D_Path_free;
    *computeError = &FrictionContact3D_compute_error;
    // *computeError = &fake_compute_error;
    frictionContact3D_Path_initialize(n, M, q, mu, iparam);
  }
  /* Fixed Point solver (Glocker Formulation) */
  else if (iparam[4] == 5)
  {
    *solve = &frictionContact3D_FixedP_solve;
    *freeSolver = &frictionContact3D_FixedP_free;
    *computeError = &fake_compute_error_nsgs;
    frictionContact3D_FixedP_initialize(n, M, q, mu, iparam);
  }
  else if (iparam[4] == 7)
  {
    *solve = &frictionContact3D_projectionOnCone_with_regularization_solve;
    *freeSolver = &frictionContact3D_projection_free;
    *computeError = &FrictionContact3D_compute_error;
    frictionContact3D_projection_initialize_with_regularization(n, M, q, mu);
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown solver number %i.\n", iparam[4]);
    exit(EXIT_FAILURE);
  }
}

void frictionContact3D_nsgs(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, Solver_Options* options)
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
  *info = checkTrivialCase(n, q, velocity, reaction, iparam, dparam);

  if (*info == 0)
    return;

  SolverPtr local_solver = NULL;
  FreeSolverPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;

  /* Connect local solver */
  initializeLocalSolver_nsgs(n, &local_solver, &freeSolver, &computeError, M, q, mu, iparam);

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */

  double localdparam[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  int localiparam[5] = {0, 0, 0, 0, 0};

  localiparam[0] = itermax / 10;
  localdparam[0] = dparam[2];


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
        (*local_solver)(contact, n, reaction, localiparam, localdparam);
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
      (*computeError)(problem, reaction , velocity, tolerance, &absolute_error);
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
        (*local_solver)(contact, n, reaction, iparam, dparam);

      /* **** Criterium convergence **** */
      (*computeError)(problem, reaction , velocity, tolerance, &error);

      if (verbose > 0)
        printf("----------------------------------- FC3D - NSGS - Iteration %i Error = %14.7e\n", iter, error);

      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }
  dparam[0] = tolerance;
  dparam[1] = error;


  /***** Free memory *****/
  (*freeSolver)();
}

