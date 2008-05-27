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
#include "FrictionContact3D_Newton.h"
#include "FrictionContact3D_projection.h"
#include "FrictionContact3D_Solvers.h"
#include "NCP_Solvers.h"
#include "LA.h"
#include <stdio.h>
#include <stdlib.h>

void fake_compute_error(int n, double* velocity, double* reaction, double* error)
{}

void initializeLocalSolver(int n, SolverPtr* solve, FreeSolverPtr* freeSolver, ComputeErrorPtr* computeError, const double* const M, const double* const q, const double* const mu, int* iparam)
{
  /** Connect to local solver */
  /* Projection */
  if (iparam[4] == 0)
  {
    *solve = &frictionContact3D_projection_solve;
    *freeSolver = &frictionContact3D_projection_free;
    frictionContact3D_projection_initialize(n, M, q, mu);
  }
  /* Newton solver (Alart-Curnier or Fischer-Burmeister)*/
  else if (iparam[4] == 1 || iparam[4] == 2)
  {
    *solve = &frictionContact3D_Newton_solve;
    *freeSolver = &frictionContact3D_Newton_free;
    frictionContact3D_Newton_initialize(n, M, q, mu, iparam);
  }
  /* Path solver */
  else if (iparam[4] == 3)
  {

  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown solver.\n");
    exit(EXIT_FAILURE);
  }
  *computeError = &fake_compute_error;
}

void frictionContact3D_nsgs(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, Solver_Options* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  double* M = problem->M->matrix0;
  double* mu = problem->mu;
  /* Dimension of the problem */
  int n = 3 * nc;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];

  /* Check for trivial case */
  *info = checkTrivialCase(n, q, velocity, reaction, iparam, dparam);

  int incx = 1;
  double  qs = DNRM2(n , q , incx);
  double den = 1. / qs;

  if (*info == 0)
    return;

  SolverPtr local_solver = NULL;
  FreeSolverPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;

  /* Connect local solver */
  initializeLocalSolver(n, &local_solver, &freeSolver, &computeError, M, q, mu, iparam);

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;


  double * W    = (double*)malloc(n * sizeof(double));

  int incy = 1;
  double num;
  int contact; /* Number of the current row of blocks in M */
  while ((iter < itermax) && (hasNotConverged > 0))
  {
    DCOPY(n , velocity , incx , W , incy);
    DCOPY(n , q , incx , velocity , incy);
    ++iter;
    /* Loop through the contact points */
    for (contact = 0 ; contact < nc ; ++contact)
      (*local_solver)(contact, n, reaction, iparam, dparam);

    /* **** Criterium convergence **** */
    //   (*computeError)(n,velocity,reaction,&error);

    //NCP_compute_error( n , M , q , reaction , verbose , velocity, &error );
    //     lcp_compute_error( n , M , q , reaction, verbose, velocity , &error);
    DGEMV(LA_NOTRANS , n , n , 1.0 , M , n , reaction , incx , 1.0 , velocity , incy);
    qs = -1.0;
    DAXPY(n , qs , velocity , incx , W , incy);
    num = DNRM2(n, W , incx);
    error = num * den;
    if (verbose > 0)
      printf("-----------------------------------Iteration %i Erreur = %14.7e\n", iter, error);

    if (error < tolerance) hasNotConverged = 0;

    *info = hasNotConverged;
  }

  /***** Free memory *****/
  (*freeSolver)();
}

