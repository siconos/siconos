/* Siconos-Numerics version 2.0.1, Copyright INRIA 2005-2007.
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
/*!\file FrictionContact3D_nsgs.c
 *
 * This subroutine allows the resolution of contact problems with friction, using a non-smooth Gauss Seidel method.\n
 *
 *   Try \f$(reaction,velocity)\f$ such that:\n
 *   \f$
 *    \left\lbrace
 *     \begin{array}{l}
 *      M reaction + q = velocity\\
 *      0 \le reaction_n \perp velocity_n \ge 0\\
 *      -velocity_t \in \partial\psi_{[-\mu reaction_n, \mu reaction_n]}(reaction_t)\\
 *     \end{array}
 *    \right.
 *   \f$
 *
 * here M is an n by n  matrix, q, reaction and velocity some n-dimensional vectors.
 *
 * \fn  frictionContact3D_nsgs( int nc , double *M , double *q , double *reaction , double *velocity , double *mu,
 *                         int *info\n, int *iparam , double *dparam)
 *
 * with:\n
 *
 * \param nc      Unchanged parameter which represents the number of contacts. The dimension of the system is 3*nc.
 * \param M       Unchanged parameter which contains the components of the matrix with a fortran storage (column major).
 * \param q       Unchanged parameter which contains the components of the right hand side vector.
 * \param reaction       Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param velocity       Modified parameter which returns the solution of the problem.
 * \param mu   the list of friction coefficients. mu[i] corresponds to contact number i.
 * \param info    Modified parameter which returns the termination value\n
 *                0 - convergence\n
 *                1 - iter = itermax, ie the simulation reached the maximum number of iterations allowed\n
 *                2 - negative diagonal term(s) in M.\n
 *
 * and: \n
 *
 * \param iparam[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
 * \param iparam[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
 *                       0 - no output\n
 *                       1 - active screen output\n
 * \param iparam[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
 * \param iparam[3] = local iter_max
 * \param iparam[4] = iter local (output)
 * \param iparam[5] = local formulation (0: Alart-Curnier, 1: Fischer-Burmeister)
 * \param iparam[6] = local solver (0: projection, 1: Newton). Projection only for AC case.
 *
 * \param dparam[0] = tol     Input unchanged parameter which represents the tolerance required.
 * \param dparam[1] = error   Output modified parameter which returns the final error value.
 * \param dparam[2] = local tolerance
 * \param dparam[3] = Output modified parameter which returns the local error
 *
 *
 * \author Houari Khenous and Franck Perignon - Creation: 12/11/2007 - Last modification 03/12/2007.
 *
 * Steps:
 *  1 - Check for trivial case (if norm(q) = 0 -> obvious sol. and return)
 *  2 - If no trivial solution
 *    a - connect to local solver according to iparam[6]
 *    b - iterations until error<tolerance or max iteration number reached
 *             -> loop over through the contact points, call to the local solver.
 *             -> compute error
 */

#include "FrictionContact3D_Newton.h"
#include "FrictionContact3D_projection.h"
#include "FrictionContact3DSolvers.h"
#include <stdio.h>
#include <stdlib.h>

void fake_compute_error(int n, double* velocity, double* reaction, double* error)
{}

void initializeLocalSolver(int n, SolverPtr* solve, FreeSolverPtr* freeSolver, ComputeErrorPtr* computeError, const double* const M, const double* const q, const double* const mu, int* iparam)
{
  /** Connect to local solver */
  /* Projection */
  if (iparam[6] == 0)
  {
    *solve = &frictionContact3D_projection_solve;
    *freeSolver = &frictionContact3D_projection_free;
    frictionContact3D_projection_initialize(n, M, q, mu);
  }
  /* Newton solver */
  else if (iparam[6] == 1)
  {
    *solve = &frictionContact3D_Newton_solve;
    *freeSolver = &frictionContact3D_Newton_free;
    frictionContact3D_Newton_initialize(n, M, q, mu, iparam);
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown solver.\n");
    exit(EXIT_FAILURE);
  }
  *computeError = &fake_compute_error;
}

void frictionContact3D_nsgs(int nc, double *M, double *q, double *reaction, double *velocity, double *mu, int *info, int *iparam, double *dparam)
{
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
  initializeLocalSolver(n, &local_solver, &freeSolver, &computeError, M, q, mu, iparam);

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */
  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    /* Loop through the contact points */
    for (contact = 0 ; contact < nc ; ++contact)
      (*local_solver)(contact, n, reaction, iparam, dparam);

    /* **** Criterium convergence **** */
    //   (*computeError)(n,velocity,reaction,&error);

    NCP_compute_error(n , M , q , reaction , iparam[1] , velocity, &error);

    if (error < tolerance) hasNotConverged = 0;

    *info = hasNotConverged;
  }

  /***** Free memory *****/
  (*freeSolver)();
}

