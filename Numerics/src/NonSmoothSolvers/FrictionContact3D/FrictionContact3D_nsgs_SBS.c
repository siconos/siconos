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

#include "FrictionContact3D_Solvers.h"
#include "NCP_Solvers.h"
#include <stdio.h>
#include <stdlib.h>

/* Local function */
void initializeLocalSolver_SBS(int n, SolverPtr* solve, FreeSolverPtr* freeSolver, ComputeErrorPtr* computeError, const SparseBlockStructuredMatrix* const M, const double* const q, const double* const mu, int* iparam)
{
  /** Connect to local solver */
  /* Projection */
  if (iparam[6] == 0)
  {
    *solve = &frictionContact3D_projectionWithDiagonalization_solve;
    *freeSolver = &frictionContact3D_projection_free;
    frictionContact3D_projection_initialize_SBS(n, M, q, mu);
    /*       computeError = & */
  }
  /* Newton solver */
  else if (iparam[6] == 1)
  {
    *solve = &frictionContact3D_Newton_solve;
    *freeSolver = &frictionContact3D_Newton_free;
    frictionContact3D_Newton_initialize_SBS(n, M, q, mu, iparam);
    /*       computeError = & */
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown solver.\n");
    exit(EXIT_FAILURE);
  }
}

/* NSGS solver */
void frictionContact3D_nsgs_SBS(int nc, SparseBlockStructuredMatrix *M, double *q, double *reaction, double *velocity, double *mu, int *info, int *iparam, double *dparam)
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
  initializeLocalSolver_SBS(n, &local_solver, &freeSolver, &computeError, M, q, mu, iparam);

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

    /* Criterium convergence
    computes velocity and error
    */
    NCP_block_compute_error(n , M , q , reaction , iparam[1] , velocity, &error);
    //     (*computeError)(n, velocity, reaction, &error);

    if (error < tolerance) hasNotConverged = 0;

    *info = hasNotConverged;
  }

  /***** Free memory *****/
  (*freeSolver)();
}

