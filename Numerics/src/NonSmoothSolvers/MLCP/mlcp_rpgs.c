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
#include "LA.h"
#include "MLCP_Solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#define EPSDIAG DBL_EPSILON

int mixedLinearComplementarity_rpgs_setDefaultSolverOptions(MixedLinearComplementarity_Problem* problem, Solver_Options* pSolver)
{
  pSolver->dparam[2] = 0.5; /*rho*/
  mixedLinearComplementarity_default_setDefaultSolverOptions(problem, pSolver);
  return 0;
}


/*
 *
 * double *z : size n+m
 * double *w : size n+m
 */
void mlcp_rpgs(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options)
{
  double* A = problem->A;
  double* B = problem->B;
  double* C = problem->C;
  double* D = problem->D;
  double* a = problem->a;
  double* b = problem->b;
  int n = problem->n;
  int m = problem->m;
  double *u = &z[0];
  double *v = &z[n];

  int incx, incy, incAx, incAy, incBx, incBy;
  int i, iter;
  int itermax, verbose;
  int incxn;
  double err, vi, viprev, uiprev;
  double tol, rho;
  double *diagA, *diagB;
  verbose = 0;
  incx = 1;
  incy = 1;
  incxn = n;
  /* Recup input */

  itermax = options->iparam[0];
  tol   = options->dparam[0];
  rho   = options->dparam[2];

  /* Initialize output */

  options->iparam[1] = 0;
  options->dparam[1] = 0.0;

  /* Allocation */

  diagA = (double*)malloc(n * sizeof(double));
  diagB = (double*)malloc(m * sizeof(double));



  incx = 1;
  incy = 1;

  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < n ; ++i)
  {
    if (A[i * n + i] < -EPSDIAG)
    {

      if (verbose > 0)
      {
        printf(" Negative diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diagA);
      free(diagB);

      return;
    }
    else
    {
      diagA[i] = 1.0 / (A[i * n + i] + rho);

    }
  }
  for (i = 0 ; i < m ; ++i)
  {
    if (B[i * m + i] < -EPSDIAG)
    {

      if (verbose > 0)
      {
        printf(" Negative diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diagA);
      free(diagB);

      return;
    }
    else
    {
      diagB[i] = 1.0 / (B[i * m + i] + rho);

    }
  }
  /*start iterations*/

  iter = 0;
  err  = 1.;

  incx = 1;
  incy = 1;
  incAx = n;
  incAy = 1;
  incBx = m;
  incBy = 1;


  mlcp_compute_error(problem, z, w, tol, &err);
  printf("Error = %12.8e\n", err);

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    incx = 1;
    incy = 1;


    for (i = 0 ; i < n ; ++i)
    {
      uiprev = u[i];
      u[i] = 0.0;
      //zi = -( q[i] + DDOT( n , &vec[i] , incx , z , incy ))*diag[i];
      //u[i] = -( a[i]  - (rho*uiprev) +DDOT( n , &A[i] , n , u , 1 )   + DDOT( m , &C[i] , n , v , 1 )         )*diagA[i];
      u[i] = -(a[i]   - (rho * uiprev) + DDOT(n , &A[i] , incAx , u , incy)   + DDOT(m , &C[i] , incAx , v , incy)) * diagA[i];
    }

    for (i = 0 ; i < m ; ++i)
    {
      viprev = v[i];
      v[i] = 0.0;
      //zi = -( q[i] + DDOT( n , &vec[i] , incx , z , incy ))*diag[i];
      //v[i] = -( b[i] -(rho*viprev) + DDOT( n , &D[i] , m , u , 1 )   + DDOT( m , &B[i] , m , v , 1 )         )*diagB[i];
      vi = -(b[i] - (rho * viprev) + DDOT(n , &D[i] , incBx , u , incy)   + DDOT(m , &B[i] , incBx , v , incy)) * diagB[i];
      if (vi > 0)
        v[i] = vi;
    }



    /* **** Criterium convergence compliant with filter_result_MLCP **** */
    mlcp_compute_error(problem, z, w, tol, &err);


    if (verbose == 2)
    {
      printf(" # i%d -- %g : ", iter, err);
      for (i = 0 ; i < n ; ++i) printf(" %g", u[i]);
      for (i = 0 ; i < m ; ++i) printf(" %g", v[i]);
      for (i = 0 ; i < m ; ++i) printf(" %g", w[i]);
      printf("\n");
    }

    /* **** ********************* **** */

  }
  options->iparam[1] = iter;
  options->dparam[1] = err;

  if (err > tol)
  {
    printf("Siconos/Numerics: mlcp_rpgs: No convergence of RPGS after %d iterations\n" , iter);
    printf("Siconos/Numerics: mlcp_rpgs: The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if (verbose > 0)
    {
      printf("Siconos/Numerics: mlcp_rpgs: Convergence of RPGS after %d iterations\n" , iter);
      printf("Siconos/Numerics: mlcp_rpgs: The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(diagA);
  free(diagB);
  return;
}
