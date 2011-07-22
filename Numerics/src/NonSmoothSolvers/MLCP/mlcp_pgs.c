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

#include "MLCP_Solvers.h"
#include "LA.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
/*
 *
 * double *z : size n+m
 * double *w : size n+m
 */

int mixedLinearComplementarity_pgs_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver)
{

  mixedLinearComplementarity_default_setDefaultSolverOptions(problem, pSolver);
  pSolver->iparam[2] = 0; //implicit
  return 0;
}
void mlcp_pgs(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
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
  double *Buf;

  int incx, incy, incAx, incAy, incBx, incBy;
  int i, iter;
  int itermax, verbose;
  int incxn;
  int pgsExplicit;
  double err, vi;
  double tol;
  double prev;
  double *diagA, *diagB;
  verbose = 0;

  incx = 1;
  incy = 1;
  incxn = n;
  /* Recup input */

  itermax = options->iparam[0];
  pgsExplicit = options->iparam[2];
  tol   = options->dparam[0];

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
    if ((fabs(A[i * n + i]) < DBL_EPSILON))
    {

      if (verbose > 0)
      {
        printf(" Vanishing diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diagA);
      free(diagB);
      *info = 1;
      return;
    }
    else
    {
      diagA[i] = 1.0 / A[i * n + i];

    }
  }
  for (i = 0 ; i < m ; ++i)
  {
    if ((fabs(B[i * m + i]) < DBL_EPSILON))
    {

      if (verbose > 0)
      {
        printf(" Vanishing diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diagA);
      free(diagB);

      return;
    }
    else
    {
      diagB[i] = 1.0 / B[i * m + i];

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

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    incx = 1;
    incy = 1;

    if (pgsExplicit)
    {
      /*Use w like a buffer*/
      DCOPY(n , w , incx , u , incy);  //w <- q
      Buf = w;

      for (i = 0 ; i < n ; ++i)
      {
        prev = Buf[i];
        Buf[i] = 0;
        //zi = -( q[i] + DDOT( n , &vec[i] , incx , z , incy ))*diag[i];
        u[i] =  - (a[i] + DDOT(n , &A[i] , incAx , Buf , incAy)   + DDOT(m , &C[i] , incAx , v , incBy)) * diagA[i];
        Buf[i] = prev;
      }
      for (i = 0 ; i < m ; ++i)
      {
        v[i] = 0.0;
        //zi = -( q[i] + DDOT( n , &vec[i] , incx , z , incy ))*diag[i];
        vi = -(b[i] + DDOT(n , &D[i] , incBx , u , incAy)   + DDOT(m , &B[i] , incBx , v , incBy)) * diagB[i];

        if (vi < 0) v[i] = 0.0;
        else v[i] = vi;
      }
    }
    else
    {

      for (i = 0 ; i < n ; ++i)
      {
        u[i] = 0.0;

        //zi = -( q[i] + DDOT( n , &vec[i] , incx , z , incy ))*diag[i];
        u[i] =  - (a[i] + DDOT(n , &A[i] , incAx , u , incAy)   + DDOT(m , &C[i] , incAx , v , incBy)) * diagA[i];
      }

      for (i = 0 ; i < m ; ++i)
      {
        v[i] = 0.0;
        //zi = -( q[i] + DDOT( n , &vec[i] , incx , z , incy ))*diag[i];
        vi = -(b[i] + DDOT(n , &D[i] , incBx , u , incAy)   + DDOT(m , &B[i] , incBx , v , incBy)) * diagB[i];

        if (vi < 0) v[i] = 0.0;
        else v[i] = vi;
      }
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
    printf("Siconos/Numerics: mlcp_pgs: No convergence of PGS after %d iterations\n" , iter);
    printf("Siconos/Numerics: mlcp_pgs: The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if (verbose > 0)
    {
      printf("Siconos/Numerics: mlcp_pgs: Convergence of PGS after %d iterations\n" , iter);
      printf("Siconos/Numerics: mlcp_pgs: The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(diagA);
  free(diagB);
  return;
}
