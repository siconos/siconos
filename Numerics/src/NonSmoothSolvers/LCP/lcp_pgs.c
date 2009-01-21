/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "LA.h"
#include "LCP_Solvers.h"

void lcp_pgs(LinearComplementarity_Problem* problem, double *z, double *w, int *info , Solver_Options* options)
{
  /* matrix M/vector q of the lcp */
  double * M = problem->M->matrix0;
  double * q = problem->q;
  int n = problem->size;

  int incx = 1, incy = 1;
  int i;
  double qs, den, zi;

  /* Solver parameters */
  int itermax = options->iparam[0];
  double tol = options->dparam[0];
  /* Initialize output */
  options->iparam[1] = 0;
  options->dparam[1] = 0.0;

  /* Check for non trivial case */
  qs = DNRM2(n , q , incx);
  if (verbose > 0)
  {
    printf("===== Starting of LCP solving with Projected Gauss Seidel algorithm.\n");
    printf("\n ||q||= %g \n", qs);
  }
  den = 1.0 / qs;

  /* Preparation of the diagonal of the inverse matrix */
  double * diag = (double*)malloc(n * sizeof(double));
  for (i = 0 ; i < n ; ++i)
  {
    if (fabs(M[i * n + i]) < DBL_EPSILON)
    {
      if (verbose > 0)
      {
        printf("Numerics::lcp_pgs, error: vanishing diagonal term \n");
        printf(" The problem cannot be solved with this method \n");
      }

      *info = 2;
      free(diag);
      return;
    }
    else diag[i] = 1.0 / M[i * n + i];
  }

  /* Iterations*/
  int iter = 0;
  double err  = 1.;
  while ((iter < itermax) && (err > tol))
  {

    ++iter;
    incx = 1;
    incy = 1;
    /* Initialization of w with q */
    DCOPY(n , q , incx , w , incy);

    incx = n;
    for (i = 0 ; i < n ; ++i)
    {
      z[i] = 0.0;
      zi = -(q[i] + DDOT(n , &M[i] , incx , z , incy)) * diag[i];
      if (zi < 0) z[i] = 0.0;
      else z[i] = zi;
      /* z[i]=fmax(0.0,-( q[i] + ddot_( (integer *)&n , &M[i] , (integer *)&incxn , z , (integer *)&incy ))*diag[i]);*/
    }
    /* **** Criterium convergence **** */
    lcp_compute_error(problem, z, w, tol, &err);

    if (verbose == 2)
    {
      printf(" # i%d -- %g : ", iter, err);
      for (i = 0 ; i < n ; ++i) printf(" %g", z[i]);
      for (i = 0 ; i < n ; ++i) printf(" %g", w[i]);
      printf("\n");
    }
  }

  options->iparam[1] = iter;
  options->dparam[1] = err;

  if (err > tol)
  {
    printf("Siconos/Numerics: lcp_pgs: No convergence of PGS after %d iterations\n" , iter);
    printf("The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if (verbose > 0)
    {
      printf("Siconos/Numerics: lcp_pgs: Convergence of PGS after %d iterations\n" , iter);
      printf("The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(diag);
}
