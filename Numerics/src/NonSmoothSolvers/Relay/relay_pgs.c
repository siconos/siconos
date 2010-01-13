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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "LA.h"
#include "Relay_Solvers.h"
#include <assert.h>

void relay_pgs(Relay_Problem* problem, double *z, double *w, int *info, Solver_Options* options)
{


  double* M = problem->M->matrix0;
  double* q = problem->q;
  int n = problem -> size;
  double *a = problem->lb;
  double *b = problem->ub;

  assert(M);
  assert(q);
  assert(n);
  assert(a);
  assert(b);



  int itermax = options->iparam[0];
  double tol = options->dparam[0];


  int i;
  double zi;
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
  *info = 1;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    /* Initialization of w with q */
    DCOPY(n , q , 1 , w , 1);

    for (i = 0 ; i < n ; ++i)
    {
      z[i] = 0.0;
      zi = -(q[i] + DDOT(n , &M[i] , n , z , 1)) * diag[i];
      z[i] = zi;
      if (zi < a[i]) z[i] = a[i];
      else if (zi > b[i]) z[i] = b[i];
    }
    /* **** Criterium convergence **** */
    *info = relay_compute_error(problem, z, w, tol, &err);

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
    printf("Siconos/Numerics: relay_pgs: No convergence of PGS after %d iterations\n" , iter);
    printf("The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if (verbose > 0)
    {
      printf("Siconos/Numerics: relay_pgs: Convergence of PGS after %d iterations\n" , iter);
      printf("The residue is : %g \n", err);
    }
    *info = 0;
  }




  free(diag);



}
