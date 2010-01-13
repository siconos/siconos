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
#include "LCP_Solvers.h"

void lcp_cpg(LinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options)
{
  /* matrix M/vector q of the lcp */
  double * M = problem->M->matrix0;

  double * q = problem->q;

  /* size of the LCP */
  int n = problem->size;

  int incx, incy;
  int i, iter;
  int itermax = options->iparam[0];


  double err, a1, b1 , qs;

  double alpha, beta, rp, pMp;
  double den, num;
  double tol = options->dparam[0];

  int *status;
  double *zz , *pp , *rr, *ww, *Mp;

  *info = 1;
  incx  = 1;

  /*output*/

  options->iparam[1] = 0;
  options->dparam[1] = 0.0;

  qs = DNRM2(n , q , incx);

  /*printf( " Norm: %g \n", qs );*/

  den = 1.0 / qs;

  /* Allocations */

  status = (int*)malloc(n * sizeof(int));

  ww = (double*)malloc(n * sizeof(double));
  rr = (double*)malloc(n * sizeof(double));
  pp = (double*)malloc(n * sizeof(double));
  zz = (double*)malloc(n * sizeof(double));

  Mp = (double*)malloc(n * sizeof(double));

  incx = 1;

  for (i = 0; i < n; ++i)
  {

    status[i] = 0;

    ww[i] = 0.;
    rr[i] = 0.;
    pp[i] = 0.;
    zz[i] = 0.;

    Mp[i] = 0.;

  }

  /* rr = -Wz + q */

  incx = 1;
  incy = 1;

  DCOPY(n, q, incx, rr, incy);

  a1 = -1.;
  b1 = -1.;

  DGEMV(LA_NOTRANS, n, n, a1, M, n, z, incx, b1, rr, incy);

  /* Initialization of gradients */
  /* rr -> p and rr -> w */

  DCOPY(n, rr, incx, ww, incy);
  DCOPY(n, rr, incx, pp, incy);

  iter = 0.0;
  err  = 1.0 ;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    /* Compute initial pMp */

    incx = 1;
    incy = 1;

    DCOPY(n, pp, incx, Mp, incy);

    a1 = 1.0;
    b1 = 0.0;

    DGEMV(LA_NOTRANS, n, n, a1, M, n, Mp, incx, b1, w, incy);

    pMp = DDOT(n, pp, incx, w, incy);

    if (fabs(pMp) < DBL_EPSILON)
    {

      if (verbose > 0)
      {
        printf(" Operation no conform at the iteration %d \n", iter);
        printf(" Alpha can be obtained with pWp = %10.4g  \n", pMp);
      }

      free(Mp);
      free(ww);
      free(rr);
      free(pp);
      free(zz);
      free(status);

      options->iparam[1] = iter;
      options->dparam[1] = err;
      *info = 3;
      return;
    }

    rp  = DDOT(n, pp, incx, rr, incy);

    alpha = rp / pMp;

    /*
     * Iterate prediction:
     * z' = z + alpha*p
     *
     */

    DAXPY(n, alpha, pp, incx, z, incy);

    /* Iterate projection*/

    for (i = 0; i < n; ++i)
    {
      if (z[i] > 0.0)
      {
        status[i] = 1;
      }
      else
      {
        z[i] = 0.0;
        status[i] = 0;
      }
    }

    /* rr = -Wz + q */

    DCOPY(n, rr, incx, w , incy);
    DCOPY(n, q , incx, rr, incy);

    a1 = -1.;
    b1 = -1.;

    DGEMV(LA_NOTRANS, n, n, a1, M, n, z, incx, b1, rr, incy);

    /* Gradients projection
     * rr --> ww
     * pp --> zz
     */

    for (i = 0; i < n; ++i)
    {

      if (status[i])
      {
        ww[i] = rr[i];
        zz[i] = pp[i];
      }
      else
      {
        if (rr[i] < 0)
        {
          ww[i] = 0.0;
          zz[i] = 0.0;
        }
        else
        {
          ww[i] = rr[i];
          if (pp[i] < 0) zz[i] = 0.0;
          else zz[i] = pp[i];
        }
      }
    }

    /*   beta = -w.Mp / pMp  */

    rp = DDOT(n , ww, incx, w, incy);

    beta = -rp / pMp;

    DCOPY(n, ww, incx, pp, incy);
    DAXPY(n, beta, zz, incx, pp, incy);

    /* **** Criterium convergence **** */

    qs   = -1.0;
    DAXPY(n, qs, rr, incx, w, incy);
    num = DNRM2(n, w, incx);
    err = num * den;

  }

  options->iparam[1] = iter;
  options->dparam[1] = err;

  DCOPY(n, rr, incx, w, incy);

  qs   = -1.0;
  DSCAL(n, qs, w, incx);


  *info = 1;
  if (verbose > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
    }
    else
    {
      printf(" Convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else if (err <= tol) *info = 0;

  free(Mp);
  free(status);

  free(ww);
  free(rr);
  free(pp);
  free(zz);

}
