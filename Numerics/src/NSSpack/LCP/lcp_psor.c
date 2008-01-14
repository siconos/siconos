/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
#include "LA.h"

int lcp_compute_error(int n, double *vec , double *q , double *z , int verbose, double *w, double *err);


void lcp_psor(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP)
{

  int n;
  int incx = 1, incy = 1;
  int incxn;
  int i, iter;
  int itermax, verbose;

  double qs, err, den;
  double tol, omega;
  double *ww, *diag;

  n = *nn;
  incxn = n;
  /* Recup input */

  itermax = iparamLCP[0];
  verbose  = iparamLCP[1];

  tol   = dparamLCP[0];
  omega = dparamLCP[1];

  /* Initialize output */

  iparamLCP[2] = 0;
  dparamLCP[2] = 0.0;

  /* Allocation */

  ww   = (double*)malloc(n * sizeof(double));
  diag = (double*)malloc(n * sizeof(double));

  /* Check for non trivial case */

  qs = DNRM2(n , q , incx);

  if (verbose > 0) printf("\n ||q||= %g \n", qs);

  if (qs > 1e-16) den = 1.0 / qs;
  else
  {
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }

    free(ww);
    free(diag);

    *info = 0;
    return;
  }

  for (i = 0 ; i < n ; ++i)
  {
    ww[i] = 0.;
  }

  DCOPY(n , q , incx , w , incy);
  /* Intialization of w and z */
  /*if(initmethod == 0) {*/
  /* dcopy_( (integer *)&n , q , (integer *)&incx , w , (integer *)&incy );*/
  /*dcopy_( (integer *)&n , ww , (integer *)&incx , z , (integer *)&incy );*/
  /*}    else if(initmethod == 0) { */


  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < n ; ++i)
  {
    if (fabs(vec[i * n + i]) < 1e-16)
    {

      if (verbose > 0)
      {
        printf(" Warning negative diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diag);
      free(ww);

      return;
    }
    else diag[i] = 1.0 / vec[i * n + i];
  }

  /*start iterations*/

  iter = 0;
  err  = 1.;

  qs   = -1.0;

  DCOPY(n , q , incx , w , incy);

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    DCOPY(n , w , incx , ww , incy);   /* w --> ww */
    DCOPY(n , q , incx , w , incy);    /* q --> w */

    for (i = 0 ; i < n ; ++i)
    {

      z[i] = 0.0;

      /*      zi = -( q[i] + ddot_( (integer *)&n , &vec[i] , (integer *)&incxn , z , (integer *)&incy ))*diag[i]; */

      /*       if( zi < 0 ) z[i] = 0.0;  */
      /*       else z[i] = zi; */

      z[i] = fmax(0.0, -(q[i] + DDOT(n , &vec[i] , incxn , z , incy)) * diag[i]);

    }

    /* **** Criterium convergence **** */

    lcp_compute_error(n, vec, q, z, verbose, w, &err);

    /* **** ********************* **** */
  }

  iparamLCP[2] = iter;
  dparamLCP[2] = err;


  if (err > tol)
  {
    printf("Siconos/Numerics: lcp_psor: No convergence of PSOR after %d iterations\n" , iter);
    printf("Siconos/Numerics: lcp_psor: The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if (verbose > 0)
    {
      printf("Siconos/Numerics: lcp_psor: Convergence of PSOR after %d iterations\n" , iter);
      printf("Siconos/Numerics: lcp_psor: The residue is : %g \n", err);
    }
    *info = 0;
  }


  free(ww);
  free(diag);

  return;
}
