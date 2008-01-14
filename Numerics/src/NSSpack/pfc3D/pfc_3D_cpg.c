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
#include "pfc_3D_Solvers.h"

void pfc_3D_cpg(int nc , double *vec , double *q , double *z , double *w , double *mu, int *info,
                int *iparamLCP , double *dparamLCP)
{

  FILE *f101;
  int n;
  int incx, incy;
  int i, iter, itermax, ispeak;

  double err, a1, b1, qs;

  double alpha, beta, rp, pMp;
  double den, num, tol;

  int *status;
  double *zz , *pp , *rr , *ww , *Mp;

  printf("pfc_3D_* Algorithms are not reliable yet, report to siconos.gforge.inria.fr if you need it soon \n");
  return ;

  incx = 1;
  incy = 1;
  n    = 3 * nc;

  /* Recup input */

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  tol = dparamLCP[0];

  /* Initialize output */

  iparamLCP[2] = 0;
  dparamLCP[1] = 0.0;

  if (ispeak == 2) f101 = fopen("pfc_3D_cpg.log" , "w+");

  qs = DNRM2(n , q , incx);

  if (ispeak > 0) printf("\n ||q||= %g \n" , qs);

  if (qs > 1e-16) den = 1.0 / qs;
  else
  {
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }
    *info = 9;
    if (ispeak == 2) fclose(f101);
    return;
  }

  /* Allocations */

  status = (int*)malloc(nc * sizeof(int));

  ww = (double*)malloc(n * sizeof(double));
  rr = (double*)malloc(n * sizeof(double));
  pp = (double*)malloc(n * sizeof(double));
  zz = (double*)malloc(n * sizeof(double));

  Mp = (double*)malloc(n * sizeof(double));

  incx = 1;

  for (i = 0; i < n; ++i)
  {

    ww[i] = 0.;
    rr[i] = 0.;
    pp[i] = 0.;
    zz[i] = 0.;

    Mp[i] = 0.;

  }

  for (i = 0; i < nc ; ++i) status[i] = 0;


  /* rr = -Wz + q */

  incy = 1;

  DCOPY(n , q , incx , rr , incy);

  a1 = -1.;
  b1 = -1.;

  DGEMV(LA_NOTRANS , n , n , a1 , vec , n , z , incx , b1 , rr , incy);

  /* Initialization of gradients */
  /* rr -> p and rr -> w */

  DCOPY(n , rr , incx , ww , incy);
  DCOPY(n , rr , incx , pp , incy);

  iter = 0;
  err  = 1.0 ;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    /* Compute initial pMp */

    incx = 1;
    incy = 1;

    DCOPY(n , pp , incx , Mp , incy);

    a1 = 1.0;
    b1 = 0.0;

    DGEMV(LA_NOTRANS , n , n , a1 , vec , n , Mp , incx , b1 , w , incy);

    pMp = DDOT(n , pp , incx , w  , incy);
    //printf( " pWp = %10.4g  \n", pMp );
    if (fabs(pMp) < 1e-16)
    {

      if (ispeak > 0)
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

      iparamLCP[2] = iter;
      dparamLCP[1] = err;

      if (ispeak == 2) fclose(f101);
      *info = 3;
      return;
    }

    rp  = DDOT(n , pp , incx , rr , incy);
    //printf( " rp = %10.4g  \n", rp );
    alpha = rp / pMp;

    /*
     * Iterate prediction:
     * z' = z + alpha*p
     *
     */

    DAXPY(n , alpha , pp , incx , z , incy);

    /* Iterate projection*/

    pfc_3D_projc(nc , mu , z , pp , status);

    /* rr = -Wz + q */

    DCOPY(n , rr , incx , w  , incy);
    DCOPY(n , q  , incx , rr , incy);

    a1 = -1.;
    b1 = -1.;

    DGEMV(LA_NOTRANS , n , n , a1 , vec , n , z , incx , b1 , rr , incy);

    /* Gradients projection
     * rr --> ww
     * pp --> zz
     */

    pfc_3D_projf(nc , ww , zz , rr , pp , status);

    /*   beta = -w.Mp / pMp  */

    rp = DDOT(n , ww , incx, w , incy);

    beta = -rp / pMp;

    DCOPY(n , ww , incx , pp , incy);
    DAXPY(n, beta , zz , incx , pp , incy);

    /* **** Criterium convergence **** */

    qs   = -1.0;
    DAXPY(n , qs , rr , incx , w , incy);

    /*
     * Note: the function dnrm2_ leads to the generation of core
     */
    num = sqrt(DDOT(n , w , incx, w , incy));

    err = num * den;

    if (ispeak == 2) for (i = 0 ; i < n ; ++i) fprintf(f101, "%d  %d  %14.7e\n", iter , i , z[i]);

  }

  iparamLCP[2] = iter;
  dparamLCP[1] = err;

  DCOPY(n , rr , incx , w , incy);

  qs   = -1.0;
  DSCAL(n , qs , w , incx);

  if (ispeak > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 1;
    }
    else
    {
      printf(" Convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else
  {
    if (err > tol) *info = 1;
    else *info = 0;
  }


  free(Mp);

  free(ww);
  free(rr);
  free(pp);
  free(zz);
  free(status);

  if (ispeak == 2) fclose(f101);

}
