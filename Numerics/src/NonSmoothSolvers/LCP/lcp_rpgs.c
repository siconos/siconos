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

#include "LA.h"
#include <math.h>
#include "LCP_Solvers.h"

#define EPSDIAG 1e-16

void lcp_rpgs(int *nn , double *M , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP)
{

  int n, incx, incy;
  int i, iter;
  int itermax, ispeak;
  double qs, err, den, zi;
  double tol, omega, rho;
  double *diag;
  double Mii, ziprev;

  /*  double *buffer_errors;
    FILE *ficbuffer_errors;*/

  n = *nn;
  /* Recup input */

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  /*  buffer_errors = malloc( itermax*sizeof( double ) );*/

  tol   = dparamLCP[0];
  rho   = dparamLCP[1];
  omega = dparamLCP[2];

  /* Initialize output */

  iparamLCP[2] = 0;
  dparamLCP[3] = 0.0;

  /* Allocation */

  /*  ww   = ( double* )malloc( n*sizeof( double ) );*/
  /*  zprev = ( double* )malloc( n*sizeof( double ) );*/
  diag = (double*)malloc(n * sizeof(double));
  /*  diagprev = ( int* )malloc( n*sizeof( int ) );*/

  /*  qs = 0.;*/
  incx = 1;
  qs = DNRM2(n , q , incx);
  if (ispeak > 0) printf("\n ||q||= %g \n", qs);
  den = 1.0 / qs;

  /* Initialization of z & w */
  /*
    for( i = 0 ; i < n ; i++ ) z[i] = 0;

    incx = 1;
    incy = 1;
    dcopy_( (integer *)&n , q , (integer *)&incx , w , (integer *)&incy );
  */
  /* Preparation of the diagonal of the inverse matrix */
  /*  maxdiag = 0.0;
    for( i = 0 ; i < n ; i++ ){
      if (M[i*(n+1)] > maxdiag) maxdiag = M[i*(n+1)];
    }

    if (maxdiag > EPSDIAG) rho = maxdiag;
    invrho = 1.0/rho;*/

  for (i = 0 ; i < n ; ++i)
  {
    Mii = M[i * (n + 1)];
    /* if(abs(Mii+rho)<EPSDIAG){   */ /* Version of Pascal */
    if (Mii < -EPSDIAG)
    {
      if (ispeak > 0)
      {
        printf(" RPGS : Warning negative diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diag);
      /*      free(zprev);*/
      /*      free(ww);*/
      /*      free(diagprev);*/
      /*      free(buffer_errors);*/
      return;
    }
    else
    {
      diag[i] = 1.0 / (Mii + rho);
      /*        qs += pow(q[i]*diag[i],2);*/
      /*        if (Mii < EPSDIAG ){
                  diag[i] = invrho;
                  diagprev[i] = 1;
              }
              else {
                  diag[i] = 1.0/Mii;
                  diagprev[i] = 0;
              }
      */
    }
  }
  /*  den = 1.0/sqrt(qs);*/

  /*start iterations*/

  iter = 0;
  err  = 1.;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    incx = n;
    incy = 1;
    for (i = 0 ; i < n ; ++i)
    {
      ziprev = z[i];
      z[i] = 0.0;

      zi = -(q[i] - (rho * ziprev) + DDOT(n , &M[i] , incx , z , incy)) * diag[i];

      if (zi > 0) z[i] = zi;

    }
    /* **** Criterium convergence **** */

    /*     incx =  1;
         incy =  1;
    */
    /*     a1 = 1.0; */
    /*     b1 = 1.0; */

    /*     dgemv_( &NOTRANS , (integer *)&n , (integer *)&n , &a1 , M , (integer *)&n , z , (integer *)&incx , &b1 , w , (integer *)&incy ); */

    /*     a1   = -1.0;
         daxpy_( (integer *)&n , &a1 , z , (integer *)&incx , zprev , (integer *)&incy );
         num = dnrm2_( (integer *)&n, zprev , (integer *)&incx );*/
    /*     daxpy_( (integer *)&n , &a1 , w , (integer *)&incx , ww , (integer *)&incy ); */
    /*     num = dnrm2_( (integer *)&n, ww , (integer *)&incx ); */
    /*     err = sqrt(num)*den;*/

    /* **** Criterium convergence compliant with filter_result_LCP **** */
    lcp_compute_error(n, M, q, z, ispeak, w, &err);


    /*    buffer_errors[iter-1] = err;*/

    if (ispeak == 2)
    {
      printf(" # i%d -- %g : ", iter, err);
      for (i = 0 ; i < n ; ++i) printf(" %g", z[i]);
      for (i = 0 ; i < n ; ++i) printf(" %g", w[i]);
      printf("\n");
    }

    /* **** ********************* **** */

  }

  iparamLCP[2] = iter;
  dparamLCP[3] = err;

  if (ispeak > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of RPGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 1;
    }
    else
    {
      printf(" Convergence of RPGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else
  {
    if (err > tol) *info = 1;
    else *info = 0;
  }

  /*  if (iter == itermax)
    {
          ficbuffer_errors = fopen("errors.dat","w");
          if (ficbuffer_errors == NULL)
          {
              printf("Impossible d'ouvrir errors.dat !\n");
          } else {
              for(i = 0 ; i < itermax ; i++)
              {
                  fprintf(ficbuffer_errors,"%g\n",buffer_errors[i]);
              }
              fclose(ficbuffer_errors);
          }
    }
  */
  /*  free(ww);*/
  free(diag);
  /*  free(zprev);*/
  /*   free(diagprev);*/
  /*  free(buffer_errors);*/

  return;
}
