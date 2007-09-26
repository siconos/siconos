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
/*!\file pfc_3D_Fischer_.c\n
 *
 * \fn pfc_3D_Fischer( int *nn , double *vec , double *q , double *z , double *w , double coef)\n
 *
 *
 *
 *
 * \author Houari Khenous (24/09/2007)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include <time.h>
#include <NSSpack.h>
#include <pfc_3D_Fischer_Burmeister.h>




/* Compute function G */
void Compute_G_FB(int m, double *G, double *C, double *y , double *x , double *b, double *Ip, double *IP, double *I3, double coef21 , double r, double coef)
{

  double coef2;
  coef2 = coef * coef;
  double *A;
  A     = (double*)malloc(m * sizeof(double));

  /* Ip = matrix_inv2(Ip`)*/
  /* IP = matrix_inv2(Ip`)*mup */
  /* I3 = matrix_inv2(Ip)*e3 */


  A[0] = (C[0 * 3 + 0] + coef * (C[1 * 3 + 0] * IP[0] + C[2 * 3 + 0] * IP[1])) * y[0]
         - (C[1 * 3 + 0] * Ip[0 * 2 + 0] + C[2 * 3 + 0] * Ip[0 * 2 + 1]) * y[1]
         - (C[1 * 3 + 0] * Ip[1 * 2 + 0] + C[2 * 3 + 0] * Ip[1 * 2 + 1]) * y[2]
         + b[0]
         + x[0];

  A[1] = (-(C[0 * 3 + 1] * Ip[0 * 2 + 0] + C[0 * 3 + 2] * Ip[0 * 2 + 1]) + C[1 * 3 + 1] * IP[0] + C[2 * 3 + 1] * IP[1]) * y[0]
         + (C[1 * 3 + 1] * Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + C[2 * 3 + 1] * Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + C[1 * 3 + 2] * Ip[0 * 2 + 0] * Ip[0 * 2 + 1] + C[2 * 3 + 2] * Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * y[1]
         + (C[1 * 3 + 1] * Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + C[2 * 3 + 1] * Ip[0 * 2 + 0] * Ip[1 * 2 + 1] + C[1 * 3 + 2] * Ip[1 * 2 + 0] * Ip[0 * 2 + 1] + C[2 * 3 + 2] * Ip[1 * 2 + 1] * Ip[1 * 2 + 0]) * y[2]
         - I3[0] * y[3]
         - 2 * ((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2])) * y[4]
         - Ip[0 * 2 + 0] * b[1] - Ip[0 * 2 + 1] * b[2]
         + x[1];


  A[2] = (-(C[0 * 3 + 1] * Ip[1 * 2 + 0] + C[0 * 3 + 2] * Ip[1 * 2 + 1]) + C[1 * 3 + 2] * IP[0] + C[2 * 3 + 2] * IP[1]) * y[0]
         + (C[1 * 3 + 1] * Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + C[2 * 3 + 1] * Ip[0 * 2 + 1] * Ip[1 * 2 + 0] + C[1 * 3 + 2] * Ip[0 * 2 + 0] * Ip[1 * 2 + 1] + C[2 * 3 + 2] * Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * y[1]
         + (C[1 * 3 + 1] * Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + C[2 * 3 + 1] * Ip[1 * 2 + 0] * Ip[1 * 2 + 1] + C[1 * 3 + 2] * Ip[1 * 2 + 0] * Ip[1 * 2 + 1] + C[2 * 3 + 2] * Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * y[2]
         - I3[1] * y[3]
         - 2 * ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2])) * y[4]
         - Ip[1 * 2 + 0] * b[1] - Ip[1 * 2 + 1] * b[2]
         + x[2];

  A[3] = coef * (1 - I3[0] - I3[1]) * y[0] + I3[0] * y[1] + I3[1] * y[2]
         + x[3];

  A[4] = coef2 * y[0] * y[0] - sqrt((coef * y[0] - y[1]) * (coef * y[0] - y[1]) + (coef * y[0] - y[2]) * (coef * y[0] - y[2]))
         + x[4];
  int incx = 1;
  r = DNRM2(5, A , incx);

  for (int i; i < m; ++i)
    G[i] = A[i] - r;
  free(A);
}



/* Compute Jacobian of function G */
void Compute_JacG_FB(int m, double *JacG , double *C , double *y , double *x , double *Ip , double *IP, double *I3, double coef3 , double coef21, double coef)
{

  double *A, *B;
  double coef2;
  int mm;
  mm = m * m;

  A     = (double*)malloc(mm * sizeof(double));
  B     = (double*)malloc(mm * sizeof(double));

  coef2 = coef * coef;

  A[0 * 5 + 0] = C[0 * 3 + 0] + coef * (C[1 * 3 + 0] * IP[0] + C[2 * 3 + 0] * IP[1]);
  A[1 * 5 + 0] = C[1 * 3 + 0] * Ip[0 * 2 + 0] + C[2 * 3 + 0] * Ip[0 * 2 + 1];
  A[2 * 5 + 0] = C[1 * 3 + 0] * Ip[1 * 2 + 0] + C[2 * 3 + 0] * Ip[1 * 2 + 1];

  A[0 * 5 + 1] = -(C[0 * 3 + 1] * Ip[0 * 2 + 0] + C[0 * 3 + 2] * Ip[0 * 2 + 1]) + C[1 * 3 + 1] * IP[0] + C[2 * 3 + 1] * IP[1];
  A[1 * 5 + 1] = C[1 * 3 + 1] * Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + C[2 * 3 + 1] * Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + C[1 * 3 + 2] * Ip[0 * 2 + 0] * Ip[0 * 2 + 1] + C[2 * 3 + 2] * Ip[0 * 2 + 1] * Ip[0 * 2 + 1];
  A[2 * 5 + 1] = C[1 * 3 + 1] * Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + C[2 * 3 + 1] * Ip[0 * 2 + 0] * Ip[1 * 2 + 1] + C[1 * 3 + 2] * Ip[1 * 2 + 0] * Ip[0 * 2 + 1] + C[2 * 3 + 2] * Ip[1 * 2 + 1] * Ip[1 * 2 + 0];
  A[3 * 5 + 1] = -I3[0];

  A[0 * 5 + 2] = -(C[0 * 3 + 1] * Ip[1 * 2 + 0] + C[0 * 3 + 2] * Ip[1 * 2 + 1]) + C[1 * 3 + 2] * IP[0] + C[2 * 3 + 2] * IP[1];
  A[1 * 5 + 2] = C[1 * 3 + 1] * Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + C[2 * 3 + 1] * Ip[0 * 2 + 1] * Ip[1 * 2 + 0] + C[1 * 3 + 2] * Ip[0 * 2 + 0] * Ip[1 * 2 + 1] + C[2 * 3 + 2] * Ip[0 * 2 + 1] * Ip[1 * 2 + 1];
  A[2 * 5 + 2] = C[1 * 3 + 1] * Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + C[2 * 3 + 1] * Ip[1 * 2 + 0] * Ip[1 * 2 + 1] + C[1 * 3 + 2] * Ip[1 * 2 + 0] * Ip[1 * 2 + 1] + C[2 * 3 + 2] * Ip[1 * 2 + 1] * Ip[1 * 2 + 1];
  A[3 * 5 + 2] = -I3[1];

  A[0 * 5 + 3] = coef * (1 - I3[0] - I3[1]);
  A[1 * 5 + 3] = I3[0];
  A[2 * 5 + 3] = I3[1];

  B[0 * 5 + 1] = -2 * ((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1])) * y[4];
  B[0 * 5 + 2] = -2 * ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1])) * y[4];
  B[1 * 5 + 1] = -2 * (Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * y[4];
  B[1 * 5 + 2] = -2 * (Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * y[4];
  B[2 * 5 + 1] = -2 * (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * y[4];
  B[2 * 5 + 2] = -2 * (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * y[4];
  B[4 * 5 + 1] = -2 * ((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));
  B[4 * 5 + 2] = -2 * ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));
  B[0 * 5 + 4] = 2 * coef2 * y[0] - 2 * coef * (((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]))
                 + ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2])));
  B[1 * 5 + 4] = 2 * ((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));
  B[2 * 5 + 4] = 2 * ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));
  for (int i; i < m; ++i)
  {
    for (int j; j < m; ++j)
    {
      JacG[j * m + i] = A[j * m + i] + B[j * m + i];
    }
  }

  free(A);
  free(B);
}

//_/_/   Inverse Matrix 3x3  _/_//
void matrix_inv2(double *a, double *b)
{
  double det;
  det = a[0 * 2 + 0] * a[1 * 2 + 1] - a[1 * 2 + 0] * a[0 * 2 + 1];

  b[0 * 2 + 0] =  a[1 * 2 + 1] / det;
  b[0 * 2 + 1] = -a[0 * 2 + 1] / det;
  b[1 * 2 + 0] = -a[1 * 2 + 0] / det;
  b[1 * 2 + 1] =  a[0 * 2 + 0] / det;

}


void Linesearch_FB(int n, double *Ip, double *zz, double *ww, double *www, double *b, double *C, double *IP, double *I3, double coef21, double coef3, double mu, double err1)
{

  double alpha, err2;
  double *G, *zzzz;
  double qs, a1;
  int i, incx, incy;

  zzzz = (double*)malloc(n * sizeof(double));
  G    = (double*)malloc(n * sizeof(double));

  for (i = 0 ; i < n ; ++i)
    zzzz[i] = G[i] = 0.;

  incx =  1;
  incy =  1;
  a1 = -1.;
  qs = 1.;

  alpha = 1.;
  while (alpha > 0.05)
  {

    DCOPY(n , zz , incx , zzzz , incy);
    DAXPY(n , alpha , www , incx , zzzz , incy);

    Compute_G_FB(n , G , C , zz , ww , b , Ip , IP , I3 , coef21 , coef3 , mu);
    err2 = DNRM2(n, G , incx);

    if (err2 < err1) break;
    alpha = alpha * 0.5;
  }
  err1 = err2;

  DCOPY(n , zzzz , incx , zz , incy);

  free(zzzz);
  free(G);

}
