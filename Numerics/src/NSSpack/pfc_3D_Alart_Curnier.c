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
/*!\file pfc_3D_Alart_Curnier.c\n
 *
 * \fn pfc_3D_Alart_Curnier( int *nn , double *vec , double *q , double *z , double *w , double coef)\n
 *
 *
 * This subroutine allows the computation of the equation to be solved (We try to solve G = 0) with Newton method. \n
 *
 *            |          w - M z - q             |\n
 *            |                                  |\n
 *       G =  | 1/rn*[wn - pos_part(wn - rn*zn)] |;\n
 *            |                                  |\n
 *            |   1/rt*[wT - proj(wT - rt*zT)]   |\n
 *
 *
 *
 * where M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an  * n-dimensional vector.\n
 *
 * Inside, we compute G, JacG and the inverse of JacG. \n
 *
 * \fn  void Compute_G( int m, double *G, double *x , double *y , double rn, double rt, double coef )\n
 *
 * \fn  void Compute_JacG( int m, double *A, double *B, double *JacG, double *C, double *x , double *y , double rn, double rt, double coef )\n
 *
 * \fn  void matrix_inv3(double *a, double *b)\n
 *
 *
 * Generic pfc_3D parameters:\n
 *
 * \param y       Modified parameter which contains the initial solution and returns the solution of the problem.\n
 * \param x       Modified parameter which returns the solution of the problem.\n
 * \param rn, rt  Modified parameter which contains the augmentation parameters values.\n
 * \param coef    Modified parameter which contains the friction coefficient value.\n
 *
 *
 * \author Houari Khenous (20/09/2007)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include <time.h>
#include <NSSpack.h>
#include <pfc_3D_Alart_Curnier.h>




/* Compute function G */
void Compute_G_AC(int m, double *G, double *C, double *x , double *y , double *b, double *param2, double *param3, double *param4, double rn, double rt, double coef)
{
  int incx, incy;
  double qs, zn , zt, zs, num, mrn, coef2, a2, b2, ab;
  coef2 = coef * coef;
  incx =  1;
  incy =  1;
  qs = 1.;


  DCOPY(m , b , incx , y , incy);

  DGEMV(LA_NOTRANS , m , m , qs , C , m , x , incx , qs , y , incy);

  /* Projection on [0, +infty[ and on D(0, mu*zn) */
  zn = x[0] - rn * y[0];
  if (zn > 0)
  {
    G[0] = y[0];
  }
  else
  {
    G[0] = x[0] / rn;
  }
  zt = x[1] - rt * y[1];
  zs = x[2] - rt * y[2];

  a2 = zt * zt;
  b2 = zs * zs;
  ab = zt * zs;
  mrn = zt * zt + zs * zs;

  // if the radius is negative or the vector is null projection on the disk = 0
  if (mrn < 1e-16)
  {
    G[1] = x[1] / rt;
    G[2] = x[2] / rt;
  }

  // if the radius is positive and the vector is non null, we compute projection on the disk
  else
  {
    if (mrn <= coef2 * x[0]*x[0])
    {
      G[1] = y[1];
      G[2] = y[2];
    }
    else
    {
      num  = coef / sqrt(mrn);
      G[1] = (x[1] - zt * x[0] * num) / rt;
      G[2] = (x[2] - zs * x[0] * num) / rt;
    }
  }
}



/* Compute Jacobian of function G */
void Compute_JacG_AC(int m, double *JacG , double *C , double *x , double *y , double *b , double *A , double *B, double rn, double rt, double coef)
{

  double a1, qs, zn , zt, zs, mrn, mrn3, coef2, num, a2, b2, ab;
  double *AC;
  int mm;
  int incx, incy;
  mm = m * m;
  AC   = (double*)malloc(mm * sizeof(double));

  incx =  1;
  incy =  1;
  qs = 1.;

  coef2 = coef * coef;

  DCOPY(m , b , incx , y , incy);

  DGEMV(LA_NOTRANS , m , m , qs , C , m , x , incx , qs , y , incy);


  /* Projection on [0, +infty[ and on D(0, mu*zn) */
  zn = x[0] - rn * y[0];
  if (zn > 0)
  {
    A[0 * m + 0] = 1.;
  }
  else
  {
    B[0 * m + 0] = 1. / rn;
  }
  zt = x[1] - rt * y[1];
  zs = x[2] - rt * y[2];

  a2 = zt * zt;
  b2 = zs * zs;
  ab = zt * zs;
  mrn = zt * zt + zs * zs;

  // if the radius is negative or the vector is null projection on the disk = 0
  if (mrn < 1e-16)
    B[1 * m + 1] = B[2 * m + 2] = 1. / rt;
  // if the radius is positive and the vector is non null, we compute projection on the disk
  else
  {
    if (mrn <= coef2 * x[0]*x[0])
      A[1 * m + 1] = A[2 * m + 2] = 1.;
    else
    {
      num  = coef / sqrt(mrn);
      B[0 * m + 1] = - num * zt / rt;
      B[0 * m + 2] = - num * zs / rt;
      mrn3 = sqrt(mrn) * sqrt(mrn) * sqrt(mrn);
      A[1 * m + 1] =  coef * x[0] * b2 / mrn3;
      A[2 * m + 1] = -coef * x[0] * ab / mrn3;
      A[1 * m + 2] = -coef * x[0] * ab / mrn3;
      A[2 * m + 2] =  coef * x[0] * a2 / mrn3;
      B[1 * m + 1] = (1. - rt * coef * x[0] * b2 / mrn3) / rt;
      B[2 * m + 1] =  coef * x[0] * ab / mrn3;
      B[1 * m + 2] =  coef * x[0] * ab / mrn3;
      B[2 * m + 2] = (1. - rt * coef * x[0] * a2 / mrn3) / rt;
    }
  }

  DCOPY(mm , B , incx , JacG , incy);

  incx =  3;
  a1 = 1.0;
  DGEMM(LA_NOTRANS , LA_NOTRANS , m , m , m , a1 , A , m , C , incx , a1 , JacG , incx);

  free(AC);
}

//_/_/   Inverse Matrix 3x3  _/_//
void matrix_inv3(double *a, double *b)
{
  double det;
  det = a[0 * 3 + 0] * (a[1 * 3 + 1] * a[2 * 3 + 2] - a[1 * 3 + 2] * a[2 * 3 + 1]) + a[0 * 3 + 1] * (a[2 * 3 + 0] * a[1 * 3 + 2] - a[2 * 3 + 2] * a[1 * 3 + 0]) + a[0 * 3 + 2] * (a[1 * 3 + 0] * a[2 * 3 + 1] - a[1 * 3 + 1] * a[2 * 3 + 0]);

  b[0 * 3 + 0] = (a[1 * 3 + 1] * a[2 * 3 + 2] - a[1 * 3 + 2] * a[2 * 3 + 1]) / det;
  b[1 * 3 + 0] = (a[2 * 3 + 0] * a[1 * 3 + 2] - a[2 * 3 + 2] * a[1 * 3 + 0]) / det;
  b[2 * 3 + 0] = (a[1 * 3 + 0] * a[2 * 3 + 1] - a[1 * 3 + 1] * a[2 * 3 + 0]) / det;

  b[0 * 3 + 1] = (a[2 * 3 + 1] * a[0 * 3 + 2] - a[2 * 3 + 2] * a[0 * 3 + 1]) / det;
  b[1 * 3 + 1] = (a[0 * 3 + 0] * a[2 * 3 + 2] - a[0 * 3 + 2] * a[2 * 3 + 0]) / det;
  b[2 * 3 + 1] = (a[2 * 3 + 0] * a[0 * 3 + 1] - a[2 * 3 + 1] * a[0 * 3 + 0]) / det;

  b[0 * 3 + 2] = (a[0 * 3 + 1] * a[1 * 3 + 2] - a[1 * 3 + 1] * a[0 * 3 + 2]) / det;
  b[1 * 3 + 2] = (a[1 * 3 + 0] * a[0 * 3 + 2] - a[1 * 3 + 2] * a[0 * 3 + 0]) / det;
  b[2 * 3 + 2] = (a[0 * 3 + 0] * a[1 * 3 + 1] - a[0 * 3 + 1] * a[1 * 3 + 0]) / det;


}

void Linesearch_AC(int n, double *G, double *zz, double *ww, double *www, double *b, double *C, double *zzzz, double *wwww, double an, double at, double mu, double err1)
{
  double err2, alpha, qs, a1;
  int mm, incx, incy;
  double *A, *B, *AC;

  mm = n * n;

  AC   = (double*)malloc(mm * sizeof(double));
  A    = (double*)malloc(mm * sizeof(double));
  B    = (double*)malloc(mm * sizeof(double));

  incx =  1;
  incy =  1;
  a1 = -1.;
  qs = 1.;

  alpha = 1.;
  while (alpha > 0.05)
  {

    DCOPY(n , zz , incx , zzzz , incy);
    DAXPY(n , alpha , www , incx , zzzz , incy);

    Compute_G_AC(n , G , C , zzzz , wwww , b , A , B , AC , an , at , mu);

    err2 = DNRM2(n, G , incx);

    if (err2 < err1) break;
    alpha = alpha * 0.5;
  }
  err1 = err2;

  DCOPY(n , zzzz , incx , zz , incy);

  DCOPY(n , wwww , incx , ww , incy);

  free(A);
  free(B);
  free(AC);
}


