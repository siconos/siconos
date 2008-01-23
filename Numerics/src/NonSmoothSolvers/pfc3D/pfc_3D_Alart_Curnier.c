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
#include <time.h>
#include "NonSmoothDrivers.h"
#include "pfc_3D_Alart_Curnier.h"

/* Compute function G */
void Compute_G_AC(int m, double *G, double *x , double *C, double *y , double *b, double *param1, double *param2, double *param3, double rn, double rt, double coef)
{

  double zn , zt, zs, num, mrn, coef2;

  coef2 = coef * coef;

  y[0] = C[0 * m + 0] * x[0] + C[1 * m + 0] * x[1] + C[2 * m + 0] * x[2] + b[0];
  y[1] = C[0 * m + 1] * x[0] + C[1 * m + 1] * x[1] + C[2 * m + 1] * x[2] + b[1];
  y[2] = C[0 * m + 2] * x[0] + C[1 * m + 2] * x[1] + C[2 * m + 2] * x[2] + b[2];

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

  mrn = zt * zt + zs * zs;

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



/* Compute Jacobian of function G */
void Compute_JacG_AC(int m, double *JacG , double *x , double *C, double *y  , double *b, double *param1, double *param2, double *param3 , double rn, double rt, double coef)
{


  double zn , zt, zs, mrn, mrn3, coef2, num, a2, b2, ab;
  int i, j;

  coef2 = coef * coef;

  y[0] = C[0 * m + 0] * x[0] + C[1 * m + 0] * x[1] + C[2 * m + 0] * x[2] + b[0];
  y[1] = C[0 * m + 1] * x[0] + C[1 * m + 1] * x[1] + C[2 * m + 1] * x[2] + b[1];
  y[2] = C[0 * m + 2] * x[0] + C[1 * m + 2] * x[1] + C[2 * m + 2] * x[2] + b[2];


  /* Projection on [0, +infty[ and on D(0, mu*zn) */
  zn = x[0] - rn * y[0];
  if (zn > 0)
  {
    for (j = 0; j < m; ++j)
      JacG[j * m + 0] = C[j * m + 0];
  }
  else
  {
    JacG[0 * m + 0] = 1. / rn;
  }
  zt = x[1] - rt * y[1];
  zs = x[2] - rt * y[2];

  a2 = zt * zt;
  b2 = zs * zs;
  ab = zt * zs;
  mrn = zt * zt + zs * zs;

  if (mrn <= coef2 * x[0]*x[0])
    for (i = 1; i < m; ++i)
      for (j = 0; j < m; ++j)
        JacG[j * m + i] = C[j * m + i];

  else
  {
    num  = 1. / sqrt(mrn);
    mrn3 = 1. / sqrt(mrn) * sqrt(mrn) * sqrt(mrn);
    double rcof = coef / rt;

    JacG[0 * m + 1] = -rcof * (num * zt + x[0] * zt * rt * mrn3 * (C[0 * m + 1] * zt + C[0 * m + 2] * zs));

    JacG[0 * m + 2] = -rcof * (num * zs + x[0] * zs * rt * mrn3 * (C[0 * m + 1] * zt + C[0 * m + 2] * zs));

    JacG[1 * m + 1] = (1 - coef * x[0] * (num * (1 - rt * C[1 * m + 1]) - zt * mrn3 * ((1 - rt * C[1 * m + 1]) * zt - rt * C[1 * m + 2] * zs))) / rt;

    JacG[2 * m + 1] =  - rcof * x[0] * ((-num * rt * C[2 * m + 1]) - zt * mrn3 * ((1 - rt * C[2 * m + 2]) * zs - rt * C[2 * m + 1] * zt));

    JacG[1 * m + 2] =  - rcof * x[0] * ((-num * rt * C[1 * m + 2]) - zs * mrn3 * ((1 - rt * C[1 * m + 1]) * zt - rt * C[1 * m + 2] * zs));

    JacG[2 * m + 2] = (1 - coef * x[0] * (num * (1 - rt * C[2 * m + 2]) - zs * mrn3 * ((1 - rt * C[2 * m + 2]) * zs - rt * C[2 * m + 1] * zt))) / rt;

  }
}

/*  Inverse Matrix 3x3  */
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

/* Linesearch */
void Linesearch_AC(int n, double *zz, double *ww, double *www, double *b, double *C, double *param1, double *param2, double *param3, double an, double at, double mu, double err1)
{
  double err2, alpha, qs, a1;
  int i, incx, incy;
  double *zzzz, *wwww, *G;
  zzzz = (double*)malloc(n * sizeof(double));
  wwww = (double*)malloc(n * sizeof(double));
  G    = (double*)malloc(n * sizeof(double));

  for (i = 0 ; i < n ; ++i)
    zzzz[i] = wwww[i] = G[i] = 0.;

  incx =  1;
  incy =  1;
  a1 = -1.;
  qs = 1.;

  alpha = 1.;
  while (alpha > 0.05)
  {

    DCOPY(n , zz , incx , zzzz , incy);
    DAXPY(n , alpha , www , incx , zzzz , incy);

    Compute_G_AC(n , G , zzzz , C , ww , b, param1, param2, param3 , an , at , mu);

    err2 = DNRM2(n, G , incx);

    if (err2 < err1) break;
    alpha = alpha * 0.5;
  }
  err1 = err2;

  DCOPY(n , zzzz , incx , zz , incy);

  DCOPY(n , wwww , incx , ww , incy);

  free(zzzz);
  free(wwww);
  free(G);


}


