/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2006.
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


/*!\file rp_gsnl.c

  This subroutine allows the primal resolution of contact problems with friction.

   Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
-w \in \partial\psi_{[-a, a]}(z)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

*/



double ddot_(int *, double [], int *, double [], int*);


/*!\fn  rp_gsnl(double vec[],double *qq,int *nn,double a[],int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)

   rp_gsnl is a specific gsnl (Gauss Seidel Non Linear)solver for relay problems.

   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param qq On enter a pointer over doubles containing the components of the double vector.
   \param nn On enter a pointer over integers, the dimension of the second member.
   \param a On enter a pointer over doubles, the bound.
   \param itermax On enter a pointer over integers, the maximum iterations required.
   \param tol On enter a pointer over doubles, the tolerance required.
   \param it_end On enter a pointer over integers, the number of iterations carried out.
   \param res On return a pointer over doubles, the error value.
   \param z On return double vector, the solution of the problem.
   \param w On return double vector, the solution of the problem.
   \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).

   \author Nineb Sheherazade.
 */


rp_gsnl(double vec[], double *q, int *nn, double a[], double b[], int * itermax, double * tol, double z[], double w[], int *it_end, double * res, int *info)
{

  int i, j, iter1, k;
  int n = *nn, incx = 1, incy = 1, itt = *itermax;
  double errmax = *tol, alpha, beta, mina;
  double err1, num, den, avn, avt, apn, apt, xn;
  double *zt, *wnum1;
  char trans = 'T';
  double M[n][n];


  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      M[i][j] = vec[i * n + j];

  wnum1 = (double*) malloc(n * sizeof(double));
  zt = (double*) malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    w[i] = 0.;
    z[i] = 0.;
    zt[i] = 0.;
  }

  iter1 = 1;
  err1 = 1.;

  while ((iter1 < itt) && (err1 > errmax))
  {
    iter1 = iter1 + 1;
    dcopy_(&n, q, &incx, wnum1, &incy);
    alpha = 1.;
    beta = -1.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, wnum1, &incy);

    for (i = 0; i < n; i++)
    {
      avn = 0.;
      apn = 0.;
      for (j = 0; j <= i - 1; j++)
        avn = avn + M[i][j] * z[j];
      for (k = i + 1; k < n; k++)
        apn = apn + M[i][k] * z[k];
      xn = q[i] - avn - apn;
      zt[i] = 1 / M[i][i] * xn;

      if (a[i] < zt[i])
      {
        mina = a[i];
      }
      else
      {
        mina = zt[i];
      }

      if (-b[i] < mina)
      {
        z[i] = mina;
      }
      else
      {
        z[i] = -b[i];
      }
      w[i] = M[i][i] * z[i] - xn;
    }

    /* ///////// convergence criterium ///////////// */

    dcopy_(&n, w, &incx, wnum1, &incy);
    alpha = 1.;
    daxpy_(&n, &alpha, q, &incx, wnum1, &incy);
    alpha = 1.;
    beta = -1.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, wnum1, &incy);
    num = ddot_(&n, wnum1, &incx, wnum1, &incy);
    den = ddot_(&n, q, &incx, q, &incy);

    err1 = sqrt(num) / sqrt(den);
    *it_end = iter1;
    *res = err1;
  }


  if (err1 > errmax)
  {
    printf("no convergence after %d iterations, the residue is %g\n", iter1, err1);
    *info = 1;
  }
  else
  {
    printf("there is convergence after %d iterations, the residue is %g \n", iter1, err1);
    *info = 0;
  }

  free(wnum1);
  free(zt);
}
