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

/*!\file pr_nlgs.c


This subroutine allows the primal resolution of relay problems(PR).\n

   Try \f$(z,w)\f$ such that:\n
\f$
\left\lbrace
\begin{array}{l}
 w - M z = q\\
-w \in \partial\psi_{[-b, a]}(z)\\
\end{array}
\right.
\f$

 here M is an (nn \f$\times\f$nn)-matrix, q an nn-dimensional vector, z an nn-dimensional  vector and w an nn-dimensional vector.

*/


/*!\fn void  pr_nlgs(double *vec, double *q, int *nn, double *a, double *b, int * itermax, double * tol, int *chat, double *z, double *w, int *it_end, double * res, int *info)

  \n \n
   pr_nlgs is a specific nlgs(non linear Gauss-Seidel) solver for primal relay problems.\n



   \param vec       On enter, a (nn\f$\times\f$nn)-vector of doubles containing the components of the matrix with a fortran storage.
   \param qq        On enter, a nn-vector of doubles containing the components of the vector.
   \param nn        On enter, an integer, the dimension of the second member.
   \param a         On enter, a nn-vector of doubles, the upper bound.
   \param b         On enter, a nn-vector of doubles, the down bound.
   \param itermax   On enter, an integer, the maximum iterations required.
   \param tol       On enter, a double, the tolerance required.
   \param chat      On enter, an integer, the output log identifiant:\n
                     0 =  no output \n
                    >0 =  active screen output\n


   \param it_end    On return, an integer, the number of iterations carried out.
   \param res       On return, a double, the error value.
   \param z         On return, a nn-vector of doubles, the solution of the problem.
   \param w         On return, a nn-vector of double, the solution of the problem.
   \param info      On return, an integer, the termination reason:\n
                    0 = convergence,\n
        1 = no convergence,\n
        2 = Nul diagonal term\n


   \author Nineb Sheherazade.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"


void pr_nlgs(double *vec, double *q, int *nn, double *a, double *b, int * itermax, double * tol, int *chat, double *z, double *w, int *it_end, double * res, int *info)
{



  int i, j, iter1, k, ispeak = *chat;
  int n = *nn,  itt = *itermax;
  int incx = 1, incy = 1;

  double errmax = *tol, alpha, beta, mina;
  double err1, num, den, avn, apn, xn;
  double *zt, *wnum1;


  wnum1 = (double*) malloc(n * sizeof(double));
  zt    = (double*) malloc(n * sizeof(double));



  for (i = 0; i < n; i++)
  {
    w[i]     = 0.0;
    z[i]     = 0.0;
    zt[i]    = 0.0;
    wnum1[i] = 0.0;
  }


  iter1 = 0;
  err1  = 1.;





  while ((iter1 < itt) && (err1 > errmax))
  {
    iter1 = iter1 + 1;

    for (i = 0; i < n; i++)
    {
      avn = 0.;
      apn = 0.;

      for (j = 0; j <= i - 1; j++)
        avn = avn + vec[j * n + i] * z[j];

      for (k = i + 1 ; k < n ; k++)
        apn = apn + vec[k * n + i] * z[k];

      xn = -q[i] - avn - apn;


      if (fabs(vec[i * n + i]) < 1e-16)
      {
        if (ispeak > 0)
          printf("\n Warning nul diagonal term of M \n");


        free(zt);
        free(wnum1);

        *info = 2;

        return;

      }
      else
        zt[i] = 1 / vec[i * n + i] * xn;

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

      w[i] = vec[i * n + i] * z[i] - xn;

    }

    /*          Convergence criterium      */

    DCOPY(n, w, incx, wnum1, incy);

    alpha = -1.;
    DAXPY(n, alpha, q, incx, wnum1, incy);

    alpha = 1.;
    beta = -1.;
    DGEMV(LA_TRANS, n, n, alpha, vec, n, z, incx, beta, wnum1, incy);

    num = DDOT(n, wnum1, incx, wnum1, incy);

    den = DDOT(n, q, incx, q, incy);

    err1 = sqrt(num) / sqrt(den);

    *it_end = iter1;
    *res = err1;



  }


  if (err1 > errmax)
  {
    if (ispeak > 0)
      printf("No convergence after %d iterations, the residue is %g\n", iter1, err1);

    *info = 1;
  }
  else
  {
    if (ispeak > 0)
      printf("Convergence after %d iterations, the residue is %g \n", iter1, err1);

    *info = 0;
  }



  free(wnum1);
  free(zt);



}
