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
/*!\file pr_latin.c

This subroutine allows the primal resolution of relay problems (PR).\n
   Try \f$(z,w)\f$ such that:\n

   \f$
    \left\lbrace
     \begin{array}{l}
       w - M z = q \\
       -w \in \partial\psi_{[-b, a]}(z)\\
     \end{array}
    \right.
   \f$

 here M is an (\f$nn \times nn \f$)-matrix, q an nn-dimensional vector, z an nn-dimensional  vector and w an nn-dimensional vector.

*/

/*!\fn void pr_latin(double *vec, double *qq, int *nn, double * k_latin, double *a, double *b, int * itermax, double * tol, int *chat, double *z, double *w, int *it_end, double * res, int *info)

   pr_latin is a specific latin solver for primal relay problems.


   \param vec       On enter, a (\f$nn \times nn\f$)-vector of doubles containing the components of the matrix with a fortran storage.
   \param qq        On enter, a nn-vector of doubles containing the components of the vector.
   \param nn        On enter, an integer which represents the dimension of the second member.
   \param k_latin   On enter, a double, the latin coefficient (strictly non negative).
   \param a         On enter, a nn-vector of doubles, the upper bound.
   \param b         On enter, a nn-vector of doubles, the down bound.
   \param itermax   On enter, an integer which represents the maximum iterations required.
   \param tol       On enter, a double which contains the tolerance required.
   \param chat      On enter, an integer the output log identifiant:\n                              0 : no output\n
                    >0 : active screen output
   \param it_end    On return, an integer which represents the number of iterations carried out.
   \param res       On return, a double, the error value.
   \param z         On return,a nn-vector of doubles wich contains the solution of the problem.
   \param w         On return, a nn-vector of doubles which contains the solution of the problem.
   \param info      On return, an integer which represents the termination reason: \n
                    0 = convergence,\n
        1 = no convergence,\n
        2 = Cholesky factorization failed,\n
        3 = Nul diagonal term\n


   \author Nineb Sheherazade.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"



void pr_latin(double *vec, double *qq, int *nn, double * k_latin, double *a, double *b, int * itermax, double * tol, int *chat, double *z, double *w, int *it_end, double * res, int *info)
{


  int i, j, iter1, nrhs, info2, ispeak = *chat;
  int n = *nn,  itt = *itermax;
  integer incx = 1, incy = 1;
  double errmax = *tol, alpha, beta, mina;
  double  err1, num11, err0;
  double den11, den22;
  double *wc, *zc, *wnum1, *znum1;
  double *zt;
  double *k, *kinv, *DPO;

  char trans = 'T', uplo = 'U', notrans = 'N', diag = 'N';




  /*             Allocations                           */

  k     = (double*)malloc(n * n * sizeof(double));
  DPO   = (double*)malloc(n * n * sizeof(double));
  kinv  = (double*)malloc(n * n * sizeof(double));
  wc    = (double*) malloc(n * sizeof(double));
  zc    = (double*) malloc(n * sizeof(double));
  znum1 = (double*) malloc(n * sizeof(double));
  wnum1 = (double*) malloc(n * sizeof(double));
  zt    = (double*) malloc(n * sizeof(double));


  /*             Initialization                          */


  for (i = 0; i < n * n; i++)
  {
    k[i]    = 0.0;
    kinv[i] = 0.0;
    DPO[i]  = 0.0;
    if (i < n)
    {
      wc[i]    = 0.0;
      zc[i]    = 0.0;
      z[i]     = 0.0;
      w[i]     = 0.0;
      znum1[i] = 0.0;
      wnum1[i] = 0.0;
      zt[i]    = 0.0;
    }
  }


  for (i = 0; i < n; i++)
  {
    k[i + n * i] =  *k_latin * vec[i * n + i];

    if (fabs(k[i + n * i]) < 1e-12)
    {

      if (ispeak > 0)
        printf("\n Warning nul diagonal term in k matrix \n ");

      free(k);
      free(kinv);
      free(DPO);
      free(wc);
      free(zc);
      free(znum1);
      free(wnum1);
      free(zt);

      *info = 3;
      return;

    }
    else

      kinv[i + n * i] = 1 / k[i + n * i];
  }






  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      DPO[i + n * j] =  vec[j * n + i] + k[i + n * j];
    }


  /*             Cholesky                           */

  dpotrf_(&uplo, (integer *)&n, DPO , (integer *)&n, (integer *)&info2);



  if (info2 != 0)
  {
    printf("\n Matter with Cholesky factorization \n");

    free(k);
    free(kinv);
    free(DPO);
    free(wc);
    free(zc);
    free(znum1);
    free(wnum1);
    free(zt);

    *info = 2;
    return;
  }


  /*             End of Cholesky                        */



  /*              Iteration loops                         */

  iter1 = 0;
  err1 = 1.;

  while ((iter1 < itt) && (err1 > errmax))
  {

    /*          Linear stage (zc,wc) -> (z,w)            */

    alpha = 1.;
    beta  = 1.;
    dgemv_(&trans, (integer *)&n, (integer *)&n, &alpha, k, (integer *)&n, zc, &incx, &beta, wc, &incy);

    dcopy_((integer *)&n, qq, &incx, znum1, &incy);


    alpha = -1.;
    dscal_((integer *)&n , &alpha , znum1 , &incx);

    alpha = 1.;
    daxpy_((integer *)&n, &alpha, wc, &incx, znum1, &incy);


    nrhs = 1;
    dtrtrs_(&uplo, &trans, &diag, (integer *)&n, (integer *)&nrhs, DPO, (integer *)&n, znum1, (integer *)&n, (integer*)&info2);

    dtrtrs_(&uplo, &notrans, &diag, (integer *)&n, (integer *)&nrhs, DPO, (integer *)&n, znum1, (integer *)&n, (integer*)&info2);

    dcopy_((integer *)&n, znum1, &incx, z, &incy);

    dcopy_((integer *)&n, wc, &incx, w, &incy);

    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, (integer *)&n, (integer *)&n, &alpha, k, (integer *)&n, z, &incx, &beta, w, &incy);

    /*            Local stage (z,w)->(zc,wc)            */


    dcopy_((integer *)&n, z, &incx, zt, &incy);

    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, (integer *)&n, (integer *)&n, &alpha, kinv, (integer *)&n, w, &incx, &beta, zt, &incy);

    for (i = 0; i < n; i++)
    {
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
        zc[i] = mina;
      }
      else
      {
        zc[i] = -b[i];
      }
    }

    dcopy_((integer *)&n, w, &incx, wc, &incy);

    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, (integer *)&n, (integer *)&n, &alpha, k, (integer *)&n, z, &incx, &beta, wc, &incy);

    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, (integer *)&n, (integer *)&n, &alpha, k, (integer *)&n, zc, &incx, &beta, wc, &incy);

    /*             Convergence criterium                     */

    dcopy_((integer *)&n, w, &incx, wnum1, &incy);
    alpha = -1.;
    daxpy_((integer *)&n, &alpha, wc, &incx, wnum1, &incy);
    dcopy_((integer *)&n, z, &incx, znum1, &incy);
    daxpy_((integer *)&n, &alpha, zc, &incx, znum1, &incy);
    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, (integer *)&n, (integer *)&n, &alpha, k, (integer *)&n, znum1, &incx, &beta, wnum1, &incy);

    /*    num1(:) =(w(:)-wc(:))+matmul( k(:,:),(z(:)-zc(:)))   */

    num11 = 0.;
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, (integer *)&n, (integer *)&n, &alpha, kinv, (integer *)&n, wnum1, &incx, &beta, znum1, &incy);

    num11 = ddot_((integer *)&n, wnum1, &incx, znum1, &incy);

    dcopy_((integer *)&n, z, &incx, znum1, &incy);

    daxpy_((integer *)&n, &alpha, zc, &incx, znum1, &incy);

    beta = 0.;
    alpha = 1.;
    dgemv_(&trans, (integer *)&n, (integer *)&n, &alpha, k, (integer *)&n, znum1, &incx, &beta, wnum1, &incy);

    den22 = ddot_((integer *)&n, znum1, &incx, wnum1, &incy);

    dcopy_((integer *)&n, w, &incx, wnum1, &incy);

    alpha = 1.;
    daxpy_((integer *)&n, &alpha, wc, &incx, wnum1, &incy);

    beta = 0.;
    alpha = 1.;
    dgemv_(&trans, (integer *)&n, (integer *)&n, &alpha, kinv, (integer *)&n, wnum1, &incx, &beta, znum1, &incy);

    den11 = ddot_((integer *)&n, wnum1, &incx, znum1, &incy);

    err0 = num11 / (den11 + den22);

    err1 = sqrt(err0);

    iter1 = iter1 + 1;

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



  free(wc);
  free(zc);
  free(znum1);
  free(wnum1);
  free(zt);
  free(DPO);
  free(k);
  free(kinv);



}
