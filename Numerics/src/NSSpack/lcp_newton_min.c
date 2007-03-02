/* Siconos-Numerics version 2.0.1, Copyright INRIA 2005-2006.
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

/*!\file lcp_newton_min.c


  This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
  Try \f$(z,w)\f$ such that:\n

  \f$
   \left\lbrace
    \begin{array}{l}
    0 \le z \perp w - M z  = q \ge 0\\
    \end{array}
  \right.
  \f$

  M is an (\f$nn \times nn\f$)  matrix , q , w and z nn-vectors.
*/
/*!\fn void lcp_newton_min( int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP )

 lcp_newton_min use a nonsmooth newton method based on the min formulation  (or max formulation) of the LCP

 \f$
   0 \le z \perp w \ge 0 \Longrightarrow \min(w,\rho z)=0 \Longrightarrow w = \max(0,w - \rho z)
 \f$

 \f$
   H(z) = H(\left[ \begin{array}{c} z \\ w \end{array}\right])= \left[ \begin{array}{c} w-Mz-q \\ min(w,\rho z) \end{array}\right] =0\\
 \f$


 References: Alart & Curnier 1990, Pang 1990


 \param nn      On enter, an integer which represents the dimension of the system.
 \param vec     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
 \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
 \param z       On return, a nn-vector of doubles which contains the solution of the problem.
 \param w       On return, a nn-vector of doubles which contains the solution of the problem.
 \param info    On return, an integer which returns the termination value:\n
                0 : convergence\n
                1 : iter = itermax\n
                2 : Problem in resolution in DGESV\n

 \param iparamLCP On enter/return a vectr of integers:\n
               - iparamLCP[0] = itermax On enter, the maximum number of iterations allowed.
               - iparamLCP[1] = verbose  On enter, the output log identifiant:\n
                       0 : no output\n
                       >0 : active screen output\n
               - iparamLCP[2] = it_end  On return, the number of iterations performed by the algorithm.

 \param dparamLCP On enter/return, a vector of doubles:\n
                 - dparamLCP[0] = tol     On enter the tolerance required.
                 - dparamLCP[1] = res     On return the final error value.

 \author Vincent Acary

 \todo Optimizing the memory allocation (Try to avoid the copy of JacH into A)
 \todo Add rules for the computation of the penalization rho
 \todo Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

void lcp_newton_min(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP)
{


  int i, j, iter;
  int n = *nn, m, mm, k;
  int itermax, verbose;

  integer  incx, incy;
  char NOTRANS = 'N';
  double err, tol, a1, b1;
  double alpha;
  int infoDGESV;

  int *ipiv;

  double *JacH, *H, *A;

  double *rho;

  incx = 1;
  incy = 1;
  /*input*/

  itermax = iparamLCP[0];
  verbose  = iparamLCP[1];

  tol   = dparamLCP[0];

  /*output*/

  iparamLCP[2] = 0;
  dparamLCP[1] = 0.0;

  for (i = 0; i < n; i++)
  {
    z[i] = 1.0;
    w[i] = 1.0;
  }

  /* rho*/
  rho = (double *)malloc(n * sizeof(double));
  for (i = 0; i < n; i++) rho[i] = 1.0 / vec[i * n + i] ;
  /* /for (i=0;i<n;i++) rho[i]=1.0/n ;
  // Sizw of the problem*/
  m = 2 * n;
  mm = m * m;
  /* / Creation of the gradient of the function H*/

  JacH = (double *)malloc(m * m * sizeof(double));
  A   = (double *)malloc(m * m * sizeof(double));

  for (j = 0; j < n; j++)
  {
    for (i = 0; i < n; i++) JacH[j * m + i] = -vec[j * n + i]; /* / should be replaced by a tricky use of BLAS*/
  }
  for (j = n; j < m; j++)
  {
    for (i = 0; i < n; i++) JacH[j * m + i] = 0.0;
    JacH[j * m + j - n] = 1.0;
  }
  for (j = 0; j < m; j++)
  {
    for (i = n; i < m; i++) JacH[j * m + i] = 0.0;
  }


  /* / Creation of the RHS H, */
  H = (double *)malloc(m * sizeof(double));
  /* / Construction of the RHS*/
  a1 = -1.;
  b1 = -1.;
  /* / q --> H*/
  dcopy_((integer *)&n , q , &incx , H , &incy);
  /* / -Mz-q --> H*/
  dgemv_(&NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , z , &incx , &b1 , H , &incy);
  /* / w+H --> H*/
  alpha = 1.0;
  daxpy_((integer *)&n , &alpha , w , &incx , H , &incy);     /* / c'est faux*/


  for (i = n; i < m; i++)
  {
    if (w[i - n] > rho[i - n]*z[i - n]) H[i] = rho[i - n] * z[i - n];
    else H[i] = w[i - n];
  }


  ipiv = (int *)malloc(m * sizeof(int));



  iter = 0;
  err  = 1.;



  while ((iter < itermax) && (err > tol))
  {
    ++iter;
    /* / Construction of the directional derivatives of H, JacH*/
    for (i = 0; i < n; i++)
    {
      if (w[i] > rho[i]*z[i])
      {
        JacH[i * m + i + n] = rho[i];
        JacH[(i + n)*m + (i + n)] = 0.0;
      }
      else
      {
        JacH[i * m + i + n] = 0.0;
        JacH[(i + n)*m + (i + n)] = 1.0;
      }
    }


    /* / Computation of the element of the subgradient.*/

    dcopy_((integer *)&mm , JacH , &incx , A , &incy);
    k = 1;
    F77NAME(dgesv)((integer *)&m, (integer *)&k, A, (integer *)&m, (integer *)ipiv, H, (integer *)&m, (integer *)&infoDGESV);

    if (infoDGESV)
    {
      if (verbose > 0)
      {
        printf("Problem in DGESV\n");
      }
      iparamLCP[2] = iter;
      dparamLCP[1] = err;

      free(H);
      free(A);
      free(JacH);
      free(ipiv);
      free(rho);
      *info = 2;
      return;

    }


    /* / iteration*/
    alpha = -1.0;
    daxpy_((integer *)&n , &alpha , H , &incx , z , &incy);      /* /  z-H --> z*/
    daxpy_((integer *)&n , &alpha , &H[n] , &incx , w , &incy);  /* /  w-H --> w*/

    /* / Construction of the RHS for the next iterate and for the error evalutaion*/
    a1 = 1.;
    b1 = 1.;
    dcopy_((integer *)&n , q , &incx , H , &incy);                                         /* / q --> H*/
    dgemv_(&NOTRANS , (integer *)&n , (integer *)&n , &a1 , vec , (integer *)&n , z , &incx , &b1 , H , &incy);  /* / Mz+q --> H*/
    alpha = -1.0;
    daxpy_((integer *)&n , &alpha , w , &incx , H , &incy);                                /* / w-Mz-q --> H*/

    for (i = n; i < m; i++)
    {
      if (w[i - n] > rho[i - n]*z[i - n]) H[i] = rho[i - n] * z[i - n];
      else H[i] = w[i - n];
    }

    /* / Error Evaluation*/

    // err = dnrm2_( (integer *)&m , H , &incx ); // necessary ?
    lcp_compute_error(n, vec, q, z, verbose, w, &err);


  }

  iparamLCP[2] = iter;
  dparamLCP[1] = err;


  if (err > tol)
  {
    printf("Siconos/Numerics: lcp_newton_min: No convergence of NEWTON_MIN after %d iterations\n" , iter);
    printf("Siconos/Numerics: lcp_newton_min: The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if (verbose > 0)
    {
      printf("Siconos/Numerics: lcp_newton_min: Convergence of NEWTON_MIN after %d iterations\n" , iter);
      printf("Siconos/Numerics: lcp_newton_min: The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(H);
  free(A);
  free(JacH);
  free(ipiv);
  free(rho);

}
