#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

/*!\file lcp_newton.c
 *
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).
 * Try \f$(z,w)\f$ such that:
 *
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *   0 \le z \perp M z + b = w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 *
 * M is an (n x n)  matrix , q , w and z n-vectors.
 *
 *
 *!\fn lcp_newtonmin(  double *vec , double *q , int *nn , int *itermax , double *tol , int *ispeak , double *z , double *w , int *it_end , double *res , int *info )
 *
 * lcp_newton_min use a nonsmooth newton method based on the min formulation  (or max formulation) of the LCP
 *
 * \f$
 *   0 \le z \perp w \ge 0 \Longrightarrow \min(w,\rho z)=0 \Longrightarrow w = \max(0,w - \rho z)
 * \f$

 * \f$
 *   H(z) = H(\left[ \begin{array}{c} z \\ w \end{array}\right])= \left[ \begin{array}{c} w-Mz-q \\ min(w,\rho z) \end{array}\right] =0\\
 * \f$
 *
 *
 * References: Alart & Curnier 1990, Pang 1990
 *
 * Generic lcp parameters:\n
 *
 * \param nn      Unchanged parameter which represents the dimension of the system.
 * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
 * \param q       Unchanged parameter which contains the components of the right hand side vector.
 * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param w       Modified parameter which returns the solution of the problem.
 * \param info    Modified parameter which returns the termination value\n
 *                0 - convergence\n
 *                1 - iter = itermax\n
 *                2 - negative diagonal term\n
 *
 * Specific NLGS parameters:\n
 *
 * \param iparamLCP[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
 * \param iparamLCP[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
 *                       0 - no output\n
 *                       0 < active screen output\n
 * \param iparamLCP[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
 *
 * \param dparamLCP[0] = tol     Input unchanged parameter which represents the tolerance required.
 * \param dparamLCP[1] = res     Output modified parameter which returns the final error value.
 *
 * \author Vincent Acary
 *
 * \todo Optimizing the memory allocation (Try to avoid the copy of JacH into A)
 * \todo Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002
 */

void lcp_newton_min(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                    int *iparamLCP , double *dparamLCP)
{


  int i, j, iter;
  int n = *nn, m, mm, k;
  int itermax, ispeak;

  int incx, incy;
  char NOTRANS = 'N';
  double err, tol, a1, b1;
  double alpha;
  int infoDGESV;

  int *ipiv;

  double *JacH, *H, *A;

  double rho;

  incx = 1;
  incy = 1;
  /*input*/

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  tol   = dparamLCP[0];

  /*output*/

  iparamLCP[2] = 0;
  dparamLCP[1] = 0.0;

  for (i = 0; i < n; i++)
  {
    z[i] = 1.0;
    w[i] = 1.0;
  }

  // rho
  rho = 1.0 / n;
  // Sizw of the problem
  m = 2 * n;
  mm = m * m;
  // Creation of the gradient of the function H

  JacH = (double *)malloc(m * m * sizeof(double));
  A   = (double *)malloc(m * m * sizeof(double));

  for (j = 0; j < n; j++)
  {
    for (i = 0; i < n; i++) JacH[j * m + i] = -vec[j * n + i]; // should be replaced by a tricky use of BLAS
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


  // Creation of the RHS H,
  H = (double *)malloc(m * sizeof(double));
  // Construction of the RHS
  a1 = -1.;
  b1 = -1.;
  // q --> H
  dcopy_(&n , q , &incx , H , &incy);
  // -Mz-q --> H
  dgemv_(&NOTRANS , &n , &n , &a1 , vec , &n , z , &incx , &b1 , H , &incy);
  // w+H --> H
  alpha = 1.0;
  daxpy_(&n , &alpha , w , &incx , H , &incy);     // c'est faux


  for (i = n; i < m; i++)
  {
    if (w[i - n] > rho * z[i - n]) H[i] = rho * z[i - n];
    else H[i] = w[i - n];
  }


  ipiv = (int *)malloc(m * sizeof(int));

  //

  iter = 0;
  err  = 1.;



  while ((iter < itermax) && (err > tol))
  {
    ++iter;
    // Construction of the directional derivatives of H, JacH
    for (i = 0; i < n; i++)
    {
      if (w[i] > rho * z[i])
      {
        JacH[i * m + i + n] = rho;
        JacH[(i + n)*m + (i + n)] = 0.0;
      }
      else
      {
        JacH[i * m + i + n] = 0.0;
        JacH[(i + n)*m + (i + n)] = 1.0;
      }
    }


    // Computation of the element of the subgradient.

    dcopy_(&mm , JacH , &incx , A , &incy);
    k = 1;
    F77NAME(dgesv)(&m, &k, A, &m, ipiv, H, &m, &infoDGESV);

    if (infoDGESV)
    {
      if (ispeak > 0)
      {
        printf("Problem in DGESV\n");
      }
      iparamLCP[2] = iter;
      dparamLCP[1] = err;

      free(H);
      free(A);
      free(JacH);
      free(ipiv);

      return (*info = 2);

    }


    // iteration
    alpha = -1.0;
    daxpy_(&n , &alpha , H , &incx , z , &incy);     //  z-H --> z
    daxpy_(&n , &alpha , &H[n] , &incx , w , &incy);  //  w-H --> w

    // Construction of the RHS for the next iterate and for the error evalutaion
    a1 = 1.;
    b1 = 1.;
    dcopy_(&n , q , &incx , H , &incy);                                         // q --> H
    dgemv_(&NOTRANS , &n , &n , &a1 , vec , &n , z , &incx , &b1 , H , &incy);  // Mz+q --> H
    alpha = -1.0;
    daxpy_(&n , &alpha , w , &incx , H , &incy);                               // w-Mz-q --> H

    for (i = n; i < m; i++)
    {
      if (w[i - n] > rho * z[i - n]) H[i] = rho * z[i - n];
      else H[i] = w[i - n];
    }

    // Error Evaluation

    err = dnrm2_(&m , H , &incx);

  }

  iparamLCP[2] = iter;
  dparamLCP[1] = err;

  if (ispeak > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of NLGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 1;
    }
    else
    {
      printf(" Convergence of NLGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else
  {
    if (err > tol) *info = 1;
    else *info = 0;
  }
  free(H);
  free(A);
  free(JacH);
  free(ipiv);

}
