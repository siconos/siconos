
/*!\file dr_gsnl.c


This subroutine allows the resolution of DR (Dual Relay) problem.
Try \f$(z,w)\f$ such that:

\f$
\left\lbrace
\begin{array}{l}
M z + q = w\\
-z \in \partial\psi_{[-b,a]}(w)\\
\end{array}
\right.
\f$

here M is an n by n  matrix, q an n-dimensional vector, w, z, a and b are n-dimensional vectors.

*/

/*!\fn  dr_gsnl(double vec[], double *q, int *nn, double a[], double b[], int * itermax, double * tol, int *chat, double z[], double w[], int *it_end, double * res, int *info)

   dr_gsnl is a specific nlgs (Non Linear Gauss Seidel) solver for dual relay problems.

   \param vec      On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param q        On enter a pointer over doubles containing the components of the double vector.
   \param nn       On enter a pointer over integers, the dimension of the second member.
   \param a        On enter a pointer over doubles, the upper bound.
   \param b        On enter a pointer over doubles, the lower bound.
   \param itermax  On enter a pointer over integers, the maximum iterations required.
   \param tol      On enter a pointer over doubles, the tolerance required.
   \param chat      On enter a pointer over integer, the output log identifiant
                    0 > =  no output
                    0 < =  active screen output


   \param it_end   On return a pointer over integers, the number of iterations carried out.
   \param res      On return a pointer over doubles, the error value.
   \param z        On return double vector, the solution of the problem.
   \param w        On return double vector, the solution of the problem.
   \param info     On return a pointer over integers, the termination reason
                    0 = convergence
        1 = no convergence
        2 = nul diagonal term


   \author Nineb Sheherazade.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"



void dr_nlgs(double vec[], double *q, int *nn, double a[], double b[], int * itermax, double * tol, int *chat, double z[], double w[], int *it_end, double * res, int *info)
{



  int i, j, iter1, k, ispeak = *chat;
  int n = *nn, incx = 1, incy = 1, itt = *itermax;

  double errmax = *tol, alpha, beta, mina;
  double err1, num, den, avn, xn, apn;
  double *zt, *wnum1;

  char trans = 'N';





  wnum1    = (double*) malloc(n * sizeof(double));
  zt       = (double*) malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    w[i]     = 0.;
    z[i]     = 0.;
    zt[i]    = 0.;
    wnum1[i] = 0.;
  }


  iter1 = 1;
  err1 = 1.;


  while ((iter1 < itt) && (err1 > errmax))
  {

    iter1 = iter1 + 1;

    for (i = 0; i < n; i++)
    {
      avn = 0.;
      apn = 0.;

      for (j = 0; j <= i - 1; j++)
        avn = avn + vec[j * n + i] * z[j];

      for (k = i + 1; k < n; k++)
        apn = apn + vec[k * n + i] * z[k];

      xn = -q[i] - avn - apn;

      zt[i] = -xn;

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
        w[i] = mina;
      }
      else
      {
        w[i] = -b[i];
      }

      if (fabs(vec[i * n + i]) < 1e-16)
      {
        printf("\n Warning nul diagonal term of M \n");

        free(zt);
        free(wnum1);

        *info = 2;

        return;

      }
      else
        z[i] = 1 / vec[i * n + i] * (w[i] + xn);


    }

    /*              Convergence criterium              */

    dcopy_(&n, w, &incx, wnum1, &incy);

    alpha = -1.;
    daxpy_(&n, &alpha, q, &incx, wnum1, &incy);

    alpha = 1.;
    beta = -1.;
    dgemv_(&trans, &n, &n, &alpha, vec, &n, z, &incx, &beta, wnum1, &incy);

    num = ddot_(&n, wnum1, &incx, wnum1, &incy);

    den = ddot_(&n, q, &incx, q, &incy);

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
