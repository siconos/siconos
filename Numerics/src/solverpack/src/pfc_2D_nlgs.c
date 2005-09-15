/*!\file pfc_2D_nlgs.c

  This subroutine allows the primal resolution of contact problems with friction.

   Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
0 \le z_n \perp w_n \ge 0\\
-w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

\fn  pfc_2D_nlgs(double vec[],double *qq,int *nn,double *mumu,int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)

   cfp_gsnl is a specific gsnl (Gauss Seidel Non Linear) solver for primal contact problem with friction.

   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param qq On enter a pointer over doubles containing the components of the double vector.
   \param nn On enter a pointer over integers, the dimension of the second member.
   \param mumu On enter a pointer over doubles, the friction coefficient.
   \param itermax On enter a pointer over integers, the maximum iterations required.
   \param tol On enter a pointer over doubles, the tolerance required.
   \param it_end On enter a pointer over integers, the number of iterations carried out.
   \param res On return a pointer over doubles, the error value.
   \param z On return double vector, the solution of the problem.
   \param w On return double vector, the solution of the problem.
   \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).

   \author Nineb Sheherazade.

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

pfc_2D_nlgs(double vec[], double *q, int *nn, double *mumu, int * itermax, double * tol, double z[], double w[], int *it_end, double * res, int *info)
{

  FILE *f101;
  double errmax = *tol, alpha, beta, mu = *mumu;
  int n = *nn, incx = 1, incy = 1, nc = n / 2, i, j, k, iter, maxit = *itermax;
  double /*M[n][n],*/ *y;
  /*  double fric1[nc], fric[nc], ww[n];*/
  double *fric1, *fric, *ww;
  double normr, eps, avn, avt, apn, apt, zn , zt, den1, num1;
  char trans = 'T';
  double(*M)[n];


  M = malloc(n * n * sizeof(double));



  f101 = fopen("resultat_gsnl.dat", "w+");

  iter = 0;
  eps = 1.e-08;


  y = (double*) malloc(n * sizeof(double));
  fric1 = (double*) malloc(nc * sizeof(double));
  fric = (double*) malloc(nc * sizeof(double));
  ww = (double*) malloc(n * sizeof(double));

  /*  for (i = 0; i < n; i++)
  {
      z[i] = 0.0;
      w[i] = 0.0;
      ww[i] = 0.0;
      }*/

  for (i = 0; i < nc; i++)
  {
    fric1[i] = 1.0;
    fric[i] = mu * fric1[i];
  }


  for (i = 0; i < n; i++)
  {
    z[i] = 0.0;
    w[i] = 0.0;
    ww[i] = 0.0;
    for (j = 0; j < n; j++)
      M[i][j] = vec[i * n + j];
  }

  normr = 1.;

  while ((iter < maxit) && (normr > errmax))
  {

    iter = iter + 1;

    /*//  loop over contacts*/

    for (i = 0; i < nc; i++)
    {
      avn = 0.;
      avt = 0.;
      apn = 0.;
      apt = 0.;

      for (j = 0; j <= 2 * i - 1; j++)
      {
        avn += M[2 * i][j] * z[j];
        avt += M[2 * i + 1][j] * z[j];
      }

      for (k = 2 * i + 2; k < n; k++)
      {
        apn = apn + M[2 * i][k] * z[k];
        apt = apt + M[2 * i + 1][k] * z[k];
      }


      zn = q[2 * i] - avn - apn;
      zt = q[2 * i + 1] - avt - apt;

      if (zn > eps)
      {
        z[2 * i] = zn / M[2 * i][2 * i];
        w[2 * i] = 0.;
        z[2 * i + 1] = zt / M[2 * i + 1][2 * i + 1];
        w[2 * i + 1] = 0.;
        if (z[2 * i + 1] > fric[i]*z[2 * i])
        {
          z[2 * i + 1] = fric[i] * z[2 * i];
          w[2 * i + 1] = -zt + M[2 * i + 1][2 * i + 1] * z[2 * i + 1];
        }
        else if (z[2 * i + 1] < -fric[i]*z[2 * i])
        {
          z[2 * i + 1] = -fric[i] * z[2 * i];
          w[2 * i + 1] = -zt + M[2 * i + 1][2 * i + 1] * z[2 * i + 1];
        }
      }
      else
      {
        z[2 * i] = 0.0;
        w[2 * i] = -zn;
        z[2 * i + 1] = 0.;
        w[2 * i + 1] = -zt;

      }
    }


    dcopy_(&n, q, &incx, y, &incy);
    alpha = 1.;
    beta = -1.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, ww, &incy);


    /*    // convergence criterium*/

    alpha = -1.;
    daxpy_(&n, &alpha, w, &incx, y, &incy);


    num1 = ddot_(&n, y, &incx, y, &incy);
    /*   den1 = ddot_(&n, w, &incx, w, &incy);*/
    den1 = ddot_(&n, q, &incx, q, &incy);

    normr = sqrt(num1 / den1);
    *it_end = iter;
    *res = normr;
    for (i = 0; i < n; i++)
    {
      /*result_gs[i][iter1-1] = z[i]; */
      fprintf(f101, "%d  %d  %14.7e\n", iter - 1, i, z[i]);
    }


  }





  if (normr > errmax)
  {
    printf("no convergence after %d iterations, the residue is %g\n", iter, normr);
    *info = 1;
  }
  else
  {
    printf("there is convergence after %d iterations, the residue is %g \n", iter, normr);
    *info = 0;
  }

  free(y);
  free(M);
  free(fric1);
  free(fric);
  free(ww);


  fclose(f101);
}
