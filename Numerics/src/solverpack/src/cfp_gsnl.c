#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*!\file cfp_gsnl.c

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

*/



double ddot_(int *, double [], int *, double [], int*);


/*!\fn int cfp_gsnl(double vec[],double *qq,int *nn,double *mumu,int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)

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
\return int : \todo tell whiwh result is good
   \author Nineb Sheherazade.

 */
int cfp_gsnl(double vec[], double *qq, int *nn, double *mumu, int * itermax, double * tol, double z[], double w[], int *it_end, double * res, int *info)
{
  FILE *f13;
  double errmax = *tol, alpha, beta, mu = *mumu;
  double num1;
  int maxit = *itermax;
  double den1;
  double *q, *y;
  int n = *nn, incx, incy, nc = n / 2, i, j, k;
  double M[n][n];
  double fric1[nc], fric[nc], zn, zt, wt[nc];
  double yy[n], Az[n], Azz[n], detM[nc], ww[n], av_z[n], num, den;
  double zz[n][maxit], resvecgs[maxit + 1], normr, normM, epsM, eps;
  int status[nc], iter, ss[nc][maxit];
  char trans;
  double avn, avt, apn, apt, equi;


  iter = 0;
  eps = 1.e-08;

  q = (double*)malloc(n * sizeof(double));
  y = (double*)malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    z[i] = 0.0;
    w[i] = 0.0;
    Az[i] = 0.0;
    ww[i] = w[i];
    av_z[i] = z[i];
    q[i] = qq[i];
  }

  for (i = 0; i < nc; i++)
  {
    fric1[i] = 1.0;
    status[i] = 0;
    fric[i] = mu * fric1[i];
  }


  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      M[i][j] = vec[i * n + j];


  incx = 1;
  incy = 1;
  dcopy_(&n, w, &incx, y, &incy);

  trans = 'T';
  alpha = -1.;
  incx = 1;
  incy = 1;
  daxpy_(&n, &alpha, Az, &incx, y, &incy);
  alpha = 1.;
  daxpy_(&n, &alpha, q, &incx, y, &incy);

  num = ddot_(&n, y, &incx, y, &incy);
  den = ddot_(&n, q, &incx, q, &incy);

  normr = sqrt(num / den);

  resvecgs[1] = normr;

  //  normr=1.;
  //  normM=1.;
  //  espM=0.0001;



  while ((iter < maxit) && (normr > errmax))
  {

    iter = iter + 1;

    printf("iteration numbers %d and error evaluation %e \n ", iter, normr);


    for (i = 0; i < n; i++)
    {
      ww[i] = w[i];
      av_z[i] = z[i];
    }

    //  loop over contacts

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


      detM[i] = M[2 * i][2 * i] * M[2 * i + 1][2 * i + 1] - M[2 * i + 1][2 * i] * M[2 * i][2 * i + 1];



      z[2 * i] = 0.;
      w[2 * i] = -zn;
      z[2 * i + 1] = 0.;
      w[2 * i + 1] = -zt;

      if (w[2 * i] > eps) status[i] = 1; // no contact
      else
      {
        w[2 * i] = 0.0;
        w[2 * i + 1] = 0.0;
        z[2 * i] = (M[2 * i + 1][2 * i + 1] * zn - M[2 * i][2 * i + 1] * zt) / detM[i];
        z[2 * i + 1] = (-M[2 * i + 1][2 * i] * zn + M[2 * i][2 * i] * zt) / detM[i];


        if ((z[2 * i] > eps) && (fabs(z[2 * i + 1]) - mu * z[2 * i]) < eps)
        {
          status[i] = 2; // adherent status
        }
        else
        {
          z[2 * i] = zn / (M[2 * i][2 * i] + mu * M[2 * i][2 * i + 1]);
          z[2 * i + 1] = mu * z[2 * i];
          w[2 * i] = 0.;
          w[2 * i + 1] = -zt + ((M[2 * i][2 * i + 1] + mu * M[2 * i + 1][2 * i + 1]) / (M[2 * i][2 * i] + mu * M[2 * i][2 * i + 1])) * zn;

          if (z[2 * i] > eps && w[2 * i + 1] > eps)
          {
            status[i] = 3; // G+ status
          }
          else
          {
            z[2 * i] = zn / (M[2 * i][2 * i] - mu * M[2 * i][2 * i + 1]);
            w[2 * i] = 0.;
            z[2 * i + 1] = -mu * z[2 * i];
            w[2 * i + 1] = -zt + zn * ((M[2 * i][2 * i + 1] - mu * M[2 * i + 1][2 * i + 1]) / (M[2 * i][2 * i] + mu * M[2 * i][2 * i + 1]));
            status[i] = 4; // G- status
          }
        }
      }

    }


    // convergence criterium
    incx = 1;
    incy = 1;
    dcopy_(&n, z, &incx, y, &incy);
    trans = 'T';
    alpha = -1.;
    daxpy_(&n, &alpha, av_z, &incx, y, &incy);
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, y, &incx, &beta, y, &incy);
    dcopy_(&n, y, &incx, Azz, &incy);


    dcopy_(&n, ww, &incx, y, &incy);
    trans = 'T';
    alpha = -1.;
    daxpy_(&n, &alpha, w, &incx, y, &incy);


    num1 = ddot_(&n, y, &incx, y, &incy);
    den1 = ddot_(&n, w, &incx, w, &incy);

    normr = sqrt(num1 / den1);

    dcopy_(&n, z, &incx, y, &incy);
    trans = 'T';
    alpha = -1.;
    daxpy_(&n, &alpha, av_z, &incx, y, &incy);

    num1 = ddot_(&n, y, &incx, Azz, &incy);

    dcopy_(&n, z, &incx, y, &incy);
    trans = 'T';
    alpha = 1.;
    daxpy_(&n, &alpha, av_z, &incx, y, &incy);

    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, y, &incx, &beta, y, &incy);

    den1 = ddot_(&n, y, &incx, y, &incy);

    normM = sqrt(num1 / den1);

    resvecgs[iter] = normr;

    for (i = 0; i < n; i++)
      zz[i][iter] = z[i];

    for (i = 0; i < nc; i++)
      ss[i][iter] = status[i];
  }


  if (normr > errmax)
  {
    printf("no convergence after %d iterations, the residue is %g\n", iter, normr);
  }


  if (normr < errmax)
  {
    printf("there is convergence after %d iterations, the residue is %g \n", iter, normr);
  }

  it_end = &iter;

  /*  verification equilibre  */
  dcopy_(&n, q, &incx, y, &incy);
  trans = 'T';
  alpha = -1.;
  beta = 1.;
  dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, y, &incy);

  alpha = 1.;
  daxpy_(&n, &alpha, w, &incx, y, &incy);


  num1 = ddot_(&n, y, &incx, y, &incy);
  den1 = ddot_(&n, q, &incx, q, &incy);

  equi = sqrt(num1 / den1);

  printf("equilibrium is %g \n", equi);
  /* fin verification equilibre */

  if (normr < errmax) *info = 0;
  else *info = 1;

  free(q);
  free(y);
  return *info ;
}
