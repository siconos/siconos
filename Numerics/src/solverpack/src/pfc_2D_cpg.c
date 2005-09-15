/*!\file pfc_2D_cpg.c
 *
 *This subroutine allows the primal resolution of contact problems with friction.
 *
 * Try \f$(z,w)\f$ such that:
 *\f$
 *\left\lbrace
 *\begin{array}{l}
 *M z- w=q\\
 *0 \le z_n \perp w_n \ge 0\\
 *-w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
 *\end{array}
 *\right.
 *\f$
 *
 *here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.
 *
 *\fn  pfc_2D_cpg(double vec[],double *qq,int *nn,double *mu3,int * itermax, double * tol,double xout[],
 *                double rout[],int *it_end,double * res,int *info)
 *
 * pfc_2D_cpg is a specific gcp (gradient conjugated projected) solver for primal contact problem with friction.
 *
 * \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
 * \param qq On enter a pointer over doubles containing the components of the double vector.
 * \param nn On enter a pointer over integers, the dimension of the second member.
 * \param mu3 On enter a pointer over doubles, the friction coefficient.
 * \param itermax On enter a pointer over integers, the maximum iterations required.
 * \param tol On enter a pointer over doubles, the tolerance required.
 * \param it_end On enter a pointer over integers, the number of iterations carried out.
 * \param res On return a pointer over doubles, the error value.
 * \param xout On return double vector, the solution of the problem.
 * \param rout On return double vector, the solution of the problem.
 * \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).
 *
 * \author Nineb Sheherazade.
 * \coauthor
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"


pfc_2D_cpg(double vec[], double *b, int *nn, double *mu3, int * itermax, double * tol, double x[], double rout[], int *it_end, double * res, int *info)
{

  FILE *f101;
  int n = *nn, maxit = *itermax;
  double mu = *mu3, eps = 1.e-08;
  int nc = n / 2, i, j, iter, k, ii, incx = 1, incy = 1;
  double pAp, alpha, beta, wAp, rp, normr;
  char trans = 'T';
  double /*A[n][n],*/ alphaf, betaf, den, num;
  int *stat, *statusi;
  double *p, *fric;
  double *fric1, *v, *w, *Ap, *xi, *z;
  double /*r[n]*/ *r;
  double(*A)[n];

  A = malloc(n * n * sizeof(double));

  r = (double*) malloc(n * sizeof(double));
  fric = (double*) malloc(n * sizeof(double));
  p = (double*) malloc(n * sizeof(double));
  v = (double*) malloc(n * sizeof(double));
  w = (double*) malloc(n * sizeof(double));
  Ap = (double*) malloc(n * sizeof(double));
  xi = (double*) malloc(n * sizeof(double));
  z = (double*) malloc(n * sizeof(double));
  fric1 = (double*) malloc(n * sizeof(double));
  stat = (int*) malloc(nc * sizeof(int));
  statusi = (int*) malloc(nc * sizeof(int));

  f101 = fopen("resultat_gcp.dat", "w+");

  for (i = 0; i < n; i++)
  {
    x[i] = 0.;
    xi[i] = 0.;
    r[i] = 0.;
    v[i] = 0.;
    p[i] = 0.;
    w[i] = 0.;
    Ap[i] = 0.;
    z[i] = 0.;
    fric1[i] = 1.;
    fric[i] = mu * fric1[i];
  }

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = vec[i * n + j];

  for (i = 0; i < nc; i++)
  {
    stat[i] = 0;
    statusi[i] = 0;
  }


  dcopy_(&n, b, &incx, r, &incy);
  alphaf = -1.;
  betaf = 1.;
  dgemv_(&trans, &n, &n, &alphaf, A, &n, x, &incx, &betaf, r, &incy);

  iter = 1;


  /*//  !Check for initial status*/
  for (i = 0; i < nc; i++)
  {
    mu = fric[i];
    if (x[2 * i] <= eps)
    {
      /*// !no contact*/
      stat[i] = 0;
    }
    else if (x[2 * i + 1] <= -mu * x[2 * i])
    {
      /*/ !slide backward*/
      stat[i] = 1;
    }
    else if (x[2 * i + 1] >= mu * x[2 * i])
    {
      /*/  !slide forward*/
      stat[i] = 3;
    }
    else
    {
      /*/ !stick contact*/
      stat[i] = 2;
    }
  }

  for (ii = 0; ii < maxit; ii++)
  {
    for (i = 0; i < nc; i++)
      statusi[i] = stat[i];


    dcopy_(&n, r, &incx, v, &incy);

    if (ii == 0)
    {
      dcopy_(&n, r, &incx, w, &incy);
      dcopy_(&n, w, &incx, p, &incy);
      alphaf = 1.;
      betaf = 0.;
      dgemv_(&trans, &n, &n, &alphaf, A, &n, p, &incx, &betaf, Ap, &incy);
      pAp = ddot_(&n, p, &incx, Ap, &incy);
    }
    else
    {
      alphaf = 1.;
      betaf = 0.;
      dgemv_(&trans, &n, &n, &alphaf, A, &n, p, &incx, &betaf, Ap, &incy);
      pAp = ddot_(&n, p, &incx, Ap, &incy);

      if (pAp == 0)
      {
        printf("operation non conform alpha at the iteration %d \n", iter);
        break;
      }

    }

    rp = ddot_(&n, r, &incx, p, &incy);
    alpha = rp / pAp;


    dcopy_(&n, x, &incx, xi, &incy);
    alphaf = alpha;
    daxpy_(&n, &alphaf, p, &incx, xi, &incy);

    projc(xi, &n, statusi, p, fric, x, stat);


    /*/   r(:)=b(:)-matmul(A,x)*/

    dcopy_(&n, b, &incx, r, &incy);
    alphaf = -1.;
    betaf = 1.;
    dgemv_(&trans, &n, &n, &alphaf, A, &n, x, &incx, &betaf, r, &incy);

    projf(statusi, &n, r, fric, w);
    projf(statusi, &n, p, fric, z);


    wAp = ddot_(&n, w, &incx, Ap, &incy);
    beta = - wAp / pAp;

    dcopy_(&n, w, &incx, p, &incy);
    alphaf = beta;
    betaf = 1.;
    daxpy_(&n, &alphaf, z, &incx, p, &incy);


    alphaf = 1.;
    betaf = 0.;
    dgemv_(&trans, &n, &n, &alphaf, A, &n, p, &incx, &betaf, Ap, &incy);

    pAp = ddot_(&n, p, &incx, Ap, &incy);

    dcopy_(&n, r, &incx, xi, &incy);
    alphaf = -1.;
    daxpy_(&n, &alphaf, v, &incx, xi, &incy);

    num = ddot_(&n, xi, &incx, xi, &incy);
    den = ddot_(&n, v, &incx, v, &incy);

    normr = sqrt(num / den);

    for (i = 0; i < n; i++)
    {
      /*result_gs[i][iter1-1] = z[i]; */
      fprintf(f101, "%d  %d  %14.7e\n", ii, i, x[i]);
    }

    if (normr < *tol)
    {
      iter = ii;
      /*printf("convergence after %d iterations with a residual %g\n",iter,normr);*/
      *info = 0;
      break;
    }
    iter = iter + 1;



  }




  if (normr < *tol)
  {
    *info = 0;
    printf("convergence after %d iterations with a residual %g\n", iter, normr);
  }
  else
  {
    *info = 1;
    printf("no convergence after %d iterations with a residual %g\n", iter, normr);
  }


  *it_end = iter;
  *res = normr;

  for (i = 0; i < n; i++)
    rout[i] = -r[i];


  free(fric);
  free(p);
  free(v);
  free(w);
  free(Ap);
  free(xi);
  free(z);
  free(fric1);
  free(stat);
  free(statusi);
  free(r);
  free(A);

  fclose(f101);

}
