/*!\file rd_nlgs.c
 *
 *\fn  rd_nlgs(double vec[],double *qq,int *nn,double a[],int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)
 *
 * rd_nlgs is a specific nlgs (Gauss Seidel Non Linear)solver for relay problems.
 *
 * \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
 * \param qq On enter a pointer over doubles containing the components of the double vector.
 * \param nn On enter a pointer over integers, the dimension of the second member.
 * \param a On enter a pointer over doubles, the bound.
 * \param itermax On enter a pointer over integers, the maximum iterations required.
 * \param tol On enter a pointer over doubles, the tolerance required.
 * \param it_end On enter a pointer over integers, the number of iterations carried out.
 * \param res On return a pointer over doubles, the error value.
 * \param z On return double vector, the solution of the problem.
 * \param w On return double vector, the solution of the problem.
 * \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).
 *
 * \author Nineb Sheherazade & Mathieu Renouf.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

void dr_nlgs(double *vec , double *q , int *nn , double *a , double *b , int *itermax , double *tol ,
             double *z , double *w , int *it_end , double *res , int *info)
{

  FILE *f101;
  int i, j, iter, k, ispeak;
  int n = *nn, incx = 1, incy = 1;
  double alpha, beta, mina, a1, b1, zi;
  double err, num, den, avn, avt, apn, apt, xn, qs;
  double *zt, *diag, *ww;

  char NOTRANS = 'N';

  zt   = (double*) malloc(n * sizeof(double));
  ww   = (double*)malloc(n * sizeof(double));
  diag = (double*)malloc(n * sizeof(double));

  if (0) f101 = fopen("resultat_nlgs.dat", "w+");

  for (i = 0 ; i < n ; ++i)
  {
    w[i]  = 0.;
    zt[i] = 0.;
  }

  /* Preparation of the diagonal of the inverse matrix */

  for (i = 0 ; i < n ; ++i)
  {
    if (fabs(vec[i * n + i]) < 1e-16)
    {

      if (ispeak > 0)
      {
        printf(" Warning negative diagonal term \n");
        printf(" The local problem can be solved \n");
      }

      *info = 2;
      free(diag);
      free(ww);

      return;
    }
    else diag[i] = 1.0 / vec[i * n + i];
  }

  iter = 0;
  err  = 1.;
  dcopy_(&n , q , &incx , w , &incy);

  while ((iter < *itermax) && (err > *tol))
  {

    ++iter;

    dcopy_(&n , w , &incx , ww , &incy);
    dcopy_(&n , q , &incx , w , &incy);

    for (i = 0 ; i < n ; ++i)
    {

      incx = n;
      incy = 1;

      z[i] = 0.0;

      zi = -(q[i] + ddot_(&n , &vec[i] , &incx , z , &incy));

      zt[i] = -zi;

      mina = min(a[i] , z[i]);
      w[i] = max(mina , -b[i]);

      z[i] = diag[i] * (w[i] + zi);

    }

    /* **** Criterium convergence **** */

    incx =  1;
    incy =  1;

    a1 = 1.0;
    b1 = 1.0;

    dgemv_(&NOTRANS , &n , &n , &a1 , vec , &n , z , &incx , &b1 , w , &incy);

    qs   = -1.0;
    daxpy_(&n , &qs , w , &incx , ww , &incy);

    num = dnrm2_(&n, ww , &incx);
    err = num * den;

    if (0)
    {
      for (i = 0; i < n; i++) fprintf(f101, "%d  %d  %14.7e %14.7e\n", iter - 1, i, z[i], w[i]);
    }

  }


  if (err > *tol)
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

  free(ww);
  free(diag);
  free(zt);

  if (0) fclose(f101);

}
