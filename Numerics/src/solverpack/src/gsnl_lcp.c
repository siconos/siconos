#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*!\file gsnl_lcp.c
 *
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).
 * Try \f$(z,w)\f$ such that:
 * \f$
 * \left\lbrace
 * \begin{array}{l}
 * M z- w= -q\\
 * 0 \le z \perp w \ge 0\\
 * \end{array}
 * \right.
 * \f$
 *
 * where M is an (n x n)-matrix, q , w and z n-vectors.
 *
 *
 * \fn  gsnl_lcp(double vec[],double *qq,int *nn,int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)
 *
 * gsnl_lcp is a solver for LCP based on the principle of splitting method.
 * (Non Linear Gauss-Seidel)
 *
 * \param vec On enter a pointer over doubles containing the components of the double matrix with a fortran90 allocation.
 * \param qq On enter a pointer over doubles containing the components of the double vector.
 * \param nn On enter a pointer over integers, the dimension of the second member.
 * \param itermax On enter a pointer over integers, the maximum iterations required.
 * \param tol On enter a pointer over doubles, the tolerance required.
 * \param it_end On enter a pointer over integers, the number of iterations carried out.
 * \param res On return a pointer over doubles, the error value.
 * \param z On return double vector, the solution of the problem.
 * \param w On return double vector, the solution of the problem.
 * \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).
 *
 * \author Nineb Sheherazade.
 * \author Last Modif: Mathieu Renouf
 *
 * ===========================================================================
 * Prototypes for level 1 BLAS functions
 * ===========================================================================
 */

double dnrm2_(int* , double* , int*);
double ddot_(int* , double* , int* , double* , int*);

/*
 */

void gsnl_lcp(double *vec , double *qq , int *nn , int *itermax , double *tol , double *z ,
              double *w , int *it_end , double *res , int *info)
{


  int n = *nn;
  int itt = *itermax;
  double errmax = *tol;

  int i , j , iter;

  int incx, incy;
  double qs , err , num, den, zi;

  double *ww;

  /* Note:
   * vec has a F90 storage
   */

  ww = (double*)malloc(n * sizeof(double));

  /* Check for trivial case */
  /* Usefull or not ?*/
  /* Note:
   * If M is PSD then the trivial solution if one of
   * the admissible ones. In this case the test can
   * be avoid to obtain an other admissible one.
   */

  /*   incx = 1; */
  /*   qs = idamax_( &n , qq , &incx ); */

  /*   if( qs <= 0.0 ){ */
  /*     for( i = 0 ; i < n ; ++i ){ */
  /*       w[i] = 0.; */
  /*       z[i] = 0.; */
  /*     } */
  /*     info = 0; */
  /*     return; */
  /*   } */

  /* Non trivial case */

  incx = 1;
  qs = dnrm2_(&n , &qq[0] , &incx);

  if (qs > 1e-16) den = 1.0 / qs;
  else
  {
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }
    info = 0;
    return;
  }

  for (i = 0 ; i < n ; ++i)
  {
    ww[i] = 0.;
    w[i] = 0.;
    /*z[i] = 0.; Start from an initial vector */
  }

  iter = 0;
  err  = 1.;

  while ((iter < itt) && (err > errmax))
  {

    ++iter;

    incx = 1;
    incy = 1;

    dcopy_(&n , w , &incx , ww , &incy);

    for (i = 0 ; i < n ; ++i)
    {

      incx = n;
      incy = 1;

      z[i] = 0.0;

      zi = (qq[i] - ddot_(&n , &vec[i] , &incx , z , &incy)) / vec[i * n + i];

      if (zi < 0) z[i] = 0.0;
      else z[i] = zi;

      w[i] = -zi +  z[i] * vec[i * n + i];

    }

    /* **** Criterium convergence **** */

    qs   = -1;
    incx =  1;
    incy =  1;

    daxpy_(&n , &qs , w , &incx , ww , &incy);

    num = dnrm2_(&n, ww , &incx);
    err = num * den;

    /* **** ********************* **** */

  }

  *it_end = iter;
  *res    = err;

  /* Compute w */

  if (err > errmax)
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

}
