#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*!\file gcp_lcp.c


   This subroutine allows the resolution of LCP (Linear Complementary Problem).
   Try \f$(z,w)\f$ such that:

   \f$
    \left\lbrace
     \begin{array}{l}
      M z- w=q\\
      0 \le z \perp w \ge 0\\
     \end{array}
    \right.
   \f$

  here M is an n by n  matrix, q an n-dimensional vector, w an n-dimensional  vector and z an n-dimensional vector.
*/

/*!\fn  gcp_lcp(double vec[],double *q,int *nn,int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)

   gcp_lcp is a basic gcp (gradient conjugated projected) solver for LCP.

   \param vec On enter a pointer over doubles containing the components of the double matrix with a fortran90 allocation.
   \param q On enter a pointer over doubles containing the components of the double vector.
   \param nn On enter a pointer over integers, the dimension of the second member.
   \param itermax On enter a pointer over integers, the maximum iterations required.
   \param tol On enter a pointer over doubles, the tolerance required.
   \param it_end On enter a pointer over integers, the number of iterations carried out.
   \param res On return a pointer over doubles, the error value.
   \param z On return double vector, the solution of the problem.
   \param w On return double vector, the solution of the problem.
   \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).

 * \author Nineb Sheherazade.
 * \author Last Modif: Mathieu Renouf
 *
 * ===========================================================================
 * Prototypes for level 1 BLAS functions
 * ===========================================================================
 */

double dnrm2_(int* , double* , int*);
double  ddot_(int* , double* , int* , double* , int*);
void   dcopy_(int* , double* , int* , double* , int*);
void   daxpy_(int* , double* , double* , int* , double* , int*);
void   dgemv_(char* , int* , int* , double* , double* , int* , double* , int* , double* , double* , int*);

void gcp_lcp(double *vec , double *q , int *nn , int *itermax , double *tol , double *z ,
             double *vv , int *it_end , double * res , int *info)
{


  int n = *nn, incx , incy, itt = *itermax;
  double errmax = *tol;

  int i , iter;
  double err, a1, b1 , qs ;

  double alpha , beta , rp , pMp;
  double den, num;

  char NOTRANS = 'N';

  int *status;
  double *zz , *dz , *pp , *rr, *ww, *Mp;

  status = (int*)malloc(n * sizeof(int));

  dz = (double*)malloc(n * sizeof(double));

  ww = (double*)malloc(n * sizeof(double));
  rr = (double*)malloc(n * sizeof(double));
  pp = (double*)malloc(n * sizeof(double));
  zz = (double*)malloc(n * sizeof(double));

  Mp = (double*)malloc(n * sizeof(double));

  incx = 1;
  qs = dnrm2_(&n , &q[0] , &incx);

  //printf( " Norm: %g \n", qs );

  if (qs > 1e-16) den = 1.0 / qs;
  else
  {
    for (i = 0 ; i < n ; ++i)
    {
      vv[i] = 0.;
      z[i] = 0.;
    }
    info = 0;
    return;
  }

  for (i = 0 ; i < n ; ++i)
  {
    ww[i] = 0.;
    /*z[i] = 0.; Start from an initial vector */
  }

  for (i = 0; i < n; ++i)
  {

    status[i] = 0;

    ww[i] = 0.;
    rr[i] = 0.;
    pp[i] = 0.;
    zz[i] = 0.;

    Mp[i] = 0.;

  }

  /* rr = -Wz + q */
  incx = 1;
  incy = 1;

  dcopy_(&n , q , &incx , rr , &incy);

  a1 = -1.;
  b1 =  1.;

  dgemv_(&NOTRANS , &n , &n , &a1 , vec , &n , z , &incx , &b1 , rr , &incy);

  /* Initialization of gradients */
  /* rr -> p and rr -> w */

  dcopy_(&n , rr , &incx , ww , &incy);
  dcopy_(&n , rr , &incx , pp , &incy);

  iter = 0.0;
  err  = 1.0 ;

  while ((iter < itt) && (err > errmax))
  {

    ++iter;

    /* Compute initial pMp */

    incx = 1;
    incy = 1;

    dcopy_(&n , pp , &incx , Mp , &incy);

    a1 = 1.0;
    b1 = 0.0;

    dgemv_(&NOTRANS , &n , &n , &a1 , vec , &n , Mp , &incx , &b1 , vv , &incy);

    pMp = ddot_(&n , pp , &incx , vv  , &incy);

    if (fabs(pMp) < 1e-16)
    {
      printf(" Operation no conform at the iteration %d \n", iter);
      printf(" Alpha can be obtained with pWp = %10.4g  \n", pMp);
      return (*info = 3);
    }

    rp  = ddot_(&n , pp , &incx , rr , &incy);

    alpha = rp / pMp;

    /*
     * Iterate prediction:
     * z' = z + alpha*p
     *
     */

    dcopy_(&n , z , &incx , dz , &incy);

    daxpy_(&n , &alpha , pp , &incx , z , &incy);

    /* Iterate projection*/

    for (i = 0; i < n; ++i)
    {
      if (z[i] > 0.0)
      {
        status[i] = 1;
      }
      else
      {
        z[i] = 0.0;
        status[i] = 0;
      }
    }

    /* rr = -Wz + q */

    dcopy_(&n , q , &incx , rr , &incy);

    a1 = -1.;
    b1 =  1.;

    dgemv_(&NOTRANS , &n , &n , &a1 , vec , &n , z , &incx , &b1 , rr , &incy);

    /* Gradients projection
     * rr --> ww
     * pp --> zz
     */

    for (i = 0; i < n; ++i)
    {

      if (status[i])
      {
        ww[i] = rr[i];
        zz[i] = pp[i];
      }
      else
      {
        if (rr[i] < 0)
        {
          ww[i] = 0.0;
          zz[i] = 0.0;
        }
        else
        {
          ww[i] = rr[i];
          if (pp[i] < 0) zz[i] = 0.0;
          else zz[i] = pp[i];
        }
      }
    }

    /*   beta = -w.Mp / pMp  */

    rp = ddot_(&n , ww , &incx, vv , &incy);

    beta = -rp / pMp;

    dcopy_(&n , ww , &incx , pp , &incy);
    daxpy_(&n, &beta , zz , &incx , pp , &incy);

    a1 = -1;
    daxpy_(&n , &a1 , z , &incx , dz , &incy);
    num = dnrm2_(&n , dz , &incx);
    err = num * den;

  }

  *it_end = iter;
  *res    = err;

  /*    printf("iteration numbers %d and error evaluation %g \n ",iter,err);   */

  if (err > errmax)
  {
    printf(" No convergence of CPG after %d iterations\n" , iter);
    printf(" The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    printf(" Convergence of CPG after %d iterations\n" , iter);
    printf(" The residue is : %g \n", err);
    *info = 0;
  }

  free(Mp);

  free(ww);
  free(rr);
  free(pp);
  free(zz);

  free(dz);

}
