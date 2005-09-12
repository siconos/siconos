#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*!\file latin_lcp.c


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


/*!\fn  latin_lcp(double vec[],double *qq,int *nn, double * k_latin,int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)

   latin_lcp is a basic latin solver for LCP.


   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param qq On enter a pointer over doubles containing the components of the double vector.
   \param nn On enter a pointer over integers, the dimension of the second member.
   \param k_latin On enter a pointer over doubles, the k_latin coefficient (positive).
   \param itermax On enter a pointer over integers, the maximum iterations required.
   \param tol On enter a pointer over doubles, the tolerance required.
   \param it_end On enter a pointer over integers, the number of iterations carried out.
   \param res On return a pointer over doubles, the error value.
   \param z On return double vector, the solution of the problem.
   \param w On return double vector, the solution of the problem.
   \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).
   \author Nineb Sheherazade.
*/

void lcp_latin(int *nn , double *vec , double *qq , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP)
{

  FILE *f101;

  int i, j, kk, iter, info2, nrhs;
  int n, iout, itermax, tol, k_latin;

  int incx = 1, incy = 1;

  double alpha, beta;
  double rr, rrr, r1, r2, r3, invR0, invRT0, err, z0, num11, err0, invRTinvR0;
  double den11, den22, vv;

  double  *wc, *zc, *kinvden1, *kinvden2, *wt, *maxwt, *wnum1, *znum1, *ww, *zz;
  double *num1, *kinvnum1, *den1, *den2, *wden1, *zden1;

  char trans = 'T', uplo = 'U', diag = 'N', trans2 = 'N', uplo2 = 'U';
  double  *kinvwden1, *kzden1;

  double  **k, **kinv, **DPO;

  /*input*/

  itermax = iparamLCP[0];
  iout    = iparamLCP[1];
  tol     = dparamLCP[0];
  k_latin = dparamLCP[1];

  /* *** */

  n = *nn;

  /* ********* */

  if (iout > 0) f101 = fopen("resultat_latin.dat", "w+");


  /*** Allocations ***/

  ww = (double*) malloc(n * sizeof(double));
  zz = (double*) malloc(n * sizeof(double));
  wc = (double*) malloc(n * sizeof(double));
  zc = (double*) malloc(n * sizeof(double));
  wt = (double*) malloc(n * sizeof(double));

  znum1 = (double*) malloc(n * sizeof(double));
  wnum1 = (double*) malloc(n * sizeof(double));
  maxwt = (double*) malloc(n * sizeof(double));

  kinvden1 = (double*) malloc(n * sizeof(double));
  kinvden2 = (double*) malloc(n * sizeof(double));

  num1 = (double*) malloc(n * sizeof(double));

  kinvnum1  = (double*) malloc(n * sizeof(double));
  kinvwden1 = (double*) malloc(n * sizeof(double));

  den1 = (double*) malloc(n * sizeof(double));
  den2 = (double*) malloc(n * sizeof(double));
  wden1 = (double*) malloc(n * sizeof(double));
  zden1 = (double*) malloc(n * sizeof(double));

  kzden1 = (double*) malloc(n * sizeof(double));

  DPO = (double**) malloc(n * sizeof(double*));
  for (i = 0 ; i < n ; ++i) DPO[i] = (double *)malloc(n * sizeof(double));

  k = (double**) malloc(n * sizeof(double*));
  for (i = 0 ; i < n ; ++i) k[i] = (double *)malloc(n * sizeof(double));

  kinv = (double**) malloc(n * sizeof(double*));
  for (i = 0 ; i < n ; ++i) kinv[i] = (double *)malloc(n * sizeof(double));

  /* Initialization */

  for (i = 0; i < n; i++)
  {

    wc[i] = 0.0;
    zc[i] = 0.;
    z[i] = 0.;
    w[i] = 0.;
    znum1[i] = 0.;
    wnum1[i] = 0.;
    kinvden1[i] = 0.;
    kinvden2[i] = 0.;
    wt[i] = 0.;
    maxwt[i] = 0.;
    num1[i] = 0.;
    kinvnum1[i] = 0.;
    den1[i] = 0.;
    den2[i] = 0.;


    for (j = 0 ; j < n ; j++)
    {
      k[i][j] = 0.;
      kinv[i][j] = 0.;
      DPO[i][j] = 0.;
    }

    k[i][i] = k_latin * vec[i * n + i];

  }

  /* get C storage format for matrix*/

  for (i = 0 ; i < n ; i++)
  {

    for (j = 0 ; j < n ; j++) DPO[i][j] = vec[j * n + i] + k[i][j];

    kinv[i][i] = 1 / k[i][i];

  }

  /* **** Call Cholesky ***** */

  dpotrf_(&uplo, &n, DPO , &n, &info2);

  /*  printf("success of cholesky? %d\n",info2);*/

  if (info2 != 0)
  {
    printf("\nmatter with cholesky info2 %d\n", info2);
    return (*info = 2);
  }

  /* **** End of Cholesky ***** */

  /*    //  !iteration loops*/

  iter = 0;
  err  = 1.;

  alpha = -1.;
  dscal_(&n, &alpha , qq, &incx);

  while ((iter < itermax) && (err > tol))
  {

    /*      //   !linear stage (zc,wc) -> (z,w)*/

    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, zc, &incx, &beta, /*wt*/ wc, &incy);

    dcopy_(&n, qq, &incx, znum1, &incy);

    alpha = 1.;
    daxpy_(&n, &alpha, wc, &incx, znum1, &incy);

    nrhs = 1;
    dtrtrs_(&uplo, &trans, &diag, &n, &nrhs, DPO, &n, znum1, &n, &info2);

    trans2 = 'N';
    uplo2 = 'U';
    dtrtrs_(&uplo2, &trans2, &diag, &n, &nrhs, DPO, &n, znum1, &n, &info2);
    dcopy_(&n, znum1, &incx, z, &incy);

    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, z, &incx, &beta, wc, &incy);

    dcopy_(&n, wc, &incx, w, &incy);

    /*Local Stage*/

    dcopy_(&n, w, &incx, wt, &incy);
    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, z, &incx, &beta, wt, &incy);

    for (i = 0; i < n; i++)
    {
      if (wt[i] > 0.0) wc[i] = wt[i];
      else wc[i] = 0.0;
    }


    for (i = 0; i < n; i++)
    {
      if (-wt[i] < 0.0) zc[i] = 0.;
      else zc[i] = -wt[i] / k[i][i];
    }


    /*      // convergence criterium */


    dcopy_(&n, w, &incx, wnum1, &incy);
    alpha = -1.;
    daxpy_(&n, &alpha, wc, &incx, wnum1, &incy);

    dcopy_(&n, z, &incx, znum1, &incy);
    daxpy_(&n, &alpha, zc, &incx, znum1, &incy);

    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, znum1, &incx, &beta, wnum1, &incy);

    /*      //    wnum1(:) =(w(:)-wc(:))+matmul( k(:,:),(z(:)-zc(:)))*/

    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, kinv, &n, wnum1, &incx, &beta, kinvnum1, &incy);

    num11 = ddot_(&n, wnum1, &incx, kinvnum1, &incy);

    dcopy_(&n, z, &incx, zz, &incy);
    dcopy_(&n, w, &incx, ww, &incy);

    alpha = 1.;
    daxpy_(&n, &alpha, wc, &incx, ww, &incy);

    daxpy_(&n, &alpha, zc, &incx, zz, &incy);

    beta = 0.;
    alpha = 1.;
    dgemv_(&trans, &n, &n, &alpha, k, &n, zz, &incx, &beta, kzden1, &incy);

    den22 = ddot_(&n, zz, &incx, kzden1, &incy);

    beta = 0.;
    alpha = 1.;
    dgemv_(&trans, &n, &n, &alpha, kinv, &n, ww, &incx, &beta, kinvwden1, &incy);

    den11 = ddot_(&n, ww, &incx, kinvwden1, &incy);

    err0 = num11 / (den11 + den22);
    err = sqrt(err0);
    ++iter;

    /*for(i=0;i<n;i++)
      {

      printf("conv i %d wc %14.7e zc %14.7e z %14.7e w %14.7e \n",i,wc[i],zc[i],z[i],w[i]);
      }*/

    if (iout > 0)
    {
      for (i = 0 ; i < n ; ++i)
      {
        /*result_gs[i][iter-1] = z[i]; */
        fprintf(f101, "%d  %d  %14.7e \n", iter - 1, i, z[i]);
      }

    }
  }

  iparamLCP[2] = iter;
  dparamLCP[2] = err;

  if (iout > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of Latin after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 1;
    }
    else
    {
      printf(" Convergence of Latin after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else if (err < tol) *info = 0;


  free(wc);
  free(zz);
  free(ww);
  free(zc);
  free(znum1);
  free(wnum1);
  free(kinvden1);
  free(kinvden2);
  free(wt);
  free(maxwt);
  free(num1);
  free(kinvnum1);
  free(den1);
  free(den2);
  free(wden1);
  free(zden1);
  free(kinvwden1);
  free(kzden1);

  for (i = 0 ; i < n; ++i) free(kinv[i]);
  for (i = 0 ; i < n; ++i) free(DPO[i]);
  for (i = 0 ; i < n; ++i) free(k[i]);

  free(kinv);
  free(k);
  free(DPO);

  if (iout > 0) fclose(f101);

}

