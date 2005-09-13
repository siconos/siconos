#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"

/*!\file lcp_latin.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
 * Try \f$(z,w)\f$ such that:\n
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *    M z + q= w\\
 *    0 \le z \perp w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 *
 * where M is an (n x n)-matrix, q , w and z n-vectors.
 *
 * \fn  lcp_latin( int *nn , double *vec , double *q , double *z , int *info ,
 *                int *iparamLCP , double *dparamLCP )
 *
 * lcp_latin (LArge Time INcrements) is a solver for LCP based on the principle of splitting method\n
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
 * Specific LATIN parameters:\n
 *
 * \param iparamLCP[0] = itermax  Input unchanged parameter which represents the maximum number of iterations allowed.
 * \param iparamLCP[1] = iout     Input unchanged parameter which represents the output log identifiant\n
 *                       0 - no output\n
 *                       0 < active screen output\n
 * \param iparamLCP[2] = it_end   Output modified parameter which returns the number of iterations performed by the algorithm.
 *
 * \param dparamLCP[0] = tol      Input unchanged parameter which represents the tolerance required.
 * \param dparamLCP[1] = k_latin  Input unchanged parameter which represents the latin parameter.
 * \param dparamLCP[2] = res      Output modified parameter which returns the final error value.
 *
 * \author Nineb Sheherazade.
 * Last modifications: Mathieu Renouf
 */

void lcp_latin(int *nn , double *vec , double *qq , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP)
{

  FILE *f101;

  int i, j, iter, info2, nrhs;
  int n, n2, iout, itermax, tol;

  int incx = 1, incy = 1;

  double alpha, beta, k_latin;
  double err, num11, err0;
  double den11, den22;

  double  *wc, *zc, *kinvden1, *kinvden2, *wt, *maxwt, *wnum1, *znum1, *ww, *zz;
  double *num1, *kinvnum1, *den1, *den2, *wden1, *zden1;

  char TRANS = 'T', UPLO = 'U', DIAG = 'N', NOTRANS = 'N';

  double  *kinvwden1, *kzden1;

  double  *k, *kinv, *DPO;

  /*input*/

  itermax = iparamLCP[0];
  iout    = iparamLCP[1];
  tol     = dparamLCP[0];
  k_latin = dparamLCP[1];

  /* *** */

  n  = *nn;
  n2 = n * n;

  /* ********* */

  if (iout > 0) f101 = fopen("resultat_latin.dat", "w+");

  /* Preparation of the diagonal of the inverse matrix */

  k    = (double*) malloc(n * sizeof(double));
  kinv = (double*) malloc(n * sizeof(double));

  for (i = 0 ; i < n ; ++i)
  {

    k[i] = k_latin * vec[i * n + i];

    if (fabs(k[i]) < 1e-16)
    {

      if (iout > 0)
      {
        printf(" Warning negative diagonal term \n");
        printf(" The local problem can be solved \n");
      }

      free(k);
      free(kinv);
      *info = 4;

      return;
    }
    else kinv[i] = 1.0 / k[i];
  }

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

  DPO  = (double*) malloc(n2 * sizeof(double));
  k    = (double*) malloc(n * sizeof(double));
  kinv = (double*) malloc(n * sizeof(double));

  /* Initialization */

  for (i = 0; i < n; i++)
  {

    z[i] = 0.;
    w[i] = 0.;
    wc[i] = 0.;
    wt[i] = 0.;
    zc[i] = 0.;
    num1[i] = 0.;
    den1[i] = 0.;
    den2[i] = 0.;
    znum1[i] = 0.;
    wnum1[i] = 0.;
    maxwt[i] = 0.;
    kinvden1[i] = 0.;
    kinvden2[i] = 0.;
    kinvnum1[i] = 0.;

  }

  for (i = 0; i < n2 ; ++i) DPO[i] = 0.;

  for (i = 0 ; i < n ; ++i)
  {
    for (j = 0 ; j < n ; ++j) DPO[n * j + i] = vec[j * n + i];
    DPO[n * i + i] += k[i];
  }

  /* **     Call Cholesky     ** */
  /* ** Only for PSD matrices ** */

  dpotrf_(&UPLO , &n , DPO , &n, &info2);

  if (info2 != 0)
  {
    printf("\n Cholesky Factorization failed \n");
    printf(" Minor %d non PSD\n", info2);
    return (*info = 2);
  }

  /* **** End of Cholesky ***** */

  /* start  iteration loops */

  iter = 0;
  err  = 1.;

  /*Check formulation ?*/

  alpha = -1.;
  dscal_(&n , &alpha , qq , &incx);

  while ((iter < itermax) && (err > tol))
  {

    /* linear stage (zc,wc) -> (z,w)*/

    for (i = 0 ; i < n ; ++i) wc[i] +=  k[i] * z[i];

    dcopy_(&n , qq , &incx , znum1 , &incy);

    alpha = 1.;
    daxpy_(&n , &alpha , wc , &incx , znum1, &incy);

    nrhs = 1;

    dtrtrs_(&UPLO , &NOTRANS , &DIAG , &n , &nrhs , DPO , &n , znum1 , &n , &info2);

    dtrtrs_(&UPLO ,   &TRANS , &DIAG , &n , &nrhs , DPO , &n , znum1 , &n , &info2);

    dcopy_(&n , znum1 , &incx , z , &incy);

    alpha = -1.;
    beta  =  1.;

    for (i = 0 ; i < n ; ++i) wc[i] +=  k[i] * z[i];

    dcopy_(&n , wc , &incx , w , &incy);

    /* Local Stage */

    dcopy_(&n , w , &incx , wt , &incy);

    alpha = -1.;
    beta  =  1.;

    for (i = 0 ; i < n ; ++i)
    {

      wt[i] +=  k[i] * zc[i];

      if (wt[i] > 0.0)
      {
        wc[i] = wt[i];
        zc[i] = 0.;
      }
      else
      {
        wc[i] = 0.0;
        zc[i] = -wt[i] * kinv[i];
      }
    }

    /* convergence criterium */

    dcopy_(&n , w , &incx , wnum1 , &incy);
    alpha = -1.;
    daxpy_(&n , &alpha , wc , &incx , wnum1 , &incy);

    dcopy_(&n , z      , &incx , znum1 , &incy);
    daxpy_(&n , &alpha , zc    , &incx , znum1 , &incy);

    alpha = 1.;
    beta  = 1.;

    for (i = 0 ; i < n ; ++i) wnum1[i] +=  k[i] * znum1[i];

    /*      //    wnum1(:) =(w(:)-wc(:))+matmul( k(:,:),(z(:)-zc(:)))*/

    alpha = 1.;
    beta  = 0.;

    for (i = 0 ; i < n ; ++i) kinvnum1[i] =  kinv[i] * wnum1[i];

    num11 = ddot_(&n, wnum1, &incx, kinvnum1, &incy);

    dcopy_(&n, z, &incx, zz, &incy);
    dcopy_(&n, w, &incx, ww, &incy);

    alpha = 1.;

    daxpy_(&n, &alpha, wc, &incx, ww, &incy);

    daxpy_(&n, &alpha, zc, &incx, zz, &incy);

    for (i = 0 ; i < n ; ++i) kzden1[i] =  k[i] * zz[i];

    den22 = ddot_(&n, zz, &incx, kzden1, &incy);

    for (i = 0 ; i < n ; ++i) kinvwden1[i] =  kinv[i] * ww[i];

    den11 = ddot_(&n, ww, &incx, kinvwden1, &incy);

    err0 = num11 / (den11 + den22);
    err = sqrt(err0);

    ++iter;

    /*for(i=0;i<n;i++) printf("conv i %d wc %14.7e zc %14.7e z %14.7e w %14.7e \n",i,wc[i],zc[i],z[i],w[i]); */

    if (iout > 0)
    {
      /*Print result result_gs[i][iter-1] = z[i]; */
      for (i = 0 ; i < n ; ++i)
      {
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

  free(k);
  free(wc);
  free(zz);
  free(ww);
  free(zc);
  free(wt);
  free(DPO);
  free(num1);
  free(den1);
  free(den2);
  free(kinv);
  free(maxwt);
  free(wden1);
  free(zden1);
  free(znum1);
  free(wnum1);
  free(kzden1);
  free(kinvnum1);
  free(kinvden1);
  free(kinvden2);
  free(kinvwden1);

  if (iout > 0) fclose(f101);

}

