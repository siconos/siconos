#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*!\file cfp_latin.c

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


/*!\fn cfp_latin(double vec[],double *qq,int *nn, double * k_latin,double *mumu,int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)

   cfp_latin  is a specific latin solver for primal contact problem with friction.

   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param qq On enter a pointer over doubles containing the components of the double vector.
   \param nn On enter a pointer over integers, the dimension of the second member.
   \param k_latin On enter a pointer over doubles, the latin coefficient (positive).
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

cfp_latin(double vec[], double *qq, int *nn, double * k_latin, double *mumu, int * itermax, double * tol, double z[], double w[], int *it_end, double * res, int *info)
{

  int i, j, kk, iter1, ino, ddl;
  int n = *nn, incx = 1, incy = 1, nc = n / 2, idim, jdim, nbno, taille, taillet, taillen, itt = *itermax;
  double errmax = *tol, alpha, beta, maxa, bb, cc, zw, aa, nt, wn, tc, zc0, mu = *mumu;
  double rr, rrr, r1, r2, r3, invR0, invRT0, err1, z0, num11, err0, invRTinvR0;
  double den11, den22, vv, knz0, ktz0, ktz[nc], wf[nc];
  double  *wc, *zc, *wt, *maxwt, *wnum1, *znum1;
  double *zt, *maxzt, kn[nc][nc], kt[nc][nc];
  char trans = 'T';
  double k[n][n], A[n][n], R[n][n], RT[n][n], invRT[n][n], invR[n][n], kf[n][n], kninv[nc][nc];
  double invRTinvR[n][n], kinvwden1[n], kzden1[n], kfinv[n][n], knz[nc], wtnc[nc];
  int ddln[nc];
  int ddlt[nc], vectnt[n];


  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      k[i][j] = 0.;
      kf[i][j] = 0.;
      kfinv[i][j] = 0.;
    }

  for (i = 0; i < n; i++)
  {
    k[i][i] = *k_latin / vec[i * n + i];
    vectnt[i] = i + 1;
  }

  for (i = 0; i < nc; i++)
  {
    ddln[i] = vectnt[2 * i];
    if (i != 0)
    {
      ddlt[i] = vectnt[2 * i - 1];
    }
    else
    {
      ddlt[i] = 0;
    }
  }


  for (i = 0; i < nc; i++)
    for (j = 0; j < nc; j++)
    {
      kn[i][j] = 0.;
      kt[i][j] = 0.;
      kninv[i][j] = 0.;
    }

  for (i = 0; i < nc; i++)
  {
    kn[i][i] = k[ddln[i]][ddln[i]];
    kt[i][i] = k[ddlt[i]][ddlt[i]];
  }

  taillen = sizeof(ddln) / sizeof(ddln[0]);
  taillet = sizeof(ddlt) / sizeof(ddlt[0]);

  idim = 1 +  taillen / taillet;
  jdim = idim - 1;


  if (idim != 2)
  {
    printf("case not yet treated\n");
    return (*info = 3);
  }



  taille = 0;
  for (i = 0; i < n; i++)
    taille = sizeof(qq[i]) + taille;

  taille = taille / sizeof(qq[0]);
  nbno = taille / idim;

  for (i = 0; i < nc; i++)
  {
    kf[ddln[i]][ddln[i]] = kn[i][i];
    kf[ddlt[i]][ddlt[i]] = kt[i][i];
  }

  for (i = 0; i < n; i++)
    kfinv[i][i] = 1 / kf[i][i];

  for (i = 0; i < nc; i++)
    kninv[i][i] = 1 / kt[i][i];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      A[i][j] = vec[i * n + j] + kfinv[i][j];
      R[i][j] = 0.;
    }


  /*/   !!!!!!!!!!!!!!!!!!!!!Cholesky!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  R[0][0] = sqrt(A[0][0]);

  for (i = 1; i < n; i++)
  {
    rr = 0.0;
    rrr = 0.0;
    for (j = 0; j <= i; j++)
    {
      r1 = 0.0;
      r2 = 0.0;
      for (kk = 0; kk <= j - 1; kk++)
      {
        rr = R[i][kk] * R[j][kk] + r1;
        r1 = rr;
      }

      if (fabs(R[j][j]) <= 1.e-10)
      {
        printf("nul pivot %d ,and R(%d,%d) %g \n", j, j, j, R[j][j]);
        break;
      }
      R[i][j] = (1 / R[j][j]) * (A[i][j] - rr);
      r3 = 0.0;
      for (kk = 0; kk <= i - 1; kk++)
      {
        rrr = R[i][kk] * R[i][kk] + r3;
        r3 = rrr;
      }
      R[i][i] = sqrt(A[i][i] - rrr);
    }
  }


  /*    //  !!!!!end of cholesky!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*   //  !determination of the R tranposeted*/
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      RT[i][j] = R[j][i];

  /*   //  !inverse of R and RT*/

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      invRT[i][j] = 0.;
      invR[i][j] = 0.;
    }


  /*   //  !!!!!!!!!inversion of the inf triangular matrix!!!!!!!!!!!!!*/
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      if (i == j)
      {
        invR[i][j] = 1 / R[i][j];
      }
      else
      {
        invR0 = 0.;
        for (kk = j; kk <= (i - 1); kk++)
        {
          invR[i][j] = (-1 / R[i][i]) * R[i][kk] * invR[kk][j] + invR0;
          invR0 = invR[i][j];
        }
      }
    }

  /*   //  !!!!!!!!!!!!!!!!!!!end of inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*  //  !!!!!!!!!!!!!!!!!!!!!!inversion of the sup triangular matrix!!!!!!!*/
  for (i = 0; i < n; i++)
    invRT[i][i] = 1 / RT[i][i];

  for (i = n - 2; i >= 0; i--)
    for (j = n - 1; j >= 0; j--)
    {
      invRT0 = 0.;
      for (kk = i + 1; kk <= j; kk++)
      {
        invRT[i][j] = (-1 / RT[i][i]) * RT[i][kk] * invRT[kk][j] + invRT0;
        invRT0 = invRT[i][j];
      }
    }

  /*  //  !!!!!!!!!!!!!!!!!!!end of inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


  /*  //  ! initialisation  s^ = 0*/

  wc = (double*) malloc(n * sizeof(double));
  zc = (double*) malloc(n * sizeof(double));
  znum1 = (double*) malloc(n * sizeof(double));
  wnum1 = (double*) malloc(n * sizeof(double));
  wt = (double*) malloc(n * sizeof(double));
  maxwt = (double*) malloc(nc * sizeof(double));
  zt = (double*) malloc(nc * sizeof(double));
  maxzt = (double*) malloc(n * sizeof(double));



  for (i = 0; i < n; i++)
  {
    wc[i] = 0.0;
    zc[i] = 0.;
    z[i] = 0.;
    w[i] = 0.;
    znum1[i] = 0.;
    wnum1[i] = 0.;
    wt[i] = 0.;
    maxzt[i] = 0.;

  }


  for (i = 0; i < nc; i++)
  {
    maxwt[i] = 0.;
    zt[i] = 0.;
    knz[i] = 0.;
    ktz[i] = 0.;
    wf[i] = 0.;
    wtnc[i] = 0.;
  }

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      invRTinvR[i][j] = 0.;
    }

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      invRTinvR0 = 0.;
      for (kk = 0; kk < n; kk++)
      {
        invRTinvR[i][j] = invRT[i][kk] * invR[kk][j] + invRTinvR0;
        invRTinvR0 = invRTinvR[i][j];
      }
    }



  /*    //  !iteration loops*/

  iter1 = 0;
  err1 = 1.;

  while ((iter1 < itt) && (err1 > errmax))
  {

    /*   !linear stage (zc,wc) -> (z,w)*/


    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, kfinv, &n, zc, &incx, &beta, wc, &incy);
    dcopy_(&n, qq, &incx, znum1, &incy);
    alpha = 1.;
    daxpy_(&n, &alpha, wc, &incx, znum1, &incy);
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, invRTinvR, &n, znum1, &incx, &beta, z, &incy);
    alpha = -1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, kfinv, &n, z, &incx, &beta, wc, &incy);
    dcopy_(&n, wc, &incx, w, &incy);

    /*      // Local stage (z,w)->(zc,wc)*/

    for (i = 0; i < n; i++)
    {
      zc[i] = 0.;
      wc[i] = 0.0;
    }

    /*      // normal party*/
    for (i = 0; i < nc; i++)
    {
      knz0 = 0.;
      for (kk = 0; kk < nc; kk++)
      {
        knz[i] = kt[i][kk] * w[ddlt[kk]] + knz0;
        knz0 = knz[i];
      }

      zt[i] = z[ddlt[i]] - knz[i];

      if (zt[i] > 0.0)
      {
        zc[ddlt[i]] = zt[i];
        maxzt[i] = 0.0;
      }
      else
      {
        zc[ddlt[i]] = 0.0;
        maxzt[i] = -zt[i];
      }
    }


    for (i = 0; i < nc; i++)
    {
      zc0 = 0.;
      ktz0 = 0.;
      for (j = 0; j < nc; j++)
      {
        wc[ddlt[i]] = kninv[i][j] * maxzt[j] + zc0;
        zc0 = wc[ddlt[i]];
        ktz[i] = kn[i][j] * w[ddln[j]] + ktz0;
        ktz0 =  ktz[i];
      }
      wf[i] = z[ddln[i]] - ktz[i];
    }


    /*   loop other nodes*/

    for (ino = 0; ino < nbno; ino++)
    {
      ddl = ddln[ino];
      nt = fabs(wf[ino]);
      /*  tangential vector*/
      if (nt < 1.e-8) tc = 0.;
      else tc = wf[ino] / nt;
      /*   composante selon le vecteur tangentiel*/
      wn = zc[ddlt[ino]];
      aa = nt - mu * wn;
      if (aa > 0.0)
      {
        maxa = aa;
      }
      else
      {
        maxa = 0.0;
      }
      wc[ddl] = (maxa / (-1 * kn[ino][ino])) * tc;
      aa = -nt + mu * wn;
      if (aa > 0.0)
      {
        maxa = aa;
      }
      else
      {
        maxa = 0.0;
      }
      zc[ddl] = (mu * wn - maxa) * tc;
    }

    /*  convergence criterium */

    dcopy_(&n, z, &incx, znum1, &incy);
    alpha = -1.;
    daxpy_(&n, &alpha, zc, &incx, znum1, &incy);
    dcopy_(&n, w, &incx, wnum1, &incy);
    daxpy_(&n, &alpha, wc, &incx, wnum1, &incy);
    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, kf, &n, wnum1, &incx, &beta, znum1, &incy);
    num11 = 0.;
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, kfinv, &n, znum1, &incx, &beta, wnum1, &incy);
    num11 = ddot_(&n, wnum1, &incx, znum1, &incy);
    dcopy_(&n, z, &incx, znum1, &incy);
    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, kf, &n, w, &incx, &beta, znum1, &incy);
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, kfinv, &n, znum1, &incx, &beta, wnum1, &incy);
    den11 = ddot_(&n, wnum1, &incx, znum1, &incy);

    dcopy_(&n, zc, &incx, znum1, &incy);
    alpha = 1.;
    beta = 1.;
    dgemv_(&trans, &n, &n, &alpha, kf, &n, wc, &incx, &beta, znum1, &incy);
    alpha = 1.;
    beta = 0.;
    dgemv_(&trans, &n, &n, &alpha, kfinv, &n, znum1, &incx, &beta, wnum1, &incy);
    den22 = ddot_(&n, znum1, &incx, wnum1, &incy);
    err0 = num11 / (den11 + den22);
    err1 = sqrt(err0);
    iter1 = iter1 + 1;
    *it_end = iter1;
    *res = err1;
  }


  if (err1 > errmax)
  {
    printf("no convergence after %d iterations, the residue is %g\n", iter1, err1);
    *info = 1;
  }
  else
  {
    printf("there is convergence after %d iterations, the residue is %g \n", iter1, err1);
    *info = 0;
  }



  free(wc);
  free(zc);
  free(znum1);
  free(wnum1);
  free(wt);
  free(maxwt);
  free(zt);
  free(maxzt);

}
