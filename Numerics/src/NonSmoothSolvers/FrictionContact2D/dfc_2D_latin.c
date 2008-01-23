/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"

void dfc_2D_latin(double *vec, double *qq, int *nn, double * k_latin, double *mu, int * itermax, double * tol, int *chat, double *z, double *w, int *it_end, double * res, int *info)
{


  int       i, j, kk, iter1, ino, ddl, info2, nrhs;
  int       n = *nn,  nc = n / 2, idim, jdim, nbno, itt = *itermax;
  int       taille, taillet, taillen, ispeak = *chat;
  int       *ddln, *ddlt, *vectnt;
  int incx = 1, incy = 1;


  double    errmax = *tol, alpha, beta, maxa, aa, nt, wn, tc, zc0;
  double    err1, num11, err0;
  double    den11, den22, knz0, ktz0, *ktz, *wf;
  double    *wc, *zc, *maxwt, *wnum1, *znum1;
  double    *kn, *kt;
  double    *k, *R;
  double    *kf, *kninv;
  double    *kinvwden1, *kzden1, *kfinv, *knz, *wtnc;

  /*  char      trans='T', uplo='U', diag='N', notrans='N';*/




  ktz       = (double *) malloc(nc * sizeof(double));
  wf        = (double *) malloc(nc * sizeof(double));
  kn        = (double *) malloc(nc * nc * sizeof(double));
  kt        = (double *) malloc(nc * nc * sizeof(double));
  k         = (double *) malloc(n * n * sizeof(double));

  R         = (double *) malloc(n * n * sizeof(double));


  kf        = (double *) malloc(n * n * sizeof(double));
  kninv     = (double *) malloc(nc * nc * sizeof(double));


  kinvwden1 = (double *) malloc(n * sizeof(double));
  kzden1    = (double *) malloc(n * sizeof(double));
  kfinv     = (double *) malloc(n * n * sizeof(double));
  knz       = (double *) malloc(nc * sizeof(double));
  wtnc      = (double *) malloc(nc * sizeof(double));

  ddln      = (int *) malloc(nc * sizeof(int));
  ddlt      = (int *) malloc(nc * sizeof(int));
  vectnt    = (int *) malloc(n * sizeof(int));

  wc        = (double *) malloc(n * sizeof(double));
  zc        = (double *) malloc(n * sizeof(double));
  znum1     = (double *) malloc(n * sizeof(double));
  wnum1     = (double *) malloc(n * sizeof(double));
  maxwt     = (double *) malloc(nc * sizeof(double));





  for (i = 0; i < n; i++)
  {

    wc[i]    = 0.;
    zc[i]    = 0.;
    z[i]     = 0.;
    w[i]     = 0.;
    znum1[i] = 0.;
    wnum1[i] = 0.;

    if (i < nc)
    {
      maxwt[i] = 0.;
      knz[i]   = 0.;
      ktz[i]   = 0.;
      wf[i]    = 0.;
      wtnc[i]  = 0.;
    }


  }



  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      k[i + n * j]     = 0.;
      kf[i + n * j   ] = 0.;
      kfinv[i + n * j] = 0.;
    }

  for (i = 0; i < n; i++)
  {
    k[i + n * i]  = *k_latin * vec[i * n + i];
    vectnt[i] = i + 1;
  }

  for (i = 0; i < nc; i++)
  {
    ddln[i] = vectnt[2 * i];

    if (i != 0) ddlt[i] = vectnt[2 * i - 1];
    else ddlt[i] = 0;
  }

  for (i = 0; i < nc; i++)
    for (j = 0; j < nc; j++)
    {
      kn[i + nc * j]    = 0.;
      kt[i + nc * j]    = 0.;
      kninv[i + nc * j] = 0.;
    }

  for (i = 0; i < nc; i++)
  {
    kn[i + nc * i] = k[ddln[i] + n * ddln[i]];
    kt[i + nc * i] = k[ddlt[i] + n * ddlt[i]];
  }

  taillen = sizeof(ddln) / sizeof(ddln[0]);
  taillet = sizeof(ddlt) / sizeof(ddlt[0]);

  idim = 1 +  taillen / taillet;
  jdim = idim - 1;


  taille = 0;

  for (i = 0; i < n; i++)
    taille = sizeof(qq[i]) + taille;

  taille = taille / sizeof(qq[0]);
  nbno   = taille / idim;

  for (i = 0; i < nc; i++)
  {
    kf[ddln[i] + n * ddln[i]] = kn[i + nc * i];
    kf[ddlt[i] + n * ddlt[i]] = kt[i + nc * i];
  }





  for (i = 0; i < n; i++)
  {

    if (fabs(kf[i + n * i]) < 1e-12)
    {

      if (ispeak > 0)
        printf("\n Warning nul diagonal term \n");


      free(wc);
      free(zc);
      free(ktz);
      free(wf);
      free(kn);

      free(kt);
      free(k);
      free(R);
      free(kf);
      free(kninv);

      free(kinvwden1);
      free(kzden1);
      free(kfinv);
      free(knz);
      free(wtnc);

      free(ddln);
      free(ddlt);
      free(vectnt);
      free(znum1);
      free(wnum1);

      free(maxwt);

      *info = 3;
      return;


    }
    else

      kfinv[i + n * i] = 1 / kf[i + n * i];
  }



  for (i = 0; i < nc; i++)
  {

    if (fabs(kn[i + nc * i]) < 1e-12)
    {

      if (ispeak > 0)
        printf("\n Warning nul diagonal term\n");


      free(wc);
      free(zc);
      free(ktz);
      free(wf);
      free(kn);

      free(kt);
      free(k);
      free(R);
      free(kf);
      free(kninv);

      free(kinvwden1);
      free(kzden1);
      free(kfinv);
      free(knz);
      free(wtnc);

      free(ddln);
      free(ddlt);
      free(vectnt);
      free(znum1);
      free(wnum1);

      free(maxwt);

      *info = 4;
      return;




    }
    else
      kninv[i + nc * i] = 1 / kn[i + nc * i];
  }


  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      R[i + n * j] = vec[j * n + i] + kf[i + n * j];

    }

  /*            Cholesky              */



  DPOTRF(LA_UP, n, R , n, info2);

  if (info2 != 0)
  {

    if (ispeak > 0)
      printf(" Matter with Cholesky Factorization \n ");

    free(wc);
    free(zc);
    free(ktz);
    free(wf);
    free(kn);

    free(kt);
    free(k);
    free(R);
    free(kf);
    free(kninv);

    free(kinvwden1);
    free(kzden1);
    free(kfinv);
    free(knz);
    free(wtnc);

    free(ddln);
    free(ddlt);
    free(vectnt);
    free(znum1);
    free(wnum1);

    free(maxwt);

    *info = 2;
    return;
  }








  /*    Iteration loops     */



  iter1 = 0;
  err1  = 1.;

  while ((iter1 < itt) && (err1 > errmax))
  {

    /*   Linear stage (zc,wc) -> (z,w)          */


    alpha = 1.;
    beta  = 1.;
    DGEMV(LA_TRANS, n, n, alpha, kf, n, zc, incx, beta,  wc, incy);

    DCOPY(n, qq, incx, znum1, incy);


    alpha = -1.;
    DSCAL(n, alpha, znum1, incx);


    alpha = 1.;
    DAXPY(n, alpha, wc, incx, znum1, incy);

    nrhs = 1;
    DTRTRS(LA_UP, LA_TRANS, LA_NONUNIT, n, nrhs, R, n, znum1, n, info2);

    DTRTRS(LA_UP, LA_NOTRANS, LA_NONUNIT, n, nrhs, R, n, znum1, n, info2);

    DCOPY(n, znum1, incx, z, incy);

    alpha = -1.;
    beta  = 1.;
    DGEMV(LA_TRANS, n, n, alpha, kf, n, z, incx, beta, wc, incy);
    DCOPY(n, wc, incx, w, incy);


    /*   Local stage (z,w)->(zc,wc)   */


    for (i = 0; i < n; i++)
    {
      zc[i] = 0.;
      wc[i] = 0.0;
    }



    /*   Normal party  */


    for (i = 0; i < nc; i++)
    {
      knz0 = 0.;
      for (kk = 0; kk < nc; kk++)
      {
        knz[i] = kn[i + nc * kk] * z[ddln[kk]] + knz0;
        knz0   = knz[i];
      }

      wtnc[i] = w[ddln[i]] - knz[i];
    }


    for (i = 0; i < nc; i++)
    {
      if (wtnc[i] < 0.0)
      {
        wc[ddln[i]] = 0.0;
      }
      else
      {
        wc[ddln[i]] = wtnc[i];
      }
    }


    for (i = 0; i < nc; i++)
    {
      if (-wtnc[i] < 0.0)
      {
        maxwt[i] = 0.0;
      }
      else
      {
        maxwt[i] = -wtnc[i];
      }
    }


    for (i = 0; i < nc; i++)
    {
      zc0 = 0.;
      for (j = 0; j < nc; j++)
      {
        zc[ddln[i]] = kninv[i + nc * j] * maxwt[j] + zc0;
        zc0         = zc[ddln[i]];
      }
    }



    /*      Tangential party        */


    for (i = 0; i < nc; i++)
    {
      ktz0 = 0.;
      for (kk = 0; kk < nc; kk++)
      {
        ktz[i] = kt[i + nc * kk] * z[ddlt[kk]] + ktz0;
        ktz0   =  ktz[i];
      }
      wf[i] = w[ddlt[i]] - ktz[i];
    }


    /*     Loop other nodes        */


    for (ino = 0; ino < nbno; ino++)
    {
      ddl = ddlt[ino];
      nt  = fabs(wf[ino]);


      /*  Tangential vector*/


      if (nt < 1.e-8) tc = 0.;
      else tc = wf[ino] / nt;


      /*  Tangentiel components     */


      wn = wc[ddln[ino]];

      aa = nt - mu[ino] * wn;

      if (aa < 0.0)
      {
        maxa = 0.0;
      }
      else
      {
        maxa = aa;
      }

      zc[ddl] = (maxa / (-1 * kt[ino + nc * ino])) * tc;

      aa = -nt + mu[ino] * wn;

      if (aa < 0.0)
      {
        maxa = 0.0;
      }
      else
      {
        maxa = aa;
      }

      wc[ddl] = (mu[ino] * wn - maxa) * tc;
    }




    /*    Convergence criterium          */



    DCOPY(n, z, incx, znum1, incy);
    alpha = -1.;
    DAXPY(n, alpha, zc, incx, znum1, incy);

    DCOPY(n, w, incx, wnum1, incy);

    DAXPY(n, alpha, wc, incx, wnum1, incy);

    alpha = 1.;
    beta = 1.;
    DGEMV(LA_TRANS, n, n, alpha, kf, n, znum1, incx, beta, wnum1, incy);

    num11 = 0.;

    alpha = 1.;
    beta = 0.;
    DGEMV(LA_TRANS, n, n, alpha, kfinv, n, wnum1, incx, beta, znum1, incy);

    num11 = DDOT(n, wnum1, incx, znum1, incy);

    DCOPY(n, z, incx, znum1, incy);

    alpha = 1.;
    DAXPY(n, alpha, zc, incx, znum1, incy);

    beta  = 0.;
    alpha = 1.;
    DGEMV(LA_TRANS, n, n, alpha, kf, n, znum1, incx, beta, wnum1, incy);

    den11 = DDOT(n, znum1, incx, wnum1, incy);


    DCOPY(n, w, incx, wnum1, incy);

    alpha = 1.;
    DAXPY(n, alpha, wc, incx, wnum1, incy);

    beta  = 0.;
    alpha = 1.;
    DGEMV(LA_TRANS, n, n, alpha, kfinv, n, wnum1, incx, beta, znum1, incy);

    den22 = DDOT(n, wnum1, incx, znum1, incy);


    err0    = num11 / (den11 + den22);
    err1    = sqrt(err0);
    *res    = err1;
    iter1   = iter1 + 1;
    *it_end = iter1;





  }


  if (err1 > errmax)
  {
    if (ispeak > 0)
      printf("No convergence after %d iterations, the residue is %g\n", iter1, err1);

    *info = 1;
  }
  else if (err1 < errmax)
  {

    if (ispeak > 0)
      printf("Convergence after %d iterations, the residue is %g \n", iter1, err1);

    *info = 0;
  }

  free(wc);
  free(zc);
  free(ktz);
  free(wf);
  free(kn);

  free(kt);
  free(k);
  free(R);
  free(kf);
  free(kninv);


  free(kinvwden1);
  free(kzden1);
  free(kfinv);
  free(knz);
  free(wtnc);


  free(ddln);
  free(ddlt);
  free(vectnt);
  free(znum1);
  free(wnum1);

  free(maxwt);



}

