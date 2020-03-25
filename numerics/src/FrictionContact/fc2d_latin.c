/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <assert.h>                  // for assert
#include <float.h>                   // for DBL_EPSILON
#include <math.h>                    // for fabs, sqrt
#include <stdio.h>                   // for printf, size_t, NULL
#include <stdlib.h>                  // for free, calloc
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "Friction_cst.h"            // for SICONOS_FRICTION_2D_K_LATIN
#include "NumericsFwd.h"             // for SolverOptions, FrictionContactPr...
#include "NumericsMatrix.h"          // for NumericsMatrix, RawNumericsMatrix
#include "SolverOptions.h"           // for SolverOptions, solver_options_nu...
#include "fc2d_Solvers.h"            // for fc2d_latin, fc2d_latin_setDefaul...
#include "numerics_verbose.h"        // for verbose
#include "SiconosBlas.h"             // for cblas_dgemv, cblas_dcopy, cblas_...
#include "SiconosLapack.h"           // for DTRTRS, DPOTRF, lapack_int, LA_UP, LA_NONUNIT, LA_NOTRANS


void fc2d_latin(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
{
  int nc = problem->numberOfContacts;
  assert(nc>0);
  double * vec = problem->M->matrix0;
  double *qq = problem->q;
  double * mu = problem->mu;



  lapack_int info77 = 0;
  int i, j, kk, ddl, nrhs;
  lapack_int info2 = 0;
  int n = 2 * nc;
  size_t idim, nbno, ino;
  int incx = 1, incy = 1;
  int itt, iter1;
  size_t taille, taillet, taillen;
  int *ddln;
  int *ddlt, *vectnt;
  assert(n>0);

  double  errmax, alpha, beta, maxa, k_latin;
  double  aa, nt, wn, tc, zc0;
  double  err1, num11, err0;
  double  den11, den22, knz0, ktz0, *ktz, *wf;
  double  *wc, *zc, *wt, *maxwt, *wnum1, *znum1;
  double  *zt, *maxzt;

  double  *kn, *kt;

  // char    trans='T', diag='N';
  // char    uplo='U', notrans='N';



  double  *k, *DPO, *kf, *kninv;
  double  *kinvwden1, *kzden1, *kfinv, *knz, *wtnc;



  /*                Recup input                    */


  itt     = options->iparam[SICONOS_IPARAM_MAX_ITER];
  errmax  = options->dparam[SICONOS_DPARAM_TOL];
  k_latin = options->dparam[SICONOS_FRICTION_2D_K_LATIN];

  /*               Initialize output                */


  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;


  /*               Allocations                      */

  k         = (double*) calloc(n * n, sizeof(double));
  DPO       = (double*) calloc(n * n, sizeof(double));
  kf        = (double*) calloc(n * n, sizeof(double));
  kfinv     = (double*) calloc(n * n, sizeof(double));

  kninv     = (double*) calloc(nc * nc, sizeof(double));
  kn        = (double*) calloc(nc * nc, sizeof(double));
  kt        = (double*) calloc(nc * nc, sizeof(double));

  kinvwden1 = (double*) calloc(n, sizeof(double));
  kzden1    = (double*) calloc(n, sizeof(double));
  wc        = (double*) calloc(n, sizeof(double));
  zc        = (double*) calloc(n, sizeof(double));
  znum1     = (double*) calloc(n, sizeof(double));
  wnum1     = (double*) calloc(n, sizeof(double));
  wt        = (double*) calloc(n, sizeof(double));
  maxzt     = (double*) calloc(n, sizeof(double));



  knz       = (double*) calloc(nc, sizeof(double));
  wtnc      = (double*) calloc(nc, sizeof(double));
  ktz       = (double*) calloc(nc, sizeof(double));
  wf        = (double*) calloc(nc, sizeof(double));
  maxwt     = (double*) calloc(nc, sizeof(double));
  zt        = (double*) calloc(nc, sizeof(double));


  vectnt    = (int*) calloc(n, sizeof(int));

  ddln      = (int*) calloc(nc, sizeof(int));
  ddlt      = (int*) calloc(nc, sizeof(int));

  /*                    Initialization                   */



  for(i = 0; i < n * n; i++)
  {
    if(i < nc * nc)
    {
      if(i < n)
      {
        reaction[i]     = 0.;
        velocity[i]     = 0.;
      }
    }
  }


  for(i = 0; i < n; i++)
  {

    if(fabs(vec[i * n + i]) < DBL_EPSILON)
    {

      if(verbose > 0)
        printf("\n Warning nul diagonal term in M matrix \n");

      free(k);
      free(DPO);
      free(kf);
      free(kfinv);
      free(kninv);
      free(kn);
      free(kt);
      free(kinvwden1);
      free(kzden1);
      free(wc);
      free(zc);
      free(znum1);
      free(wnum1);
      free(wt);
      free(maxzt);
      free(knz);
      free(wtnc);
      free(ktz);
      free(wf);
      free(maxwt);
      free(zt);
      free(vectnt);
      free(ddln);
      free(ddlt);

      *info = 3;

      return;


    }
    else
    {

      k[i + n * i] = k_latin / vec[i * n + i];
      vectnt[i] = i + 1;

    }

  }


  for(i = 0; i < nc; i++)
  {
    ddln[i] = vectnt[2 * i];
    if(i != 0) ddlt[i] = vectnt[2 * i - 1];
    else ddlt[i] = 0;

  }


  for(i = 0; i < nc; i++)
  {
    kn[i + nc * i] = k[ddln[i] + n * ddln[i]];
    kt[i + nc * i] = k[ddlt[i] + n * ddlt[i]];
  }




  taillen = nc; //sizeof(ddln) / sizeof(ddln[0]);
  taillet = nc; //sizeof(ddlt) / sizeof(ddlt[0]);

  idim = 1 +  taillen / taillet;

  taille = 0;
  for(i = 0; i < n; i++)
    taille = sizeof(qq[i]) + taille;

  taille = taille / sizeof(qq[0]);
  nbno = taille / idim;


  for(i = 0; i < nc; i++)
  {
    kf[ddln[i] + n * ddln[i]] = kn[i + nc * i];
    kf[ddlt[i] + n * ddlt[i]] = kt[i + nc * i];
  }


  for(i = 0; i < n; i++)
  {
    kfinv[i + n * i] = 1. / kf[i + n * i];

    if(i < nc)
      kninv[i + nc * i] = 1. / kt[i + nc * i];

  }


  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      DPO[i + n * j] = vec[j * n + i] + kfinv[i + n * j];



  DPOTRF(LA_UP, n, DPO, n, &info2);

  if(info2 != 0)
  {
    if(verbose > 0)
      printf("\n Matter with Cholesky factorization \n");

    free(k);
    free(DPO);
    free(kf);
    free(kfinv);
    free(kninv);
    free(kn);
    free(kt);
    free(kinvwden1);
    free(kzden1);
    free(wc);
    free(zc);
    free(znum1);
    free(wnum1);
    free(wt);
    free(maxzt);
    free(knz);
    free(wtnc);
    free(ktz);
    free(wf);
    free(maxwt);
    free(zt);
    free(vectnt);
    free(ddln);
    free(ddlt);

    *info = 2;
    return;
  }

  /*                Iteration loops                  */


  iter1 = 0;
  err1  = 1.;

  while((iter1 < itt) && (err1 > errmax))
  {

    /*               Linear stage (zc,wc) -> (z,w)         */

    alpha  = 1.;
    beta   = 1.;
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alpha, kfinv, n, zc, incx, beta, wc, incy);

    cblas_dcopy(n, qq, incx, znum1, incy);

    alpha = -1.;
    cblas_dscal(n, alpha, znum1, incx);

    alpha = 1.;
    cblas_daxpy(n, alpha, wc, incx, znum1, incy);

    nrhs = 1;
    DTRTRS(LA_UP, LA_TRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, &info77);

    DTRTRS(LA_UP, LA_NOTRANS, LA_NONUNIT, n, nrhs, DPO, n, znum1, n, &info77);

    cblas_dcopy(n, znum1, incx, reaction, incy);

    alpha = -1.;
    beta = 1.;
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alpha, kfinv, n, reaction, incx, beta, wc, incy);

    cblas_dcopy(n, wc, incx, velocity, incy);



    /*               Local stage (z,w)->(zc,wc)          */


    for(i = 0; i < n; i++)
    {
      zc[i] = 0.;
      wc[i] = 0.0;
    }


    /*          Normal party                           */



    for(i = 0; i < nc; i++)
    {
      knz0 = 0.;
      for(kk = 0; kk < nc; kk++)
      {
        knz[i] = kt[i + nc * kk] * velocity[ddlt[kk]] + knz0;
        knz0 = knz[i];
      }

      zt[i] = reaction[ddlt[i]] - knz[i];

      if(zt[i] > 0.0)
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

    for(i = 0; i < nc; i++)
    {
      zc0 = 0.;
      ktz0 = 0.;
      for(j = 0; j < nc; j++)
      {
        wc[ddlt[i]] = kninv[i + nc * j] * maxzt[j] + zc0;
        zc0 = wc[ddlt[i]];
        ktz[i] = kn[i + nc * j] * velocity[ddln[j]] + ktz0;
        ktz0 =  ktz[i];
      }
      wf[i] = reaction[ddln[i]] - ktz[i];
    }


    /*             Loop other nodes              */


    for(ino = 0; ino < nbno; ino++)
    {
      ddl  = ddln[ino];
      nt   = fabs(wf[ino]);


      /*          Tangential vector              */



      if(nt < 1.e-8) tc = 0.;
      else tc = wf[ino] / nt;



      /*               Tangentiel component             */


      wn = zc[ddlt[ino]];

      aa = nt - mu[ino] * wn;

      if(aa > 0.0)
      {
        maxa = aa;
      }
      else
      {
        maxa = 0.0;
      }

      wc[ddl] = (maxa / (-1 * kn[ino + nc * ino])) * tc;

      aa = -nt + mu[ino] * wn;

      if(aa > 0.0)
      {
        maxa = aa;
      }
      else
      {
        maxa = 0.0;
      }

      zc[ddl] = (mu[ino] * wn - maxa) * tc;

    }

    /*               Convergence criterium                */



    cblas_dcopy(n, reaction, incx, znum1, incy);

    alpha = -1.;
    cblas_daxpy(n, alpha, zc, incx, znum1, incy);

    cblas_dcopy(n, velocity, incx, wnum1, incy);

    cblas_daxpy(n, alpha, wc, incx, wnum1, incy);

    alpha  = 1.;
    beta   = 1.;
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alpha, kf, n, wnum1, incx, beta, znum1, incy);

    num11  = 0.;
    alpha  = 1.;
    beta = 0.;
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alpha, kfinv, n, znum1, incx, beta, wnum1, incy);

    num11 = cblas_ddot(n, wnum1, incx, znum1, incy);

    cblas_dcopy(n, reaction, incx, znum1, incy);

    alpha  = 1.;
    beta   = 1.;
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alpha, kf, n, velocity, incx, beta, znum1, incy);

    alpha  = 1.;
    beta   = 0.;
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alpha, kfinv, n, znum1, incx, beta, wnum1, incy);

    den11  = cblas_ddot(n, wnum1, incx, znum1, incy);

    cblas_dcopy(n, zc, incx, znum1, incy);

    alpha  = 1.;
    beta   = 1.;
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alpha, kf, n, wc, incx, beta, znum1, incy);

    alpha  = 1.;
    beta   = 0.;
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alpha, kfinv, n, znum1, incx, beta, wnum1, incy);

    den22  = cblas_ddot(n, znum1, incx, wnum1, incy);

    err0   = num11 / (den11 + den22);

    err1   = sqrt(err0);

    iter1   = iter1 + 1;
  }
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter1;
  options->dparam[SICONOS_DPARAM_RESIDU] = err1;

  if(err1 > errmax)
  {

    if(verbose > 0)
      printf("No convergence after %d iterations, the residue is %g\n", iter1, err1);

    *info = 1;
  }
  else
  {

    if(verbose > 0)
      printf("Convergence after %d iterations, the residue is %g \n", iter1, err1);

    *info = 0;
  }

  free(k);
  free(DPO);
  free(kf);
  free(kfinv);
  free(kninv);
  free(kn);
  free(kt);
  free(kinvwden1);
  free(kzden1);
  free(wc);
  free(zc);
  free(znum1);
  free(wnum1);
  free(wt);
  free(maxzt);
  free(knz);
  free(wtnc);
  free(ktz);
  free(wf);
  free(maxwt);
  free(zt);
  free(vectnt);
  free(ddln);
  free(ddlt);



}

void fc2d_latin_set_default(SolverOptions *options)
{
  options->dparam[SICONOS_FRICTION_2D_K_LATIN] = 0.3;
}
