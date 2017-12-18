/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "SiconosBlas.h"
#include "LinearComplementarityProblem.h"
#include "LCP_Solvers.h"
#include "lcp_cst.h"
#include "SolverOptions.h"
#include "NumericsMatrix.h"

#include "numerics_verbose.h"

void lcp_cpg(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  /* matrix M/vector q of the lcp */
  double * M = problem->M->matrix0;

  double * q = problem->q;

  /* size of the LCP */
  int n = problem->size;

  int incx, incy;
  int i, iter;
  int itermax = options->iparam[0];


  double err, a1, b1 , qs;

  double alpha, beta, rp, pMp;
  double tol = options->dparam[0];

  int *status;
  double *zz , *pp , *rr, *ww, *Mp;

  *info = 1;
  incx  = 1;

  /*output*/

  options->iparam[1] = 0;
  options->dparam[1] = 0.0;

  qs = cblas_dnrm2(n , q , incx);

  /*printf( " Norm: %g \n", qs );*/

  /* Allocations */

  status = (int*)malloc(n * sizeof(int));

  ww = (double*)malloc(n * sizeof(double));
  rr = (double*)malloc(n * sizeof(double));
  pp = (double*)malloc(n * sizeof(double));
  zz = (double*)malloc(n * sizeof(double));

  Mp = (double*)malloc(n * sizeof(double));

  incx = 1;

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

  cblas_dcopy(n, q, incx, rr, incy);

  a1 = -1.;
  b1 = -1.;

  cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, a1, M, n, z, incx, b1, rr, incy);

  /* Initialization of gradients */
  /* rr -> p and rr -> w */

  cblas_dcopy(n, rr, incx, ww, incy);
  cblas_dcopy(n, rr, incx, pp, incy);

  iter = 0;
  err  = 1.0 ;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;

    /* Compute initial pMp */

    incx = 1;
    incy = 1;

    cblas_dcopy(n, pp, incx, Mp, incy);

    a1 = 1.0;
    b1 = 0.0;

    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, a1, M, n, Mp, incx, b1, w, incy);

    pMp = cblas_ddot(n, pp, incx, w, incy);

    if (fabs(pMp) < DBL_EPSILON)
    {

      if (verbose > 0)
      {
        printf(" Operation not conform at the iteration %d \n", iter);
        printf(" Alpha can be obtained with pWp = %10.4g  \n", pMp);
        printf(" The residue is : %g \n", err);
      }

      free(Mp);
      free(ww);
      free(rr);
      free(pp);
      free(zz);
      free(status);

      options->iparam[1] = iter;
      options->dparam[1] = err;
      *info = 3;
      return;
    }

    rp  = cblas_ddot(n, pp, incx, rr, incy);

    alpha = rp / pMp;

    /*
     * Iterate prediction:
     * z' = z + alpha*p
     *
     */

    cblas_daxpy(n, alpha, pp, incx, z, incy);

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

    cblas_dcopy(n, rr, incx, w , incy);
    cblas_dcopy(n, q , incx, rr, incy);

    a1 = -1.;
    b1 = -1.;

    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, a1, M, n, z, incx, b1, rr, incy);

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

    /*  beta = -w.Mp / pMp  */

    rp = cblas_ddot(n , ww, incx, w, incy);

    beta = -rp / pMp;

    cblas_dcopy(n, ww, incx, pp, incy);
    cblas_daxpy(n, beta, zz, incx, pp, incy);

    /* **** Criterium convergence **** */

    cblas_dcopy(n, rr, incx, w, incy);
    qs   = -1.0;
    cblas_dscal(n, qs, w, incx);

    lcp_compute_error(problem, z, w, tol,  &err);

  }

  options->iparam[1] = iter;
  options->dparam[1] = err;

  cblas_dcopy(n, rr, incx, w, incy);

  qs   = -1.0;
  cblas_dscal(n, qs, w, incx);


  *info = 1;
  if (verbose > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
    }
    else
    {
      printf(" Convergence of CPG after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else if (err <= tol) *info = 0;

  free(Mp);
  free(status);

  free(ww);
  free(rr);
  free(pp);
  free(zz);

}
int linearComplementarity_cpg_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the CPG Solver\n");
  }
  options->solverId = SICONOS_LCP_CPG;


  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-6;


  return 0;
}
