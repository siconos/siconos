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
#include <math.h>                    // for sqrt
#include <stdio.h>                   // for printf
#include <stdlib.h>                  // for free, malloc, calloc, NULL
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "Friction_cst.h"            // for SICONOS_FRICTION_2D_CPG
#include "NumericsFwd.h"             // for SolverOptions, FrictionContactPr...
#include "NumericsMatrix.h"          // for NumericsMatrix, RawNumericsMatrix
#include "SolverOptions.h"           // for SolverOptions, solver_options_nu...
#include "fc2d_Solvers.h"            // for fc2d_projf, fc2d_projc, fc2d_cpg
#include "numerics_verbose.h"        // for verbose
#include "SiconosBlas.h"                   // for cblas_dcopy, cblas_ddot, cblas_d...

void fc2d_cpg(FrictionContactProblem* problem , double *reaction , double *velocity , int *info, SolverOptions* options)
{
  int nc = problem->numberOfContacts;
  assert(nc>0);
  double * vec = problem->M->matrix0;
  double * mu = problem->mu;

  int       n = 2 * nc, i, iter;
  assert(n>0);
  int       incx = 1, incy = 1;
  int       *stat, *statusi, it_end;


  double    eps = 1.e-12;
  double    pAp, alpha, beta, wAp, rp, normr;
  double    alphaf, betaf, den, num, res;

  double    *p, *fric, *r;
  double    *fric1, *v, *w, *Ap, *xi, *z;


  int maxit = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tol = options->dparam[SICONOS_DPARAM_TOL];
  options->iparam[SICONOS_IPARAM_ITER_DONE]  = 0;
  options->dparam[SICONOS_DPARAM_RESIDU]  = 0.0;


  r       = (double*) malloc(n * sizeof(double));
  p       = (double*) malloc(n * sizeof(double));
  v       = (double*) malloc(n * sizeof(double));
  w       = (double*) malloc(n * sizeof(double));
  Ap      = (double*) malloc(n * sizeof(double));
  xi      = (double*) malloc(n * sizeof(double));
  z       = (double*) malloc(n * sizeof(double));
  fric1   = (double*) malloc(n * sizeof(double));

  fric    = (double*) malloc(nc * sizeof(double));
  stat    = (int*)    malloc(nc * sizeof(int));
  statusi = (int*)    malloc(nc * sizeof(int));





  for (i = 0; i < n ; i++)
  {
    reaction[i]     = 0.0;
    xi[i]    = 0.0;
    r[i]     = 0.0;
    v[i]     = 0.0;
    p[i]     = 0.0;
    w[i]     = 0.0;
    Ap[i]    = 0.0;
    z[i]     = 0.0;
    fric1[i] = 1.0;

    if (i < nc)
    {
      fric[i]  = mu[i] * fric1[i];
      stat[i]    = 0;
      statusi[i] = 0;

    }

  }



  cblas_dcopy(n, problem->q, incx, r, incy);

  alphaf = -1.;
  betaf  = -1.;

  cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alphaf, vec, n, reaction, incx, betaf, r, incy);





  /*             Check for initial status             */


  for (i = 0; i < nc; i++)
  {
    mu[i] = fric[i];
    if (reaction[2 * i] <= eps)
    {
      /*       No contact            */
      stat[i] = 0;
    }
    else if (reaction[2 * i + 1] <=  -mu[i]*reaction[2 * i])
    {
      /*     Slide backward         */
      stat[i] = 1;
    }
    else if (reaction[2 * i + 1] >=  mu[i]*reaction[2 * i])
    {
      /*   Slide forward          */
      stat[i] = 3;
    }
    else
    {
      /*     Stick contact        */
      stat[i] = 2;
    }
  }


  iter  = 0;
  normr = 1.0;



  while ((iter < maxit) && (normr > tol))
  {



    for (i = 0 ; i < nc ; i++)
      statusi[i] = stat[i];


    cblas_dcopy(n, r, incx, v, incy);

    if (iter == 0)
    {
      cblas_dcopy(n, r, incx, w, incy);

      cblas_dcopy(n, w, incx, p, incy);
    }

    alphaf = 1.0;
    betaf  = 0.0;
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alphaf, vec, n, p, incx, betaf, Ap, incy);

    pAp    = cblas_ddot(n, p, incx, Ap, incy);

    /*}
    else
    {
    alphaf = 1.0;
    betaf  = 0.0;
    dgemv_( &notrans, (integer *)&n, (integer *)&n, &alphaf, vec, (integer *)&n, p, &incx, &betaf, Ap, &incy );

    pAp    = ddot_( (integer *)&n, p, &incx, Ap, &incy );*/

    if (pAp == 0)
    {
      if (verbose > 0)
        printf("\n Operation non conform alpha at the iteration %d \n", iter);

      free(r);
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

      *info = 2;

      return;
    }

    /*} */

    rp     = cblas_ddot(n, r, incx, p, incy);

    alpha  = rp / pAp;

    cblas_dcopy(n, reaction, incx, xi, incy);

    alphaf = alpha;
    cblas_daxpy(n, alphaf, p, incx, xi, incy);

    fc2d_projc(xi, &n, statusi, p, fric, reaction, stat);


    /*         r(:)=b(:)-matmul(A,x)          */

    cblas_dcopy(n, problem->q, incx, r, incy);

    alphaf = -1.;
    betaf  = -1.;
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, alphaf, vec, n, reaction, incx, betaf, r, incy);

    fc2d_projf(statusi, &n, r, fric, w);

    fc2d_projf(statusi, &n, p, fric, z);


    wAp    = cblas_ddot(n, w, incx, Ap, incy);

    beta   = - wAp / pAp;

    cblas_dcopy(n, w, incx, p, incy);

    alphaf  = beta;
    cblas_daxpy(n, alphaf, z, incx, p, incy);


    /*  alphaf  = 1.;
    betaf   = 0.;
    dgemv_( &notrans, (integer *)&n, (integer *)&n, &alphaf, vec , (integer *)&n, p, &incx, &betaf, Ap, &incy );

    pAp     = ddot_( (integer *)&n, p, &incx, Ap, &incy );*/

    cblas_dcopy(n, r, incx, xi, incy);

    alphaf  = -1.;
    cblas_daxpy(n, alphaf, v, incx, xi, incy);

    num     = cblas_ddot(n, xi, incx, xi, incy);

    den     = cblas_ddot(n, v, incx, v, incy);

    normr   = sqrt(num / den);

    it_end  = iter;
    res     = normr;


    options->iparam[SICONOS_IPARAM_ITER_DONE] = it_end;
    options->dparam[SICONOS_DPARAM_RESIDU] = res;


    iter = iter + 1;

  }




  if (normr < tol)
  {

    if (verbose > 0)
      printf("convergence after %d iterations with a residual %g\n", iter - 1, normr);

    *info = 0;


  }
  else
  {
    if (verbose > 0)
      printf("no convergence after %d iterations with a residual %g\n", iter - 1, normr);

    *info = 1;
  }


  alpha = -1.;
  cblas_dscal(n , alpha , r , incx);

  cblas_dcopy(n, r, incx, velocity, incy);



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




}

