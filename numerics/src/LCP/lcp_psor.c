/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include <float.h>                         // for DBL_EPSILON
#include <math.h>                          // for fabs, fmax
#include <stdio.h>                         // for printf, NULL
#include <stdlib.h>                        // for free, malloc
#include "LCP_Solvers.h"                   // for lcp_compute_error, lcp_psor
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsFwd.h"                   // for SolverOptions, LinearCompl...
#include "NumericsMatrix.h"                // for NumericsMatrix
#include "SolverOptions.h"                 // for SolverOptions, solver_opti...
#include "lcp_cst.h"                       // for SICONOS_LCP_DPARAM_RHO
#include "numerics_verbose.h"              // for verbose
#include "SiconosBlas.h"                   // for cblas_dcopy, cblas_ddot


/*\warning omega is not explicitely used. must be completed    */
void lcp_psor(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  /* matrix M/vector q of the lcp */
  double * M = problem->M->matrix0;

  double * q = problem->q;

  /* size of the LCP */
  int n = problem->size;
  int incx = 1, incy = 1;
  int incxn;
  int i, iter;

  double qs, err;
  double *ww, *diag;
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tol = options->dparam[SICONOS_DPARAM_TOL];
  double omega = options->dparam[SICONOS_LCP_DPARAM_RHO]; // Not yet used
  printf("Warning : omega %f is not used !!!!!\n", omega);

  incxn = n;

  /* Initialize output */

  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;

  /* Allocation */

  ww   = (double*)malloc(n * sizeof(double));
  diag = (double*)malloc(n * sizeof(double));

  /* Check for non trivial case */

  qs = cblas_dnrm2(n, q, incx);

  if(verbose > 0) printf("\n ||q||= %g \n", qs);

  // Note FP : den never used ...
  //if (qs > DBL_EPSILON) den = 1.0 / qs;
  //else
  if(qs <= DBL_EPSILON)
  {
    for(i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }

    free(ww);
    free(diag);

    *info = 0;
    return;
  }

  for(i = 0 ; i < n ; ++i)
  {
    ww[i] = 0.;
  }

  cblas_dcopy(n, q, incx, w, incy);
  /* Intialization of w and z */
  /*if(initmethod == 0) {*/
  /* dcopy_( (integer *)&n , q , (integer *)&incx , w , (integer *)&incy );*/
  /*dcopy_( (integer *)&n , ww , (integer *)&incx , z , (integer *)&incy );*/
  /*}    else if(initmethod == 0) { */


  /* Preparation of the diagonal of the inverse matrix */

  for(i = 0 ; i < n ; ++i)
  {
    if(fabs(M[i * n + i]) < DBL_EPSILON)
    {

      if(verbose > 0)
      {
        printf(" Warning negative diagonal term \n");
        printf(" The local problem cannot be solved \n");
      }

      *info = 2;
      free(diag);
      free(ww);

      return;
    }
    else diag[i] = 1.0 / M[i * n + i];
  }

  /*start iterations*/

  iter = 0;
  err  = 1.;

  qs   = -1.0;

  cblas_dcopy(n, q, incx, w, incy);

  while((iter < itermax) && (err > tol))
  {

    ++iter;

    cblas_dcopy(n, w, incx, ww, incy);       /* w --> ww */
    cblas_dcopy(n, q, incx, w, incy);        /* q --> w */

    for(i = 0 ; i < n ; ++i)
    {

      z[i] = 0.0;

      /*      zi = -( q[i] + ddot_( (integer *)&n , &M[i] , (integer *)&incxn , z , (integer *)&incy ))*diag[i]; */

      /*       if( zi < 0 ) z[i] = 0.0;  */
      /*       else z[i] = zi; */

      z[i] = fmax(0.0, -(q[i] + cblas_ddot(n, &M[i], incxn, z, incy)) * diag[i]);

    }

    /* **** Criterium convergence **** */
    lcp_compute_error(problem, z, w, tol, &err);

    /* **** ********************* **** */
  }

  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  options->dparam[SICONOS_DPARAM_RESIDU] = err;

  if(err > tol)
  {
    printf("Siconos/Numerics: lcp_psor: No convergence of PSOR after %d iterations\n", iter);
    printf("Siconos/Numerics: lcp_psor: The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if(verbose > 0)
    {
      printf("Siconos/Numerics: lcp_psor: Convergence of PSOR after %d iterations\n", iter);
      printf("Siconos/Numerics: lcp_psor: The residue is : %g \n", err);
    }
    *info = 0;
  }


  free(ww);
  free(diag);

  return;
}
void lcp_psor_set_default(SolverOptions* options)
{
  options->dparam[SICONOS_LCP_DPARAM_RHO] = 0.1;
}
