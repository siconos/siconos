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

#include <assert.h>                        // for assert
#include <float.h>                         // for DBL_EPSILON
#include <math.h>                          // for fabs
#ifndef __cplusplus
#include <stdbool.h>                       // for false
#endif
#include <stdio.h>                         // for printf
#include <stdlib.h>                        // for free, malloc
#include "LCP_Solvers.h"                   // for lcp_compute_error, lcp_pgs
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsFwd.h"                   // for SolverOptions, LinearCompl...
#include "NumericsMatrix.h"                // for NM_get_value, NM_row_prod_...
#include "SiconosBlas.h"                   // for cblas_dcopy, cblas_dnrm2
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_DPA...
#include "debug.h"                         // for DEBUG_PRINTF
#include "numerics_verbose.h"              // for verbose

void lcp_pgs(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  /* matrix M/vector q of the lcp */
  NumericsMatrix * M = problem->M;
  double * q = problem->q;

  assert(M);
  assert(q);
  assert(options->iSize > 1);
  assert(options->dSize > 1);

  int n = problem->size;

  int i;
  double qs, zi;

  /* Solver parameters */
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tol = options->dparam[SICONOS_DPARAM_TOL];
  /* Initialize output */

  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;

  if(verbose > 0)
  {
    qs = cblas_dnrm2(n, q, 1);
    printf("===== Starting of LCP solving with Projected Gauss-Seidel algorithm.\n");
    printf("\n ||q||= %g \n", qs);
    printf("itermax = %d \n", itermax);
    printf("tolerance = %g \n", tol);
  }

  /* Preparation of the diagonal of the inverse matrix */
  double * diag = (double*)malloc(n * sizeof(double));
  double diag_i = 0.0;
  for(i = 0 ; i < n ; ++i)
  {
    diag_i = NM_get_value(M,i,i);
    if(fabs(diag_i) < DBL_EPSILON)
    {
      if(verbose > 0)
      {
        printf("Numerics::lcp_pgs, error: vanishing diagonal term \n");
        printf(" The problem cannot be solved with this method \n");
      }

      *info = 2;
      free(diag);
      return;
    }
    else diag[i] = 1.0 / diag_i;
  }

  /* Iterations*/
  int iter = 0;
  double err  = 1.;

  while((iter < itermax) && (err > tol))
  {

    ++iter;

    /* Initialization of w with q */
    cblas_dcopy(n, q, 1, w, 1);

    for(i = 0 ; i < n ; ++i)
    {
      //z[i] = 0.0;

      zi =q[i];
      DEBUG_PRINTF("zi = %e\n", zi);
      NM_row_prod_no_diag1x1(n, i, i, M, z, &zi, false);
      DEBUG_PRINTF("diag[i] = %e\t zi = %e\n",diag[i], zi);
      zi = -(zi) * diag[i];


      if(zi < 0) z[i] = 0.0;
      else z[i] = zi;
      /* z[i]=fmax(0.0,-( q[i] + ddot_( (integer *)&n , &M[i] , (integer *)&incxn , z , (integer *)&incy ))*diag[i]);*/
    }
    /* **** Criterium convergence **** */
    lcp_compute_error(problem, z, w, tol, &err);

    if(verbose == 2)
    {
      printf(" # i%d -- %g : ", iter, err);
      for(i = 0 ; i < n ; ++i) printf(" %g", z[i]);
      for(i = 0 ; i < n ; ++i) printf(" %g", w[i]);
      printf("\n");
    }
  }

  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  options->dparam[SICONOS_DPARAM_RESIDU] = err;

  if(err > tol)
  {
    if(verbose > 0)
    {
      printf("Siconos/Numerics: lcp_pgs: No convergence of PGS after %d iterations\n", iter);
      printf("The residue is : %g \n", err);
    }
    *info = 1;
  }
  else
  {
    if(verbose > 0)
    {
      printf("Siconos/Numerics: lcp_pgs: Convergence of PGS after %d iterations\n", iter);
      printf("The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(diag);
}

