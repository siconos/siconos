/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include <assert.h>            // for assert
#include <float.h>             // for DBL_EPSILON
#include <math.h>              // for fabs
#include <stdio.h>             // for printf
#include <stdlib.h>            // for free, malloc
#include "NumericsFwd.h"       // for RelayProblem, SolverOptions, NumericsM...
#include "NumericsMatrix.h"    // for NumericsMatrix
#include "RelayProblem.h"      // for RelayProblem
#include "Relay_Solvers.h"     // for relay_compute_error, relay_pgs
#include "SiconosBlas.h"       // for cblas_dcopy, cblas_ddot
#include "SolverOptions.h"     // for SolverOptions, SICONOS_DPARAM_RESIDU
#include "numerics_verbose.h"  // for verbose

void relay_pgs(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{


  double* M = problem->M->matrix0;
  double* q = problem->q;
  int n = problem -> size;
  double *a = problem->lb;
  double *b = problem->ub;

  assert(M);
  assert(q);
  assert(n);
  assert(a);
  assert(b);



  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tol = options->dparam[SICONOS_DPARAM_TOL];


  int i;
  double zi;
  double * diag = (double*)malloc(n * sizeof(double));


  for(i = 0 ; i < n ; ++i)
  {
    if(fabs(M[i * n + i]) < DBL_EPSILON)
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
    else diag[i] = 1.0 / M[i * n + i];
  }

  /* Iterations*/
  int iter = 0;
  double err  = 1.;
  *info = 1;

  while((iter < itermax) && (err > tol))
  {

    ++iter;

    /* Initialization of w with q */
    cblas_dcopy(n, q, 1, w, 1);

    for(i = 0 ; i < n ; ++i)
    {
      z[i] = 0.0;
      zi = -(q[i] + cblas_ddot(n, &M[i], n, z, 1)) * diag[i];
      z[i] = zi;
      if(zi < a[i]) z[i] = a[i];
      else if(zi > b[i]) z[i] = b[i];
    }
    /* **** Criterium convergence **** */
    *info = relay_compute_error(problem, z, w, tol, &err);

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
    printf("Siconos/Numerics: relay_pgs: No convergence of PGS after %d iterations\n", iter);
    printf("The residue is : %g \n", err);
    *info = 1;
  }
  else
  {
    if(verbose > 0)
    {
      printf("Siconos/Numerics: relay_pgs: Convergence of PGS after %d iterations\n", iter);
      printf("The residue is : %g \n", err);
    }
    *info = 0;
  }




  free(diag);



}

