/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include "NCP_FixedP.h"

#include <stdlib.h>  // for free, malloc

#include "SiconosBlas.h"  // for cblas_dnrm2
#include "SolverOptions.h"
#include "fc3d_2NCP_Glocker.h"  // for compute_Z_GlockerFixedP
#include "numerics_verbose.h"   // for verbose
#include "stdio.h"              // for printf

/*============================ Fixed point Solver ==================================*/

int Fixe(int n, double* z, int* iparam, double* dparam) {
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];  // maximum number of iterations allowed
  int niter = 0;                                  // current iteration number
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  if (verbose > 0) {
    printf(" ============= Starting of Newton process =============\n");
    printf(" - tolerance: %f\n - maximum number of iterations: %i\n", tolerance, itermax);
  }

  int i;

  /* Connect F and its jacobian to input functions */
  //  setFuncEval(F);

  /* Memory allocation for phi and its jacobian */
  double* www = (double*)malloc(sizeof(double) * n);

  double terminationCriterion = 1;

  /** Iterations ... */
  while ((niter < itermax) && (terminationCriterion > tolerance)) {
    ++niter;

    // printf(" ============= Fixed Point Iteration ============= %i\n",niter);
    for (i = 0; i < n; ++i) compute_Z_GlockerFixedP(i, www);

    terminationCriterion = cblas_dnrm2(n, www, 1);
    // printf(" error = %14.7e\n", terminationCriterion);
    if (verbose > 0) {
      printf("Non Smooth Newton, iteration number %i, error equal to %14.7e .\n", niter,
             terminationCriterion);
      printf(" -----------");
    }
  }

  /* Total number of iterations */
  iparam[SICONOS_IPARAM_ITER_DONE] = niter;
  /* Final error */
  dparam[SICONOS_DPARAM_RESIDU] = terminationCriterion;

  /** Free memory*/
  free(www);

  if (verbose > 0) {
    if (dparam[SICONOS_DPARAM_RESIDU] > tolerance)
      printf("Non Smooth Newton warning: no convergence after %i iterations\n", niter);

    else
      printf("Non Smooth Newton: convergence after %i iterations\n", niter);
    printf(" The residue is : %e \n", dparam[SICONOS_DPARAM_RESIDU]);
  }

  if (dparam[SICONOS_DPARAM_RESIDU] > tolerance)
    return 1;
  else
    return 0;
}
