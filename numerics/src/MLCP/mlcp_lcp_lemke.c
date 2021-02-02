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

//#include <float.h>                            // for DBL_EPSILON
//#include <math.h>                             // for fabs


#include <stdio.h>                              // for printf
#include <stdlib.h>                             // for free, malloc, exit
#include <string.h>                             // for free, malloc, exit
#include "MLCP_Solvers.h"                       // for mlcp_compute_error
#include "LCP_Solvers.h"                        // for mlcp_compute_error
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "mlcp_cst.h"                           // for SICONOS_MLCP_LCP_LEMKE
#include "lcp_cst.h"                            // for SICONOS_LCP_LEMKE
#include "mlcp_to_lcp.h"                        // for mlcp_to_lcp
#include "NumericsMatrix.h"                     // for storageType
#include "numerics_verbose.h"                   // for numerics_printf


/* #define DEBUG_MESSAGES */
#include "debug.h"

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"                     // for NV_display
#endif

void mlcp_lcp_lemke(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  DEBUG_BEGIN("mlcp_lcp_lemke(...)\n");
  assert(problem);
  assert(z);
  assert(w);

  DEBUG_EXPR(mixedLinearComplementarity_display(problem););
  LinearComplementarityProblem* lcp =  mlcp_to_lcp(problem);

  if(!lcp)
  {
    numerics_printf_verbose(0, "mlcp_lcp_lemke: mlcp_to_lcp conversion failed");
    *info = 1;
    DEBUG_PRINT("mlcp_lcp_lemke: mlcp_to_lcp conversion failed\n");
    DEBUG_END("mlcp_lcp_lemke(...)\n");
    return;
  }
  DEBUG_EXPR(linearComplementarity_display(lcp););

  int n = problem->n;
  int m = problem->m;

  double * z_lcp = z+ n;
  double * w_lcp = w+ n;

  options->solverId = SICONOS_LCP_LEMKE;
  lcp_lexicolemke(lcp, z_lcp, w_lcp, info, options);

  DEBUG_EXPR(
    double *  lcp_error;
    lcp_compute_error(lcp, z_lcp, w_lcp,
                      options->dparam[SICONOS_DPARAM_TOL], lcp_error);
    printf("lcp_error = %12.8e\n",  *lcp_error);
    );
  if(*info)
  {
    numerics_printf_verbose(0, "mlcp_lcp_lemke: lcp_lexicolemke failed");
    freeLinearComplementarityProblem(lcp);
    return;
  }

  DEBUG_EXPR(NV_display(z_lcp,m));
  DEBUG_EXPR(NV_display(w_lcp,m));


  /* compute u  = -A^{1}(Cv + a)*/
  double *  u = z;
  memcpy(u,problem->a,problem->n*sizeof(double));
  NumericsMatrix * C= NM_new();
  NM_fill(C, NM_DENSE, n, m, problem->C);
  NM_gemv(-1.0, C, z_lcp, -1.0, u);
  double * A_copy = (double* )calloc(n*n,sizeof(double));
  memcpy(A_copy, problem->A, n*n*sizeof(double)); /* we preserve the original A */
  NumericsMatrix * A = NM_new();
  NM_fill(A, NM_DENSE, n, n, A_copy);
  NM_gesv_expert(A, u, NM_KEEP_FACTORS);




  /* compute error (this also recompute w !!) */

  double tol   = options->dparam[SICONOS_DPARAM_TOL];
  double error;
  mlcp_compute_error(problem, z, w, tol, &error);

  if(error > tol)
  {
    numerics_printf_verbose(1,"---- MLCP - LCP(LEMKE)  - error = %14.7e ", error);
    *info = 1;
  }
  else
  {
    numerics_printf_verbose(1,"---- MLCP - LCP(LEMKE)  - error = %14.7e ", error);
    *info = 0;
  }

  options->dparam[SICONOS_DPARAM_RESIDU] = error;

  options->solverId = SICONOS_MLCP_LCP_LEMKE;


  freeLinearComplementarityProblem(lcp);
  NM_clear(A);
  C->matrix0=NULL;
  NM_clear(C);
  DEBUG_END("mlcp_lcp_lemke(...)\n");
  return;
}

void mlcp_lcp_lemke_default(SolverOptions* options)
{
  options->filterOn = false;

}
