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

#include <stdio.h>   // for printf
#include <stdlib.h>  // for free

#include "NumericsFwd.h"                                        // for Varia...
#include "SOCLCP_Solvers.h"                                     // for soclc...
#include "SecondOrderConeLinearComplementarityProblem.h"        // for Secon...
#include "SecondOrderConeLinearComplementarityProblem_as_VI.h"  // for Secon...
#include "SiconosBlas.h"                                        // for cblas...
#include "SolverOptions.h"                                      // for Solve...
#include "VariationalInequality.h"                              // for Varia...
#include "VariationalInequality_Solvers.h"                      // for varia...
#include "numerics_verbose.h"
#include "soclcp_compute_error.h"

/* Solver registration system */
#include "utils/solver_registry.h"
#include "utils/numerics_errors.h"                               // for soclc...

void soclcp_VI_ExtraGradient(SecondOrderConeLinearComplementarityProblem *problem,
                             double *reaction, double *velocity, int *info,
                             SolverOptions *options) {
  /* Dimension of the problem */
  int n = problem->n;

  VariationalInequality *vi = (VariationalInequality *)malloc(sizeof(VariationalInequality));

  // vi.self = &vi;
  vi->F = &Function_VI_SOCLCP;
  vi->ProjectionOnX = &Projection_VI_SOCLCP;

  double error = 1e24;

  SecondOrderConeLinearComplementarityProblem_as_VI *soclcp_as_vi =
      (SecondOrderConeLinearComplementarityProblem_as_VI *)malloc(
          sizeof(SecondOrderConeLinearComplementarityProblem_as_VI));
  vi->env = soclcp_as_vi;
  vi->size = n;

  /*Set the norm of the VI to the norm of problem->q  */
  vi->normVI = cblas_dnrm2(n, problem->q, 1);
  vi->istheNormVIset = 1;

  soclcp_as_vi->vi = vi;
  soclcp_as_vi->soclcp = problem;
  /* soclcp_display(fc3d_as_vi->fc3d); */

  variationalInequality_ExtraGradient(vi, reaction, velocity, info, options);

  /* **** Criterium convergence **** */
  soclcp_compute_error(problem, reaction, velocity, options->dparam[0], options, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);
   * printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  if (verbose > 0) {
    printf(
        "--------------- SOCLCP - VI Extra Gradient (VI_EG) - #Iteration %i Final Residual = "
        "%14.7e\n",
        options->iparam[SICONOS_IPARAM_MAX_ITER], options->dparam[SICONOS_DPARAM_RESIDU]);
  }
  free(vi);

  free(soclcp_as_vi);
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers SICONOS_SOCLCP_VI_EG in the global solver registry.
 */

static void soclcp_vi_eg_set_default(SolverOptions* options) {
  /* Call VI_EG set_default for proper initialization */
  variationalInequality_ExtraGradient_set_default(options);
}

static int soclcp_vi_eg_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  /* set_default already called by solver_options_create */
  return NUMERICS_OK;
}

static int soclcp_vi_eg_solve_wrap(void* problem, double* reaction,
                                   double* velocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  soclcp_VI_ExtraGradient(
      (SecondOrderConeLinearComplementarityProblem*)problem, reaction, velocity, &info, options);
  return info;
}

static void soclcp_vi_eg_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_SOCLCP_VI_EG, "SOCLCP_VI_EG",
                "Variational Inequality Extra Gradient for SOCLCP",
                soclcp_vi_eg_init_wrap,
                soclcp_vi_eg_solve_wrap,
                soclcp_vi_eg_free_wrap,
                NULL,  /* error function */
                soclcp_vi_eg_set_default,  /* set_default */
                1000,  /* default_max_iter */
                1e-4,  /* default_tol */
                0      /* is_local_solver */);
