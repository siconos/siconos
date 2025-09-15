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
#include "numerics_verbose.h"                                   // for verbose
#include "soclcp_compute_error.h"                               // for soclc...

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
