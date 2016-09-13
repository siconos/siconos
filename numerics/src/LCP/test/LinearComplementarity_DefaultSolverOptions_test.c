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
#include "SolverOptions.h"
#include "NumericsFwd.h"
#include "LCP_Solvers.h"
#include "lcp_cst.h"
#include "LinearComplementarityProblem.h"
#include "assert.h"

int main(void)
{
  printf("\n Start of test on Default SolverOptions\n");
  int info = 0 ;
  SolverOptions * options = (SolverOptions *)malloc(sizeof(SolverOptions));

  FILE * finput  =  fopen("./data/lcp_mmc.dat", "r");
  LinearComplementarityProblem* problem = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));

  info = linearComplementarity_newFromFile(problem, finput);

  fclose(finput);


  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_NSGS_SBM);
  assert(options->internalSolvers);
  solver_options_set(options->internalSolvers, SICONOS_LCP_LEMKE);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_PGS);
  solver_options_print(options);
  solver_options_delete(options);


  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_RPGS);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_QP);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_NSQP);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_CPG);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_PSOR);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_LATIN);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_LATIN_W);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_LEMKE);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_PATH);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_ENUM);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_NEWTONMIN);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_AVI_CAOFERRIS);
  solver_options_print(options);
  solver_options_delete(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_PIVOT);
  solver_options_print(options);
  solver_options_delete(options);


  freeLinearComplementarityProblem(problem);
  free(options);


  printf("\n End of test on Default SolverOptions\n");
  return info;
}
