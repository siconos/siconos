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
#include "NonSmoothDrivers.h"

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
  set_SolverOptions(options->internalSolvers, SICONOS_LCP_LEMKE);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_PGS);
  printSolverOptions(options);
  deleteSolverOptions(options);


  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_RPGS);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_QP);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_NSQP);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_CPG);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_PSOR);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_LATIN);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_LATIN_W);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_LEMKE);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_PATH);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_ENUM);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_NEWTONMIN);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_AVI_CAOFERRIS);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_PIVOT);
  printSolverOptions(options);
  deleteSolverOptions(options);


  freeLinearComplementarityProblem(problem);
  free(options);


  printf("\n End of test on Default SolverOptions\n");
  return info;
}
