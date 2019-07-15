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

#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "MCP_cst.h"
#include "MixedComplementarityProblem.h"
#include "SolverOptions.h"
#include "MCP_Solvers.h"


void testF(int size, double *z, double * F);
void testF(int size, double *z, double * F)
{
  printf("call to MCP function F(z) ...\n");
}

void testNablaF(int size, double *z, double *F);
void testNablaF(int size, double *z, double *F)
{
  printf("call to MCP function nablaF(z) ...\n");
}

int main(void)
{
  printf(" Start tests for MCP solvers.\n");

  int info = 0 ;

  /* Set solver options */
  SolverOptions options;

  /* FB solver */
  options.solverId = SICONOS_MCP_OLD_FB;

  /* Create a MixedComplementarityProblem */
  MixedComplementarityProblem_old* problem = (MixedComplementarityProblem_old *)malloc(sizeof(MixedComplementarityProblem_old));

  problem->sizeEqualities = 2;
  problem->sizeInequalities = 3;
  problem->computeFmcp = &testF ;
  problem->computeNablaFmcp = &testNablaF ;
  problem->Fmcp = NULL;
  problem->nablaFmcp = NULL;

  mcp_old_setDefaultSolverOptions(problem, &options);

  int size = 5;
  double z[4];
  double F[4];
  double nablaF[16];
  problem->computeFmcp(size, z, F);
  problem->computeNablaFmcp(size, z, nablaF);

  /* Initialize the solver */
  mcp_old_driver_init(problem, &options) ;

  /// TODO : write a real test ... ////

  printf("End of MCP solvers test. \n");
  mcp_old_driver_reset(problem, &options);
  mixedComplementarityProblem_old_free(problem);
  solver_options_delete(&options);

  return info;
}
