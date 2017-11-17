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
#include <string.h>
#include <time.h>
#include <float.h>
#include "MLCP_Solvers.h"
#include "SiconosCompat.h"
#include "NonSmoothDrivers.h"
#include "numerics_verbose.h"


void  mixedLinearComplementarity_default_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pOptions)
{
  pOptions->isSet = 0;
  pOptions->iSize = 10;
  pOptions->iparam = 0;
  pOptions->dSize = 10;
  pOptions->dparam = 0;
  pOptions->filterOn = 0;
  pOptions->dWork = 0;
  pOptions->iWork = 0;
  pOptions->iparam = (int*)malloc(10 * sizeof(int));
  pOptions->dparam = (double*)malloc(10 * sizeof(double));
  pOptions->numberOfInternalSolvers = 0;
  solver_options_nullify(pOptions);


  pOptions->dparam[0] = 10 - 7;
  /*default number of it*/
  pOptions->iparam[0] = 1000;
  /*enum case : do not use dgels*/
  pOptions->iparam[4] = 0;
  pOptions->iparam[5] = 3; /*Number of registered configurations*/
  pOptions->iparam[8] = 0; /*Prb nedd a update*/
  pOptions->dparam[5] = 1e-12; /*tol used by direct solver to check complementarity*/
  pOptions->dparam[6] = 1e-12; /*tol for direct solver to determinate if a value is positive*/

  int sizeOfIwork = mlcp_driver_get_iwork(problem, pOptions);
  if (sizeOfIwork)
    pOptions->iWork = (int*)malloc(sizeOfIwork * sizeof(int));
  int sizeOfDwork = mlcp_driver_get_dwork(problem, pOptions);
  if (sizeOfDwork)
    pOptions->dWork = (double*)malloc(sizeOfDwork * sizeof(double));



}

void  mixedLinearComplementarity_deleteDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pOptions)
{
  if (pOptions->iparam)
    free(pOptions->iparam);
  if (pOptions->dparam)
    free(pOptions->dparam);
  if (pOptions->iWork)
    free(pOptions->iWork);
  if (pOptions->dWork)
    free(pOptions->dWork);
  // FP : I comment the lines below, that results in failures
  // in some test-GMP-REDUCED3_3D_QUARTIC-GMP. 
  // Todo : fix this ...
  //if (pOptions->callback)
  //  free(pOptions->callback);
  solver_options_nullify(pOptions);

}

int mixedLinearComplementarity_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pOptions)
{
  solver_options_nullify(pOptions);
  int info = -1;

  switch (pOptions->solverId)
  {
  case SICONOS_MLCP_DIRECT_ENUM:
  {
    info =    mixedLinearComplementarity_directEnum_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_PATH_ENUM:
  {
    info =    mixedLinearComplementarity_pathEnum_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case  SICONOS_MLCP_DIRECT_PATH_ENUM:
  {
    info =    mixedLinearComplementarity_directPathEnum_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_DIRECT_SIMPLEX:
  {
    info =    mixedLinearComplementarity_directSimplex_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_DIRECT_PATH:
  {
    info =    mixedLinearComplementarity_directPath_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_DIRECT_FB:
  {
    info =    mixedLinearComplementarity_directFB_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_SIMPLEX:
  {
    info =    mixedLinearComplementarity_simplex_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_PGS:
  {
    info =    mixedLinearComplementarity_pgs_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_PGS_SBM:
  {
    info =    mixedLinearComplementarity_pgs_SBM_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_RPGS:
  {
    info =    mixedLinearComplementarity_rpgs_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_RPSOR:
  {
    info =    mixedLinearComplementarity_rpsor_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_PATH:
  {
    info =    mixedLinearComplementarity_path_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_ENUM:
  {
    info =    mixedLinearComplementarity_enum_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_FB:
  {
    info =    mixedLinearComplementarity_fb_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_PSOR:
  {
    info = mixedLinearComplementarity_psor_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  default:
  {
    numerics_error("mixedLinearComplementarity_setDefaultSolverOptions", "Unknown Solver");

  }
  }
  return info;
}

