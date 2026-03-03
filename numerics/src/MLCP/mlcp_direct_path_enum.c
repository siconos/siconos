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

/*
|A C| |u| |a| |0|
|   |*| |+| |=| |
|D B| |v| |b| |w|
0<z*v>0
dim(u)=mm
dim(v)=nn

**************************************************************************/
#include "mlcp_direct_path_enum.h"

/* Solver registration system */
#include "utils/solver_registry.h"
#include "utils/numerics_errors.h"

#include <stdio.h>  // for printf

#include "MLCP_Solvers.h"                       // for mixedLinearComplement...
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "SolverOptions.h"                      // for SolverOptions
#include "mlcp_direct.h"                        // for mlcp_direct_getNbDWork
#include "mlcp_path_enum.h"                     // for mlcp_path_enum, mlcp_...
static int sN;
static int sM;

static int* siWorkPathEnum = 0;
static int* siWorkDirect = 0;
static double* sdWorkPathEnum = 0;
static double* sdWorkDirect = 0;

void mlcp_direct_path_enum_init(MixedLinearComplementarityProblem* problem,
                                SolverOptions* options) {
  sN = problem->n;
  sM = problem->m;
  int iOffset = mlcp_direct_getNbIWork(problem, options);
  int dOffset = mlcp_direct_getNbDWork(problem, options);
  siWorkPathEnum = options->iWork + iOffset;
  siWorkDirect = options->iWork;
  sdWorkPathEnum = options->dWork + dOffset;
  sdWorkDirect = options->dWork;
  mlcp_direct_init(problem, options);
  options->dWork = sdWorkPathEnum;
  options->iWork = siWorkPathEnum;
  mlcp_path_enum_init(problem, options);
}
void mlcp_direct_path_enum_reset() {
  mlcp_direct_reset();
  mlcp_path_enum_reset();
  siWorkPathEnum = 0;
  siWorkDirect = 0;
  sdWorkPathEnum = 0;
  sdWorkDirect = 0;
}

void mlcp_direct_path_enum(MixedLinearComplementarityProblem* problem, double* z, double* w,
                           int* info, SolverOptions* options) {
  if (!siWorkPathEnum) {
    *info = 1;
    printf(
        "MLCP_DIRECT_PATH_ENUM error, call a non initialised method!!!!!!!!!!!!!!!!!!!!!\n");
    return;
  }
  /*First, try direct solver*/
  options->dWork = sdWorkDirect;
  options->iWork = siWorkDirect;
  mlcp_direct(problem, z, w, info, options);
  if (*info) {
    options->dWork = sdWorkPathEnum;
    options->iWork = siWorkPathEnum;
    /*solver direct failed, so run the enum solver.*/
    mlcp_path_enum(problem, z, w, info, options);
    if (!(*info)) {
      mlcp_direct_addConfigFromWSolution(problem, w + sN);
    }
  }
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers SICONOS_MLCP_DIRECT_PATH_ENUM in the global solver registry.
 */

static void mlcp_direct_path_enum_set_default(SolverOptions* options) {
  /* Use mlcp_direct defaults as the base (for NUMBER_OF_CONFIGURATIONS, etc.) */
  mlcp_direct_set_default(options);
}

static int mlcp_direct_path_enum_init_wrap(void* problem, SolverOptions* options) {
  mlcp_direct_path_enum_init((MixedLinearComplementarityProblem*)problem, options);
  return NUMERICS_OK;
}

static int mlcp_direct_path_enum_solve_wrap(void* problem, double* z, double* w, SolverOptions* options) {
  int info = NUMERICS_OK;
  mlcp_direct_path_enum((MixedLinearComplementarityProblem*)problem, z, w, &info, options);
  return info;
}

static void mlcp_direct_path_enum_free_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  mlcp_direct_path_enum_reset();
}

REGISTER_SOLVER(SICONOS_MLCP_DIRECT_PATH_ENUM, "MLCP_DIRECT_PATH_ENUM",
                "Direct-PATH-Enum hybrid solver for Mixed Linear Complementarity Problems",
                mlcp_direct_path_enum_init_wrap,
                mlcp_direct_path_enum_solve_wrap,
                mlcp_direct_path_enum_free_wrap,
                NULL,  /* error function */
                mlcp_direct_path_enum_set_default,
                1000,  /* default_max_iter */
                1e-12, /* default_tol */
                0      /* is_local_solver */);
