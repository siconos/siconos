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

#include "mlcp_path_enum.h"

/* Solver registration system */
#include <stdio.h>  // for printf

#include "MLCP_Solvers.h"  // for mixedLinearComplement...
#include "mlcp_cst.h"
#include "mlcp_enum.h"  // for mlcp_enum_getNbDWork
#include "numerics_errors.h"
#include "solver_registry.h"
static int sN;
static int sM;

static int* siWorkEnum = 0;
static int* siWorkPath = 0;
static double* sdWorkEnum = 0;
static double* sdWorkPath = 0;

int mlcp_path_enum_getNbIWork(MixedLinearComplementarityProblem* problem,
                              SolverOptions* options) {
  return mlcp_enum_getNbIWork(problem, options);
}
int mlcp_path_enum_getNbDWork(MixedLinearComplementarityProblem* problem,
                              SolverOptions* options) {
  return mlcp_enum_getNbDWork(problem, options);
}

void mlcp_path_enum_init(MixedLinearComplementarityProblem* problem, SolverOptions* options) {
  sN = problem->n;
  sM = problem->m;
  int iOffset = 0; /* mlcp_path_getNbIWork(problem,options);*/
  int dOffset = 0; /*mlcp_path_getNbDWork(problem,options);*/
  siWorkEnum = options->iWork + iOffset;
  siWorkPath = options->iWork;
  sdWorkEnum = options->dWork + dOffset;
  sdWorkPath = options->dWork;
  /*  mlcp_path_init(problem, options);*/
}
void mlcp_path_enum_reset() {
  /*mlcp_path_reset();*/
  siWorkEnum = 0;
  siWorkPath = 0;
  sdWorkEnum = 0;
  sdWorkPath = 0;
}

void mlcp_path_enum(MixedLinearComplementarityProblem* problem, double* z, double* w,
                    int* info, SolverOptions* options) {
  if (!siWorkEnum) {
    *info = 1;
    printf("MLCP_PATH_ENUM error, call a non initialised method!!!!!!!!!!!!!!!!!!!!!\n");
    return;
  }
  /*First, try direct solver*/
  //  options->dWork = sdWorkDirect;
  //  options->iWork = siWorkDirect;
  mlcp_path(problem, z, w, info, options);
  if (*info) {
    printf("MLCP_PATH_ENUM: path failed, call enum\n");
    options->dWork = sdWorkEnum;
    options->iWork = siWorkEnum;
    /*solver direct failed, so run the enum solver.*/
    mlcp_enum(problem, z, w, info, options);
  }
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers SICONOS_MLCP_PATH_ENUM in the global solver registry.
 */

static void mlcp_path_enum_set_default(SolverOptions* options) {
  /* No specific defaults for path_enum solver */
}

static int mlcp_path_enum_init_wrap(void* problem, SolverOptions* options) {
  mlcp_path_enum_init((MixedLinearComplementarityProblem*)problem, options);
  return NUMERICS_OK;
}

static int mlcp_path_enum_solve_wrap(void* problem, double* z, double* w,
                                     SolverOptions* options) {
  int info = NUMERICS_OK;
  mlcp_path_enum((MixedLinearComplementarityProblem*)problem, z, w, &info, options);
  return info;
}

static void mlcp_path_enum_free_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  mlcp_path_enum_reset();
}

REGISTER_SOLVER(SICONOS_MLCP_PATH_ENUM, "MLCP_PATH_ENUM",
                "PATH-Enum hybrid solver for Mixed Linear Complementarity Problems",
                mlcp_path_enum_init_wrap, mlcp_path_enum_solve_wrap, mlcp_path_enum_free_wrap,
                NULL,                             /* error function */
                mlcp_path_enum_set_default, 1000, /* default_max_iter */
                1e-6,                             /* default_tol */
                0 /* is_local_solver */);
