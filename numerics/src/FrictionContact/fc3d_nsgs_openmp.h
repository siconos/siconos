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
#include <SiconosConfig.h>
#if defined(WITH_OPENMP) && defined(_OPENMP)
#define USE_OPENMP 1
#include <omp.h>
#endif


/** pointer to function used to update local problem */
typedef void (*Update_indexPtr)(int, FrictionContactProblem*, FrictionContactProblem*, double*, SolverOptions *,
                                unsigned int * , unsigned int);

void snPrintf(int level, SolverOptions* opts, const char *fmt, ...);


void fc3d_nsgs_openmp_redblack(FrictionContactProblem* problem, double *reaction,
                               double *velocity, int* info, SolverOptions* options);

void fc3d_nsgs_openmp_for(FrictionContactProblem* problem, double *reaction,
                               double *velocity, int* info, SolverOptions* options);

void fc3d_nsgs_openmp_ddm_naive(FrictionContactProblem* problem, double *reaction,
                               double *velocity, int* info, SolverOptions* options);

void fc3d_nsgs_error_comparison(FrictionContactProblem* problem, double *reaction,
                               double *velocity, int* info, SolverOptions* options);

void fc3d_nsgs_index_initialize_local_solver(SolverPtr* solve, Update_indexPtr* update,
                                             FreeSolverNSGSPtr* freeSolver, ComputeErrorPtr* computeError,
                                             FrictionContactProblem* problem,
                                             FrictionContactProblem* localproblem,
                                             SolverOptions * options, SolverOptions * localsolver_options);

void fc3d_nsgs_index_computeqLocal(FrictionContactProblem * problem,
                                   double *reaction, int contact,
                                   unsigned int * index, unsigned int index_size,
                                   double * qLocal);

void fc3d_nsgs_index(FrictionContactProblem* problem,
                     double *reaction, double *velocity,
                     int* info, SolverOptions* options,
                     unsigned int* index, unsigned int index_size);
