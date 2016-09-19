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
#include <math.h>
#include "VariationalInequality_Solvers.h"
#include "VariationalInequality_computeError.h"

#include "NonSmoothDrivers.h"
#include "SiconosBlas.h"
#include "SiconosSets.h"
#include "misc.h"

char *  SICONOS_VI_EG_STR = "VI_EG";
char *  SICONOS_VI_FPP_STR = "VI_FPP";
char *  SICONOS_VI_HP_STR = "VI_HP";
char *  SICONOS_VI_BOX_QI_STR = "Box VI solver based on Qi C-function";
char *  SICONOS_VI_BOX_AVI_LSA_STR = "Box VI solver based on the Newton-Josephy method";
char *  SICONOS_VI_BOX_PATH_STR = "Box VI solver based on PATH solver";

/* #define DEBUG_MESSAGES */
/* #define DEBUT_STDOUT */
#include "debug.h"


void snPrintf(int level, SolverOptions* opts, const char *fmt, ...);

int variationalInequality_driver(VariationalInequality* problem, 
                                 double *x, double *w, 
                                 SolverOptions* options)
{
  if (options == NULL)
    numericsError("variationalInequality_driver", "null input for solver and/or global options");

  assert(options->isSet);
  if (verbose > 0)
    solver_options_print(options);

  /* Solver name */
  /*char * name = options->solverName;*/

  int info = -1 ;

  /* Check for trivial case */
  info = checkTrivialCase_vi(problem, x, w, options);
  if (info == 0){
    double error;
    variationalInequality_computeError(problem, x , w, options->dparam[0], options, &error);
    printf("variationalInequality_driver. error = %8.4e\n", error);
    return info;
  }

  switch (options->solverId)
  {
    /* Non Smooth Gauss Seidel (NSGS) */
  /* Extra Gradient algorithm */
  case SICONOS_VI_EG:
  {
    snPrintf(1, options, 
             " ========================== Call ExtraGradient (EG) solver for VI problem ==========================\n");
    variationalInequality_ExtraGradient(problem, x , w , &info , options);
    break;
  }
  case SICONOS_VI_FPP:
  {
    snPrintf(1, options, 
             " ========================== Call Fixed Point Projection (FPP) solver for VI problem ==========================\n");
    variationalInequality_FixedPointProjection(problem, x , w , &info , options);
    break;
  }
  case SICONOS_VI_HP:
  {
    snPrintf(1, options, 
             " ========================== Call Hyperplane Projection (HP) solver for VI problem ==========================\n");
    variationalInequality_HyperplaneProjection(problem, x , w , &info , options);
    break;
  }
  case SICONOS_VI_BOX_QI:
  {
    variationalInequality_box_newton_QiLSA(problem, x, w, &info, options);
    break;
  }
  case SICONOS_VI_BOX_AVI_LSA:
  {
    vi_box_AVI_LSA(problem, x, w, &info, options);
    break;
  }
  case SICONOS_VI_BOX_PATH:
  {
    vi_box_path(problem, x, w, &info, options);
    break;
  }
  default:
  {
    fprintf(stderr, "Numerics, variationalInequality_driver failed. Unknown solver.\n");
    exit(EXIT_FAILURE);

  }
  }

  return info;

}

int checkTrivialCase_vi(VariationalInequality* problem, double* x, 
                        double* w, SolverOptions* options)
{
  int n = problem->size;
  if (problem->ProjectionOnX)
  {
    problem->ProjectionOnX(problem,x,w);
  }
  else
  {
    cblas_dcopy(problem->size, x, 1, w, 1);
    project_on_set(problem->size, w, problem->set);
  }
  cblas_daxpy(n, -1.0,x, 1, w , 1);
  double nnorm = cblas_dnrm2(n,w,1);
  DEBUG_PRINTF("checkTrivialCase_vi, nnorm = %6.4e\n",nnorm);
  
  if (nnorm > fmin(options->dparam[0], 1e-12))
    return 1;

  problem->F(problem,n,x,w);
  nnorm = cblas_dnrm2(n,w,1);
  DEBUG_PRINTF("checkTrivialCase_vi, nnorm = %6.4e\n",nnorm);

  if (nnorm > fmin(options->dparam[0], 1e-12))
    return 1;

  if (verbose == 1)
    printf("variationalInequality driver, trivial solution F(x) = 0, x in X.\n");
  return 0;
}
