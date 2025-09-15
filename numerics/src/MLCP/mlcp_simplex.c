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
#include "mlcp_simplex.h"  // for mlcp_simplex_init, mlcp_simplex_reset

#include "MLCP_Solvers.h"   // for mixedLinearComplementarity_simplex_setDefa...
#include "NumericsFwd.h"    // for MixedLinearComplementarityProblem, SolverO...
#include "SiconosConfig.h"  // for HAVE_MLCPSIMPLEX  // IWYU pragma: keep

#ifdef HAVE_MLCPSIMPLEX

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "SiconosCompat.h"
/*import external implementation*/
#include "external_mlcp_simplex.h"

static int sIsInitialize = 0;
#endif

void mlcp_simplex_init(MixedLinearComplementarityProblem* problem, SolverOptions* options) {
#ifdef HAVE_MLCPSIMPLEX
  int nn = problem->n;
  int mm = problem->m;
  extern_mlcp_simplex_init_with_M(&nn, &mm, problem->M->matrix0);
  sIsInitialize = 1;
#endif
}
void mlcp_simplex_reset() {
#ifdef HAVE_MLCPSIMPLEX
  extern_mlcp_simplex_stop();
  sIsInitialize = 0;
#endif
}

void mlcp_simplex(MixedLinearComplementarityProblem* problem, double* z, double* w, int* info,
                  SolverOptions* options) {
#ifdef HAVE_MLCPSIMPLEX
  //  double tol ;
  //  double * workingFloat=options->dWork;
  //  int * workingInt=options->iWork;
  //  int lin;
  //  int npm=(problem->n)+(problem->m);
  //  int npm2 = npm*npm;
  //  int NRHS=1;
  //  int one=1;
  //  int * ipiv;
  //  int check;
  //  int DGESVinfo=1;
  int nn = problem->n;
  int mm = problem->m;

  if (!sIsInitialize) extern_mlcp_simplex_init_with_M(&nn, &mm, problem->M->matrix0);

  extern_mlcp_simplex(problem->q, problem->q + nn, z, z + nn, w, info, options->iparam,
                      options->dparam);

  if (!sIsInitialize) extern_mlcp_simplex_stop();
#endif
}
