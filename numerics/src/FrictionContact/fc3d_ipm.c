/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "CSparseMatrix_internal.h"

#include "FrictionContactProblem.h"  // for FrictionContactProblem, friction...
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
#include "SiconosLapack.h"

#include "SparseBlockMatrix.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "numerics_verbose.h"
#include "NumericsVector.h"
#include "float.h"
#include "JordanAlgebra.h"

#include "NumericsSparseMatrix.h"
#include "NumericsMatrix.h"

#include "projectionOnCone.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"
#include "gfc3d_ipm.h"

const char* const   SICONOS_FRICTION_3D_IPM_STR = "FC3D IPM";


void fc3d_IPM(FrictionContactProblem* restrict problem, double* restrict reaction,
               double* restrict velocity, int* restrict info, SolverOptions* restrict options)
{
  // to complete
}
