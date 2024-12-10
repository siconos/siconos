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

#include <math.h>  // for fmax

#include "MohrCoulomb2DProblem.h"  // for MohrCoulomb2DProblem
#include "Plasticity_cst.h"            // for PLASTICITY_IPARAM_INTER...
#include "NumericsFwd.h"             // for SolverOptions, MohrCoulomb2DPr...
#include "SolverOptions.h"           // for SolverOptions
#include "mc2d_solvers.h"            // for mc2d_set_internalsolver_tolerance
#include "numerics_verbose.h"        // for numerics_printf_verbose, numeric...

void mc2d_set_internalsolver_tolerance(MohrCoulomb2DProblem* problem, SolverOptions* options,
                                       SolverOptions* internalsolver_options, double error) {
  int* iparam = options->iparam;
  if (iparam[PLASTICITY_IPARAM_INTERNAL_ERROR_STRATEGY] ==
      PLASTICITY_INTERNAL_ERROR_STRATEGY_ADAPTIVE) {
    internalsolver_options->dparam[SICONOS_DPARAM_TOL] =
        fmax(error / options->dparam[PLASTICITY_DPARAM_INTERNAL_ERROR_RATIO],
             options->dparam[SICONOS_DPARAM_TOL] / problem->numberOfCones);
    numerics_printf_verbose(2,
                            "mc2d_FixedPoint_set_internalsolver_tolerance - Internal solver "
                            "tolerance is set to %e\n",
                            internalsolver_options->dparam[SICONOS_DPARAM_TOL]);
  } else if (iparam[PLASTICITY_IPARAM_INTERNAL_ERROR_STRATEGY] ==
             PLASTICITY_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT) {
    internalsolver_options->dparam[SICONOS_DPARAM_TOL] =
        error / (options->dparam[PLASTICITY_DPARAM_INTERNAL_ERROR_RATIO] *
                 problem->numberOfCones);
    numerics_printf_verbose(2,
                            "mc2d_FixedPoint_set_internalsolver_tolerance - Internal solver "
                            "tolerance is set to %e",
                            internalsolver_options->dparam[SICONOS_DPARAM_TOL]);
  } else if (iparam[PLASTICITY_IPARAM_INTERNAL_ERROR_STRATEGY] ==
             PLASTICITY_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE) {
    // We use the user value for the error of the local solver
    numerics_printf_verbose(2,
                            "mc2d_FixedPoint_set_internalsolver_tolerance - Internal solver "
                            "tolerance is set to %e",
                            internalsolver_options->dparam[SICONOS_DPARAM_TOL]);
  } else {
    numerics_error("mc2d__set_internalsolver_tolerance",
                   "Unknown strategy for driving the tolerance");
  }
}
