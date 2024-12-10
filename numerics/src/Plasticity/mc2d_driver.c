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
#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <stdio.h>   // for NULL
#include <string.h>  // for NULL

#include "MohrCoulomb2DProblem.h"                    // for MohrCoulomb2DProblem...
#include "Plasticity_cst.h"                              // for SICONOS_FRICTI...
#include "NonSmoothDrivers.h"                          // for fc3d_driver
#include "NumericsFwd.h"                               // for SolverOptions
#include "SolverOptions.h"                             // for SolverOptions
#include "mc2d_solvers.h" 
#include "numerics_verbose.h"                          // for numerics_printf

const char* const SICONOS_MOHR_COULOMB_2D_NSGS_STR = "MC2D_NSGS";


int mc2d_driver(MohrCoulomb2DProblem* problem, double* stress, double* plastic_strain_rate,
                SolverOptions* options) {
  if (options == NULL) numerics_error("mc2d_driver", "null input for solver options");

  assert(options->isSet); /* true(1) if the SolverOptions structure has been filled in else
                             false(0) */

  if (verbose > 1) solver_options_print(options);

  int info = -1;

  if (problem->dimension != 3)
    numerics_error(
        "mc2d_driver",
        "Dimension of the problem : problem-> dimension is not compatible or is not set");

  /* Check for trivial case */
  info = mc2d_checkTrivialCase(problem, plastic_strain_rate, stress, options);
  if (info == 0) {
    /* If a trivial solution is found, we set the number of iterations to 0
       and the reached acuracy to 0.0 .
    */
    options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
    options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;
    goto exit;
  }

  switch (options->solverId) {
    /* Non Smooth Gauss Seidel (NSGS) */
    case MOHR_COULOMB_2D_NSGS: {
      numerics_printf(
          " ========================== Call NSGS solver for Mohr Coulomb 2D problem "
          "==========================\n");
      mc2d_nsgs(problem, stress, plastic_strain_rate, &info, options);
      break;
    }
    default: {
      char msg[200];
      strcpy(msg, "Unknown solver : ");
      strcat(msg, solver_options_id_to_name(options->solverId));
      strcat(msg, "\n");
      numerics_warning("mc2d_driver", msg);
      numerics_error("mc2d_driver", msg);
      info = 1;
    }
  }

exit:
  return info;
}

int mc2d_checkTrivialCase(MohrCoulomb2DProblem* problem, double* plastic_strain_rate, double* stress,
                          SolverOptions* options) {
  /* Number of contacts */
  int nc = problem->numberOfCones;
  double* q = problem->q;
  /* Dimension of the problem */
  int n = 3 * nc;
  int i = 0;
  /*take off? R=0 ?*/
  for (i = 0; i < nc; i++) {
    if (q[3 * i] < -DBL_EPSILON) return -1;
  }
  for (i = 0; i < n; ++i) {
    plastic_strain_rate[i] = q[i];
    stress[i] = 0.;
  }

  numerics_printf(
      "mc2d mc2d_checkTrivialCase,  trivial solution stress = 0, plastic_strain_rate = q.\n");
  return 0;
}
