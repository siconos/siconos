// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <math.h>  // for fmax

#include "FrictionContact_options.h"        // for SICONOS_FRICTION_3D_IPARA...
#include "NumericsFwd.h"                    // for SolverOptions, RollingFri...
#include "RollingFrictionContactProblem.h"  // for RollingFrictionContactPro...
#include "SolverOptions.h"                  // for SolverOptions
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "rolling_fc_Solvers.h"  // for rolling_friction_3d_set_internal...
#include "rolling_friction_3d_short_names.h"

int rolling_friction_2d_set_internalsolver_tolerance(RollingFrictionContactProblem* problem,
                                                     SolverOptions* options,
                                                     SolverOptions* internalsolver_options,
                                                     double error) {
  int* iparam = options->iparam;
  if (iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] ==
      SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE) {
    internalsolver_options->dparam[0] =
        fmax(error / options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO],
             options->dparam[0] / problem->numberOfContacts);
    numerics_printf_verbose(2,
                            "fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver "
                            "tolerance is set to %e\n",
                            internalsolver_options->dparam[0]);
  } else if (iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] ==
             SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT) {
    internalsolver_options->dparam[0] =
        error / (options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] *
                 problem->numberOfContacts);
    numerics_printf_verbose(2,
                            "fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver "
                            "tolerance is set to %e",
                            internalsolver_options->dparam[0]);
  } else if (iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] ==
             SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE) {
    // We use the user value for the error of the local solver
    numerics_printf_verbose(2,
                            "fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver "
                            "tolerance is set to %e",
                            internalsolver_options->dparam[0]);
  } else {
    return numerics_error("rolling_friction_3d__set_internalsolver_tolerance",
                          "Unknown strategy for driving the tolerance");
  }
  return 0;
}
