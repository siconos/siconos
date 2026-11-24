// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <math.h>  // for fmax

#include "FrictionContactProblem.h"   // for FrictionContactProblem
#include "FrictionContact_options.h"  // for SICONOS_FRICTION_3D_IPARAM_INTER...
#include "Friction_tools.h"
#include "NumericsFwd.h"    // for SolverOptions, FrictionContactPr...
#include "SolverOptions.h"  // for SolverOptions
#include "fc3d_Solvers.h"   // for fc3d_set_internalsolver_tolerance
#include "fc3d_short_names.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"

int fc3d_set_internalsolver_tolerance(int numberOfContacts, SolverOptions* options,
                                      SolverOptions* internalsolver_options, double error) {
  int* iparam = options->iparam;
  if (iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] ==
      SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE) {
    internalsolver_options->dparam[SICONOS_DPARAM_TOL] =
        fmax(error / options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO],
             options->dparam[SICONOS_DPARAM_TOL] / numberOfContacts);
    numerics_printf_verbose(2,
                            "fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver "
                            "tolerance is set to %e\n",
                            internalsolver_options->dparam[SICONOS_DPARAM_TOL]);
  } else if (iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] ==
             SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT) {
    internalsolver_options->dparam[SICONOS_DPARAM_TOL] =
        error /
        (options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] * numberOfContacts);
    numerics_printf_verbose(2,
                            "fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver "
                            "tolerance is set to %e",
                            internalsolver_options->dparam[SICONOS_DPARAM_TOL]);
  } else if (iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] ==
             SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE) {
    // We use the user value for the error of the local solver
    numerics_printf_verbose(2,
                            "fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver "
                            "tolerance is set to %e",
                            internalsolver_options->dparam[SICONOS_DPARAM_TOL]);
  } else {
    return numerics_error("fc3d__set_internalsolver_tolerance",
                          "Unknown strategy for driving the tolerance");
  }
  return 0;
}
