// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdio.h>
#include <stdlib.h>

#include "FrictionContact_options.h"
#include "NonSmoothDrivers.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"
#include "frictionContact_test_utils.h"

int main(void) {
  printf("\n Start of test on Default SolverOptions\n");
  int info = 0;
  SolverOptions* options = NULL;
  int solvers[] = {SICONOS_FRICTION_3D_NSGS, SICONOS_FRICTION_3D_NSGSV,
                   SICONOS_FRICTION_3D_PROX, SICONOS_FRICTION_3D_TFP,
                   SICONOS_FRICTION_3D_DSFP, SICONOS_FRICTION_3D_EG,
                   SICONOS_FRICTION_3D_HP};

  int n_solvers = (int)(sizeof(solvers) / sizeof(solvers[0]));

  for (int s = 0; s < n_solvers; ++s) {
    options = solver_options_create(solvers[s]);
    solver_options_print(options);
    solver_options_delete(options);
    options = NULL;
  }
  printf("\n End of test on Default SolverOptions\n");
  return info;
}
