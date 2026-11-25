// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdio.h>  // for printf, fclose, fopen, FILE

#include "NumericsFwd.h"           // for SolverOptions
#include "SOCLCP_cst.h"            // for SICONOS_SOCLCP_NSGS, SICONOS_SOCLC...
#include "SolverOptions.h"         // for SolverOptions, solver_options_delete
#include "soclcp_test_function.h"  // for soclcp_test_function

int main(void) {
  int info = 0;
  printf("Test on ./data/Example1_SOCLCP_SBM.dat\n");

  FILE* finput = fopen("./data/Example1_SOCLCP_SBM.dat", "r");
  SolverOptions* options = solver_options_create(SICONOS_SOCLCP_NSGS);
  options->dparam[SICONOS_DPARAM_TOL] = 1e-16;
  solver_options_update_internal(options, 0, SICONOS_SOCLCP_ProjectionOnCone);
  printf(" atroizeoieo zi e %d %d", options->solverId, options->internalSolvers[0]->solverId);
  info = soclcp_test_function(finput, options);

  solver_options_delete(options);
  options = NULL;

  fclose(finput);
  printf("\nEnd of test on ./data/Example1_SOCLCP_SBM.dat\n");
  return info;
}
