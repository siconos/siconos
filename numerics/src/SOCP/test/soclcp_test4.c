// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdio.h>
#include <stdlib.h>

#include "NonSmoothDrivers.h"
#include "SOCLCP_Solvers.h"
#include "SOCLCP_cst.h"
#include "SolverOptions.h"
#include "soclcp_test_function.h"

int main(void) {
  int info = 0;
  printf("Test on ./data/Capsules-i122-1617.dat \n");

  FILE* finput = fopen("./data/Capsules-i122-1617.dat", "r");
  SolverOptions* options = (SolverOptions*)malloc(sizeof(SolverOptions));
  info = soclcp_setDefaultSolverOptions(options, SICONOS_SOCLCP_VI_FPP);
  options->dparam[0] = 1e-10;
  options->iparam[0] = 2000000;
  info = soclcp_test_function(finput, options);

  solver_options_delete(options);
  options = NULL;

  fclose(finput);
  printf("\nEnd of test on ./data/Capsules-i122-1617.dat \n");
  return info;
}
