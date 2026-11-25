// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdio.h>
#include <stdlib.h>

#include "FrictionContact_options.h"
#include "NonSmoothDrivers.h"
#include "frictionContact_test_utils.h"
#include "gfc3d_Solvers.h"

int main(void) {
  int info = 0;
  const char filename[50] = "./data/LMGC_GFC3D_CubeH8.hdf5";
  printf("Test on %s\n", filename);

  TestCase current;
  current.filename = filename;
  current.options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSGS);
  current.options->dparam[SICONOS_DPARAM_TOL] = 1e-08;
  current.options->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;

  info = globalFrictionContact_test_function(&current);

  solver_options_delete(current.options);
  printf("\nEnd of test on %s\n", filename);
  return info;
}
