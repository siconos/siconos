// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#define _XOPEN_SOURCE 700
#include <math.h>    // for isfinite
#include <stdio.h>   // for printf, fclose, fopen, FILE
#include <stdlib.h>  // for calloc, free

#include "NonSmoothDrivers.h"               // for rolling_friction_3d_driver
#include "NumericsFwd.h"                    // for RollingFrictionContactPro...
#include "NumericsVerbose.h"                // for numerics_set_verbose
#include "RollingFrictionContactProblem.h"  // for rollingFrictionContactPro...
#include "SiconosBlas.h"
#include "frictionContact_test_utils.h"  // for rollingFrictionContact_te...
#include "test_utils.h"                  // for TestCase
int rollingFrictionContact_test_function(TestCase* current) {
  int k;
  RollingFrictionContactProblem* problem =
      rollingFrictionContact_new_from_filename(current->filename);
  numerics_set_verbose(1);

  FILE* foutput = fopen("checkinput.dat", "w");
  rollingFrictionContact_printInFile(problem, foutput);

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;

  int info;
  double* reaction = (double*)calloc(dim * NC, sizeof(double));
  double* velocity = (double*)calloc(dim * NC, sizeof(double));
  printf("\n rollingFrictionContact_test_function run\n");
  if (dim == 2) {
    info = 1;
  } else if (dim == 5) {
    info = rolling_friction_3d_driver(problem, reaction, velocity, current->options);
  }
  for (int k = 0; k < dim * NC; ++k) {
    info = info == 0 ? !(isfinite(velocity[k]) && isfinite(reaction[k])) : info;
  }

  printf("Norm velocity:  %12.8e\n", cblas_dnrm2(NC * dim, velocity, 1));
  printf("Norm reaction:  %12.8e\n", cblas_dnrm2(NC * dim, reaction, 1));
  int print_size = 10;
  if (dim * NC >= print_size) {
    printf("First values (%i)\n", print_size);
    for (int k = 0; k < print_size; k++) {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k,
             reaction[k]);
    }
  } else {
    for (int k = 0; k < dim * NC; k++) {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k,
             reaction[k]);
    }
  }
  printf(" ..... \n");
  if (!info) {
    printf("test succeeded\n");
  } else {
    printf("test unsuccessful\n");
  }
  free(reaction);
  free(velocity);
  fclose(foutput);

  rollingFrictionContactProblem_free(problem);

  return info;
}
