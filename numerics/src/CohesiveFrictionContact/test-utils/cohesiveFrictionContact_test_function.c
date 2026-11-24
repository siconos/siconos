// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <math.h>    // for isfinite
#include <stdio.h>   // for printf
#include <stdlib.h>  // for calloc, free
#include <time.h>    // for clock

#include "CohesiveFrictionContactProblem.h"      // for cohesiveFrictionContact_newFromFilename
#include "SiconosBlas.h"                         // for cblas_dnrm2
#include "SolverOptions.h"                       // for SolverOptions, SICONOS_DPARAM_RESIDU
#include "cohesiveFrictionContact_test_utils.h"  // for cohesiveFrictionContact_test_function
#include "cohesive_friction_3d_driver.h"         // for cohesive_friction_3d_driver
#include "test_utils.h"                          // for TestCase

int cohesiveFrictionContact_test_function(TestCase* current) {
  CohesiveFrictionContactProblem* problem =
      cohesiveFrictionContact_newFromFilename((char*)current->filename);
  if (!problem) {
    printf("test failure: unable to read %s\n", current->filename);
    return 1;
  }

  int dim = problem->dimension;
  int nc = problem->numberOfContacts;
  int ncoh = problem->numberOfCohesivePoints;
  int n = dim * (nc + ncoh);

  double* reaction = (double*)calloc(n, sizeof(double));
  double* velocity = (double*)calloc(n, sizeof(double));

  solver_options_print(current->options);

  clock_t t1 = clock();
  int info = cohesive_friction_3d_driver(problem, reaction, velocity, current->options);
  clock_t t2 = clock();

  for (int k = 0; k < n; ++k) {
    info = info == 0 ? !(isfinite(velocity[k]) && isfinite(reaction[k])) : info;
  }

  printf("Norm velocity:  %12.8e\n", cblas_dnrm2(n, velocity, 1));
  printf("Norm reaction:  %12.8e\n", cblas_dnrm2(n, reaction, 1));

  int print_size = n < 10 ? n : 10;
  printf("First values (%i)\n", print_size);
  for (int k = 0; k < print_size; k++) {
    printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k,
           reaction[k]);
  }
  printf(" ..... \n");

  printf("test %s, residual = %9.2e, info = %d, nb iter = %i\n", info ? "failure" : "success",
         current->options->dparam[SICONOS_DPARAM_RESIDU], info,
         current->options->iparam[SICONOS_IPARAM_ITER_DONE]);
  printf("\nsumry: %d  %9.2e  %5i  %10.4f", info,
         current->options->dparam[SICONOS_DPARAM_RESIDU],
         current->options->iparam[SICONOS_IPARAM_ITER_DONE],
         (double)(t2 - t1) / (double)CLOCKS_PER_SEC);
  printf("%3i %5i %5i     %s\n\n", dim, nc, ncoh, current->filename);

  free(reaction);
  free(velocity);
  cohesiveFrictionContactProblem_free(problem);

  return info;
}
