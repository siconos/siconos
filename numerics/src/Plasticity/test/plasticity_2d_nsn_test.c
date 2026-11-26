// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/* Global semi-smooth Newton (PLASTICITY_2D_NSN) on the data files of the Plasticity tests:
   Drucker-Prager footing problems (dilatancy theta > 0: natural-map merit; theta = 0: QP
   merit), Von Mises footing, small examples, FEM strip problems. Both globalizations:
   full Newton steps with watchdog (default; the one-cone footing_1 problem cycles and
   needs the watchdog) and merit function only. Checks the convergence at 1e-10 from a zero
   initial guess, and the error recomputed by plasticity_2d_compute_error. */

#include <stdio.h>
#include <stdlib.h>

#include "NonSmoothDrivers.h"
#include "NumericsVerbose.h"
#include "PlasticityProblem.h"
#include "Plasticity_options.h"
#include "SolverOptions.h"
#include "plasticity_2d_compute_error.h"

#define TOLERANCE 1e-10
#define MAX_ITER 300

static int test_file(const char* filename, int globalization) {
  PlasticityProblem* problem = plasticity_new_from_filename(filename);
  if (!problem) return 1;
  int n = 3 * problem->numberOfCones;
  double* stress = (double*)calloc((size_t)n, sizeof(double));
  double* w = (double*)calloc((size_t)n, sizeof(double));
  SolverOptions* options = solver_options_create(PLASTICITY_2D_NSN);
  options->dparam[SICONOS_DPARAM_TOL] = TOLERANCE;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = MAX_ITER;
  options->iparam[PLASTICITY_NSN_GLOBALIZATION] = globalization;
  int info = plasticity_2d_driver(problem, stress, w, options);
  double error = 0.;
  plasticity_2d_compute_error(problem, stress, w, TOLERANCE, options, 0., &error);
  printf("%-50s %-8s %4d cones: info %d, %3d iterations, error %.2e\n", filename,
         globalization == PLASTICITY_NSN_GLOBALIZATION_MERIT ? "merit" : "watchdog",
         problem->numberOfCones, info, options->iparam[SICONOS_IPARAM_ITER_DONE], error);
  int failed = info != 0 || error > TOLERANCE;
  solver_options_delete(options);
  free(stress);
  free(w);
  plasticityProblem_free(problem);
  return failed;
}

int main(void) {
  numerics_set_verbose(0);
  int failed = 0;
  const char* files[] = {"./data/plasticity_2d_example1.dat",
                         "./data/plasticity_2d_example1_theta0.dat",
                         "./data/plasticity_2d_footing_1.dat",
                         "./data/plasticity_2d_footing_100.dat",
                         "./data/plasticity_2d_footing_100_theta0.dat",
                         "./data/plasticity_2d_footing_100_theta0.05.dat",
                         "./data/plasticity_2d_vonmises_example1.dat",
                         "./data/plasticity_2d_vonmises_footing_100.dat",
                         "./data/plasticity_2d_vonmises_strip_4.dat",
                         "./data/plasticity_2d_vonmises_strip_10.dat",
                         "./data/plasticity_2d_vonmises_strip_100.dat"};
  for (size_t i = 0; i < sizeof(files) / sizeof(files[0]); i++) {
    failed += test_file(files[i], PLASTICITY_NSN_GLOBALIZATION_WATCHDOG);
    failed += test_file(files[i], PLASTICITY_NSN_GLOBALIZATION_MERIT);
  }
  printf("failures: %d\n", failed);
  return failed ? 1 : 0;
}
