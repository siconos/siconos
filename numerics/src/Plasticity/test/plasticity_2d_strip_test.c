// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/* NSGS on Von Mises plasticity problems dumped from a FEM simulation.

   The data files plasticity_2d_vonmises_strip_<n>.dat are the problems of time step n of
   the strip with a circular hole (siconos-tutorials, mechanics/FEM/Plasticity/StripWithHole:
   60 T3 elements, plane stress, Voigt coordinates, physical yield stress 210), written by
   fem::nonsmooth_formulations::Plasticity::setProblemDump:
   - step 1: elastic (the solution is inside the yield set in every element);
   - step 4: first plastic element;
   - step 10: 18 plastic elements;
   - step 100: plastic plateau (23 plastic elements).
   They come from a quasi-static regime where M (the condensed iteration matrix) is
   ill-conditioned: NSGS needs 1e4-1e5 sweeps (the motivation for global solvers).

   Checks, for PLASTICITY_2D_NSGS and PLASTICITY_2D_NSGS_GENERIC: convergence at the
   required tolerance, error recomputed by plasticity_2d_compute_error, Von Mises stress
   <= yield stress, and the two solvers give the same stresses. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "NonSmoothDrivers.h"
#include "NumericsVerbose.h"
#include "PlasticityProblem.h"
#include "Plasticity_options.h"
#include "SolverOptions.h"
#include "plasticity_2d_compute_error.h"

#define TOLERANCE 1e-8
#define MAX_ITER 200000

static int solve(PlasticityProblem* problem, int solver_id, double* stress) {
  int n = 3 * problem->numberOfCones;
  double* strainrate = (double*)calloc((size_t)n, sizeof(double));
  for (int i = 0; i < n; i++) stress[i] = 0.;
  SolverOptions* options = solver_options_create(solver_id);
  options->dparam[SICONOS_DPARAM_TOL] = TOLERANCE;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = MAX_ITER;
  int info = plasticity_2d_driver(problem, stress, strainrate, options);
  int iter = options->iparam[SICONOS_IPARAM_ITER_DONE];
  double residual = options->dparam[SICONOS_DPARAM_RESIDU];

  double error = 0.;
  double* w = (double*)calloc((size_t)n, sizeof(double));
  plasticity_2d_compute_error(problem, stress, w, TOLERANCE, options, 0., &error);

  /* Von Mises stress (plane stress, Voigt) <= sigma_y */
  double max_ratio = 0.;
  int n_plastic = 0;
  for (int c = 0; c < problem->numberOfCones; c++) {
    double* s = &stress[3 * c];
    double vm = sqrt(s[0] * s[0] + s[1] * s[1] - s[0] * s[1] + 3. * s[2] * s[2]);
    double ratio = vm / problem->model.von_mises->sigma_y[c];
    if (ratio > max_ratio) max_ratio = ratio;
    if (ratio > 1. - 1e-6) n_plastic++;
  }
  printf(
      "  %-28s info %d, %6d iterations, residual %.2e, error %.2e, plastic %2d, max "
      "vm/sy - 1 = %.1e\n",
      solver_options_id_to_name(solver_id), info, iter, residual, error, n_plastic,
      max_ratio - 1.);
  int failed = info != 0 || error > TOLERANCE || max_ratio > 1. + 1e-10;
  free(w);
  free(strainrate);
  solver_options_delete(options);
  return failed;
}

static int test_file(const char* filename) {
  printf("%s\n", filename);
  PlasticityProblem* problem = plasticity_new_from_filename(filename);
  if (!problem || problem->model_type != PLASTICITY_MODEL_VON_MISES ||
      problem->model.von_mises->hypothesis != PLASTICITY_VONMISES_PLANE_STRESS ||
      problem->stress_coordinates != PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT) {
    printf("  unexpected problem (model, hypothesis or coordinates)\n");
    if (problem) plasticityProblem_free(problem);
    return 1;
  }
  int n = 3 * problem->numberOfCones;
  double* s1 = (double*)malloc((size_t)n * sizeof(double));
  double* s2 = (double*)malloc((size_t)n * sizeof(double));
  int failed = solve(problem, PLASTICITY_2D_NSGS, s1);
  failed += solve(problem, PLASTICITY_2D_NSGS_GENERIC, s2);

  /* same solution (unique: M positive definite), up to the tolerance */
  double diff = 0., norm = 0.;
  for (int i = 0; i < n; i++) {
    diff = fmax(diff, fabs(s1[i] - s2[i]));
    norm = fmax(norm, fabs(s1[i]));
  }
  printf("  max |sigma_NSGS - sigma_GENERIC| / max |sigma| = %.1e\n", diff / norm);
  if (diff > 1e3 * TOLERANCE * norm) failed++;
  free(s1);
  free(s2);
  plasticityProblem_free(problem);
  return failed;
}

int main(void) {
  numerics_set_verbose(0);
  int failed = 0;
  failed += test_file("./data/plasticity_2d_vonmises_strip_1.dat");
  failed += test_file("./data/plasticity_2d_vonmises_strip_4.dat");
  failed += test_file("./data/plasticity_2d_vonmises_strip_10.dat");
  failed += test_file("./data/plasticity_2d_vonmises_strip_100.dat");
  printf("failures: %d\n", failed);
  return failed ? 1 : 0;
}
