// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file plasticity_2d_driver.c
 * \brief PLASTICITY_2D driver using the solver registration system
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "NumericsMatrix.h"
#include "PlasticityProblem.h"
#include "SiconosBlas.h"
#include "naming_conventions.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "plasticity_2d_compute_error.h"
#include "plasticity_2d_solvers.h"
#include "plasticity_2d_stress_coordinates.h"
#include "solver_registry.h"

/* ===========================================================================
 * Trivial Case Check
 * =========================================================================== */

int plasticity_2d_checkTrivialCase(PlasticityProblem* problem, double* strainrate,
                                   double* stress, SolverOptions* options) {
  (void)options;
  int nc = problem->numberOfCones;
  double* q = problem->q;
  int n = 3 * nc;

  /* sigma = 0 must solve every cone: Drucker-Prager, q + theta |q_T| e_N in the dual cone
   * (q_N >= 0 is not enough when theta < eta); Von Mises, 0 is inside the yield set: w = q
   * must vanish */
  for (int i = 0; i < nc; i++) {
    if (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER) {
      if (!plasticity_2d_DruckerPrager_zero_is_solution(
              &q[3 * i], problem->model.drucker_prager->eta[i],
              problem->model.drucker_prager->theta[i]))
        return -1;
    } else {
      if (fabs(q[3 * i]) > DBL_EPSILON || fabs(q[3 * i + 1]) > DBL_EPSILON ||
          fabs(q[3 * i + 2]) > DBL_EPSILON)
        return -1;
    }
  }
  for (int i = 0; i < n; ++i) {
    strainrate[i] = q[i];
    stress[i] = 0.;
  }

  numerics_printf("mc2d: trivial solution, stress = 0, strainrate = q.");
  return 0;
}

/* ===========================================================================
 * Registration-Based Driver
 * =========================================================================== */

static int plasticity_2d_driver_without_cohesion(PlasticityProblem* problem, double* stress,
                                                 double* strainrate, SolverOptions* options) {
  /* Input validation */
  if (!problem || !stress || !strainrate || !options) {
    numerics_error_log("plasticity_2d_driver", "null input argument");
    return NUMERICS_ERR_NULL_POINTER;
  }

  /* Check dimension */
  CHECK_DIMENSION(problem->dimension, 3);

  /* Initialize output */
  SET_SOLVER_ITER_DONE(options, 0);
  SET_SOLVER_RESIDUAL(options, 0.0);

  /* Check for trivial case */
  int trivial_status = plasticity_2d_checkTrivialCase(problem, strainrate, stress, options);
  if (trivial_status == 0) {
    return NUMERICS_OK;
  }

  /* Lookup solver in registry */
  const SolverEntry* solver = solver_registry_lookup(options->solverId);

  CHECK_COND(solver != NULL, NUMERICS_ERR_INVALID_SOLVER, "Solver not found");

  numerics_printf_verbose(1, "plasticity_2d_driver: using solver '%s' (%s)", solver->name,
                          solver->description);

  /* Validate solver is appropriate for this problem type */
  if (solver->is_local_solver) {
    numerics_printf(
        "plasticity_2d_driver: solver '%s' is a local solver, cannot be used as main solver",
        solver->name);
    return NUMERICS_ERR_INVALID_SOLVER;
  }

  /* Check solve function exists */
  CHECK_COND(solver->solve != NULL, NUMERICS_ERR_INVALID_SOLVER,
             "Solver has no solve function");

  /* Initialize solver if init function provided */
  if (solver->init) {
    int init_status = solver->init(problem, options);
    if (init_status != NUMERICS_OK) {
      numerics_printf("plasticity_2d_driver: solver initialization failed with error %d",
                      init_status);
      return init_status;
    }
  }

  /* Call the solver */
  numerics_printf_verbose(1, "plasticity_2d_driver: calling solver...");
  int solve_status = solver->solve(problem, stress, strainrate, options);

  /* Cleanup if needed */
  if (solver->free) {
    solver->free(problem, options);
  }

  /* Log result */
  if (solve_status == NUMERICS_OK) {
    numerics_printf_verbose(1, "plasticity_2d_driver: solver converged successfully");
  } else {
    numerics_printf_verbose(1, "plasticity_2d_driver: solver returned status %d",
                            solve_status);
  }

  return solve_status;
}

/* Cartesian stresses: the problem is written once for all in cone coordinates,
 * sigma_c = P s, M_c = P^{-T} M P^{-1}, q_c = P^{-T} q (the model parameters are the same),
 * solved, and the solution mapped back, s = P^{-1} sigma_c, w_s = P^T w_c. */
static int plasticity_2d_driver_cartesian(PlasticityProblem* problem, double* stress,
                                          double* strainrate, SolverOptions* options) {
  PlasticityStressCoordinates coordinates = problem->stress_coordinates;
  int nc = problem->numberOfCones;
  PlasticityProblem cone = *problem; /* shares the model parameters */
  cone.stress_coordinates = PLASTICITY_STRESS_COORDINATES_CONE;
  cone.M = plasticity_2d_matrix_to_cone(coordinates, problem->M);
  cone.q = (double*)malloc((size_t)(3 * nc) * sizeof(double));
  plasticity_2d_strainrate_to_cone(coordinates, nc, problem->q, cone.q);
  plasticity_2d_stress_to_cone(coordinates, nc, stress, stress);

  int info = plasticity_2d_driver(&cone, stress, strainrate, options);

  plasticity_2d_stress_from_cone(coordinates, nc, stress, stress);
  plasticity_2d_strainrate_from_cone(coordinates, nc, strainrate, strainrate);
  NM_free(cone.M);
  free(cone.q);
  return info;
}

/* Drucker-Prager with cohesion c: with a = c/eta e_N for each cone, sigma~ = sigma + a is
 * solution of the problem without cohesion with q~ = q - M a (w = M sigma + q is unchanged).
 * The problem is translated, solved without cohesion, and translated back. */
int plasticity_2d_driver(PlasticityProblem* problem, double* stress, double* strainrate,
                         SolverOptions* options) {
  if (!problem || !stress || !strainrate || !options) {
    numerics_error_log("plasticity_2d_driver", "null input argument");
    return NUMERICS_ERR_NULL_POINTER;
  }
  if (problem->stress_coordinates != PLASTICITY_STRESS_COORDINATES_CONE)
    return plasticity_2d_driver_cartesian(problem, stress, strainrate, options);
  double* shift = NULL;
  int err = plasticity_DruckerPrager_cohesion_shift(problem, &shift);
  if (err) return err;
  if (!shift)
    return plasticity_2d_driver_without_cohesion(problem, stress, strainrate, options);

  int n = 3 * problem->numberOfCones;
  double* q = problem->q;
  double* cohesion = problem->model.drucker_prager->cohesion;
  double* q_shifted = (double*)malloc((size_t)n * sizeof(double));
  /* q~ = q - M a, sigma~ = sigma + a */
  cblas_dcopy(n, q, 1, q_shifted, 1);
  NM_gemv(-1.0, problem->M, shift, 1.0, q_shifted);
  cblas_daxpy(n, 1.0, shift, 1, stress, 1);
  problem->q = q_shifted;
  problem->model.drucker_prager->cohesion = NULL; /* the solvers see the problem without
                                                     cohesion (error computation included) */

  int info = plasticity_2d_driver_without_cohesion(problem, stress, strainrate, options);

  problem->q = q;
  problem->model.drucker_prager->cohesion = cohesion;
  cblas_daxpy(n, -1.0, shift, 1, stress, 1);
  free(q_shifted);
  free(shift);
  return info;
}
