/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "NumericsMatrix.h"
#include "SolverOptions.h"
#include "fc2d_Solvers.h"
#include "fc3d_Solvers.h"

/* Independent reference implementations used to cross-check the public helpers. */

static int reference_status_3d_bnrt(const double velocity[3], const double reaction[3]) {
  double somme[3], dur[3];
  for (int k = 0; k < 3; ++k) {
    somme[k] = velocity[k] + reaction[k];
    dur[k] = velocity[k] - reaction[k];
  }
  const double cesp = sqrt(DBL_EPSILON);
  const double somme_norm =
      sqrt(somme[0] * somme[0] + somme[1] * somme[1] + somme[2] * somme[2]);
  const double somme_tang_norm = sqrt(somme[1] * somme[1] + somme[2] * somme[2]);
  const double dur_tang_norm = sqrt(dur[1] * dur[1] + dur[2] * dur[2]);

  const double ns = somme[0] - somme_tang_norm;
  if (ns <= cesp * somme_norm) return SICONOS_FRICTION_CONTACT_STATUS_T;
  if (dur[0] >= dur_tang_norm) return SICONOS_FRICTION_CONTACT_STATUS_N;
  if (-dur[0] >= dur_tang_norm) return SICONOS_FRICTION_CONTACT_STATUS_B;
  return SICONOS_FRICTION_CONTACT_STATUS_R;
}

static int reference_status_3d(const double velocity[3], const double reaction[3], double mu) {
  const double tol = 1e-10;
  const double rN = reaction[0];
  if (rN <= tol) return 0;
  const double rT_norm = hypot(reaction[1], reaction[2]);
  const double uT_norm = hypot(velocity[1], velocity[2]);
  const double radius = mu * rN;
  if (rT_norm < radius - tol && uT_norm <= tol) return 1;
  return 2;
}

static int reference_status_2d(const double velocity[2], const double reaction[2], double mu) {
  const double tol = 1e-10;
  const double rN = reaction[0];
  const double rT = reaction[1];
  const double uT = velocity[1];
  const double radius = mu * rN;
  if (rN <= tol) return 0;
  if (fabs(rT) < radius - tol && fabs(uT) <= tol) return 1;
  if (uT > tol) return 2;
  if (uT < -tol) return 3;
  return (rT < 0.0) ? 2 : 3;
}

static int check_status_consistency_3d(FrictionContactProblem* problem, const double* velocity,
                                       const double* reaction) {
  int nc = problem->numberOfContacts;
  int* statuses = (int*)malloc(nc * sizeof(int));
  int nB = -1, nN = -1, nR = -1, nT = -1;

  FrictionContactProblem_get_contact_statuses_3d_bnrt(velocity, reaction, (size_t)nc,
                                                      statuses);
  FrictionContactProblem_count_contact_statuses_3d_bnrt(velocity, reaction, (size_t)nc, &nB,
                                                        &nN, &nR, &nT);

  int count_B = 0, count_N = 0, count_R = 0, count_T = 0;
  int failures = 0;
  for (int i = 0; i < nc; ++i) {
    const double* v = &velocity[3 * i];
    const double* r = &reaction[3 * i];
    const double mu = problem->mu[i];

    int status_bnrt = FrictionContactProblem_compute_contact_status_3d_bnrt(v, r);
    int status_vi = FrictionContactProblem_compute_contact_status_3d(v, r, mu);
    int ref_bnrt = reference_status_3d_bnrt(v, r);
    int ref_vi = reference_status_3d(v, r, mu);

    if (status_bnrt != ref_bnrt) {
      printf("3D contact %d: BNRT helper %d != reference %d\n", i, status_bnrt, ref_bnrt);
      failures++;
    }
    if (status_vi != ref_vi) {
      printf("3D contact %d: VI helper %d != reference %d\n", i, status_vi, ref_vi);
      failures++;
    }
    if (status_bnrt != statuses[i]) {
      printf("3D contact %d: batch BNRT status %d != per-contact %d\n", i, statuses[i],
             status_bnrt);
      failures++;
    }

    /* Cross-check the two conventions against each other on a physical solution. */
    if (status_vi == 0 && status_bnrt != SICONOS_FRICTION_CONTACT_STATUS_N) {
      printf("3D contact %d: VI no-contact but BNRT status %d (expected N)\n", i, status_bnrt);
      failures++;
    }
    if (status_vi == 1 && status_bnrt != SICONOS_FRICTION_CONTACT_STATUS_B) {
      printf("3D contact %d: VI stick but BNRT status %d (expected B)\n", i, status_bnrt);
      failures++;
    }
    if (status_vi == 2 && status_bnrt != SICONOS_FRICTION_CONTACT_STATUS_R) {
      printf("3D contact %d: VI slip but BNRT status %d (expected R)\n", i, status_bnrt);
      failures++;
    }

    switch (status_bnrt) {
      case SICONOS_FRICTION_CONTACT_STATUS_B:
        count_B++;
        break;
      case SICONOS_FRICTION_CONTACT_STATUS_N:
        count_N++;
        break;
      case SICONOS_FRICTION_CONTACT_STATUS_R:
        count_R++;
        break;
      case SICONOS_FRICTION_CONTACT_STATUS_T:
        count_T++;
        break;
    }
  }

  if (count_B != nB || count_N != nN || count_R != nR || count_T != nT) {
    printf(
        "3D BNRT counts mismatch: per-contact (%d,%d,%d,%d) vs count helper (%d,%d,%d,%d)\n",
        count_B, count_N, count_R, count_T, nB, nN, nR, nT);
    failures++;
  }

  printf("3D status summary: nc=%d, B=%d, N=%d, R=%d, T=%d (failures=%d)\n", nc, nB, nN, nR,
         nT, failures);

  free(statuses);
  return failures;
}

static int check_status_consistency_2d(FrictionContactProblem* problem, const double* velocity,
                                       const double* reaction) {
  int nc = problem->numberOfContacts;
  int failures = 0;
  int count_no = 0, count_stick = 0, count_slip_p = 0, count_slip_m = 0;

  for (int i = 0; i < nc; ++i) {
    const double* v = &velocity[2 * i];
    const double* r = &reaction[2 * i];
    const double mu = problem->mu[i];

    int status = FrictionContactProblem_compute_contact_status_2d(v, r, mu);
    int ref = reference_status_2d(v, r, mu);

    if (status != ref) {
      printf("2D contact %d: helper %d != reference %d\n", i, status, ref);
      failures++;
    }

    switch (status) {
      case 0:
        count_no++;
        break;
      case 1:
        count_stick++;
        break;
      case 2:
        count_slip_p++;
        break;
      case 3:
        count_slip_m++;
        break;
      default:
        printf("2D contact %d: unexpected status %d\n", i, status);
        failures++;
    }
  }

  printf(
      "2D status summary: nc=%d, no-contact=%d, stick=%d, slip+=%d, slip-=%d (failures=%d)\n",
      nc, count_no, count_stick, count_slip_p, count_slip_m, failures);

  return failures;
}

static int test_3d() {
  printf("\n=== 3D friction contact status consistency ===\n");

  /* Build a tiny 3D problem with one no-contact, one stick and one slip contact.
   * M is the identity, so the exact solution is known by construction. */
  const int nc = 3;
  const int n = nc * 3;
  NumericsMatrix* M = NM_eye((size_t)n);
  double q_data[] = {/* no-contact: r = 0, v = q  */
                     1.0, 0.0, 0.0,
                     /* stick:      r = -q, v = 0 */
                     -1.0, 0.0, 0.0,
                     /* slip:       r = (1, -0.5, 0), v = (0, 1, 0) */
                     -1.0, 1.5, 0.0};
  double mu_data[] = {0.5, 0.5, 0.5};
  double* q = (double*)malloc(n * sizeof(double));
  double* mu = (double*)malloc(nc * sizeof(double));
  memcpy(q, q_data, n * sizeof(double));
  memcpy(mu, mu_data, nc * sizeof(double));

  FrictionContactProblem* problem =
      frictionContactProblem_new_with_data(3, (size_t)nc, M, q, mu);
  if (!problem) {
    printf("Failed to create 3D problem\n");
    return 1;
  }

  double* reaction = (double*)calloc(n, sizeof(double));
  double* velocity = (double*)calloc(n, sizeof(double));

  SolverOptions* options = solver_options_create(SICONOS_FRICTION_3D_NSGS);
  options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  solver_options_update_internal(options, 0, SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);

  int info = -1;
  fc3d_nsgs(problem, reaction, velocity, &info, options);

  if (info != 0) {
    printf("3D solver did not converge (info=%d, residual=%e)\n", info,
           options->dparam[SICONOS_DPARAM_RESIDU]);
  }

  int failures = check_status_consistency_3d(problem, velocity, reaction);

  solver_options_delete(options);
  free(reaction);
  free(velocity);
  frictionContactProblem_free(problem);
  return failures;
}

static int test_2d() {
  printf("\n=== 2D friction contact status consistency ===\n");

  /* Build a tiny 2D problem with one no-contact, one stick and one slip+ contact.
   * M is the identity, so the exact solution is known by construction. */
  const int nc = 3;
  const int n = nc * 2;
  NumericsMatrix* M = NM_eye((size_t)n);
  double q_data[] = {/* no-contact: r = 0, v = q  */
                     1.0, 0.0,
                     /* stick:      r = -q, v = 0 */
                     -1.0, 0.0,
                     /* slip+:      r = (1, -0.3), v = (0, 0.5) */
                     -1.0, 0.8};
  double mu_data[] = {0.3, 0.3, 0.3};
  double* q = (double*)malloc(n * sizeof(double));
  double* mu = (double*)malloc(nc * sizeof(double));
  memcpy(q, q_data, n * sizeof(double));
  memcpy(mu, mu_data, nc * sizeof(double));

  FrictionContactProblem* problem =
      frictionContactProblem_new_with_data(2, (size_t)nc, M, q, mu);
  if (!problem) {
    printf("Failed to create 2D problem\n");
    return 1;
  }

  double* reaction = (double*)calloc(n, sizeof(double));
  double* velocity = (double*)calloc(n, sizeof(double));

  SolverOptions* options = solver_options_create(SICONOS_FRICTION_2D_NSGS);
  options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;

  int info = -1;
  fc2d_nsgs(problem, reaction, velocity, &info, options);

  if (info != 0) {
    printf("2D solver did not converge (info=%d, residual=%e)\n", info,
           options->dparam[SICONOS_DPARAM_RESIDU]);
  }

  int failures = check_status_consistency_2d(problem, velocity, reaction);

  solver_options_delete(options);
  free(reaction);
  free(velocity);
  frictionContactProblem_free(problem);
  return failures;
}

int main() {
  int failures = 0;
  failures += test_3d();
  failures += test_2d();

  if (failures == 0) {
    printf("\nAll contact-status consistency checks passed.\n");
    return 0;
  } else {
    printf("\n%d contact-status consistency checks failed.\n", failures);
    return 1;
  }
}
