/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/* Several contacts of 2D rolling friction: rolling_friction_2d_driver
   (NSGS) with the default local solver and with ONECONTACT_ENUM on random
   chains of disks (rfc2d_test_tools.h), the law checked independently of
   rolling_friction_2d_compute_error; for mu_r = 0, the normal and tangent
   reactions against fc2d_driver (NSGS) where both converge.  NSGS may not
   converge on a few random chains (fc2d neither): at most 1% failures. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "NonSmoothDrivers.h"
#include "NumericsMatrix.h"
#include "RollingFrictionContactProblem.h"
#include "RollingFrictionContact_options.h"
#include "SolverOptions.h"
#include "rfc2d_test_tools.h"

int main(void) {
  const int sizes[] = {2, 5, 20};
  const int np = 500;
  int failed = 0;
  for (int use_enum = 0; use_enum <= 1; ++use_enum)
    for (int s = 0; s < 3; ++s) {
      int nc = sizes[s], n = 3 * nc, m2 = 2 * nc;
      srand(1);
      int not_converged = 0, law_failures = 0, fc2d_differs = 0;
      double worst = 0., worst_diff = 0.;
      double* W = malloc((size_t)n * n * sizeof(double));
      double* Wc = malloc((size_t)n * n * sizeof(double));
      double* W2 = malloc((size_t)m2 * m2 * sizeof(double));
      double *q = malloc(n * sizeof(double)), *q2 = malloc(m2 * sizeof(double));
      double *mu = malloc(nc * sizeof(double)), *mur = malloc(nc * sizeof(double));
      double *r = malloc(n * sizeof(double)), *u = malloc(n * sizeof(double));
      double *r2 = malloc(m2 * sizeof(double)), *u2 = malloc(m2 * sizeof(double));
      for (int k = 0; k < np; ++k) {
        int with_mur = k % 2;
        rfc2d_random_chain(nc, with_mur, W, q, mu, mur);
        for (int i = 0; i < n * n; ++i) Wc[i] = W[i];
        for (int i = 0; i < n; ++i) r[i] = u[i] = 0.;
        NumericsMatrix* M = NM_create_with_data(NM_DENSE, n, n, Wc);
        RollingFrictionContactProblem* P =
            rollingFrictionContactProblem_new_with_data(3, nc, M, q, mu, mur);
        SolverOptions* so = solver_options_create(SICONOS_ROLLING_FRICTION_2D_NSGS);
        so->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
        so->dparam[SICONOS_DPARAM_TOL] = 1e-12;
        if (use_enum) {
          solver_options_delete(so->internalSolvers[0]);
          so->internalSolvers[0] =
              solver_options_create(SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ENUM);
        }
        int info = rolling_friction_2d_driver(P, r, u, so);
        solver_options_delete(so);
        free(P);
        free(M);
        if (info) {
          not_converged++;
          continue;
        }
        double v = rfc2d_law_violation(nc, W, q, mu, mur, r);
        worst = fmax(worst, v);
        if (v > 1e-8) law_failures++;

        if (!with_mur) {
          for (int c = 0; c < nc; ++c)
            for (int i = 0; i < 2; ++i) {
              q2[2 * c + i] = q[3 * c + i];
              r2[2 * c + i] = u2[2 * c + i] = 0.;
              for (int d = 0; d < nc; ++d)
                for (int j = 0; j < 2; ++j)
                  W2[(2 * c + i) + (2 * d + j) * m2] = W[(3 * c + i) + (3 * d + j) * n];
            }
          NumericsMatrix* M2 = NM_create_with_data(NM_DENSE, m2, m2, W2);
          FrictionContactProblem* P2 = frictionContactProblem_new_with_data(2, nc, M2, q2, mu);
          SolverOptions* so2 = solver_options_create(SICONOS_FRICTION_2D_NSGS);
          so2->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
          so2->dparam[SICONOS_DPARAM_TOL] = 1e-12;
          int info2 = fc2d_driver(P2, r2, u2, so2);
          solver_options_delete(so2);
          free(P2);
          free(M2);
          if (!info2) {
            double diff = 0., norm = 1.;
            for (int c = 0; c < nc; ++c)
              for (int i = 0; i < 2; ++i) {
                diff = fmax(diff, fabs(r2[2 * c + i] - r[3 * c + i]));
                norm = fmax(norm, fabs(r2[2 * c + i]));
              }
            diff /= norm;
            worst_diff = fmax(worst_diff, diff);
            if (diff > 1e-6) fc2d_differs++;
          }
        }
      }
      int bad = not_converged + law_failures + fc2d_differs;
      int ok = bad <= np / 100;
      failed += !ok;
      printf(
          "%-8s nc=%2d: %d problems, not converged %d, law failures %d (worst %.2e), "
          "mu_r=0 against fc2d > 1e-6: %d (worst %.2e) %s\n",
          use_enum ? "enum" : "default", nc, np, not_converged, law_failures, worst,
          fc2d_differs, worst_diff, ok ? "ok" : "FAILED");
      free(W);
      free(Wc);
      free(W2);
      free(q);
      free(q2);
      free(mu);
      free(mur);
      free(r);
      free(u);
      free(r2);
      free(u2);
    }
  return failed ? 1 : 0;
}
