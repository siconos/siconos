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
/* One contact of 2D rolling friction, rolling_friction_2d_driver (NSGS)
   with the default local solver and with ONECONTACT_ENUM.
   1. Known answers: a disk (m = 1, R = 1, J = 1/2) on the line y = 0, one
      time step h of gravity, rows normal (0 1 0), tangent (1 0 R),
      rolling (0 0 1), mu = 0.5, mu_r = 0.1:
      A. rolling without slip (v = (v0, 0, -v0 / R)): r_n = h g, the
         rolling moment saturated r_r = mu_r r_n, the tangent sticking
         r_t = -(2/3) mu_r r_n;
      B. separating (q_n > 0): r = 0;
      C. sliding without rotation (v = (v0, 0, 0)): r_n = h g,
         r_t = -mu r_n, r_r = mu_r r_n.
   2. Random problems: the law checked by rfc2d_law_violation, no failure
      with ENUM, at most 1% with the default local solver. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "NonSmoothDrivers.h"
#include "NumericsMatrix.h"
#include "RollingFrictionContactProblem.h"
#include "RollingFrictionContact_options.h"
#include "SolverOptions.h"
#include "rfc2d_test_tools.h"

static int solve(double* W, double* q, double mu, double mur, int use_enum, double r[3]) {
  double u[3] = {0., 0., 0.};
  r[0] = r[1] = r[2] = 0.;
  NumericsMatrix* M = NM_create_with_data(NM_DENSE, 3, 3, W);
  RollingFrictionContactProblem* P =
      rollingFrictionContactProblem_new_with_data(3, 1, M, q, &mu, &mur);
  SolverOptions* so = solver_options_create(SICONOS_ROLLING_FRICTION_2D_NSGS);
  so->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  so->dparam[SICONOS_DPARAM_TOL] = 1e-12;
  if (use_enum) {
    solver_options_delete(so->internalSolvers[0]);
    so->internalSolvers[0] =
        solver_options_create(SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ENUM);
  }
  int info = rolling_friction_2d_driver(P, r, u, so);
  solver_options_delete(so);
  free(P);
  free(M); /* W and q belong to the caller */
  return info;
}

int main(void) {
  const double m = 1., R = 1., J = 0.5, h = 0.01, g = 9.81, v0 = 2., mu = 0.5, mur = 0.1;
  /* H rows n, t, r; M^-1 = diag(1/m, 1/m, 1/I) */
  const double H[3][3] = {{0., 1., 0.}, {1., 0., R}, {0., 0., 1.}};
  const double Minv[3] = {1. / m, 1. / m, 1. / J};
  double W[9];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      double s = 0.;
      for (int l = 0; l < 3; ++l) s += H[i][l] * Minv[l] * H[j][l];
      W[i + 3 * j] = s;
    }
  /* q = H v_free, v_free = v + h M^-1 f */
  double vA[3] = {v0, -h * g, -v0 / R}, vB[3] = {v0, 1., -v0 / R}, vC[3] = {v0, -h * g, 0.};
  double* vs[3] = {vA, vB, vC};
  double rn = m * h * g;
  double expected[3][3] = {
      {rn, -2. / 3. * mur * rn, mur * rn}, {0., 0., 0.}, {rn, -mu * rn, mur * rn}};
  const char* names[3] = {"A rolling", "B separating", "C sliding"};
  int failures = 0;
  for (int use_enum = 0; use_enum <= 1; ++use_enum)
    for (int k = 0; k < 3; ++k) {
      double q[3], r[3], Wc[9];
      for (int i = 0; i < 9; ++i) Wc[i] = W[i];
      for (int i = 0; i < 3; ++i)
        q[i] = H[i][0] * vs[k][0] + H[i][1] * vs[k][1] + H[i][2] * vs[k][2];
      int info = solve(Wc, q, mu, mur, use_enum, r);
      double err = 0.;
      for (int i = 0; i < 3; ++i) err = fmax(err, fabs(r[i] - expected[k][i]));
      err /= rn;
      int ok = info == 0 && err < 1e-10;
      failures += !ok;
      printf("%-14s %-8s r = (% .10e, % .10e, % .10e) relative error %.2e %s\n", names[k],
             use_enum ? "enum" : "default", r[0], r[1], r[2], err, ok ? "ok" : "FAILED");
    }

  /* random one-contact problems */
  srand(1);
  for (int use_enum = 0; use_enum <= 1; ++use_enum) {
    double worst = 0.;
    int bad = 0;
    for (int k = 0; k < 2000; ++k) {
      double Wr[9], q[3], mu_k, mur_k, r[3];
      rfc2d_random_chain(1, k % 2, Wr, q, &mu_k, &mur_k);
      double Wc[9];
      for (int i = 0; i < 9; ++i) Wc[i] = Wr[i];
      int info = solve(Wc, q, mu_k, mur_k, use_enum, r);
      double v = rfc2d_law_violation(1, Wr, q, &mu_k, &mur_k, r);
      worst = fmax(worst, v);
      if (info || v > 1e-8) bad++;
    }
    /* ENUM is exact; the default local solver iterates (projection with
       local iterations) and may stop short on a few problems: at most 1% */
    int ok = use_enum ? bad == 0 : bad <= 20;
    printf("random one contact, %-8s: 2000 problems, %d failures, worst violation %.2e %s\n",
           use_enum ? "enum" : "default", bad, worst, ok ? "ok" : "FAILED");
    failures += !ok;
  }
  return failures ? 1 : 0;
}
