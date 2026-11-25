// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
         r_t = -mu r_n, r_r = mu_r r_n;
      D. without rolling row (zero rolling row of H, mu_r = 0), sliding and
         sticking: r_r = 0.
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

  /* D. a contact without rolling row (zero third row and column of W,
     mu_r = 0), as assembled for a body without rotational resistance:
     sliding (v = (v0, 0, 0)) and sticking (v = (vs, 0, 0), |r_t| < mu r_n) */
  {
    const double H0[3][3] = {{0., 1., 0.}, {1., 0., R}, {0., 0., 0.}};
    double W0[9];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        double s = 0.;
        for (int l = 0; l < 3; ++l) s += H0[i][l] * Minv[l] * H0[j][l];
        W0[i + 3 * j] = s;
      }
    const double vs = 0.03; /* r_t = -vs / W_tt = -vs / 3 */
    double vD[2][3] = {{v0, -h * g, 0.}, {vs, -h * g, 0.}};
    double expectedD[2][3] = {{rn, -mu * rn, 0.}, {rn, -vs / W0[1 + 3 * 1], 0.}};
    const char* namesD[2] = {"D sliding", "D sticking"};
    for (int use_enum = 0; use_enum <= 1; ++use_enum)
      for (int k = 0; k < 2; ++k) {
        double q[3], r[3], Wc[9];
        for (int i = 0; i < 9; ++i) Wc[i] = W0[i];
        for (int i = 0; i < 3; ++i)
          q[i] = H0[i][0] * vD[k][0] + H0[i][1] * vD[k][1] + H0[i][2] * vD[k][2];
        int info = solve(Wc, q, mu, 0., use_enum, r);
        double err = 0.;
        for (int i = 0; i < 3; ++i) err = fmax(err, fabs(r[i] - expectedD[k][i]));
        err /= rn;
        /* fmax drops a NaN: the reaction is checked finite on its own */
        int ok = info == 0 && err < 1e-10 && isfinite(r[0] + r[1] + r[2]);
        failures += !ok;
        printf("%-14s %-8s r = (% .10e, % .10e, % .10e) relative error %.2e %s\n", namesD[k],
               use_enum ? "enum" : "default", r[0], r[1], r[2], err, ok ? "ok" : "FAILED");
      }
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

  /* random one-contact problems without rolling row (third row and column
     of W and q_r zeroed, mu_r random or 0): r_r = 0 */
  for (int use_enum = 0; use_enum <= 1; ++use_enum) {
    double worst = 0.;
    int bad = 0;
    for (int k = 0; k < 2000; ++k) {
      double Wr[9], q[3], mu_k, mur_k, r[3];
      rfc2d_random_chain(1, k % 2, Wr, q, &mu_k, &mur_k);
      Wr[2] = Wr[5] = Wr[6] = Wr[7] = Wr[8] = 0.;
      q[2] = 0.;
      double Wc[9];
      for (int i = 0; i < 9; ++i) Wc[i] = Wr[i];
      int info = solve(Wc, q, mu_k, mur_k, use_enum, r);
      double v = rfc2d_law_violation(1, Wr, q, &mu_k, &mur_k, r);
      worst = fmax(worst, v);
      if (info || v > 1e-8 || !isfinite(r[0] + r[1] + r[2])) bad++;
    }
    int ok = use_enum ? bad == 0 : bad <= 20;
    printf(
        "random without rolling row, %-8s: 2000 problems, %d failures, worst violation %.2e "
        "%s\n",
        use_enum ? "enum" : "default", bad, worst, ok ? "ok" : "FAILED");
    failures += !ok;
  }
  return failures ? 1 : 0;
}
