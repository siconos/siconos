// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "rolling_friction_2d_enum.h"

#include <float.h>   // for DBL_MAX
#include <math.h>    // for fabs, fmax
#include <stdlib.h>  // for realloc
#include <string.h>  // for memcpy

#include "FrictionContact_options.h"        // for SICONOS_FRICTION_3D_CURRENT_CON...
#include "NumericsMatrix.h"                 // for NumericsMatrix
#include "RollingFrictionContactProblem.h"  // for RollingFrictionContactProblem
#include "SolverOptions.h"                  // for SolverOptions, SICONOS_DPARAM_R...
#include "numerics_errors.h"
#include "op3x3.h"                                    // for solve_3x3_gepp
#include "rolling_friction_2d_local_problem_tools.h"  // for rolling_friction_2d_local_pro...
#include "rolling_friction_3d_short_names.h"          // for RFC2D_OC_ENUM
#include "solver_registry.h"

/* a case is accepted when its relative violation is below this value */
#define RFC2D_ENUM_TOLERANCE 1e-12

/* a block without rolling row (a contact whose relation has none: zero
   third row and column of W, W being SPD): its rolling reaction does not
   act, it is taken 0 */
static inline int no_rolling_row(const double* W) {
  const double wmax = fmax(W[0], W[4]);
  return !(W[8] > 1e-14 * wmax);
}

/* the reaction of a closed case with at least one row on the boundary of
   the cone (st or sr not 0, see closed_case): r_n from u_n = 0, and the
   other unknown, if any, from the u = 0 of its sticking row, in closed form;
   0 if the system is (numerically) singular */
static int boundary_case_reaction(const double* W, const double* q, double mu, double mur,
                                  int st, int sr, double r[3]) {
#define W_(i, j) W[(i) + 3 * (j)]
  const double wmax = fmax(fmax(W_(0, 0), W_(1, 1)), W_(2, 2)); /* W is SPD */
  if (no_rolling_row(W)) {
    /* u_r = 0 whatever r_r: r_r = 0, the cases with the rolling row on the
       boundary at mu_r = 0 */
    mur = 0.;
    if (!sr) sr = 1;
  }
  if (st && sr) {
    /* one unknown: r_n */
    double d = W_(0, 0) + st * mu * W_(0, 1) + sr * mur * W_(0, 2);
    if (!(fabs(d) > 1e-14 * wmax)) return 0;
    r[0] = -q[0] / d;
    r[1] = st * mu * r[0];
    r[2] = sr * mur * r[0];
  } else if (st) {
    /* the tangent row slides, the rolling row sticks: r_n, r_r from the
       rows 0 and 2 */
    double a00 = W_(0, 0) + st * mu * W_(0, 1), a01 = W_(0, 2);
    double a10 = W_(2, 0) + st * mu * W_(2, 1), a11 = W_(2, 2);
    double det = a00 * a11 - a01 * a10;
    if (!(fabs(det) > 1e-14 * wmax * wmax)) return 0;
    r[0] = (-q[0] * a11 + a01 * q[2]) / det;
    r[2] = (-a00 * q[2] + a10 * q[0]) / det;
    r[1] = st * mu * r[0];
  } else {
    /* the rolling row rolls, the tangent row sticks: r_n, r_t from the
       rows 0 and 1 */
    double a00 = W_(0, 0) + sr * mur * W_(0, 2), a01 = W_(0, 1);
    double a10 = W_(1, 0) + sr * mur * W_(1, 2), a11 = W_(1, 1);
    double det = a00 * a11 - a01 * a10;
    if (!(fabs(det) > 1e-14 * wmax * wmax)) return 0;
    r[0] = (-q[0] * a11 + a01 * q[1]) / det;
    r[1] = (-a00 * q[1] + a10 * q[0]) / det;
    r[2] = sr * mur * r[0];
  }
  return 1;
#undef W_
}

/* the closed contact case (st, sr): st, sr = 0 for a sticking row, +1 or -1
   for a row on the boundary r_t = st mu r_n, r_r = sr mu_r r_n.  r is the
   reaction of the case, the return value its relative violation of the
   law (DBL_MAX if its system is singular) */
static double closed_case(const double* W, const double* q, double mu, double mur, int st,
                          int sr, double r[3]) {
#define W_(i, j) W[(i) + 3 * (j)]
  if (st == 0 && sr == 0) {
    /* stick-stick: W r = -q */
    r[0] = -q[0];
    r[1] = -q[1];
    r[2] = -q[2];
    if (solve_3x3_gepp(W, r)) return DBL_MAX;
  } else if (!boundary_case_reaction(W, q, mu, mur, st, sr, r)) {
    return DBL_MAX;
  }
  /* solve_3x3_gepp does not test its last pivot: a block without rolling
     row (zero third row and column, mu_r = 0) gives a reaction that is not
     finite, which fmax below would not see */
  if (!(isfinite(r[0]) && isfinite(r[1]) && isfinite(r[2]))) return DBL_MAX;

  double u[3], wr = 0.;
  for (int i = 0; i < 3; ++i) {
    u[i] = q[i] + W_(i, 0) * r[0] + W_(i, 1) * r[1] + W_(i, 2) * r[2];
    wr = fmax(wr, fabs(W_(i, 0) * r[0]) + fabs(W_(i, 1) * r[1]) + fabs(W_(i, 2) * r[2]));
  }
#undef W_
  double rscale = fmax(fmax(fabs(r[0]), fabs(r[1])), fabs(r[2]));
  double uscale = fmax(fmax(fmax(fabs(q[0]), fabs(q[1])), fabs(q[2])), wr);
  if (rscale == 0.) rscale = 1.;
  if (uscale == 0.) uscale = 1.;

  double violation = fmax(0., -r[0] / rscale);
  if (st == 0)
    violation = fmax(violation, (fabs(r[1]) - mu * r[0]) / rscale);
  else /* the velocity opposed to the reaction */
    violation = fmax(violation, st * u[1] / uscale);
  if (sr == 0)
    violation = fmax(violation, (fabs(r[2]) - mur * r[0]) / rscale);
  else
    violation = fmax(violation, sr * u[2] / uscale);
  return violation;
}

/* all the closed cases in a fixed order, each checked against the whole
   law: the reaction of the first accepted case, or of the smallest
   violation; the violation is returned */
static double enumerate_all(const double* W, const double* q, double mu, double mur,
                            double best[3]) {
  static const int cases[9][2] = {{0, 0}, {1, 0},  {-1, 0}, {0, 1},  {0, -1},
                                  {1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
  double best_violation = DBL_MAX;
  for (int c = 0; c < 9; ++c) {
    double r[3];
    double violation = closed_case(W, q, mu, mur, cases[c][0], cases[c][1], r);
    if (violation < best_violation) {
      best_violation = violation;
      best[0] = r[0];
      best[1] = r[1];
      best[2] = r[2];
    }
    if (violation <= RFC2D_ENUM_TOLERANCE) break;
  }
  return best_violation;
}

/* the cases solved in closed form, each with the conditions it leaves
   open: stick-stick (from the inverse of W: r_n >= 0 and the two bounds),
   then the cases predicted by its trial (each row outside its bound slides
   or rolls with the sign of its trial reaction): the sign of the velocity
   of each sliding or rolling row, r_n >= 0, the bound of a sticking row.
   1 and r if a case is accepted, 0 otherwise. */
static int fast_cases(const double* W, const double* Winv, const double* q, double mu,
                      double mur, double r[3]) {
#define W_(i, j) W[(i) + 3 * (j)]
  const double tol = RFC2D_ENUM_TOLERANCE;
  const double wmax = fmax(fmax(W_(0, 0), W_(1, 1)), W_(2, 2)); /* W is SPD */
  const double qmax = fmax(fmax(fabs(q[0]), fabs(q[1])), fabs(q[2]));

  /* stick-stick: u = 0 */
  double rn = -(Winv[0] * q[0] + Winv[3] * q[1] + Winv[6] * q[2]);
  double rt = -(Winv[1] * q[0] + Winv[4] * q[1] + Winv[7] * q[2]);
  double rr = -(Winv[2] * q[0] + Winv[5] * q[1] + Winv[8] * q[2]);
  double rs = fmax(fmax(fabs(rn), fabs(rt)), fabs(rr));
  if (rn >= -tol * rs && fabs(rt) - mu * rn <= tol * rs && fabs(rr) - mur * rn <= tol * rs) {
    r[0] = rn;
    r[1] = rt;
    r[2] = rr;
    return 1;
  }
  const int st = fabs(rt) > mu * rn ? (rt > 0. ? 1 : -1) : 0;
  const int sr = fabs(rr) > mur * rn ? (rr > 0. ? 1 : -1) : 0;
  if (st == 0 && sr == 0) return 0;

  /* the predicted cases: both rows on the boundary first if both are out */
  for (int pass = 0; pass < 3; ++pass) {
    int cst = st, csr = sr;
    if (pass == 1) {
      if (!(st && sr)) break;
      csr = 0;
    } else if (pass == 2) {
      cst = 0;
    }
    double rc[3];
    if (!boundary_case_reaction(W, q, mu, mur, cst, csr, rc)) continue;
    rn = rc[0];
    rt = rc[1];
    rr = rc[2];
    rs = fmax(fmax(fabs(rn), fabs(rt)), fabs(rr));
    double us = fmax(qmax, wmax * rs);
    if (rn < -tol * rs) continue;
    if (cst) {
      double ut = q[1] + W_(1, 0) * rn + W_(1, 1) * rt + W_(1, 2) * rr;
      if (cst * ut > tol * us) continue;
    } else if (fabs(rt) - mu * rn > tol * rs) {
      continue;
    }
    if (csr) {
      double ur = q[2] + W_(2, 0) * rn + W_(2, 1) * rt + W_(2, 2) * rr;
      if (csr * ur > tol * us) continue;
    } else if (fabs(rr) - mur * rn > tol * rs) {
      continue;
    }
    r[0] = rn;
    r[1] = rt;
    r[2] = rr;
    return 1;
  }
  return 0;
#undef W_
}

/* the inverse of the diagonal block of contact, cached by
   rolling_friction_2d_enum_initialize, or NULL */
static const double* cached_inverse(SolverOptions* options) {
  int contact = options->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER];
  if (!options->dWork || options->dWorkSize < 18 * (size_t)(contact + 1)) return NULL;
  const double* inverse = options->dWork + 18 * contact + 9;
  return isnan(inverse[0]) ? NULL : inverse;
}

int rolling_friction_2d_enum_solve(RollingFrictionContactProblem* localproblem,
                                   double* reaction, SolverOptions* options) {
  const double* W = localproblem->M->dense;
  const double* q = localproblem->q;
  const double mu = localproblem->mu[0];
  const double mur = localproblem->mu_r[0];

  /* the opening */
  if (q[0] >= 0.) {
    reaction[0] = reaction[1] = reaction[2] = 0.;
    options->dparam[SICONOS_DPARAM_RESIDU] = 0.;
    return 0;
  }

  const double* Winv = cached_inverse(options);
  if (Winv && fast_cases(W, Winv, q, mu, mur, reaction)) {
    options->dparam[SICONOS_DPARAM_RESIDU] = 0.;
    return 0;
  }

  double best[3] = {0., 0., 0.};
  double best_violation = enumerate_all(W, q, mu, mur, best);
  if (best_violation == DBL_MAX) {
    /* every system singular: keep the reaction */
    options->dparam[SICONOS_DPARAM_RESIDU] = DBL_MAX;
    return 1;
  }
  reaction[0] = best[0];
  reaction[1] = best[1];
  reaction[2] = best[2];
  options->dparam[SICONOS_DPARAM_RESIDU] = best_violation;
  return best_violation <= RFC2D_ENUM_TOLERANCE ? 0 : 1;
}

/* the inverse of the 3 x 3 block B (column major) into inverse, NaN in
   inverse[0] if B is numerically singular */
static void invert_block(const double* B, double* inverse) {
#define B_(i, j) B[(i) + 3 * (j)]
  if (no_rolling_row(B)) {
    /* the inverse of the normal and tangent block, zero on the rolling row
       and column: the stick case gives r_r = 0 */
    double det2 = B_(0, 0) * B_(1, 1) - B_(0, 1) * B_(1, 0);
    double bmax2 =
        fmax(fmax(fabs(B_(0, 0)), fabs(B_(1, 1))), fmax(fabs(B_(0, 1)), fabs(B_(1, 0))));
    if (!(fabs(det2) > 1e-14 * bmax2 * bmax2)) {
      inverse[0] = NAN;
      return;
    }
    inverse[0 + 3 * 0] = B_(1, 1) / det2;
    inverse[1 + 3 * 0] = -B_(1, 0) / det2;
    inverse[0 + 3 * 1] = -B_(0, 1) / det2;
    inverse[1 + 3 * 1] = B_(0, 0) / det2;
    inverse[2 + 3 * 0] = inverse[2 + 3 * 1] = 0.;
    inverse[0 + 3 * 2] = inverse[1 + 3 * 2] = inverse[2 + 3 * 2] = 0.;
    return;
  }
  double c00 = B_(1, 1) * B_(2, 2) - B_(1, 2) * B_(2, 1);
  double c01 = B_(1, 2) * B_(2, 0) - B_(1, 0) * B_(2, 2);
  double c02 = B_(1, 0) * B_(2, 1) - B_(1, 1) * B_(2, 0);
  double det = B_(0, 0) * c00 + B_(0, 1) * c01 + B_(0, 2) * c02;
  double bmax = 0.;
  for (int i = 0; i < 9; ++i) bmax = fmax(bmax, fabs(B[i]));
  if (!(fabs(det) > 1e-14 * bmax * bmax * bmax)) {
    inverse[0] = NAN;
    return;
  }
  /* inverse = adjugate / det, column major: inverse(i, j) = C(j, i) / det */
  inverse[0 + 3 * 0] = c00 / det;
  inverse[1 + 3 * 0] = c01 / det;
  inverse[2 + 3 * 0] = c02 / det;
  inverse[0 + 3 * 1] = (B_(0, 2) * B_(2, 1) - B_(0, 1) * B_(2, 2)) / det;
  inverse[1 + 3 * 1] = (B_(0, 0) * B_(2, 2) - B_(0, 2) * B_(2, 0)) / det;
  inverse[2 + 3 * 1] = (B_(0, 1) * B_(2, 0) - B_(0, 0) * B_(2, 1)) / det;
  inverse[0 + 3 * 2] = (B_(0, 1) * B_(1, 2) - B_(0, 2) * B_(1, 1)) / det;
  inverse[1 + 3 * 2] = (B_(0, 2) * B_(1, 0) - B_(0, 0) * B_(1, 2)) / det;
  inverse[2 + 3 * 2] = (B_(0, 0) * B_(1, 1) - B_(0, 1) * B_(1, 0)) / det;
#undef B_
}

void rolling_friction_2d_enum_initialize(RollingFrictionContactProblem* problem,
                                         RollingFrictionContactProblem* localproblem,
                                         SolverOptions* options) {
  (void)localproblem;
  /* per contact: the diagonal block, then its inverse */
  size_t size = 18 * (size_t)problem->numberOfContacts;
  if (!options->dWork || options->dWorkSize < size) {
    options->dWork = (double*)realloc(options->dWork, size * sizeof(double));
    options->dWorkSize = size;
  }
  for (int contact = 0; contact < problem->numberOfContacts; ++contact) {
    double* cached = options->dWork + 18 * contact;
    double* block = cached;
    /* a sparse block matrix gives a pointer to its block, the other
       storages copy into block */
    NM_extract_diag_block3(problem->M, contact, &block);
    if (block != cached) memcpy(cached, block, 9 * sizeof(double));
    invert_block(cached, cached + 9);
  }
}

void rolling_friction_2d_enum_update(int contact, RollingFrictionContactProblem* problem,
                                     RollingFrictionContactProblem* localproblem,
                                     double* reaction, SolverOptions* options) {
  double* cached = options->dWork + 18 * contact;
  if (problem->M->storageType == NM_SPARSE_BLOCK)
    /* no buffer of its own (rolling_friction_2d_local_problem_allocate) */
    localproblem->M->dense = cached;
  else
    memcpy(localproblem->M->dense, cached, 9 * sizeof(double));
  /* the cached inverse of this contact for rolling_friction_2d_enum_solve */
  options->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] = contact;
  rolling_friction_2d_local_problem_compute_q(problem, localproblem, reaction, contact);
  localproblem->mu[0] = problem->mu[contact];
  localproblem->mu_r[0] = problem->mu_r[contact];
}

void rolling_friction_2d_enum_free(RollingFrictionContactProblem* problem,
                                   RollingFrictionContactProblem* localproblem,
                                   SolverOptions* options) {
  /* dWork is released with the options */
  (void)problem;
  (void)localproblem;
  (void)options;
}

void rolling_friction_2d_enum_set_default(SolverOptions* options) {
  options->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] =
      0;  // this will be set by external solver
}

/* Solver registration wrapper functions */
static int rolling_friction_2d_enum_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int rolling_friction_2d_enum_solve_wrap(void* problem, double* reaction,
                                               double* velocity, SolverOptions* options) {
  (void)velocity;
  return rolling_friction_2d_enum_solve((RollingFrictionContactProblem*)problem, reaction,
                                        options);
}

static void rolling_friction_2d_enum_free_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(RFC2D_OC_ENUM, "RFC2D_OC_ENUM",
                "Enumeration of the cases of one contact for 2D Rolling Friction Contact",
                rolling_friction_2d_enum_init_wrap, rolling_friction_2d_enum_solve_wrap,
                rolling_friction_2d_enum_free_wrap, NULL, /* error function */
                rolling_friction_2d_enum_set_default, 1,  /* default_max_iter */
                1e-12,                                    /* default_tol */
                1 /* is_local_solver */);
