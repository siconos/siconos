// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "plasticity_2d_local_newton.h"

#include <math.h>

#include "Plasticity_options.h"
#include "op3x3.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"

static double norm3(const double a[3]) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

double plasticity_2d_local_spectral_norm(const double M[9]) {
  double v[3] = {1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0)}, Mv[3], MtMv[3];
  double lambda = 0.0;
  for (int k = 0; k < 200; k++) {
    for (int i = 0; i < 3; i++) Mv[i] = M[i] * v[0] + M[i + 3] * v[1] + M[i + 6] * v[2];
    for (int j = 0; j < 3; j++)
      MtMv[j] = M[3 * j] * Mv[0] + M[3 * j + 1] * Mv[1] + M[3 * j + 2] * Mv[2];
    double nl = norm3(MtMv);
    if (nl == 0.0) return 0.0;
    for (int i = 0; i < 3; i++) v[i] = MtMv[i] / nl;
    if (fabs(nl - lambda) <= 1e-15 * nl) {
      lambda = nl;
      break;
    }
    lambda = nl;
  }
  return sqrt(lambda);
}

int plasticity_2d_local_semismooth_newton(Plasticity2DNaturalMap natural_map,
                                          Plasticity2DExactSolver exact, void* data,
                                          const double M[9], double stress[3],
                                          Plasticity2DLocalParameters* p) {
  double n2 = plasticity_2d_local_spectral_norm(M);
  double rho0 = (n2 > 0.0) ? 1.0 / n2 : 1.0;
  double rho = rho0;
  if (p->rho_strategy == PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_CONSTANT && p->rho > 0.0)
    rho = p->rho;
  int ls_max = (p->linesearch_max_iter > 0) ? p->linesearch_max_iter : 20;

  double Phi[3], H[9];
  double res = natural_map(data, rho0, stress, Phi, NULL);
  int it = 0;
  while (res > p->tol * (1.0 + norm3(stress)) && it < p->max_iter) {
    it++;
    double nPhi = natural_map(data, rho, stress, Phi, H);
    double d[3] = {-Phi[0], -Phi[1], -Phi[2]};
    if (solve_3x3_gepp(H, d)) {
      DEBUG_PRINT("local semi-smooth Newton: singular generalized Jacobian\n");
      break;
    }
    /* line search on m = 1/2 |Phi|^2; slope of m along d: -|Phi|^2 */
    double m0 = 0.5 * nPhi * nPhi, slope = -nPhi * nPhi, t = 1.0, st[3], Phit[3];
    if (p->linesearch == PLASTICITY_NSN_LINESEARCH_ARMIJO) {
      for (int ls = 0; ls < ls_max; ls++) {
        for (int i = 0; i < 3; i++) st[i] = stress[i] + t * d[i];
        double nt = natural_map(data, rho, st, Phit, NULL);
        if (0.5 * nt * nt <= m0 + 1e-4 * t * slope) break;
        t *= 0.5;
      }
    } else if (p->linesearch == PLASTICITY_NSN_LINESEARCH_GOLDSTEINPRICE) {
      /* m0 + (1-c) t slope <= m(t) <= m0 + c t slope, c = 0.1 */
      const double c = 0.1;
      double tl = 0.0, tu = INFINITY;
      for (int ls = 0; ls < ls_max; ls++) {
        for (int i = 0; i < 3; i++) st[i] = stress[i] + t * d[i];
        double nt = natural_map(data, rho, st, Phit, NULL);
        double mt = 0.5 * nt * nt;
        if (mt > m0 + c * t * slope)
          tu = t;
        else if (mt < m0 + (1.0 - c) * t * slope && t < 1.0)
          tl = t;
        else
          break;
        t = isinf(tu) ? fmin(2.0 * t, 1.0) : 0.5 * (tl + tu);
      }
    }
    for (int i = 0; i < 3; i++) stress[i] += t * d[i];
    res = natural_map(data, rho0, stress, Phi, NULL);
    DEBUG_PRINTF("local semi-smooth Newton: it %d, t = %e, residual = %e\n", it, t, res);
  }
  p->iter_done = it;
  if (res > p->tol * (1.0 + norm3(stress)) && p->hybrid && exact) {
    DEBUG_PRINT("local semi-smooth Newton: fall back to the exact solver\n");
    double s[3] = {stress[0], stress[1], stress[2]};
    if (exact(data, s) == 0 || natural_map(data, rho0, s, Phi, NULL) < res)
      for (int i = 0; i < 3; i++) stress[i] = s[i];
    res = natural_map(data, rho0, stress, Phi, NULL);
  }
  p->residual = res;
  return (res > p->tol * (1.0 + norm3(stress))) ? 1 : 0;
}
