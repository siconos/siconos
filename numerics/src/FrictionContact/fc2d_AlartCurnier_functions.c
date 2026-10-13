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

#include "fc2d_AlartCurnier_functions.h"

#include <math.h>  // for fabs
#include <stdio.h>
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/**
 * \brief Set the normal part of the 2D Alart--Curnier residual and Jacobians.
 *
 * The normal residual is
 * \f[
 *   \Phi_N = r_N - \operatorname{proj}_{\mathbb{R}^+}(r_N - \rho_N u_N)
 * \f]
 * which gives \f$\rho_N u_N\f$ when \f$r_N - \rho_N u_N \ge 0\f$ and
 * \f$r_N\f$ otherwise.
 */
static void alart_curnier_2d_normal(double rN, double uN, double rhoN, double* residual,
                                    double* jacobian_u, double* jacobian_r) {
  const double RVN = rN - rhoN * uN;
  if (RVN >= 0.0) {
    if (residual) residual[0] = rhoN * uN;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 0] = rhoN;
      jacobian_u[1 + 2 * 0] = 0.0;
      jacobian_r[0 + 2 * 0] = 0.0;
      jacobian_r[1 + 2 * 0] = 0.0;
    }
  } else {
    if (residual) residual[0] = rN;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 0] = 0.0;
      jacobian_u[1 + 2 * 0] = 0.0;
      jacobian_r[0 + 2 * 0] = 1.0;
      jacobian_r[1 + 2 * 0] = 0.0;
    }
  }
}

/**
 * \brief Zero the tangential Jacobian blocks (used in the no-contact branch).
 */
static void zero_tangential_jacobian(double* jacobian_u, double* jacobian_r) {
  if (!jacobian_u || !jacobian_r) return;
  jacobian_u[0 + 2 * 1] = 0.0;
  jacobian_u[1 + 2 * 1] = 0.0;
  jacobian_r[0 + 2 * 1] = 0.0;
  jacobian_r[1 + 2 * 1] = 0.0;
}

/**
 * \brief Compute the 2D Alart--Curnier function and its Jacobians.
 *
 * \deprecated This is the historical implementation kept for backward
 * compatibility.  Use computeAlartCurnier2D_AC() or
 * computeAlartCurnier2D_JeanMoreau() for explicit control of the friction
 * radius.
 */
void computeAlartCurnier2D(double* reaction, double* velocity, double mu, double* rho,
                           double* residual, double* jacobian_u, double* jacobian_r) {
  const double rN = reaction[0];
  const double rT = reaction[1];
  const double uN = velocity[0];
  const double uT = velocity[1];
  const double rhoN = rho[0];
  const double rhoT = rho[1];

  const double RVN = rN - rhoN * uN;
  const double RVT = rT - rhoT * uT;

  /* Normal part. */
  if (RVN >= 0.0) {
    if (residual) residual[0] = rhoN * uN;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 0] = rhoN;
      jacobian_u[1 + 2 * 0] = 0.0;
      jacobian_r[0 + 2 * 0] = 0.0;
      jacobian_r[1 + 2 * 0] = 0.0;
    }
  } else {
    if (residual) residual[0] = rN;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 0] = 0.0;
      jacobian_u[1 + 2 * 0] = 0.0;
      jacobian_r[0 + 2 * 0] = 1.0;
      jacobian_r[1 + 2 * 0] = 0.0;
    }
  }

  const double radius = mu * MAX(0, RVN);  // AC variant
  /* Tangential part: projection on [-radius, radius]. */
  if (fabs(RVT) <= radius) {
    if (residual) residual[1] = rhoT * uT;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 1] = rhoT;
      jacobian_r[0 + 2 * 1] = 0.0;
      jacobian_r[1 + 2 * 1] = 0.0;
    }
  } else {
    /* Tangential part: projection on radius */
    if (RVT > radius) {
      if (residual) residual[1] = rT - radius;
      if (jacobian_u && jacobian_r) {
        jacobian_u[0 + 2 * 1] = 0.0;
        jacobian_u[1 + 2 * 1] = 0.0;
        jacobian_r[0 + 2 * 1] = -mu;
        jacobian_r[1 + 2 * 1] = 1.0;
      }
    }
    /* Tangential part: projection on -radius */
    if (RVT < -radius) {
      if (residual) residual[1] = rT + radius;
      if (jacobian_u && jacobian_r) {
        jacobian_u[0 + 2 * 1] = 0.0;
        jacobian_u[1 + 2 * 1] = 0.0;
        jacobian_r[0 + 2 * 1] = mu;
        jacobian_r[1 + 2 * 1] = 1.0;
      }
    }
  }
}

/**
 * \brief Compute the 2D Alart--Curnier residual and Jacobians.
 *
 * Alart--Curnier variant: the friction radius depends on the normal
 * residual mapping,
 * \f[
 *   \mathrm{radius} = \mu \max(0, r_N - \rho_N u_N).
 * \f]
 * This formulation requires the cross derivative
 * \f$\partial \Phi_T / \partial u_N = \pm \mu \rho_N\f$ in the slip
 * branches.
 */
void computeAlartCurnier2D_AC(double* reaction, double* velocity, double mu, double* rho,
                              double* residual, double* jacobian_u, double* jacobian_r) {
  const double rN = reaction[0];
  const double rT = reaction[1];
  const double uN = velocity[0];
  const double uT = velocity[1];
  const double rhoN = rho[0];
  const double rhoT = rho[1];

  const double RVN = rN - rhoN * uN;
  const double RVT = rT - rhoT * uT;

  if (RVN <= 0.0) {
    if (residual) {
      residual[0] = rN;
      residual[1] = rT;
    }
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 0] = 0.0;
      jacobian_u[1 + 2 * 0] = 0.0;
      jacobian_u[0 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 1] = 0.0;

      jacobian_r[0 + 2 * 0] = 1.0;
      jacobian_r[1 + 2 * 0] = 0.0;
      jacobian_r[0 + 2 * 1] = 0.0;
      jacobian_r[1 + 2 * 1] = 1.0;
    }
    return;
  }

  /* Normal part: active contact (RVN > 0). */
  if (residual) residual[0] = rhoN * uN;
  if (jacobian_u && jacobian_r) {
    jacobian_u[0 + 2 * 0] = rhoN;
    jacobian_u[1 + 2 * 0] = 0.0;
    jacobian_r[0 + 2 * 0] = 0.0;
    jacobian_r[1 + 2 * 0] = 0.0;
  }

  const double radius = mu * RVN;

  /* Tangential part. */
  if (fabs(RVT) <= radius) {
    /* Stick. */
    if (residual) residual[1] = rhoT * uT;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 1] = rhoT;
      jacobian_r[0 + 2 * 1] = 0.0;
      jacobian_r[1 + 2 * 1] = 0.0;
    }
  } else if (RVT > radius) {
    /* Slip in the positive tangential direction. */
    if (residual) residual[1] = rT - radius;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 1] = 0.0;
      /* d/d u_N of (r_T - mu*(r_N - rho_N*u_N)) = + mu*rho_N */
      jacobian_u[1 + 2 * 0] = mu * rhoN;
      jacobian_r[0 + 2 * 1] = 0.0;
      jacobian_r[1 + 2 * 1] = 1.0;
      jacobian_r[1 + 2 * 0] = -mu;
    }
  } else {
    /* Slip in the negative tangential direction. */
    if (residual) residual[1] = rT + radius;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 1] = 0.0;
      /* d/d u_N of (r_T + mu*(r_N - rho_N*u_N)) = - mu*rho_N */
      jacobian_u[1 + 2 * 0] = -mu * rhoN;
      jacobian_r[0 + 2 * 1] = 0.0;
      jacobian_r[1 + 2 * 1] = 1.0;
      jacobian_r[1 + 2 * 0] = mu;
    }
  }
}

/**
 * \brief Compute the 2D Alart--Curnier residual and Jacobians.
 *
 * Jean--Moreau variant: the friction radius is the classical Coulomb radius
 * \f[
 *   \mathrm{radius} = \mu r_N.
 * \f]
 * With this choice the tangential residual does not depend on \f$u_N\f$ in
 * the slip branches.
 */
void computeAlartCurnier2D_JeanMoreau(double* reaction, double* velocity, double mu,
                                      double* rho, double* residual, double* jacobian_u,
                                      double* jacobian_r) {
  const double rN = reaction[0];
  const double rT = reaction[1];
  const double uN = velocity[0];
  const double uT = velocity[1];
  const double rhoN = rho[0];
  const double rhoT = rho[1];

  const double RVN = rN - rhoN * uN;
  const double RVT = rT - rhoT * uT;

  if (RVN <= 0.0) {
    if (residual) {
      residual[0] = rN;
      residual[1] = rT;
    }
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 0] = 0.0;
      jacobian_u[1 + 2 * 0] = 0.0;
      jacobian_u[0 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 1] = 0.0;

      jacobian_r[0 + 2 * 0] = 1.0;
      jacobian_r[1 + 2 * 0] = 0.0;
      jacobian_r[0 + 2 * 1] = 0.0;
      jacobian_r[1 + 2 * 1] = 1.0;
    }
    return;
  }

  /* Normal part: active contact (RVN > 0). */
  if (residual) residual[0] = rhoN * uN;
  if (jacobian_u && jacobian_r) {
    jacobian_u[0 + 2 * 0] = rhoN;
    jacobian_u[1 + 2 * 0] = 0.0;
    jacobian_r[0 + 2 * 0] = 0.0;
    jacobian_r[1 + 2 * 0] = 0.0;
  }

  const double radius = mu * rN;

  /* Tangential part. */
  if (fabs(RVT) <= radius) {
    /* Stick. */
    if (residual) residual[1] = rhoT * uT;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 1] = rhoT;
      jacobian_r[0 + 2 * 1] = 0.0;
      jacobian_r[1 + 2 * 1] = 0.0;
    }
  } else if (RVT > radius) {
    /* Slip in the positive tangential direction. */
    if (residual) residual[1] = rT - radius;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 0] = 0.0;
      jacobian_r[0 + 2 * 1] = 0.0;
      jacobian_r[1 + 2 * 1] = 1.0;
      jacobian_r[1 + 2 * 0] = -mu;
    }
  } else {
    /* Slip in the negative tangential direction. */
    if (residual) residual[1] = rT + radius;
    if (jacobian_u && jacobian_r) {
      jacobian_u[0 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 1] = 0.0;
      jacobian_u[1 + 2 * 0] = 0.0;
      jacobian_r[0 + 2 * 1] = 0.0;
      jacobian_r[1 + 2 * 1] = 1.0;
      jacobian_r[1 + 2 * 0] = mu;
    }
  }
}

/**
 * \brief Backward-compatible wrapper that defaults to the Alart--Curnier
 *        formulation with the corrected Jacobian.
 */
void computeAlartCurnier2D_new(double* reaction, double* velocity, double mu, double* rho,
                               double* residual, double* jacobian_u, double* jacobian_r) {
  computeAlartCurnier2D_AC(reaction, velocity, mu, rho, residual, jacobian_u, jacobian_r);
}
