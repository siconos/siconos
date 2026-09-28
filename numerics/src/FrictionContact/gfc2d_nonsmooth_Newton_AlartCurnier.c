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

#include "gfc2d_nonsmooth_Newton_AlartCurnier.h"

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for fabs, INFINITY



#ifndef __cplusplus
#include <stdbool.h>  // for true
#endif
#include <stdio.h>   // for NULL, printf, fprintf
#include <stdlib.h>  // for calloc, free, malloc, realloc
#include <string.h>  // for memcpy

#include "CSparseMatrix.h"                       // for CSparseMatrix, CSpar...
#include "FrictionContact_options.h"             // for SICONOS_FRICTION_3D_...
#include "GlobalFrictionContactProblem.h"        // for GlobalFrictionContactProblem
#include "NumericsMatrix.h"                      // for NM_triplet, Numerics...
#include "NumericsSparseMatrix.h"                // for NSM_new, NumericsSpa...
#include "NumericsVector.h"                      // for NV_display
#include "SiconosBlas.h"                         // for cblas_daxpy, cblas_d...
#include "SolverOptions.h"                       // for SolverOptions, SICON...
#include "fc2d_AlartCurnier_functions.h"         // for computeAlartCurnier2D
#include "fc3d_short_names.h"                    // for GFC2D_NSN_AC
#include "VI_cst.h"                              // for SICONOS_VI_DPARAM_RHO
#include "gfc2d_compute_error.h"                 // for gfc2d_compute_error
#include "gfc2d_VI_FixedPointProjection.h"       // for warm-start
#include "numerics_verbose.h"
#include "sanitizer.h"  // for cblas_dcopy_msan


/* Solver registration system */
#include "numerics_errors.h"
#include "solver_registry.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include"siconos_debug.h"


/*!
 * \file gfc2d_nonsmooth_Newton_AlartCurnier.c
 *
 * \brief Sparse global Alart--Curnier non-smooth Newton solver for 2D
 *        frictional contact problems.
 *
 * This file implements a semi-smooth Newton method for the global (sparse)
 * 2D friction contact problem
 *
 * \f[
 *   \begin{cases}
 *     M \, v = H \, r + q, \\
 *     u = H^\top v + b, \\
 *     (u_N, r_N, r_T) \text{ satisfy Signorini--Coulomb conditions},
 *   \end{cases}
 * \f]
 *
 * where \f$M\f$ is the (square) mass/dynamics matrix, \f$H\f$ is the
 * contact Jacobian, \f$q, b\f$ are data vectors, \f$v\f$ is the global
 * velocity of size \f$n\f$, \f$u\f$ is the local contact velocity of size
 * \f$m\f$, \f$r\f$ is the contact reaction of size \f$m\f$ and
 * \f$m = 2 \times \mathrm{numberOfContacts}\f$.
 *
 * The unknowns are stacked as \f$x = (v, u, r) \in \mathbb{R}^{n+2m}\f$ and
 * the semi-smooth system solved by Newton's method is
 *
 * \f[
 *   \Psi(x) =
 *   \begin{pmatrix}
 *     M v - H r - q \\
 *     u - H^\top v - b \\
 *     \Phi_{\mathrm{AC}}(u, r, \rho)
 *   \end{pmatrix} = 0,
 * \f]
 *
 * with \f$\Phi_{\mathrm{AC}}\f$ the 2D Alart--Curnier residual.
 */

/* ===========================================================================
 * Internal constants
 * =========================================================================== */

/* Goldstein--Price line-search tuning constants. */
static const double GFC2D_NSN_AC_LS_INITIAL_ALPHA_MIN = 1e-16;
static const double GFC2D_NSN_AC_LS_INITIAL_ALPHA_MAX = 1e10;
static const double GFC2D_NSN_AC_LS_SLOPE_LOWER = 0.01;
static const double GFC2D_NSN_AC_LS_SLOPE_UPPER = 0.99;

/* ===========================================================================
 * Workspace size helpers
 * =========================================================================== */

/**
 * \brief Compute the size of the persistent double workspace needed by the
 *        2D solver.
 */
static size_t compute_dwork_size(unsigned int num_global_dofs, unsigned int num_contact_dofs,
                                 unsigned int system_size) {
  return (size_t)(2 * num_global_dofs +
                  num_contact_dofs +      /* contact_ac_values */
                  2 * num_contact_dofs +  /* jac_A: 2x2 blocks -> 4 entries per contact */
                  2 * num_contact_dofs +  /* jac_B */
                  num_contact_dofs +      /* rho */
                  system_size +           /* system_residual */
                  system_size +           /* newton_direction */
                  system_size +           /* linesearch_workspace (trial point) */
                  system_size +           /* linesearch_residual */
                  system_size +           /* linear_residual_workspace */
                  system_size);           /* current_solution */
}

/**
 * \brief Compute the size of the persistent integer workspace needed by the
 *        2D solver.
 */
static size_t compute_iwork_size(unsigned int num_contact_dofs) {
  return (size_t)(2 * num_contact_dofs +  /* jac_A rows: 2 per contact */
                  2 * num_contact_dofs +  /* jac_B rows */
                  2 * num_contact_dofs +  /* jac_A cols */
                  2 * num_contact_dofs)   /* jac_B cols */
         * sizeof(CS_INT);
}

/**
 * \brief Evaluate the global residual \f$\Psi(v, u, r)\f$.
 */
static void compute_system_residual(GlobalFrictionContactProblem* problem,
                                    double* global_velocity, double* reaction, double* velocity,
                                    double* rho, double* residual) {
  assert(problem->H->size1 == problem->dimension * problem->numberOfContacts);

  const unsigned int num_contact_dofs = problem->H->size1;
  const unsigned int num_global_dofs = problem->M->size0;
  const unsigned int system_size = num_global_dofs + 2 * num_contact_dofs;

  /* residual = 0 */
  cblas_dscal(system_size, 0., residual, 1);

  /* Dynamics residual: residual[0:n] = q - M*v + H*r */
  cblas_dcopy(num_global_dofs, problem->q, 1, residual, 1);
  NM_gemv(1., problem->H, reaction, 1, residual);
  NM_gemv(-1., problem->M, global_velocity, 1, residual);

  /* Kinematic residual: residual[n:n+m] = b - u + H^T*v */
  cblas_daxpy(num_contact_dofs, -1., velocity, 1, residual + num_global_dofs, 1); // a copy of b would have been better
  cblas_daxpy(num_contact_dofs, 1., problem->b, 1, residual + num_global_dofs, 1);
  NM_tgemv(1., problem->H, global_velocity, 1, residual + num_global_dofs);

  /* Contact residual: residual[n+m:n+2m] = Phi_AC(u, r, rho) */
  for (unsigned int contact = 0; contact < problem->numberOfContacts; ++contact) {
    computeAlartCurnier2D_new(reaction + 2 * contact, velocity + 2 * contact, problem->mu[contact],
                          rho + 2 * contact,
                          residual + num_global_dofs + num_contact_dofs + 2 * contact, NULL,
                          NULL);
  }
}

/**
 * \brief Build the constant sparsity pattern of the global Jacobian \f$J\f$.
 */
static CS_INT build_jacobian_pattern(CSparseMatrix* mass_matrix, CSparseMatrix* contact_matrix,
                                     CSparseMatrix* contact_jac_A,
                                     CSparseMatrix* contact_jac_B, CSparseMatrix* jacobian) {
  /* Only triplet matrices are supported by this assembly routine. */
  assert(mass_matrix->nz >= 0);
  assert(contact_matrix->nz >= 0);
  assert(contact_jac_A->nz >= 0);
  assert(contact_jac_B->nz >= 0);

  /* M and the contact Jacobians must be square. */
  assert(mass_matrix->m == mass_matrix->n);
  assert(contact_jac_A->m == contact_jac_A->n);
  assert(contact_jac_B->m == contact_jac_B->n);

  /* A and B have the same sparsity pattern (2x2 diagonal blocks). */
  assert(contact_jac_A->nz == contact_jac_B->nz);

  /* Top-left block: -M */
  for (int e = 0; e < mass_matrix->nz; ++e) {
    CHECK_RETURN(CSparseMatrix_entry(jacobian, mass_matrix->i[e], mass_matrix->p[e],
                                     -mass_matrix->x[e]));
  }

  /* Top-right block: +H (reaction block, columns start after v and u). */
  assert(mass_matrix->n == contact_matrix->m);
  for (int e = 0; e < contact_matrix->nz; ++e) {
    const CS_INT row = contact_matrix->i[e];
    const CS_INT col = contact_matrix->p[e] + mass_matrix->n + contact_jac_A->n;
    CHECK_RETURN(CSparseMatrix_entry(jacobian, row, col, contact_matrix->x[e]));
  }

  /* Middle-left block: +H^T */
  for (int e = 0; e < contact_matrix->nz; ++e) {
    const CS_INT row = contact_matrix->p[e] + mass_matrix->m;
    const CS_INT col = contact_matrix->i[e];
    CHECK_RETURN(CSparseMatrix_entry(jacobian, row, col, contact_matrix->x[e]));
  }

  /* Middle-middle block: -I */
  for (int e = 0; e < contact_jac_A->m; ++e) {
    CHECK_RETURN(CSparseMatrix_entry(jacobian, e + mass_matrix->m, e + mass_matrix->n, -1.));
  }

  /* Remember where the contact blocks start so that only A and B can be
   * refreshed later. */
  const CS_INT contact_block_start = jacobian->nz;

  /* Bottom-middle block: A */
  for (int e = 0; e < contact_jac_A->nz; ++e) {
    const CS_INT row = contact_jac_A->i[e] + mass_matrix->m + contact_matrix->n;
    const CS_INT col = contact_jac_A->p[e] + mass_matrix->n;
    CHECK_RETURN(CSparseMatrix_entry(jacobian, row, col, contact_jac_A->x[e]));
  }

  /* Bottom-right block: B */
  for (int e = 0; e < contact_jac_B->nz; ++e) {
    const CS_INT row = contact_jac_B->i[e] + mass_matrix->m + contact_matrix->n;
    const CS_INT col = contact_jac_B->p[e] + mass_matrix->n + contact_jac_A->n;
    CHECK_RETURN(CSparseMatrix_entry(jacobian, row, col, contact_jac_B->x[e]));
  }

  return contact_block_start;
}

/**
 * \brief Refresh the contact blocks \f$A\f$ and \f$B\f$ of the Jacobian.
 */
static void update_jacobian_blocks(CSparseMatrix* mass_matrix, CSparseMatrix* contact_matrix,
                                   CSparseMatrix* contact_jac_A,
                                   CSparseMatrix* contact_jac_B, CSparseMatrix* jacobian,
                                   CS_INT contact_block_start) {
  assert(mass_matrix->nz >= 0);
  assert(contact_matrix->nz >= 0);
  assert(contact_jac_A->nz >= 0);
  assert(contact_jac_B->nz >= 0);

  assert(mass_matrix->m == mass_matrix->n);
  assert(contact_jac_A->m == contact_jac_A->n);
  assert(contact_jac_B->m == contact_jac_B->n);

  assert(jacobian->n == mass_matrix->n + contact_jac_A->m + contact_jac_B->m);
  assert(jacobian->m == mass_matrix->m + contact_jac_A->n + contact_jac_B->n);

  assert(jacobian->p);
  assert(jacobian->i);
  assert(jacobian->x);

  /* Ensure enough room for the refreshed contact blocks. */
  if ((contact_block_start + contact_jac_A->nz + contact_jac_B->nz) > jacobian->nzmax) {
    CHECK_RETURN(cs_sprealloc(jacobian, contact_block_start + contact_jac_A->nz +
                                            contact_jac_B->nz));
  }

  /* Overwrite A. */
  jacobian->nz = contact_block_start;

  for (int e = 0; e < contact_jac_A->nz; ++e) {
    if (fabs(contact_jac_A->x[e]) > DBL_EPSILON) {
      jacobian->i[jacobian->nz] = contact_jac_A->i[e] + mass_matrix->m + contact_matrix->n;
      jacobian->p[jacobian->nz] = contact_jac_A->p[e] + mass_matrix->n;
      jacobian->x[jacobian->nz] = contact_jac_A->x[e];
      jacobian->nz++;

      assert(jacobian->nz <= jacobian->nzmax);
    }
  }

  /* Overwrite B. */
  for (int e = 0; e < contact_jac_B->nz; ++e) {
    if (fabs(contact_jac_B->x[e]) > DBL_EPSILON) {
      jacobian->i[jacobian->nz] = contact_jac_B->i[e] + mass_matrix->m + contact_matrix->n;
      jacobian->p[jacobian->nz] = contact_jac_B->p[e] + mass_matrix->n + contact_jac_A->n;
      jacobian->x[jacobian->nz] = contact_jac_B->x[e];
      jacobian->nz++;

      assert(jacobian->nz <= jacobian->nzmax);
    }
  }
}

/**
 * \brief Set up triplet sparse storage for a list of 2x2 diagonal blocks.
 */
static void setup_contact_jacobian_blocks(int num_contacts, double* block_values,
                                          CSparseMatrix* block_matrix) {
  block_matrix->m = 2 * num_contacts;
  block_matrix->n = 2 * num_contacts;
  block_matrix->x = block_values;
  block_matrix->nz = 4 * num_contacts;
  block_matrix->nzmax = block_matrix->nz;

  for (int contact = 0; contact < num_contacts; ++contact) {
    for (int col = 0; col < 2; ++col) {
      for (int row = 0; row < 2; ++row) {
        const int offset = 4 * contact + row + 2 * col;
        block_matrix->i[offset] = 2 * contact + row;
        block_matrix->p[offset] = 2 * contact + col;
      }
    }
  }
}

/**
 * \brief Goldstein--Price line search on the sparse global residual.
 */
static int global_linesearch_goldstein_price(GlobalFrictionContactProblem* problem,
                                             double* current_solution, double* search_direction,
                                             double* rho, double* residual,
                                             CSparseMatrix* jacobian_csc,
                                             double* trial_point, double* trial_residual,
                                             double* step_length,
                                             unsigned int max_ls_iterations) {
  double alpha_min = GFC2D_NSN_AC_LS_INITIAL_ALPHA_MIN;
  double alpha_max = GFC2D_NSN_AC_LS_INITIAL_ALPHA_MAX;

  /* Goldstein--Price bracket constants. */
  const double slope_lower_bound = GFC2D_NSN_AC_LS_SLOPE_LOWER;
  const double slope_upper_bound = GFC2D_NSN_AC_LS_SLOPE_UPPER;

  const unsigned int num_global_dofs = (unsigned)NM_triplet(problem->M)->m;
  const unsigned int num_contact_dofs = problem->H->size1;
  const unsigned int system_size = num_global_dofs + 2 * num_contact_dofs;

  /* Merit function at the current point: q(0) = 0.5 * ||residual||^2. */
  const double merit_initial = 0.5 * cblas_ddot(system_size, residual, 1, residual, 1);

  /* Directional derivative q'(0) = residual^T * J * direction. */
  cblas_dscal(system_size, 0., trial_point, 1);
  cs_gaxpy(jacobian_csc, search_direction, trial_point);
  const double merit_derivative = cblas_ddot(system_size, residual, 1, trial_point, 1);

  DEBUG_PRINTF("global_linesearch_GP: merit_derivative=%e\n", merit_derivative);
  DEBUG_PRINTF("global_linesearch_GP: merit_initial=%e\n", merit_initial);

  for (unsigned int ls_iter = 0; ls_iter < max_ls_iterations; ++ls_iter) {
    /* Build the trial point. */
    cblas_dcopy(system_size, current_solution, 1, trial_point, 1);
    cblas_daxpy(system_size, step_length[0], search_direction, 1, trial_point, 1);

    /* Evaluate the residual at the trial point. */
    compute_system_residual(problem,
                            /* global velocity part of the trial point */
                            trial_point,
                            /* reaction part of the trial point */
                            trial_point + problem->M->size0 + problem->H->size1,
                            /* local velocity part of the trial point */
                            trial_point + problem->M->size0, rho, trial_residual);

    const double merit_trial =
        0.5 * cblas_ddot(system_size, trial_residual, 1, trial_residual, 1);
    assert(merit_trial >= 0);

    const double slope = (merit_trial - merit_initial) / step_length[0];

    const int condition_lower = (slope >= slope_upper_bound * merit_derivative);
    const int condition_upper = (slope <= slope_lower_bound * merit_derivative);

    DEBUG_PRINTF("global_linesearch_GP: lower=%i\t upper=%i\n", condition_lower,
                 condition_upper);

    if (condition_lower && condition_upper) {
      numerics_printf_verbose(2,
                              "---- GFC2D - NSN_AC - global line search success. "
                              "Number of ls iteration = %i  alpha = %.10e, q = %.10e",
                              ls_iter, step_length[0], merit_trial);
      return 0;
    } else if (!condition_lower) {
      /* Step is too small (not enough decrease): increase alpha. */
      alpha_min = step_length[0];
    } else {
      /* Step is too large (decrease too fast): reduce alpha. */
      alpha_max = step_length[0];
    }

    if (step_length[0] < GFC2D_NSN_AC_LS_INITIAL_ALPHA_MAX) {
      step_length[0] = 0.5 * (alpha_min + alpha_max);
    } else {
      step_length[0] = alpha_min;
    }
  }

  numerics_printf_verbose(2,
                          "---- GFC2D - NSN_AC - global line search unsuccessful. "
                          "Max number of ls iteration reached = %i  with alpha = %.10e",
                          max_ls_iterations, step_length[0]);

  return -1;
}

/**
 * \brief Run a VI-FPP warm start for GFC2D_NSN_AC.
 *
 * This can be used as a standalone fallback/preconditioner, or it is called
 * automatically by gfc2d_nonsmooth_Newton_AlartCurnier when the warm-start
 * option is enabled.
 */
int gfc2d_nsn_ac_warm_start(GlobalFrictionContactProblem* problem, double* reaction,
                            double* velocity, double* globalVelocity, SolverOptions* options) {
  if (options->iparam[SICONOS_FRICTION_2D_NSN_WARM_START_STRATEGY] !=
      SICONOS_FRICTION_2D_NSN_WARM_START_YES) {
    return 0;
  }

  const int warm_start_max_iter = options->iparam[SICONOS_FRICTION_2D_NSN_WARM_START_MAX_ITER];
  if (warm_start_max_iter <= 0) return 0;

  const double tolerance = options->dparam[SICONOS_DPARAM_TOL];
  numerics_printf_verbose(
      1, "---- GFC2D - NSN_AC - Warm start with GFC2D_VI_FPP (%d iterations)",
      warm_start_max_iter);

  SolverOptions* vi_options = solver_options_create(GFC2D_VI_FPP);
  if (!vi_options) {
    numerics_printf_verbose(1,
                            "---- GFC2D - NSN_AC - Warm start skipped: "
                            "could not create GFC2D_VI_FPP options");
    return NUMERICS_ERR_MEMORY;
  }

  vi_options->iparam[SICONOS_IPARAM_MAX_ITER] = warm_start_max_iter;
  vi_options->dparam[SICONOS_DPARAM_TOL] = tolerance;

  int vi_info = 0;
  gfc2d_VI_FixedPointProjection(problem, reaction, velocity, globalVelocity, &vi_info,
                                vi_options);

  /* Reuse the rho that VI_FPP computed/selected as the NSN rho. */
  const double vi_rho = fabs(vi_options->dparam[SICONOS_VI_DPARAM_RHO]);
  if (vi_rho > 0.0) {
    options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = vi_rho;
    numerics_printf_verbose(1,
                            "---- GFC2D - NSN_AC - Warm start rho = %g "
                            "(used as NSN rho)",
                            vi_rho);
  }

  solver_options_delete(vi_options);

  numerics_printf_verbose(1,
                          "---- GFC2D - NSN_AC - Warm start finished "
                          "(vi_info = %d)",
                          vi_info);
  return vi_info;
}

/**
 * \brief Solve a sparse global 2D friction contact problem with the
 *        Alart--Curnier non-smooth Newton method.
 */
void gfc2d_nonsmooth_Newton_AlartCurnier(GlobalFrictionContactProblem* problem,
                                         double* reaction, double* velocity,
                                         double* globalVelocity, int* info,
                                         SolverOptions* options) {



 
  assert(problem);
  assert(reaction);
  assert(velocity);
  assert(globalVelocity);
  assert(info);
  assert(options);

  assert(problem->dimension == 2);

  assert(options->iparam);
  assert(options->dparam);

  assert(problem->q);
  assert(problem->mu);
  assert(problem->M);
  assert(problem->H);

  /* M must be square and its row count must match H's row count. */
  assert(problem->M->size0 == problem->M->size1);
  assert(problem->M->size0 == problem->H->size0);

  /* ====================================================================== *
   * Solver options
   * ====================================================================== */
  unsigned int iteration = 0;
  const unsigned int max_iterations = options->iparam[SICONOS_IPARAM_MAX_ITER];
  assert(max_iterations > 0);

  const double tolerance = options->dparam[SICONOS_DPARAM_TOL];
  assert(tolerance > 0);

  /* ====================================================================== *
   * Sparse triplet storage for M and H
   * ====================================================================== */
  DEBUG_EXPR(NM_display(problem->M););
  DEBUG_PRINTF("M sizes = %i x %i \n", problem->M->size0, problem->M->size1);
  DEBUG_PRINTF("H sizes = %i x %i \n", problem->H->size0, problem->H->size1);

  NM_triplet(problem->M);
  NM_triplet(problem->H);

  const unsigned int num_global_dofs = (unsigned)NM_triplet(problem->M)->m;
  const unsigned int num_contact_dofs = problem->H->size1;
  const unsigned int system_size = num_global_dofs + 2 * num_contact_dofs;

  const double norm_q = cblas_dnrm2(num_global_dofs, problem->q, 1);
  const double norm_b = cblas_dnrm2(num_contact_dofs, problem->b, 1);

  numerics_printf_verbose(1, "---- GFC2D - NSN_AC - Problem information");
  numerics_printf_verbose(1, "---- GFC2D - NSN_AC - 1-norm of M = %g norm of q = %g ",
                          NM_norm_1(problem->M), norm_q);
  numerics_printf_verbose(1, "---- GFC2D - NSN_AC - inf-norm of M = %g ",
                          NM_norm_inf(problem->M));

  numerics_printf_verbose(1, "---- GFC2D - NSN_AC - 1-norm of H = %g norm of b = %g ",
                          NM_norm_1(problem->H), norm_b);
  numerics_printf_verbose(1, "---- GFC2D - NSN_AC - inf-norm of H = %g ",
                          NM_norm_inf(problem->H));

  assert((int)num_contact_dofs == problem->numberOfContacts * problem->dimension);
  assert((int)num_global_dofs == problem->H->size0);

  /* ====================================================================== *
   * Optional warm start with GFC2D_VI_FPP
   * ====================================================================== */
  gfc2d_nsn_ac_warm_start(problem, reaction, velocity, globalVelocity, options);

  /* ====================================================================== *
   * Allocate persistent solver workspace on first call
   * ====================================================================== */
  if (options->iparam[SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATED] == 0) {
    assert(options->dWork == NULL);
    assert(options->iWork == NULL);

    options->dWorkSize = compute_dwork_size(num_global_dofs, num_contact_dofs, system_size);
    options->dWork = (double*)calloc(options->dWorkSize, sizeof(double));

    options->iWorkSize = compute_iwork_size(num_contact_dofs);
    options->iWork = (int*)malloc(options->iWorkSize);

    options->iparam[SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATED] = 1;
  }

  const size_t required_dwork_size =
      compute_dwork_size(num_global_dofs, num_contact_dofs, system_size);
  const size_t required_iwork_size = compute_iwork_size(num_contact_dofs);

  if (!options->dWork || options->dWorkSize < required_dwork_size) {
    double* p = (double*)realloc(options->dWork, required_dwork_size * sizeof(double));
    if (!p) {
      numerics_error_log("gfc2d_nonsmooth_Newton_AlartCurnier", "bad alloc");
      *info = NUMERICS_ERR_MEMORY;
      return;
    }
    options->dWork = p;
    options->dWorkSize = required_dwork_size;
  }

  if (!options->iWork || options->iWorkSize < required_iwork_size) {
    int* p = (int*)realloc(options->iWork, required_iwork_size);
    if (!p) {
      numerics_error_log("gfc2d_nonsmooth_Newton_AlartCurnier", "bad alloc");
      *info = NUMERICS_ERR_MEMORY;
      return;
    }
    options->iWork = p;
    options->iWorkSize = required_iwork_size;
  }

  /* ====================================================================== *
   * Slice the workspace into named arrays
   * ====================================================================== */
  double* contact_ac_values = options->dWork + 2 * problem->M->size0;
  double* jac_A = contact_ac_values + num_contact_dofs;
  double* jac_B = jac_A + 2 * num_contact_dofs;
  double* rho = jac_B + 2 * num_contact_dofs;

  double* system_residual = rho + num_contact_dofs;
  double* newton_direction = system_residual + system_size;
  double* linesearch_workspace = newton_direction + system_size;
  double* linesearch_residual = linesearch_workspace + system_size;
  double* linear_residual_workspace = linesearch_residual + system_size;
  double* current_solution = linear_residual_workspace + system_size;

  /* The iWork buffer is reused as row/column index storage for A and B. */
  CS_INT* jac_A_rows = (CS_INT*)options->iWork;
  CS_INT* jac_B_rows = jac_A_rows + 2 * num_contact_dofs;
  CS_INT* jac_A_cols = jac_B_rows + 2 * num_contact_dofs;
  CS_INT* jac_B_cols = jac_A_cols + 2 * num_contact_dofs;

  CSparseMatrix jac_A_sparse;
  CSparseMatrix jac_B_sparse;
  CSparseMatrix* jacobian_triplet;

  jac_A_sparse.p = jac_A_cols;
  jac_B_sparse.p = jac_B_cols;
  jac_A_sparse.i = jac_A_rows;
  jac_B_sparse.i = jac_B_rows;

  setup_contact_jacobian_blocks(problem->numberOfContacts, jac_A, &jac_A_sparse);
  setup_contact_jacobian_blocks(problem->numberOfContacts, jac_B, &jac_B_sparse);

  /* ====================================================================== *
   * Allocate the global Jacobian in triplet format
   * ====================================================================== */
  jacobian_triplet =
      cs_spalloc(NM_triplet(problem->M)->n + jac_A_sparse.m + jac_B_sparse.m,
                 NM_triplet(problem->M)->n + jac_A_sparse.m + jac_B_sparse.m,
                 NM_triplet(problem->M)->nzmax + 2 * NM_triplet(problem->H)->nzmax +
                     2 * jac_A_sparse.n + jac_A_sparse.nzmax + jac_B_sparse.nzmax,
                 1, 1);

  NumericsSparseMatrix* nsm_jacobian = NSM_new();
  nsm_jacobian->triplet = jacobian_triplet;

  assert(jac_A_sparse.n == problem->H->size1);
  assert(jac_A_sparse.nz == problem->numberOfContacts * 4);
  assert(jac_B_sparse.n == problem->H->size1);
  assert(jac_B_sparse.nz == problem->numberOfContacts * 4);

  /* ====================================================================== *
   * Initialize rho and the constant Jacobian pattern
   * ====================================================================== */
  for (unsigned int i = 0; i < num_contact_dofs; ++i) {
    rho[i] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];
  }

  const CS_INT contact_block_start =
      build_jacobian_pattern(NM_triplet(problem->M), NM_triplet(problem->H), &jac_A_sparse,
                             &jac_B_sparse, jacobian_triplet);

  assert(contact_block_start > 0);
  assert(jac_A_sparse.m == jac_A_sparse.n);
  assert(jac_B_sparse.m == jac_B_sparse.n);
  assert(jac_A_sparse.m == problem->H->size1);

  /* Initialize the Newton direction to zero. */
  for (unsigned int i = 0; i < system_size; ++i) newton_direction[i] = 0.;

  NumericsMatrix* nm_jacobian = NM_new();
  assert(jacobian_triplet->m == jacobian_triplet->n);
  NM_fill(nm_jacobian, NM_SPARSE, jacobian_triplet->m, jacobian_triplet->n, nsm_jacobian);

  DEBUG_EXPR(NM_display(nm_jacobian););

  
  info[0] = 1;  /* Default to non-convergence. */

  /* ====================================================================== *
   * Initialize local velocity from the provided global velocity
   * ====================================================================== */
  cblas_dcopy(num_contact_dofs, problem->b, 1, velocity, 1);
  NM_tgemv(1., problem->H, globalVelocity, 1, velocity);

  const double norm_velocity = cblas_dnrm2(num_contact_dofs, velocity, 1);
  if (fabs(norm_velocity) < DBL_EPSILON) {
    for (unsigned int k = 0; k < num_contact_dofs; k++) velocity[k] = 0.0;
  }

  DEBUG_PRINTF("norm of velocity = %e\n", cblas_dnrm2(num_contact_dofs, velocity, 1));


  DEBUG_EXPR(NV_display(globalVelocity, num_global_dofs););
  DEBUG_EXPR(NV_display(velocity, num_contact_dofs););
  DEBUG_EXPR(NV_display(reaction, num_contact_dofs););

  /* Store the initial iterate as the current solution. */
  double* current_global_velocity = current_solution;
  double* current_velocity = &(current_solution[num_global_dofs]);
  double* current_reaction = &(current_solution[num_global_dofs + num_contact_dofs]);

  memcpy(current_global_velocity, globalVelocity, num_global_dofs * sizeof(double));
  memcpy(current_velocity, velocity, num_contact_dofs * sizeof(double));
  memcpy(current_reaction, reaction, num_contact_dofs * sizeof(double));

  numerics_printf_verbose(1, "---- GFC2D - NSN_AC - Start with tolerance = %g", tolerance);
  numerics_printf_verbose(1,
                          "| iter |   residual   | error gfc3d  |    tol      | ls step | ls info | ");
  numerics_printf_verbose(1,
                          "------------------------------------------------------------------");

  /* ====================================================================== *
   * Main Newton loop
   * ====================================================================== */
  double residual_norm = 0.0;
  while (iteration++ < max_iterations) {
    /* Evaluate the residual Psi(x_k). */
    compute_system_residual(problem, current_global_velocity, current_reaction, current_velocity,
                            rho, system_residual);

    residual_norm = cblas_dnrm2(system_size, system_residual, 1);

    /* Evaluate the contact Jacobians A(x_k) and B(x_k). */
    for (unsigned int contact = 0; contact < problem->numberOfContacts; ++contact) {
      computeAlartCurnier2D_new(current_reaction + 2 * contact, current_velocity + 2 * contact,
                            problem->mu[contact], rho + 2 * contact,
                            contact_ac_values + 2 * contact,
                            jac_A + 4 * contact, jac_B + 4 * contact);
    }

    /* Refresh the contact blocks of the global Jacobian. */
    update_jacobian_blocks(NM_triplet(problem->M), NM_triplet(problem->H), &jac_A_sparse,
                           &jac_B_sparse, jacobian_triplet, contact_block_start);

    /* Newton direction: newton_direction = -Psi(x_k). */
    cblas_dcopy(system_size, system_residual, 1, newton_direction, 1);
    cblas_dscal(system_size, -1., newton_direction, 1);

    /* Convert to compressed-column format and solve J(x_k) d = -Psi(x_k). */
    //CSparseMatrix* jacobian_csc = cs_compress(jacobian_triplet);

    NM_clearCSC(nm_jacobian);
    NM_clearCSCTranspose(nm_jacobian);
    NM_clearCSR(nm_jacobian);
    
    
    DEBUG_PRINTF("norm of jacobian = %e\n", NM_norm_1(nm_jacobian));
    DEBUG_EXPR(NM_display(nm_jacobian););
    DEBUG_EXPR(NV_display(newton_direction, num_global_dofs + num_contact_dofs*2 ););


    /* double regul=1e-6; */
    /* NM_add(regul, NM_eye(nm_jacobian->size0), 1.0, nm_jacobian); */
    
    /* NM_preserve/unpreserve ensures that the internal factorization of
     * nm_jacobian is invalidated because the underlying triplet data has
     * changed since the previous iteration. */
    const int linear_solver_info =
        NM_LU_solve(NM_preserve(NM_unpreserve(nm_jacobian)), newton_direction, 1);

    DEBUG_PRINTF("norm of newton_direction = %e\n",
                 cblas_dnrm2(system_size, newton_direction, 1));

    if (linear_solver_info != 0) {
      numerics_printf_verbose(-1, "---- GFC2D - NSN_AC - linear solver failed, info = %d",
                              linear_solver_info);
      *info = NUMERICS_ERR_LINEAR_SOLVER;
      /* // cs_spfree(jacobian_csc); */

      /* double old_residual = residual_norm;       */
      
      /* gfc2d_nsn_ac_warm_start(problem, current_reaction, current_velocity, current_global_velocity, options); */

      /* /\* Evaluate the residual Psi(x_k). *\/ */
      /* compute_system_residual(problem, current_global_velocity, current_reaction, current_velocity, */
      /*                       rho, system_residual); */

      /* residual_norm = cblas_dnrm2(system_size, system_residual, 1); */
      /* printf("residual_norm  =%e", residual_norm); */
      /* printf("old_residual  =%e", old_residual); */
      /* /\* if (residual_norm > 10 * old_residual) *\/ */
      break;

    }

    /* ====================================================================== *
     * Line search
     * ====================================================================== */
    double step_length = 1.0;
    int line_search_status = 0;

    switch (options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH]) {
      case SICONOS_FRICTION_3D_NSN_LINESEARCH_NO:
        /* Take the full Newton step. */
        break;

      case SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE:
        line_search_status = global_linesearch_goldstein_price(
            problem, current_solution, newton_direction, rho, system_residual,
            NM_csc(nm_jacobian),
            linesearch_workspace, linesearch_residual, &step_length,
            options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER]);
        break;

      default:
        *info = numerics_error("gfc2d_nonsmooth_Newton_AlartCurnier",
                               "Unknown line search option.\n");
        break;
    }

    /* Unknown line-search option: abort the Newton loop. */
    if (*info < 0) {
      break;
    }

    /* Warn if the line search could not find an acceptable step. */
    if (line_search_status != 0) {
      numerics_printf_verbose(2,
                              "---- GFC2D - NSN_AC iteration %d, line search did not "
                              "succeed; using last step length = %.10e",
                              iteration, step_length);
    }

    /* Update the iterate: x_{k+1} = x_k + alpha * d_k. */
    cblas_daxpy(system_size, step_length, newton_direction, 1, current_solution, 1);

    DEBUG_PRINTF("norm of the solution  = %e\n", cblas_dnrm2(system_size, current_solution, 1));
    DEBUG_EXPR(NV_display(current_global_velocity, num_global_dofs););
    DEBUG_EXPR(NV_display(current_velocity, num_contact_dofs););
    DEBUG_EXPR(NV_display(current_reaction, num_contact_dofs););

    /* ====================================================================== *
     * Convergence test
     * ====================================================================== */
    const int error_status =
        gfc2d_compute_error(problem, current_reaction, current_velocity, current_global_velocity,
                            tolerance, options, norm_q, norm_b,
                            &(options->dparam[SICONOS_DPARAM_RESIDU]));
    if (error_status != 0) {
      numerics_printf_verbose(-1,
                              "---- GFC2D - NSN_AC iteration %d, error computation "
                              "failed with status %d",
                              iteration, error_status);
      *info = NUMERICS_ERR_UNKNOWN;
      break;
    }

    numerics_printf_verbose(1,
                            "| %4u | %12.4e | %12.4e | %10.2e | %4.2e | %5d |",
                            iteration, residual_norm, options->dparam[SICONOS_DPARAM_RESIDU],
                            tolerance, step_length, line_search_status);

    if (options->dparam[SICONOS_DPARAM_RESIDU] < tolerance) {
      info[0] = 0;
      break;
    }

  }

  /* ====================================================================== *
   * Copy the solution back to the caller's arrays
   * ====================================================================== */
  memcpy(globalVelocity, current_solution, num_global_dofs * sizeof(double));
  memcpy(velocity, &current_solution[num_global_dofs], num_contact_dofs * sizeof(double));
  memcpy(reaction, &current_solution[num_global_dofs + num_contact_dofs],
         num_contact_dofs * sizeof(double));

  if (!info[0])
    numerics_printf_verbose(1,
                            "---- GFC2D - NSN_AC - convergence after %d iterations, "
                            "residual = %g",
                            iteration, options->dparam[SICONOS_DPARAM_RESIDU]);
  else
    numerics_printf_verbose(1,
                            "---- GFC2D - NSN_AC - no convergence after %d iterations, "
                            "residual = %g",
                            iteration, options->dparam[SICONOS_DPARAM_RESIDU]);

  options->iparam[SICONOS_IPARAM_ITER_DONE] = iteration;

  NM_clear(nm_jacobian);
  free(nm_jacobian);
}

/**
 * \brief Set default solver options for GFC2D_NSN_AC.
 */
void gfc2d_nsn_ac_set_default(SolverOptions* options) {
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 1;

  options->iparam[SICONOS_FRICTION_3D_NSN_MPI_COM] = -1;

  options->iparam[SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATED] = 0; /* > 0 memory is allocated */
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] =
      SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE;
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 10;

  options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.; /* default rho */

  /* Warm start with GFC2D_VI_FPP before launching the Newton solver. */
  options->iparam[SICONOS_FRICTION_2D_NSN_WARM_START_STRATEGY] =
    SICONOS_FRICTION_2D_NSN_WARM_START_NO;
  options->iparam[SICONOS_FRICTION_2D_NSN_WARM_START_MAX_ITER] = 50;
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers GFC2D_NSN_AC in the global solver registry, enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

static int gfc2d_nsn_ac_init_wrap(void* problem, SolverOptions* options) {
  /* set_default is now called by the registry before init_wrap */
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int gfc2d_nsn_ac_solve_wrap(void* problem, double* reaction, double* velocity,
                                   double* globalVelocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  gfc2d_nonsmooth_Newton_AlartCurnier((GlobalFrictionContactProblem*)problem, reaction, velocity,
                                      globalVelocity, &info, options);
  return info;
}

static void gfc2d_nsn_ac_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER_3VAR(GFC2D_NSN_AC, "GFC2D_NSN_AC",
                     "Non-smooth Newton Alart-Curnier for 2D Global Friction Contact",
                     gfc2d_nsn_ac_init_wrap, gfc2d_nsn_ac_solve_wrap, gfc2d_nsn_ac_free_wrap,
                     NULL,                     /* error function */
                     gfc2d_nsn_ac_set_default, /* set_default */
                     100,                      /* default_max_iter */
                     1e-6,                     /* default_tol */
                     0 /* is_local_solver */);
