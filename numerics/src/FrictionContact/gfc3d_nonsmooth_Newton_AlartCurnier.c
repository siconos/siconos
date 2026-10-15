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

#include "gfc3d_nonsmooth_Newton_AlartCurnier.h"

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for fabs, INFINITY

#include "siconos_debug.h"  // for DEBUG_PRINTF, DEBUG_EXPR, ...
#ifndef __cplusplus
#include <stdbool.h>  // for true
#endif
#include <stdio.h>   // for NULL, printf, fprintf
#include <stdlib.h>  // for calloc, free, malloc, realloc
#include <string.h>  // for memcpy

#include "AlartCurnierGenerated.h"               // for fc3d_AlartCurnierFun...
#include "CSparseMatrix.h"                       // for CSparseMatrix, CSpar...
#include "FrictionContact_options.h"             // for SICONOS_FRICTION_3D_...
#include "NumericsMatrix.h"                      // for NM_triplet, Numerics...
#include "NumericsSparseMatrix.h"                // for NSM_new, NumericsSpa...
#include "NumericsVector.h"                      // for NV_display
#include "SiconosBlas.h"                         // for cblas_daxpy, cblas_d...
#include "SolverOptions.h"                       // for SolverOptions, SICON...
#include "VI_cst.h"                              // for SICONOS_VI_DPARAM_RHO
#include "cond.h"                                // for cond
#include "fc3d_AlartCurnier_functions.h"         // for computeAlartCurnierJ...
#include "fc3d_nonsmooth_Newton_AlartCurnier.h"  // for fc3d_AlartCurnierFun...
#include "fc3d_short_names.h"                    // for GFC3D_NSN_AC
#include "gfc3d_Solvers.h"                       // for gfc3d_nonsmooth_Newt...
#include "gfc3d_compute_error.h"                 // for gfc3d_compute_error
#include "numerics_verbose.h"
#include "sanitizer.h"  // for cblas_dcopy_msan

/* Solver registration system */
#include "numerics_errors.h"
#include "solver_registry.h"

/*!
 * \file gfc3d_nonsmooth_Newton_AlartCurnier.c
 *
 * \brief Sparse global Alart--Curnier non-smooth Newton solver for 3D
 *        frictional contact problems.
 *
 * This file implements a semi-smooth Newton method for the global (sparse)
 * friction contact problem
 *
 * \f[
 *   \begin{cases}
 *     M \, v = H \, r + q, \\
 *     u = H^\top v + b, \\
 *     (u_N, r_N, r_T) \text{ satisfy Signorini--Coulomb conditions},
 *   \end{cases}
 * \f]
 *
 * where
 * \li \f$M\f$ is the (square) mass/dynamics matrix,
 * \li \f$H\f$ is the contact Jacobian mapping global velocities to local
 *     contact velocities,
 * \li \f$q\f$ and \f$b\f$ are given data vectors,
 * \li \f$v\f$ is the global velocity vector of size \f$n\f$,
 * \li \f$u\f$ is the local contact velocity vector of size \f$m\f$,
 * \li \f$r\f$ is the contact reaction vector of size \f$m\f$,
 * \li \f$m = 3 \times \mathrm{numberOfContacts}\f$.
 *
 * The unknowns are stacked as
 * \f[
 *   x = \begin{pmatrix} v \\ u \\ r \end{pmatrix} \in \mathbb{R}^{n+2m}.
 * \f]
 *
 * The Signorini--Coulomb conditions at each contact are reformulated through
 * the Alart--Curnier function \f$\Phi_{\mathrm{AC}}(u, r, \rho)\f$, where
 * \f$\rho\f$ is a vector of positive scaling parameters (one scalar per
 * contact degree of freedom).  The semi-smooth system solved by Newton's
 * method is therefore
 *
 * \f[
 *   \Psi(x) =
 *   \begin{pmatrix}
 *     M v - H r - q \\
 *     u - H^\top v - b \\
 *     \Phi_{\mathrm{AC}}(u, r, \rho)
 *   \end{pmatrix} = 0.
 * \f]
 *
 * The Jacobian of \f$\Psi\f$ has the block structure
 *
 * \f[
 *   J =
 *   \begin{pmatrix}
 *     M  &  0  & -H \\
 *     -H^\top & I & 0 \\
 *     0  & A(u, r, \rho) & B(u, r, \rho)
 *   \end{pmatrix},
 * \f]
 *
 * where the 3x3 diagonal blocks of \f$A\f$ and \f$B\f$ are the partial
 * Jacobians of \f$\Phi_{\mathrm{AC}}\f$ with respect to \f$u\f$ and
 * \f$r\f$.  At each Newton iteration the linear system
 *
 * \f[ J(x_k) \, d_k = -\Psi(x_k) \f]
 *
 * is solved for the search direction \f$d_k\f$.  The next iterate is obtained
 * by a line search (currently Goldstein--Price) applied to the merit function
 * \f$q(x) = \frac12 \|\Psi(x)\|^2\f$.
 *
 * \note The function gfc3d_sparseGlobalAlartCurnierInit() is intentionally a
 *       no-op.  All workspace allocation is performed on the first call to
 *       gfc3d_nonsmooth_Newton_AlartCurnier().
 */

/* ===========================================================================
 * Internal constants
 * =========================================================================== */

/* Goldstein--Price line-search tuning constants. */
static const double GFC3D_NSN_AC_LS_INITIAL_ALPHA_MIN = 1e-16;
static const double GFC3D_NSN_AC_LS_INITIAL_ALPHA_MAX = 1e10;
static const double GFC3D_NSN_AC_LS_SLOPE_LOWER = 0.01;
static const double GFC3D_NSN_AC_LS_SLOPE_UPPER = 0.99;

/* ===========================================================================
 * Workspace size helpers
 * =========================================================================== */

/**
 * \brief Compute the size of the persistent double workspace needed by the
 *        solver.
 *
 * \param[in] num_global_dofs number of global velocity DOFs (n)
 * \param[in] num_contact_dofs number of contact DOFs (m = 3 * n_c)
 * \param[in] system_size total size of the stacked system (n + 2*m)
 * \return required number of doubles
 */
static size_t compute_dwork_size(unsigned int num_global_dofs, unsigned int num_contact_dofs,
                                 unsigned int system_size) {
  return (size_t)(2 * num_global_dofs + num_contact_dofs + /* contact_ac_values */
                  3 * num_contact_dofs +                   /* jac_A */
                  3 * num_contact_dofs +                   /* jac_B */
                  num_contact_dofs +                       /* rho */
                  system_size +                            /* system_residual */
                  system_size +                            /* newton_direction */
                  system_size + /* linesearch_workspace (trial point) */
                  system_size + /* linesearch_residual */
                  system_size + /* linear_residual_workspace */
                  system_size); /* current_solution */
}

/**
 * \brief Compute the size of the persistent integer workspace needed by the
 *        solver.
 *
 * \param[in] num_contact_dofs number of contact DOFs (m = 3 * n_c)
 * \return required number of bytes
 */
static size_t compute_iwork_size(unsigned int num_contact_dofs) {
  return (size_t)(3 * num_contact_dofs + /* jac_A rows */
                  3 * num_contact_dofs + /* jac_B rows */
                  3 * num_contact_dofs + /* jac_A cols */
                  3 * num_contact_dofs)  /* jac_B cols */
         * sizeof(CS_INT);
}

/**
 * \brief Evaluate the global residual \f$\Psi(v, u, r)\f$.
 *
 * Compute the concatenated residual of the dynamics equation, the kinematic
 * equation and the Alart--Curnier contact equation.
 *
 * \param[in] problem the global friction contact problem
 * \param[in] compute_ac_fun_3x3 pointer to the 3x3 Alart--Curnier function
 * \param[in] global_velocity current global velocity \f$v\f$ (size: n)
 * \param[in] reaction current contact reaction \f$r\f$ (size: m)
 * \param[in] velocity current local contact velocity \f$u\f$ (size: m)
 * \param[in] rho Alart--Curnier scaling vector (size: m)
 * \param[out] residual the computed \f$\Psi\f$ vector (size: n + 2*m)
 */
static void compute_system_residual(GlobalFrictionContactProblem* problem,
                                    AlartCurnierFun3x3Ptr compute_ac_fun_3x3,
                                    double* global_velocity, double* reaction,
                                    double* velocity, double* rho, double* residual) {
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
  cblas_daxpy(num_contact_dofs, -1., velocity, 1, residual + num_global_dofs, 1);
  cblas_daxpy(num_contact_dofs, 1., problem->b, 1, residual + num_global_dofs, 1);
  NM_tgemv(1., problem->H, global_velocity, 1, residual + num_global_dofs);

  /* Contact residual: residual[n+m:n+2m] = Phi_AC(u, r, rho) */
  fc3d_AlartCurnierFunction(num_contact_dofs, compute_ac_fun_3x3, reaction, velocity,
                            problem->mu, rho, residual + num_global_dofs + num_contact_dofs,
                            NULL, NULL);
}

/**
 * \brief Build the constant sparsity pattern of the global Jacobian \f$J\f$.
 *
 * The Jacobian is assembled in triplet format.  The blocks corresponding to
 * \f$M\f$, \f$H\f$, \f$H^\top\f$ and \f$-I\f$ are inserted once.  The
 * contact blocks \f$A\f$ and \f$B\f$ are inserted with initial values; their
 * numerical content is refreshed at each Newton iteration by
 * update_jacobian_blocks().
 *
 * \param[in] mass_matrix triplet sparse matrix \f$M\f$ (size: n x n)
 * \param[in] contact_matrix triplet sparse matrix \f$H\f$ (size: n x m)
 * \param[in] contact_jac_A triplet sparse matrix \f$A\f$ (size: m x m)
 * \param[in] contact_jac_B triplet sparse matrix \f$B\f$ (size: m x m)
 * \param[out] jacobian triplet sparse matrix \f$J\f$ to be filled
 *            (size: (n+2m) x (n+2m))
 * \param[in] scaling scalar scaling applied to the \f$H\f$ and \f$B\f$ blocks
 * \return the triplet index where the \f$A\f$ block starts; this is used by
 *         update_jacobian_blocks() to overwrite only the contact blocks
 */
static CS_INT build_jacobian_pattern(CSparseMatrix* mass_matrix, CSparseMatrix* contact_matrix,
                                     CSparseMatrix* contact_jac_A,
                                     CSparseMatrix* contact_jac_B, CSparseMatrix* jacobian,
                                     double scaling) {
  /* Only triplet matrices are supported by this assembly routine. */
  assert(mass_matrix->nz >= 0);
  assert(contact_matrix->nz >= 0);
  assert(contact_jac_A->nz >= 0);
  assert(contact_jac_B->nz >= 0);

  /* M and the contact Jacobians must be square. */
  assert(mass_matrix->m == mass_matrix->n);
  assert(contact_jac_A->m == contact_jac_A->n);
  assert(contact_jac_B->m == contact_jac_B->n);

  /* A and B have the same sparsity pattern (3x3 diagonal blocks). */
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
    CHECK_RETURN(CSparseMatrix_entry(jacobian, row, col, scaling * contact_matrix->x[e]));
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
    CHECK_RETURN(CSparseMatrix_entry(jacobian, row, col, scaling * contact_jac_B->x[e]));
  }

  return contact_block_start;
}

/**
 * \brief Refresh the contact blocks \f$A\f$ and \f$B\f$ of the Jacobian.
 *
 * Only the entries with magnitude larger than DBL_EPSILON are kept, which
 * avoids storing exact zeros in the sparse pattern.
 *
 * \param[in] mass_matrix triplet sparse matrix \f$M\f$
 * \param[in] contact_matrix triplet sparse matrix \f$H\f$
 * \param[in] contact_jac_A updated contact Jacobian block \f$A\f$
 * \param[in] contact_jac_B updated contact Jacobian block \f$B\f$
 * \param[in,out] jacobian triplet sparse matrix \f$J\f$ whose contact blocks
 *                are overwritten
 * \param[in] contact_block_start triplet index returned by
 *            build_jacobian_pattern()
 */
static void update_jacobian_blocks(CSparseMatrix* mass_matrix, CSparseMatrix* contact_matrix,
                                   CSparseMatrix* contact_jac_A, CSparseMatrix* contact_jac_B,
                                   CSparseMatrix* jacobian, CS_INT contact_block_start) {
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
    CHECK_RETURN(
        cs_sprealloc(jacobian, contact_block_start + contact_jac_A->nz + contact_jac_B->nz));
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
 * \brief Set up triplet sparse storage for a list of 3x3 diagonal blocks.
 *
 * The input array \a block_values is expected to contain \f$9 \times n_c\f$
 * doubles (one dense 3x3 block per contact).  The function fills the row and
 * column index arrays so that the resulting CSparseMatrix describes a block
 * diagonal matrix with these 3x3 blocks.
 *
 * \param[in] num_contacts number of contacts \f$n_c\f$
 * \param[in,out] block_values array of length \f$9 \times n_c\f$ holding the
 *                block entries
 * \param[out] block_matrix triplet sparse matrix to initialize
 */
static void setup_contact_jacobian_blocks(int num_contacts, double* block_values,
                                          CSparseMatrix* block_matrix) {
  block_matrix->m = 3 * num_contacts;
  block_matrix->n = 3 * num_contacts;
  block_matrix->x = block_values;
  block_matrix->nz = 9 * num_contacts;
  block_matrix->nzmax = block_matrix->nz;

  for (int contact = 0; contact < num_contacts; ++contact) {
    for (int col = 0; col < 3; ++col) {
      for (int row = 0; row < 3; ++row) {
        const int offset = 9 * contact + row + 3 * col;
        block_matrix->i[offset] = 3 * contact + row;
        block_matrix->p[offset] = 3 * contact + col;
      }
    }
  }
}

/**
 * \brief Goldstein--Price line search on the sparse global residual.
 *
 * Given the current iterate \a current_solution and a descent direction
 * \a search_direction, find a step length \f$\alpha\f$ such that the merit
 * function
 *
 * \f[ q(x) = \frac12 \|\Psi(x)\|^2 \f]
 *
 * satisfies the Goldstein--Price conditions
 *
 * \f[
 *   m_2 \, q'(0) \le \frac{q(\alpha) - q(0)}{\alpha} \le m_1 \, q'(0),
 * \f]
 *
 * with \f$0 < m_1 < m_2 < 1\f$.  Here \f$q'(0) = \Psi(x)^\top J(x) d\f$.
 *
 * \param[in] problem the global friction contact problem
 * \param[in] compute_ac_fun_3x3 pointer to the 3x3 Alart--Curnier function
 * \param[in] current_solution current iterate \f$x\f$ (size: n + 2*m)
 * \param[in] search_direction Newton direction \f$d\f$ (size: n + 2*m)
 * \param[in] rho Alart--Curnier scaling vector (size: m)
 * \param[in] residual current residual \f$\Psi(x)\f$ (size: n + 2*m)
 * \param[in] jacobian_csc compressed-column Jacobian \f$J(x)\f$
 * \param[out] trial_point workspace used to build trial points (size: n + 2*m)
 * \param[out] trial_residual residual at the trial point (size: n + 2*m);
 *             kept separate from \a trial_point because
 *             compute_system_residual reads the trial point and writes the
 *             residual into distinct buffers
 * \param[in,out] step_length on input, initial trial step length;
 *                on output, accepted step length
 * \param[in] max_ls_iterations maximum number of line-search iterations
 * \return 0 if an acceptable step is found, -1 otherwise
 */

static int global_linesearch_goldstein_price(
    GlobalFrictionContactProblem* problem, AlartCurnierFun3x3Ptr compute_ac_fun_3x3,
    double* current_solution, double* search_direction, double* rho, double* residual,
    CSparseMatrix* jacobian_csc, double* trial_point, double* trial_residual,
    double* step_length, unsigned int max_ls_iterations) {
  double alpha_min = GFC3D_NSN_AC_LS_INITIAL_ALPHA_MIN;
  double alpha_max = GFC3D_NSN_AC_LS_INITIAL_ALPHA_MAX;

  /* Goldstein--Price bracket constants. */
  const double slope_lower_bound = GFC3D_NSN_AC_LS_SLOPE_LOWER;
  const double slope_upper_bound = GFC3D_NSN_AC_LS_SLOPE_UPPER;

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

    /* Evaluate the residual at the trial point.  The trial point and the
     * residual must be stored in separate buffers because
     * compute_system_residual reads the global velocity at the beginning of
     * the trial point while writing the dynamics residual at the beginning of
     * the output buffer. */
    compute_system_residual(problem, compute_ac_fun_3x3,
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
      numerics_printf_verbose(1,
                              "---- GFC3D - NSN_AC - global line search success. "
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

    if (step_length[0] < GFC3D_NSN_AC_LS_INITIAL_ALPHA_MAX) {
      step_length[0] = 0.5 * (alpha_min + alpha_max);
    } else {
      step_length[0] = alpha_min;
    }
  }

  numerics_printf_verbose(1,
                          "---- GFC3D - NSN_AC - global line search unsuccessful. "
                          "Max number of ls iteration reached = %i  with alpha = %.10e",
                          max_ls_iterations, step_length[0]);

  return -1;
}

/**
 * \brief Internal VI-FPP warm start, unconditionally executed.
 */
static int gfc3d_nsn_ac_warm_start_internal(GlobalFrictionContactProblem* problem,
                                            double* reaction, double* velocity,
                                            double* globalVelocity, SolverOptions* options) {
  const int warm_start_max_iter = options->iparam[SICONOS_FRICTION_3D_NSN_WARM_START_MAX_ITER];
  if (warm_start_max_iter <= 0) return 0;

  const double tolerance = options->dparam[SICONOS_DPARAM_TOL];
  numerics_printf_verbose(1,
                          "---- GFC3D - NSN_AC - Warm start with GFC3D_VI_FPP (%d iterations)",
                          warm_start_max_iter);

  SolverOptions* vi_options = solver_options_create(GFC3D_VI_FPP);
  if (!vi_options) {
    numerics_printf_verbose(1,
                            "---- GFC3D - NSN_AC - Warm start skipped: "
                            "could not create GFC3D_VI_FPP options");
    return NUMERICS_ERR_MEMORY;
  }

  vi_options->iparam[SICONOS_IPARAM_MAX_ITER] = warm_start_max_iter;
  vi_options->dparam[SICONOS_DPARAM_TOL] = tolerance;

  int vi_info = 0;
  gfc3d_VI_FixedPointProjection(problem, reaction, velocity, globalVelocity, &vi_info,
                                vi_options);

  /* Reuse the rho that VI_FPP computed/selected as the NSN rho. */
  const double vi_rho = fabs(vi_options->dparam[SICONOS_VI_DPARAM_RHO]);
  if (vi_rho > 0.0) {
    options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = vi_rho * 2;
    numerics_printf_verbose(1,
                            "---- GFC3D - NSN_AC - Warm start rho = %g "
                            "(used as NSN rho)",
                            vi_rho);
  }

  solver_options_delete(vi_options);

  numerics_printf_verbose(1,
                          "---- GFC3D - NSN_AC - Warm start finished "
                          "(vi_info = %d)",
                          vi_info);
  return vi_info;
}

/**
 * \brief Run a VI-FPP warm start for GFC3D_NSN_AC.
 *
 * This can be used as a standalone fallback/preconditioner, or it is called
 * automatically by gfc3d_nonsmooth_Newton_AlartCurnier when the warm-start
 * option is enabled.
 */
int gfc3d_nsn_ac_warm_start(GlobalFrictionContactProblem* problem, double* reaction,
                            double* velocity, double* globalVelocity, SolverOptions* options) {
  if (options->iparam[SICONOS_FRICTION_3D_NSN_WARM_START_STRATEGY] !=
      SICONOS_FRICTION_3D_NSN_WARM_START_YES) {
    return 0;
  }

  return gfc3d_nsn_ac_warm_start_internal(problem, reaction, velocity, globalVelocity,
                                          options);
}

/**
 * \brief Classify a single 3D contact into a discrete status.
 *
 * The tangential norms are expected to be precomputed by the caller (e.g. with
 * NV_cone_tangential_norms_3d) when the classification is performed over the
 * whole contact list.
 *
 * \return 0 = no-contact, 1 = stick, 2 = slip.
 */
static int gfc3d_contact_status_3d(double rN, double rT_norm, double uT_norm, double mu) {
  const double tol = 1e-10;
  const double radius = mu * rN;
  if (rN <= tol) {
    return 0; /* no-contact */
  }
  if (rT_norm < radius - tol && uT_norm <= tol) {
    return 1; /* stick */
  }
  return 2; /* slip */
}

/**
 * \brief Print a contact-status summary for a 3D global friction contact solution.
 */
static void gfc3d_nsn_ac_print_contact_status(GlobalFrictionContactProblem* problem,
                                              double* reaction, double* velocity) {
  const unsigned int nc = (unsigned int)problem->numberOfContacts;
  unsigned int no_contact = 0, stick = 0, slip = 0;

  double* rT_norm = (double*)malloc(nc * sizeof(double));
  double* uT_norm = (double*)malloc(nc * sizeof(double));
  NV_cone_tangential_norms_3d(reaction, nc, rT_norm);
  NV_cone_tangential_norms_3d(velocity, nc, uT_norm);

  for (unsigned int i = 0; i < nc; ++i) {
    const int status =
        gfc3d_contact_status_3d(reaction[3 * i], rT_norm[i], uT_norm[i], problem->mu[i]);
    if (status == 0)
      no_contact++;
    else if (status == 1)
      stick++;
    else
      slip++;
  }

  free(rT_norm);
  free(uT_norm);

  numerics_printf_verbose(1,
                          "---- GFC3D - NSN_AC - contact status: "
                          "no-contact=%u, stick=%u, slip=%u",
                          no_contact, stick, slip);
}

/**
 * \brief No-op initialization hook for the sparse global Alart--Curnier
 *        solver.
 *
 * \param[in,out] options the solver options structure
 */
void gfc3d_sparseGlobalAlartCurnierInit(SolverOptions* options) {
  /* All allocations are performed lazily inside
   * gfc3d_nonsmooth_Newton_AlartCurnier. */
  (void)options;
}

/**
 * \brief Solve a sparse global 3D friction contact problem with the
 *        Alart--Curnier non-smooth Newton method.
 *
 * See the file-level documentation for the mathematical formulation.
 *
 * \param[in] problem the global friction contact problem
 * \param[in,out] reaction on input, initial guess for the contact reaction;
 *                 on output, the computed reaction (size: m)
 * \param[in,out] velocity on input, initial guess for the local contact
 *                 velocity; on output, the computed velocity (size: m)
 * \param[in,out] globalVelocity on input, initial guess for the global
 *                 velocity; on output, the computed global velocity (size: n)
 * \param[out] info solver status: 0 on success, nonzero otherwise
 * \param[in,out] options solver options and parameters
 */

void gfc3d_nonsmooth_Newton_AlartCurnier(GlobalFrictionContactProblem* problem,
                                         double* reaction, double* velocity,
                                         double* globalVelocity, int* info,
                                         SolverOptions* options) {
  assert(problem);
  assert(reaction);
  assert(velocity);
  assert(globalVelocity);
  assert(info);
  assert(options);

  assert(problem->dimension == 3);

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
  unsigned int error_evaluation_frequency =
      options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY];

  if (error_evaluation_frequency == 0) {
    /* A frequency of 0 would prevent any error evaluation; force at least
     * one evaluation per iteration. */
    error_evaluation_frequency = 1;
  }

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

  /* The error function now uses the scaling selected in options->iparam
   * [SICONOS_FRICTION_3D_SCALING_NORM]. */

  numerics_printf_verbose(1, "---- GFC3D - NSN_AC - Problem information");
  numerics_printf_verbose(1, "---- GFC3D - NSN_AC - 1-norm of M = %g norm of q = %g ",
                          NM_norm_1(problem->M), norm_q);
  numerics_printf_verbose(1, "---- GFC3D - NSN_AC - inf-norm of M = %g ",
                          NM_norm_inf(problem->M));
  numerics_printf_verbose(1, "---- GFC3D - NSN_AC - 1-norm of H = %g norm of b = %g ",
                          NM_norm_1(problem->H), norm_b);
  numerics_printf_verbose(1, "---- GFC3D - NSN_AC - inf-norm of H = %g ",
                          NM_norm_inf(problem->H));

  DEBUG_EXPR(
      NumericsMatrix* Mdense = NM_create(NM_DENSE, problem->M->size0, problem->M->size1);
      NM_to_dense(problem->M, Mdense); printf(
          "conditioning = %g\n", cond(Mdense->matrix0, problem->M->size0, problem->M->size1));
      Mdense = NM_free(Mdense););

  assert((int)num_contact_dofs == problem->numberOfContacts * problem->dimension);
  assert((int)num_global_dofs == problem->H->size0);

  /* ====================================================================== *
   * Select the 3x3 Alart--Curnier function
   * ====================================================================== */
  AlartCurnierFun3x3Ptr compute_ac_fun_3x3 = NULL;

  switch (options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION]) {
    case SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD:
      compute_ac_fun_3x3 = &computeAlartCurnierSTD;
      break;

    case SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD:
      compute_ac_fun_3x3 = &computeAlartCurnierJeanMoreau;
      break;

    case SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED:
      compute_ac_fun_3x3 = &fc3d_AlartCurnierFunctionGenerated;
      break;

    case SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED:
      compute_ac_fun_3x3 = &fc3d_AlartCurnierJeanMoreauFunctionGenerated;
      break;

    default:
      *info = numerics_error("gfc3d_nonsmooth_Newton_AlartCurnier",
                             "Unknown Alart--Curnier formulation.\n");
      return;
  }

  /* ====================================================================== *
   * Optional warm start with GFC3D_VI_FPP
   * ====================================================================== */
  gfc3d_nsn_ac_warm_start(problem, reaction, velocity, globalVelocity, options);

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
      numerics_error_log("gfc3d_nonsmooth_Newton_AlartCurnier", "bad alloc");
      *info = NUMERICS_ERR_MEMORY;
      return;
    }
    options->dWork = p;
    options->dWorkSize = required_dwork_size;
  }

  if (!options->iWork || options->iWorkSize < required_iwork_size) {
    int* p = (int*)realloc(options->iWork, required_iwork_size);
    if (!p) {
      numerics_error_log("gfc3d_nonsmooth_Newton_AlartCurnier", "bad alloc");
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
  double* jac_B = jac_A + 3 * num_contact_dofs;
  double* rho = jac_B + 3 * num_contact_dofs;

  double* system_residual = rho + num_contact_dofs;
  double* newton_direction = system_residual + system_size;
  double* linesearch_workspace = newton_direction + system_size;
  double* linesearch_residual = linesearch_workspace + system_size;
  double* linear_residual_workspace = linesearch_residual + system_size;
  double* current_solution = linear_residual_workspace + system_size;

  /* The iWork buffer is reused as row/column index storage for A and B. */
  CS_INT* jac_A_rows = (CS_INT*)options->iWork;
  CS_INT* jac_B_rows = jac_A_rows + 3 * num_contact_dofs;
  CS_INT* jac_A_cols = jac_B_rows + 3 * num_contact_dofs;
  CS_INT* jac_B_cols = jac_A_cols + 3 * num_contact_dofs;

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

  DEBUG_PRINTF("NM_triplet(problem->M)->n= %lli\t,NM_triplet(problem->M)->nzmax = %lli\n",
               (long long)NM_triplet(problem->M)->n, (long long)NM_triplet(problem->M)->nzmax);
  DEBUG_PRINTF("NM_triplet(problem->M)->n + jac_A_sparse.m + jac_B_sparse.m = %lli\n",
               (long long)(NM_triplet(problem->M)->n + jac_A_sparse.m + jac_B_sparse.m));
  DEBUG_PRINTF("Jacobian nzmax= %lli\n",
               (long long)(NM_triplet(problem->M)->nzmax + 2 * NM_triplet(problem->H)->nzmax +
                           2 * jac_A_sparse.n + jac_A_sparse.nzmax + jac_B_sparse.nzmax));

  assert(jac_A_sparse.n == problem->H->size1);
  assert(jac_A_sparse.nz == problem->numberOfContacts * 9);
  assert(jac_B_sparse.n == problem->H->size1);
  assert(jac_B_sparse.nz == problem->numberOfContacts * 9);

  /* ====================================================================== *
   * Initialize rho and the constant Jacobian pattern
   * ====================================================================== */
  const CS_INT contact_block_start =
      build_jacobian_pattern(NM_triplet(problem->M), NM_triplet(problem->H), &jac_A_sparse,
                             &jac_B_sparse, jacobian_triplet, /* scaling = */ 1.0);

  assert(contact_block_start > 0);
  assert(jac_A_sparse.m == jac_A_sparse.n);
  assert(jac_B_sparse.m == jac_B_sparse.n);
  assert(jac_A_sparse.m == problem->H->size1);

  /* jacobian_triplet is square; rows = jacobian_triplet->m, cols =
   * jacobian_triplet->n. */
  assert(jacobian_triplet->m == jacobian_triplet->n);
  NumericsMatrix* nm_jacobian =
      NM_create_with_data(NM_SPARSE, jacobian_triplet->m, jacobian_triplet->n, nsm_jacobian);

  /* Optional Jacobian regularization matrix (regul * I). */
  NumericsMatrix* reg_eye = NULL;
  const int regularization = options->iparam[SICONOS_FRICTION_3D_NSN_REGULARIZATION];
  const double regul = options->dparam[SICONOS_FRICTION_3D_NSN_REGUL];
  if (regularization == SICONOS_FRICTION_3D_NSN_REGULARIZATION_YES) {
    reg_eye = NM_eye(nm_jacobian->size0);
  }

  /* Store the initial iterate as the current solution. */
  double* current_global_velocity = current_solution;
  double* current_velocity = &(current_solution[num_global_dofs]);
  double* current_reaction = &(current_solution[num_global_dofs + num_contact_dofs]);

  /* Track contact statuses across Newton iterations to report how many
   * contacts have reached a stable status. */
  int* previous_contact_status = (int*)malloc(problem->numberOfContacts * sizeof(int));
  double* rT_norm = (double*)malloc(problem->numberOfContacts * sizeof(double));
  double* uT_norm = (double*)malloc(problem->numberOfContacts * sizeof(double));

  /* Automatic warm start: only run VI-FPP if the first Newton linear solve
   * fails.  Warm start is applied at most once. */
  const int warm_start_auto = (options->iparam[SICONOS_FRICTION_3D_NSN_WARM_START_STRATEGY] ==
                               SICONOS_FRICTION_3D_NSN_WARM_START_AUTO);
  int warm_start_done = 0;

restart:
  info[0] = 1; /* Default to non-convergence. */
  iteration = 0;

  /* Initialize rho from the current NSN rho parameter. */
  for (unsigned int i = 0; i < num_contact_dofs; ++i) {
    rho[i] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];
  }

  DEBUG_EXPR(for (unsigned int i = 0; i < num_contact_dofs; ++i)
                 printf("rho[%i]=%e\t", i, rho[i]);
             printf("\n"););

  /* Initialize the Newton direction to zero. */
  for (unsigned int i = 0; i < system_size; ++i) newton_direction[i] = 0.;

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
  double linear_solver_residual = 0.0;

  DEBUG_EXPR(NV_display(globalVelocity, num_global_dofs););
  DEBUG_EXPR(NV_display(velocity, num_contact_dofs););
  DEBUG_EXPR(NV_display(reaction, num_contact_dofs););

  memcpy(current_global_velocity, globalVelocity, num_global_dofs * sizeof(double));
  memcpy(current_velocity, velocity, num_contact_dofs * sizeof(double));
  memcpy(current_reaction, reaction, num_contact_dofs * sizeof(double));

  if (warm_start_done)
    numerics_printf_verbose(1, "---- GFC3D - NSN_AC - Restarting after automatic warm start");
  else
    numerics_printf_verbose(1, "---- GFC3D - NSN_AC - Start with tolerance = %g", tolerance);

  numerics_printf_verbose(
      1,
      "| iter |   residual   | error gfc3d  |    tol      | ls step | ls info | "
      "no-c | stick | slip | stab%% |");
  numerics_printf_verbose(
      1,
      "--------------------------------------------------------------------------------"
      "---------------------------");

  for (unsigned int i = 0; i < problem->numberOfContacts; ++i) {
    previous_contact_status[i] = -1;
  }

  /* ====================================================================== *
   * Main Newton loop
   * ====================================================================== */
  while (iteration++ < max_iterations) {
    /* Evaluate the residual Psi(x_k). */
    compute_system_residual(problem, compute_ac_fun_3x3, current_global_velocity,
                            current_reaction, current_velocity, rho, system_residual);

    const double residual_norm = cblas_dnrm2(system_size, system_residual, 1);

    /* Evaluate the contact Jacobians A(x_k) and B(x_k).  The Alart--Curnier
     * function values are written to contact_ac_values; some underlying 3x3
     * functions do not tolerate a NULL result pointer even though these values
     * are not used by this solver. */
    fc3d_AlartCurnierFunction(num_contact_dofs, compute_ac_fun_3x3, current_reaction,
                              current_velocity, problem->mu, rho, contact_ac_values, jac_A,
                              jac_B);

    /* Invalidate cached CSC/CSR/transpose views because the triplet values
     * will change. */
    NM_clearCSC(nm_jacobian);
    NM_clearCSCTranspose(nm_jacobian);
    NM_clearCSR(nm_jacobian);

    /* Refresh the contact blocks of the global Jacobian. */
    update_jacobian_blocks(NM_triplet(problem->M), NM_triplet(problem->H), &jac_A_sparse,
                           &jac_B_sparse, jacobian_triplet, contact_block_start);

    /* Newton direction: newton_direction = -Psi(x_k). */
    cblas_dcopy(system_size, system_residual, 1, newton_direction, 1);
    cblas_dscal(system_size, -1., newton_direction, 1);

    /* Convert to compressed-column format.  The unregularized CSC view is
     * used by the Goldstein--Price line search. */
    CSparseMatrix* jacobian_csc = NM_csc(nm_jacobian);

    DEBUG_PRINTF("norm of jacobian = %e\n", NM_norm_1(nm_jacobian));

    /* If regularization is enabled, solve the regularized system
     * (J + regul * I) d = -Psi(x_k). */
    NumericsMatrix* J_to_solve = nm_jacobian;
    if (regularization == SICONOS_FRICTION_3D_NSN_REGULARIZATION_YES) {
      J_to_solve = NM_add(regul, reg_eye, 1.0, nm_jacobian);
    }

    /* NM_preserve/unpreserve ensures that the internal factorization of
     * the Jacobian used by the linear solver is invalidated because the
     * underlying triplet data has changed since the previous iteration. */
    const int linear_solver_info =
        NM_LU_solve(NM_preserve(NM_unpreserve(J_to_solve)), newton_direction, 1);

    if (J_to_solve != nm_jacobian) {
      NM_free(J_to_solve);
    }

    DEBUG_PRINTF("norm of newton_direction = %e\n",
                 cblas_dnrm2(system_size, newton_direction, 1));

    if (linear_solver_info != 0) {
      numerics_printf_verbose(-1, "---- GFC3D - NSN_AC - linear solver failed, info = %d",
                              linear_solver_info);

      /* Automatic warm start: if this is the first linear-solve failure and
       * the AUTO strategy is active, run VI-FPP once and restart Newton. */
      if (warm_start_auto && !warm_start_done) {
        numerics_printf_verbose(1,
                                "---- GFC3D - NSN_AC - first linear solve failed, triggering "
                                "automatic warm start");
        warm_start_done = 1;
        gfc3d_nsn_ac_warm_start_internal(problem, reaction, velocity, globalVelocity, options);
        memcpy(current_global_velocity, globalVelocity, num_global_dofs * sizeof(double));
        memcpy(current_velocity, velocity, num_contact_dofs * sizeof(double));
        memcpy(current_reaction, reaction, num_contact_dofs * sizeof(double));
        goto restart;
      }

      *info = NUMERICS_ERR_LINEAR_SOLVER;
      break;
    }

    /* Check the quality of the linear solve. */
    if (verbose > 0) {
      cblas_dcopy_msan(system_size, system_residual, 1, linear_residual_workspace, 1);
      NM_gemv(1., nm_jacobian, newton_direction, 1., linear_residual_workspace);
      linear_solver_residual = cblas_dnrm2(system_size, linear_residual_workspace, 1);
      numerics_printf_verbose(1,
                              "---- GFC3D - NSN_AC iteration %d, linear_solver_residual = %g ",
                              iteration, linear_solver_residual);
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
            problem, compute_ac_fun_3x3, current_solution, newton_direction, rho,
            system_residual, jacobian_csc, linesearch_workspace, linesearch_residual,
            &step_length, options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER]);
        break;

      default:
        *info = numerics_error("gfc3d_nonsmooth_Newton_AlartCurnier",
                               "Unknown line search option.\n");
        break;
    }

    /* Unknown formulation or line-search option: abort the Newton loop. */
    if (*info < 0) {
      break;
    }

    /* Warn if the line search could not find an acceptable step, but still
     * use the last trial step length.  This preserves the historical behavior
     * while making the failure visible. */
    if (line_search_status != 0) {
      numerics_printf_verbose(1,
                              "---- GFC3D - NSN_AC iteration %d, line search did not "
                              "succeed; using last step length = %.10e",
                              iteration, step_length);
    }

    /* Update the iterate: x_{k+1} = x_k + alpha * d_k. */
    cblas_daxpy(system_size, step_length, newton_direction, 1, current_solution, 1);

    DEBUG_PRINTF("norm of the solution  = %e\n",
                 cblas_dnrm2(system_size, current_solution, 1));
    DEBUG_EXPR(NV_display(current_global_velocity, num_global_dofs););
    DEBUG_EXPR(NV_display(current_velocity, num_contact_dofs););
    DEBUG_EXPR(NV_display(current_reaction, num_contact_dofs););

    options->dparam[SICONOS_DPARAM_RESIDU] = INFINITY;

    /* ====================================================================== *
     * Convergence test
     * ====================================================================== */
    if (!(iteration % error_evaluation_frequency)) {
      const int error_status = gfc3d_compute_error(
          problem, current_reaction, current_velocity, current_global_velocity, tolerance,
          options, norm_q, norm_b, &(options->dparam[SICONOS_DPARAM_RESIDU]));
      if (error_status < 0) {
        numerics_printf_verbose(1,
                                "---- GFC3D - NSN_AC iteration %d, error computation "
                                "failed with status %d",
                                iteration, error_status);
        *info = NUMERICS_ERR_UNKNOWN;
        break;
      }
      /* A positive status simply means the residual is still above the
       * tolerance.  Keep iterating instead of aborting the Newton loop. */
    }

    /* Classify contacts and count stabilized ones (status unchanged from
     * the previous Newton iteration). */
    unsigned int no_contact = 0, stick = 0, slip = 0;
    unsigned int stabilized = 0;
    NV_cone_tangential_norms_3d(current_reaction, problem->numberOfContacts, rT_norm);
    NV_cone_tangential_norms_3d(current_velocity, problem->numberOfContacts, uT_norm);
    for (unsigned int contact = 0; contact < problem->numberOfContacts; ++contact) {
      const double rN = current_reaction[3 * contact];
      const int status = gfc3d_contact_status_3d(rN, rT_norm[contact], uT_norm[contact],
                                                 problem->mu[contact]);
      if (status == 0)
        no_contact++;
      else if (status == 1)
        stick++;
      else
        slip++;
      if (previous_contact_status[contact] != -1 &&
          previous_contact_status[contact] == status) {
        stabilized++;
      }
      previous_contact_status[contact] = status;
    }
    const double stabilized_pct = (problem->numberOfContacts > 0)
                                      ? (100.0 * stabilized) / problem->numberOfContacts
                                      : 0.0;

    numerics_printf_verbose(1,
                            "| %4u | %12.4e | %12.4e | %10.2e | %4.2e | %5d | "
                            "%4u | %5u | %4u | %5.1f |",
                            iteration, residual_norm, options->dparam[SICONOS_DPARAM_RESIDU],
                            tolerance, step_length, line_search_status, no_contact, stick,
                            slip, stabilized_pct);

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
                            "---- GFC3D - NSN_AC - convergence after %d iterations, "
                            "residual = %g",
                            iteration, options->dparam[SICONOS_DPARAM_RESIDU]);
  else
    numerics_printf_verbose(1,
                            "---- GFC3D - NSN_AC - no convergence after %d iterations, "
                            "residual = %g",
                            iteration, options->dparam[SICONOS_DPARAM_RESIDU]);

  gfc3d_nsn_ac_print_contact_status(problem, reaction, velocity);

  options->iparam[SICONOS_IPARAM_ITER_DONE] = iteration;

#ifdef DUMP_PROBLEM
  if (info[0]) {
    static int file_counter = 0;
    char filename[64];
    printf("GLOBALAC: dumping problem\n");
    snprintf(filename, sizeof(filename), "GLOBALAC_failure%d.dat", file_counter++);
    FILE* file = fopen(filename, "w");
    frictionContact_printInFile(problem, file);
    fclose(file);
  }
#endif

  NM_clear(nm_jacobian);
  free(nm_jacobian);
  if (reg_eye) {
    NM_free(reg_eye);
  }
  free(previous_contact_status);
  free(rT_norm);
  free(uT_norm);
}

/**
 * \brief Set default solver options for GFC3D_NSN_AC.
 *
 * \param[in,out] options the solver options structure to initialize
 */

void gfc3d_nsn_ac_set_default(SolverOptions* options) {
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 1;

  options->iparam[SICONOS_FRICTION_3D_NSN_MPI_COM] = -1;

  options->iparam[SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATED] = 0; /* > 0 memory is allocated */
  options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] =
      SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD;
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] =
      SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE;
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 1000;

  options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.; /* default rho */

  /* Warm start with GFC3D_VI_FPP before launching the Newton solver.
   * By default the warm start is automatic: it is triggered only if the
   * first Newton Jacobian factorization fails. */
  options->iparam[SICONOS_FRICTION_3D_NSN_WARM_START_STRATEGY] =
      SICONOS_FRICTION_3D_NSN_WARM_START_AUTO;
  options->iparam[SICONOS_FRICTION_3D_NSN_WARM_START_MAX_ITER] = 50;

  /* Jacobian regularization: disabled by default. */
  options->iparam[SICONOS_FRICTION_3D_NSN_REGULARIZATION] =
      SICONOS_FRICTION_3D_NSN_REGULARIZATION_NO;
  options->dparam[SICONOS_FRICTION_3D_NSN_REGUL] = 1e-6;

  /* Weight for SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED (default 1.0). */
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers GFC3D_NSN_AC in the global solver registry, enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

static int gfc3d_nsn_ac_init_wrap(void* problem, SolverOptions* options) {
  /* set_default is now called by the registry before init_wrap */
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int gfc3d_nsn_ac_solve_wrap(void* problem, double* reaction, double* velocity,
                                   double* globalVelocity, SolverOptions* options) {
  int info = NUMERICS_OK;

  gfc3d_nonsmooth_Newton_AlartCurnier((GlobalFrictionContactProblem*)problem, reaction,
                                      velocity, globalVelocity, &info, options);

  return info;
}

static void gfc3d_nsn_ac_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER_3VAR(GFC3D_NSN_AC, "GFC3D_NSN_AC",
                     "Non-smooth Newton Alart-Curnier for 3D Global Friction Contact",
                     gfc3d_nsn_ac_init_wrap, gfc3d_nsn_ac_solve_wrap, gfc3d_nsn_ac_free_wrap,
                     NULL,                     /* error function */
                     gfc3d_nsn_ac_set_default, /* set_default */
                     100,                      /* default_max_iter */
                     1e-6,                     /* default_tol */
                     0 /* is_local_solver */);
