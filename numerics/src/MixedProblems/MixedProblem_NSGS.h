/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef MIXEDPROBLEM_NSGS_H
#define MIXEDPROBLEM_NSGS_H

/*!\file MixedProblem_NSGS.h
  Generic Non-Smooth Gauss-Seidel (NSGS) solver for MixedProblems.
  Handles heterogeneous combinations of friction, rolling friction,
  and plasticity models.
*/

#include "MixedProblem.h"
#include "MixedProblems_options.h"
#include "SolverOptions.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/**
 * Solve a MixedProblem w = M z + q, z_e in K_e, w_e (+ De Saxce/dilatancy term) in K_e^*
 * for each element e, using Non-Smooth Gauss-Seidel (NSGS).
 *
 * One iteration is a Gauss-Seidel sweep over the blocks, in their order:
 *  - friction 2D/3D, Drucker-Prager (Mohr-Coulomb) and Von Mises blocks: sweep over the
 *    elements (contacts, cones) of the block, directly on the global M (any storage). For
 *    element e (rows r..r+d-1) the one-element problem has the diagonal block M_ee
 *    (NM_extract_diag_block{2,3,5} or NM_extract_diag_block) and
 *    q_e + sum_{f != e} M_ef z_f (NM_row_prod_no_diag*), with the z_f already updated.
 *    Local solvers: fc2d_nsgs_local_solve (exact), fc3d_projectionOnCone_solve,
 *    plasticity_2d_projectionOnCone_solve (one projection step each),
 *    plasticity_2d_projectionOnVonMises_solve;
 *  - rolling-friction blocks: the block problem (dense diagonal block, q with the coupling to
 *    the other blocks) is solved by the rolling-friction NSGS solver.
 * With a sparse-block M, the SBM must have one diagonal block per element (block row
 * number = global element number); rolling blocks must then have one contact.
 *
 * Stopping test (iparam[SICONOS_MIXED_IPARAM_ERROR_EVALUATION]):
 *  - LIGHT: relative change of z, |z_k - z_{k-1}| / |z_k| < tol;
 *  - LIGHT_WITH_FULL_FINAL (default): when the light error is below tol, the full error
 *    (mixedProblem_compute_error) is computed and must be below tol;
 *  - FULL: full error every iparam[SICONOS_MIXED_IPARAM_ERROR_EVALUATION_FREQUENCY]
 *    iterations.
 * The initial guess is z = 0. On exit w = M z + q, dparam[SICONOS_MIXED_DPARAM_RESIDU] is the
 * last error and iparam[SICONOS_MIXED_IPARAM_ITER_DONE] the number of iterations.
 * See docs/devel/MIXED_PROBLEMS_NSGS.md.
 *
 * \param[in] problem the MixedProblem to solve
 * \param[out] z the solution vector (size: problem->total_dimension)
 * \param[out] w w = M z + q (size: problem->total_dimension)
 * \param[out] info 0 if the stopping test is satisfied, 1 otherwise
 * \param[in,out] options solver options structure
 */
void mixedProblem_NSGS(MixedProblem* problem, double* z, double* w, int* info,
                       SolverOptions* options);

/**
 * Initialize the MixedProblem NSGS solver.
 *
 * \param[in] problem the MixedProblem
 * \param[in,out] options solver options structure
 */
void mixedProblem_NSGS_init(MixedProblem* problem, SolverOptions* options);

/**
 * Free memory allocated for the MixedProblem NSGS solver.
 *
 * \param[in] problem the MixedProblem
 * \param[in,out] options solver options structure
 */
void mixedProblem_NSGS_free(MixedProblem* problem, SolverOptions* options);

/**
 * Set default solver options for MixedProblem NSGS.
 *
 * \param[in] problem the MixedProblem
 * \param[in,out] options solver options structure
 */
void mixedProblem_NSGS_set_default_solver_options(MixedProblem* problem,
                                                  SolverOptions* options);

/**
 * Full error of a MixedProblem: w = M z + q is computed, then
 *   error = sqrt(sum_e |z_e - P_{K_e}(z_e - (w_e + correction_e))|^2) / |q|
 * (no division if |q| <= DBL_EPSILON), with the unitary error functions of the
 * single-model solvers: fc2d/fc3d (De Saxce term mu |w_T|), rolling friction 2D/3D,
 * plasticity_2d (Drucker-Prager, dilatancy term theta |w_T|); for Von Mises,
 * z_e - P_C(z_e - w_e), C the Von Mises domain.
 *
 * \param[in] problem the MixedProblem
 * \param[in] z the solution vector
 * \param[out] w w = M z + q
 * \param[in] tolerance error tolerance
 * \param[in,out] options solver options (unused)
 * \param[out] error the computed error
 * \return 0 if error <= tolerance, 1 otherwise (or a numerics error code)
 */
int mixedProblem_compute_error(MixedProblem* problem, double* z, double* w, double tolerance,
                               SolverOptions* options, double* error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif /* MIXEDPROBLEM_NSGS_H */
