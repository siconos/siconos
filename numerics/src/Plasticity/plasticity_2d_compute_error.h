// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef plasticity_2d_compute_error_H
#define plasticity_2d_compute_error_H

/*!\file plasticity_2d_compute_error.h
  \brief functions related to error computation for Mohr Coulomb 2D plasticity

*/

#include "NumericsFwd.h"  // for PlasticityProblem, SolverOptions

#if defined(__cplusplus)
#undef restrict
#define restrict __restrict
extern "C" {
#endif

/**
    Error of a plasticity problem (Drucker-Prager or Von Mises), scaled as in
    fc3d_compute_error.

    For each element e, in the variables of the solvers (cone coordinates, cohesion shift,
    plane stress scaling), with r_e the stress and v_e the modified strain rate
    (w + theta ||w_T|| e_N for Drucker-Prager, w for Von Mises), the natural-map residual
    F_e = r_e - P_K(r_e - rho_e v_e) is divided by a scale s computed from r_e and rho_e v_e.
    rho balances the stress and the strain rate (options->iparam[PLASTICITY_SCALING_RHO]):
    rho_e = 1/||M_ee||_2 per element (default), 1/max_e ||M_ee||_2, or 1 (as in fc3d); the
    diagonal blocks are taken in the variables of the solvers. Then:
    - PER_CONE scope (default): error = sqrt(sum_e ||F_e||^2 / s_e^2), s_e from r_e, v_e;
    - PER_VECTOR scope: error = ||F|| / s, s from the whole vectors r, v.
    The scale combines the norms of r and v, each safeguarded by max(1, .):
    sum (PLASTICITY_SCALING_REDUCTION_MIXED, default), min of the two (MAX, the fc3d
    default), or geometric mean (MEAN); norms: Euclidean (default), cone |x_N| + ||x_T||, cone
   weighted by eta (Von Mises: eta = 1), sqrt(x_N^2 + alpha ||x_T||^2), or no scaling (absolute
   error). Options: options->iparam[PLASTICITY_SCALING_NORM / _SCOPE / _REDUCTION / _RHO],
    options->dparam[PLASTICITY_SCALING_ALPHA], set by plasticity_2d_nsgs_set_default
    (Euclidean, per cone, MIXED, rho_e); options == NULL: same defaults. Options of other
    solvers (zero-initialized) give the fc3d scaling (MAX, rho = 1).
    w is set to M z + q (in the coordinates of the problem).

    \param problem the problem
    \param z stress
    \param[out] w strain rate M z + q
    \param tolerance value for error computation
    \param options solver options (scaling), may be NULL
    \param norm unused
    \param[in,out] error value
    \return 0 if error <= tolerance, 1 otherwise
 */
/** Scale of the error (see plasticity_2d_compute_error): combination of the norms of r and v
 *  (m elements of dimension 3), each safeguarded by max(1, .).
 *  \param norm_id PLASTICITY_SCALING norm
 *  \param reduction PLASTICITY_SCALING_REDUCTION value
 *  \param r stresses (3 m)
 *  \param v modified strain rates (3 m)
 *  \param m number of elements
 *  \param eta cone coefficients (PLASTICITY_SCALING_CONE_ETA only, size m)
 *  \param alpha weight (PLASTICITY_SCALING_CONE_WEIGHTED only)
 *  \return the scale (>= 1)
 */
double plasticity_2d_scaling(int norm_id, int reduction, const double* r, const double* v,
                             int m, const double* eta, double alpha);

int plasticity_2d_compute_error(PlasticityProblem* problem, double* z, double* w,
                                double tolerance, SolverOptions* options, double norm,
                                double* error);

/**
    Error computation (using the normal map residual) for one Mohr Coulomb 2D plasticity

    \param r the reaction force
    \param u the local velocity
    \param eta coefficient of the cone
    \param theta coefficient of the dilatency
    \param worktmp work vector
    \param[in,out] error value
 */
void plasticity_2d_unitary_compute_and_add_error(double* restrict r, double* restrict u,
                                                 double eta, double theta, double* error,
                                                 double* worktmp);

/** Check whether sigma = 0 solves the one-cone Drucker-Prager problem with w = q:
 *  q + theta ||q_T|| e_N must be in the dual cone K_eta^* = {eta ||v_T|| <= v_N}, i.e.
 *  eta ||q_T|| <= q_N + theta ||q_T||. (q_N >= 0 is enough only when theta >= eta.)
 *  \param[in] q the local vector (w for sigma = 0)
 *  \param[in] eta cone coefficient
 *  \param[in] theta dilatancy coefficient
 *  \return 1 if sigma = 0 is a solution, 0 otherwise
 */
int plasticity_2d_DruckerPrager_zero_is_solution(const double q[3], double eta, double theta);

#if defined(__cplusplus)
}
#endif

#endif
