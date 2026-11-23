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
#ifndef PLASTICITY_2D_VONMISES_H
#define PLASTICITY_2D_VONMISES_H

/*!\file plasticity_2d_vonmises.h
 * \brief Von Mises (J2) plasticity in cone coordinates: yield domain, projection, error, and
 * the one-element solvers (return mapping, fixed point and semi-smooth Newton on the natural
 * map) with their SolverOptions interface used by the NSGS solvers.
 *
 * The stress of an element is given in the same coordinates as for the Drucker-Prager model,
 *   sigma = (sigma_N, sigma_T1, sigma_T2),
 * sigma_N the hydrostatic coordinate and sigma_T the deviatoric coordinates in an orthonormal
 * basis (Frobenius inner product), so that ||dev sigma|| = ||sigma_T||. With the 3D J2
 * convention sigma_eq = sqrt(3/2) ||dev sigma||, the Von Mises domain is the cylinder
 *
 *   C = { sigma : ||sigma_T|| <= R },   R = sqrt(2/3) sigma_y,
 *
 * sigma_y being the uniaxial yield stress (parameter sigma_y of Plasticity_VonMises_model).
 *
 * Problem for one element, w = M sigma + q (w: opposite of the plastic strain rate, same
 * convention as Drucker-Prager, where -(w + theta ||w_T|| e_N) is in the normal cone):
 *
 *   sigma in C,  -w in N_C(sigma)   <=>   sigma = P_C(sigma - rho w)  for any rho > 0,
 *
 * i.e. w_N = 0 (isochoric plastic flow, sigma_N is not constrained) and
 * w_T = -lambda sigma_T / R, lambda >= 0, lambda (R - ||sigma_T||) = 0.
 * See docs/devel/PLASTICITY_VON_MISES.md.
 */

#include "NumericsFwd.h"
#include "PlasticityProblem.h"
#include "SiconosConfig.h"

#if defined(__cplusplus)
extern "C" {
#endif

/** Radius of the Von Mises cylinder, R = sqrt(2/3) sigma_y.
 *  \param[in] sigma_y uniaxial yield stress
 *  \return R
 */
double plasticity_2d_vonMises_radius(double sigma_y);

/** Von Mises equivalent stress, sigma_eq = sqrt(3/2) ||sigma_T||.
 *  \param[in] stress (sigma_N, sigma_T1, sigma_T2)
 *  \return sigma_eq
 */
double plasticity_2d_vonMises_equivalent_stress(const double stress[3]);

/** Check if a stress is in the Von Mises domain.
 *  \param[in] stress (sigma_N, sigma_T1, sigma_T2)
 *  \param[in] sigma_y uniaxial yield stress
 *  \return 1 if sigma_eq <= sigma_y (elastic), 0 otherwise
 */
int plasticity_2d_vonMises_check_yield(const double stress[3], double sigma_y);

/** Euclidean projection on the Von Mises domain C (radial return of sigma_T, sigma_N
 *  unchanged).
 *  \param[in,out] stress (sigma_N, sigma_T1, sigma_T2)
 *  \param[in] sigma_y uniaxial yield stress
 *  \return 1 if the stress was in C (unchanged), 0 if it has been projected
 */
int plasticity_2d_projectionOnVonMises(double stress[3], double sigma_y);

/** Add to *error the squared natural-map residual |r - P_C(r - u)|^2 of one element.
 *  \param[in] r stress of the element
 *  \param[in] u w = M sigma + q of the element
 *  \param[in] sigma_y uniaxial yield stress
 *  \param[in,out] error accumulated squared residual
 *  \param worktmp work array of size 3
 */
void plasticity_2d_vonmises_unitary_compute_and_add_error(const double r[3], const double u[3],
                                                          double sigma_y, double* error,
                                                          double worktmp[3]);

/** Return mapping: exact solution of the one-element Von Mises problem
 *    sigma in C, -(M sigma + q) in N_C(sigma)
 *  for a 3x3 matrix M (dense, column-major, M_NN > 0 and positive definite symmetric part of
 *  the Schur complement, not necessarily symmetric). For a symmetric positive definite M
 *  (compliance), this is the closest-point projection, in the M-norm, of the trial stress
 *  -M^{-1} q on C (return mapping).
 *
 *  w_N = 0 gives sigma_N = -(q_N + M_NT sigma_T) / M_NN and the 2D problem
 *    sigma_T in D_R, -(S sigma_T + q~) in N_{D_R}(sigma_T),
 *  S = M_TT - M_TN M_NT / M_NN, q~ = q_T - M_TN q_N / M_NN, D_R the disk of radius R.
 *  Elastic case: sigma_T = -S^{-1} q~ if ||.|| <= R. Plastic case: (S + mu I) sigma_T = -q~,
 *  ||sigma_T|| = R, mu = lambda / R > 0, mu solution of 1/||sigma_T(mu)|| = 1/R, found by a
 *  safeguarded Newton / bisection method.
 *
 *  \param[in] M local matrix (3x3, column-major)
 *  \param[in] q local vector
 *  \param[in] sigma_y uniaxial yield stress
 *  \param[out] stress solution
 *  \return 0 if successful, 1 if M_NN <= 0 or S singular, 2 if the scalar equation was not
 *  solved to machine precision
 */
int plasticity_2d_vonMises_return_mapping(const double M[9], const double q[3], double sigma_y,
                                          double stress[3]);

/** Parameters and results of the iterative one-element solvers (fixed point and
 *  semi-smooth Newton on the natural map). */
typedef struct {
  int max_iter;            /**< maximum number of iterations */
  double tol;              /**< tolerance on the natural-map residual (rho = 1 / ||M||_2) */
  int rho_strategy;        /**< PLASTICITY_NSN_RHO_STRATEGY values (constant, spectral norm,
                                adaptive) */
  double rho;              /**< rho for the constant strategy (<= 0: 1 / ||M||_2) */
  int linesearch;          /**< PLASTICITY_NSN_LINESEARCH values (Newton) */
  int linesearch_max_iter; /**< maximum number of line-search / rho-adaptation steps */
  int hybrid;              /**< Newton: if not 0, fall back to the return mapping on failure */
  int iter_done;           /**< (out) number of iterations */
  double residual;         /**< (out) natural-map residual, rho = 1 / ||M||_2 */
} Plasticity2DVonMisesLocalParameters;

/** Natural-map residual |sigma - P_C(sigma - rho (M sigma + q))| of one element.
 *  \param[in] M local matrix (3x3, column-major)
 *  \param[in] q local vector
 *  \param[in] sigma_y uniaxial yield stress
 *  \param[in] rho rho > 0
 *  \param[in] stress stress
 *  \return the residual
 */
double plasticity_2d_vonMises_natural_map_residual(const double M[9], const double q[3],
                                                   double sigma_y, double rho,
                                                   const double stress[3]);

/** Fixed point on the natural map for one element:
 *    sigma <- P_C(sigma - rho (M sigma + q)),
 *  with rho constant (PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_CONSTANT), equal to 1/||M||_2
 *  (SPECTRAL_NORM strategies) or adapted at each iteration
 *  (PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE, default: rho is decreased until
 *  rho ||w(sigma+) - w(sigma)|| <= 0.9 ||sigma+ - sigma||, increased when the ratio is
 *  below 0.3, as in the friction solvers). Converges when M is positive definite (any
 *  symmetric part) and rho is small enough.
 *  \param[in] M local matrix (3x3, column-major)
 *  \param[in] q local vector
 *  \param[in] sigma_y uniaxial yield stress
 *  \param[in,out] stress initial guess / solution
 *  \param[in,out] p parameters and results
 *  \return 0 if the residual is below p->tol, 1 otherwise
 */
int plasticity_2d_vonMises_fixed_point(const double M[9], const double q[3], double sigma_y,
                                       double stress[3],
                                       Plasticity2DVonMisesLocalParameters* p);

/** Semi-smooth Newton method on the natural map of one element,
 *    Phi(sigma) = sigma - P_C(sigma - rho (M sigma + q)) = 0,
 *  with the generalized Jacobian H = I - J_P(y) (I - rho M), y = sigma - rho w, J_P the
 *  Jacobian of the projection (identity inside C, (R/||y_T||)(I - n n^T) on the deviatoric
 *  part outside, n = y_T/||y_T||), and a line search on 1/2 ||Phi||^2 (Armijo,
 *  Goldstein-Price or none). rho = 1/||M||_2 unless the constant strategy is chosen. If the
 *  method fails and p->hybrid != 0, the return mapping is used.
 *  \param[in] M local matrix (3x3, column-major)
 *  \param[in] q local vector
 *  \param[in] sigma_y uniaxial yield stress
 *  \param[in,out] stress initial guess / solution
 *  \param[in,out] p parameters and results
 *  \return 0 if the residual is below p->tol, 1 otherwise
 */
int plasticity_2d_vonMises_semismooth_newton(const double M[9], const double q[3],
                                             double sigma_y, double stress[3],
                                             Plasticity2DVonMisesLocalParameters* p);

/* ============================================================================
 * Local solvers (SolverOptions interface, used by the NSGS solvers)
 * ============================================================================ */

/** Von Mises (J2) one-element problem (cone coordinates, see
 * plasticity_2d_vonmises.h): one projection step
 *   sigma <- P_C(sigma - rho (M sigma + q)),  rho = 1 / M_NN.
 * \param localproblem the local problem (one element, model Von Mises)
 * \param stress_local the local stress (input/output)
 * \param options solver options (unused)
 * \return 0
 */
int plasticity_2d_projectionOnVonMises_solve(PlasticityProblem* localproblem,
                                             double* stress_local, SolverOptions* options);

/** Von Mises (J2) one-element problem: return mapping (exact solution,
 * plasticity_2d_vonMises_return_mapping).
 * \param localproblem the local problem (one element, model Von Mises)
 * \param stress_local the local stress (output)
 * \param options solver options (may be NULL)
 * \return 0 if successful, 1 if the local matrix is not suitable (M_NN <= 0, S singular)
 */
int plasticity_2d_vonMises_return_mapping_solve(PlasticityProblem* localproblem,
                                                double* stress_local, SolverOptions* options);

/** Von Mises (J2) one-element problem: fixed point on the natural map
 * (plasticity_2d_vonMises_fixed_point). Options: iparam[SICONOS_IPARAM_MAX_ITER],
 * dparam[SICONOS_DPARAM_TOL], iparam[PLASTICITY_NSN_RHO_STRATEGY] (constant, spectral norm,
 * adaptive (default)), dparam[PLASTICITY_NSN_RHO] (constant rho, 0: 1/||M||_2),
 * iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER] (rho adaptation steps).
 * \param localproblem the local problem (one element, model Von Mises)
 * \param stress_local the local stress (initial guess / solution)
 * \param options solver options
 * \return 0 if converged
 */
int plasticity_2d_vonMises_fixed_point_solve(PlasticityProblem* localproblem,
                                             double* stress_local, SolverOptions* options);

/** Von Mises (J2) one-element problem: semi-smooth Newton on the natural map
 * (plasticity_2d_vonMises_semismooth_newton). Options: iparam[SICONOS_IPARAM_MAX_ITER],
 * dparam[SICONOS_DPARAM_TOL], iparam[PLASTICITY_NSN_RHO_STRATEGY] (constant or spectral norm
 * (default)), dparam[PLASTICITY_NSN_RHO], iparam[PLASTICITY_NSN_LINESEARCH] (Armijo
 * (default), Goldstein-Price, none), iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER],
 * iparam[PLASTICITY_NSN_HYBRID_STRATEGY] (not NO: fall back to the return mapping on
 * failure, default).
 * \param localproblem the local problem (one element, model Von Mises)
 * \param stress_local the local stress (initial guess / solution)
 * \param options solver options
 * \return 0 if converged
 */
int plasticity_2d_vonMises_nsn_solve(PlasticityProblem* localproblem, double* stress_local,
                                     SolverOptions* options);

/** Default options of the Von Mises fixed point local solver
 * \param options the solver options to set
 */
void plasticity_2d_vonmises_fp_set_default(SolverOptions* options);

/** Default options of the Von Mises semi-smooth Newton local solver
 * \param options the solver options to set
 */
void plasticity_2d_vonmises_nsn_set_default(SolverOptions* options);

/** Set default solver options for Von Mises projection
 * \param options the solver options to set
 */
void plasticity_2d_vonmises_set_default(SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif /* PLASTICITY_2D_VONMISES_H */
