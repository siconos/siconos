// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef PLASTICITY_2D_LOCAL_NEWTON_H
#define PLASTICITY_2D_LOCAL_NEWTON_H

/*!\file plasticity_2d_local_newton.h
 * \brief Semi-smooth Newton method on the natural map of one element (3 unknowns), shared by
 * the Von Mises and Drucker-Prager local solvers (see docs/devel/PLASTICITY.md).
 *
 * The model provides the natural map Phi(sigma) = sigma - P_K(sigma - rho v(sigma)) and an
 * element of its generalized Jacobian H = I - J_P(y)(I - rho J_v), y = sigma - rho v.
 * The Newton direction d solves H d = -Phi; the step is chosen by a line search on
 * 1/2 ||Phi||^2 (Armijo, Goldstein-Price or none). The stopping test is
 * ||Phi_{rho0}(sigma)|| <= tol (1 + ||sigma||), rho0 = 1/||M||_2.
 */

#if defined(__cplusplus)
extern "C" {
#endif

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
  int hybrid;              /**< Newton: if not 0, fall back to the exact (return mapping)
                                solver on failure */
  int ball;        /**< Von Mises yield set: 0 cylinder ||sigma_T|| <= sqrt(2/3) sigma_y (plane
                        strain), 1 ball ||sigma|| <= sigma_y (plane stress, scaled variables) */
  int iter_done;   /**< (out) number of iterations */
  double residual; /**< (out) natural-map residual, rho = 1 / ||M||_2 */
} Plasticity2DLocalParameters;

/** Natural map of one element: returns ||Phi||, Phi = sigma - P_K(sigma - rho v(sigma)),
 *  and, if H is not NULL, an element of the generalized Jacobian of Phi (3x3,
 *  column-major).
 *  \param data model data
 *  \param rho rho > 0
 *  \param stress sigma
 *  \param[out] Phi natural map
 *  \param[out] H generalized Jacobian, or NULL
 */
typedef double (*Plasticity2DNaturalMap)(void* data, double rho, const double stress[3],
                                         double Phi[3], double H[9]);

/** Exact solver used as a fallback (return mapping); returns 0 if successful. */
typedef int (*Plasticity2DExactSolver)(void* data, double stress[3]);

/** ||M||_2 (largest singular value) of a 3x3 column-major matrix, by power iteration on
 *  M^T M. */
double plasticity_2d_local_spectral_norm(const double M[9]);

/** Semi-smooth Newton method on the natural map of one element.
 *  rho = 1/||M||_2, or p->rho with the constant strategy (p->rho > 0). If the method does
 *  not converge and p->hybrid != 0, the exact solver (if not NULL) is called.
 *  \param natural_map natural map and generalized Jacobian of the model
 *  \param exact exact solver (fallback), may be NULL
 *  \param data model data passed to natural_map and exact
 *  \param M local matrix (for rho)
 *  \param[in,out] stress initial guess / solution
 *  \param[in,out] p parameters and results
 *  \return 0 if the residual is below p->tol, 1 otherwise
 */
int plasticity_2d_local_semismooth_newton(Plasticity2DNaturalMap natural_map,
                                          Plasticity2DExactSolver exact, void* data,
                                          const double M[9], double stress[3],
                                          Plasticity2DLocalParameters* p);

#if defined(__cplusplus)
}
#endif

#endif
