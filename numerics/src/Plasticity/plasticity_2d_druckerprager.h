// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef PLASTICITY_2D_DRUCKERPRAGER_H
#define PLASTICITY_2D_DRUCKERPRAGER_H

/*!\file plasticity_2d_druckerprager.h
  \brief Projection local solvers for the Drucker-Prager (Mohr-Coulomb) model (one cone).
  The Von Mises local solvers are in plasticity_2d_vonmises.h.

  Each solver must have 4 functions in its interface:
  - initialize: link global static variables to the considered problem (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for
  a specific contact
  - solve: solve the local problem
  - free
  We consider a "global" (ie for several contacts) problem, used to initialize the static
  global variables. Then a "local" (ie for one contact => size = 3) problem is built (update
  function) and solved (solve function).

  Two different storages are available for M: dense and sparse block.

*/
#include "NumericsFwd.h"  // for PlasticityProblem, SolverOptions
#include "plasticity_2d_local_newton.h"

#if defined(__cplusplus)
extern "C" {
#endif

/** Initialize Mohr Coulomb 2D projection
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 */
void plasticity_2d_projection_initialize(PlasticityProblem* problem,
                                         PlasticityProblem* localproblem);

/** Initialize Mohr Coulomb 2D projection with regularization
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 */
void plasticity_2d_projection_initialize_with_regularization(PlasticityProblem* problem,
                                                             PlasticityProblem* localproblem);

/** Update Mohr Coulomb 2D projection solver: formalize local problem for one contact.
 * \param number (position in global matrix) of the considered contact
 * \param problem :  the global problem to solve
 * \param localproblem :  the local problem to initialize
 * \param reaction (only the block corresponding to the current contact will be modified,
 *   the rest is used to formalize the local problem)
 * \param options
 */
void plasticity_2d_projection_update(int number, PlasticityProblem* problem,
                                     PlasticityProblem* localproblem, double* reaction,
                                     SolverOptions* options);

/** solve Mohr Coulomb 2D problem with projection on the Cone
 * \param localproblem :  the local problem to initialize
 * \param reaction
 * \param options
 * \return 0 if successfull
 */
int plasticity_2d_projectionOnCone_solve(PlasticityProblem* localproblem, double* reaction,
                                         SolverOptions* options);

/** solve Mohr Coulomb 2D problem with projection on the Cone with local
 *   iteration up to convergence of the local problem
 * \param localproblem :  the local problem to initialize
 * \param reaction
 * \param options
 * \return 0 if successfull
 */
int plasticity_2d_projectionOnConeWithLocalIteration_solve(PlasticityProblem* localproblem,
                                                           double* reaction,
                                                           SolverOptions* options);
void plasticity_2d_projectionOnConeWithLocalIteration_free(PlasticityProblem* problem,
                                                           PlasticityProblem* localproblem,
                                                           SolverOptions* localsolver_options);
void plasticity_2d_projectionOnConeWithLocalIteration_initialize(
    PlasticityProblem* problem, PlasticityProblem* localproblem,
    SolverOptions* localsolver_options);
/** free memory for Mohr Coulomb 2D projection solver
 * \param problem :  the  problem to free
 * \param localproblem :  the  problem to free
 * \param localsolver_options
 */
void plasticity_2d_projection_free(PlasticityProblem* problem, PlasticityProblem* localproblem,
                                   SolverOptions* localsolver_options);

/* ===========================================================================
 * One cone: return mapping and semi-smooth Newton on the natural map
 *   sigma in K, -(w + theta ||w_T|| e_N) in N_K(sigma), w = M sigma + q,
 *   K = {||sigma_T|| <= eta sigma_N} (without cohesion: see plasticity_2d_driver)
 * =========================================================================== */

/** Natural-map residual ||sigma - P_K(sigma - rho (w + theta ||w_T|| e_N))|| of one cone.
 *  \param M local matrix (3x3, column-major)
 *  \param q local vector
 *  \param eta cone coefficient
 *  \param theta dilatancy coefficient
 *  \param rho rho > 0
 *  \param stress stress
 *  \return the residual
 */
double plasticity_2d_DruckerPrager_natural_map_residual(const double M[9], const double q[3],
                                                        double eta, double theta, double rho,
                                                        const double stress[3]);

/** Generalized Jacobian (3x3, column-major) of the projection on the Drucker-Prager cone
 *  ||y_T|| <= eta y_N at y (same cases as projectionOnCone: polar cone, inside, boundary).
 *  Used by the local and global semi-smooth Newton methods.
 *  \param y the point
 *  \param eta cone coefficient
 *  \param[out] J the Jacobian
 */
void plasticity_2d_DruckerPrager_projection_jacobian(const double y[3], double eta,
                                                     double J[9]);

/** Return mapping of one cone: exact solution by enumeration of the cases, for any 3x3 M
 *  (not necessarily symmetric) and any dilatancy theta:
 *  - elastic: w = 0, sigma = -M^{-1} q (if in K);
 *  - apex: sigma = 0 (if q + theta ||q_T|| e_N is in the dual cone);
 *  - sliding: sigma = sigma_N (1, eta t), w = lambda (eta - theta, -t), t = (cos phi, sin
 *    phi), sigma_N > 0, lambda >= 0: sigma_N M a(phi) - lambda b(phi) = -q. Solvability
 *    det[M a, b, q] = 0 is a trigonometric polynomial of degree 2 in phi, i.e. a quartic in
 *    tan(phi/2) (solved by BIQUADROOTS); each root is polished by Newton's method on
 *    (sigma_N, lambda, phi).
 *  The candidates with a zero natural-map residual (up to round-off) are solutions. For a
 *  symmetric positive definite M and theta = eta (associated flow) the solution is unique
 *  (closest-point projection of the trial stress -M^{-1} q on K in the M-norm); otherwise
 *  the solution closest to the input stress is returned (continuity in the Gauss-Seidel
 *  iterations).
 *  \param M local matrix (3x3, column-major)
 *  \param q local vector
 *  \param eta cone coefficient
 *  \param theta dilatancy coefficient
 *  \param[in,out] stress input: reference stress for the choice among several solutions;
 *  output: solution
 *  \return 0 if a solution has been found, 1 otherwise (the candidate with the smallest
 *  residual is returned)
 */
int plasticity_2d_DruckerPrager_return_mapping(const double M[9], const double q[3],
                                               double eta, double theta, double stress[3]);

/** Semi-smooth Newton method on the natural map of one cone,
 *    Phi(sigma) = sigma - P_K(sigma - rho v(sigma)),  v = w + theta ||w_T|| e_N,
 *  with H = I - J_P(y)(I - rho J_v), J_v = (I + theta e_N [0, n^T]) M, n = w_T/||w_T||, and
 *  a line search (see plasticity_2d_local_semismooth_newton). If p->hybrid != 0, falls back
 *  to the return mapping when the method does not converge.
 *  \param M local matrix (3x3, column-major)
 *  \param q local vector
 *  \param eta cone coefficient
 *  \param theta dilatancy coefficient
 *  \param[in,out] stress initial guess / solution
 *  \param[in,out] p parameters and results (p->ball unused)
 *  \return 0 if the residual is below p->tol, 1 otherwise
 */
int plasticity_2d_DruckerPrager_semismooth_newton(const double M[9], const double q[3],
                                                  double eta, double theta, double stress[3],
                                                  Plasticity2DLocalParameters* p);

/** Local solver PLASTICITY_2D_ONECONE_DRUCKERPRAGER_RETURN_MAPPING (SolverOptions
 *  interface). \return 0 if successful */
int plasticity_2d_DruckerPrager_return_mapping_solve(PlasticityProblem* localproblem,
                                                     double* stress_local,
                                                     SolverOptions* options);

/** Local solver PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN (SolverOptions interface): same
 *  options as PLASTICITY_2D_ONECONE_VONMISES_NSN. \return 0 if converged */
int plasticity_2d_DruckerPrager_nsn_solve(PlasticityProblem* localproblem,
                                          double* stress_local, SolverOptions* options);

/** Default options of PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN (rho = 1/||M||_2, Armijo line
 *  search, 20 steps, fallback to the return mapping) */
void plasticity_2d_druckerprager_nsn_set_default(SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
