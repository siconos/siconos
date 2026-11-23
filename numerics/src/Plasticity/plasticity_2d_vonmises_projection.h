// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef PLASTICITY_2D_VONMISES_PROJECTION_H
#define PLASTICITY_2D_VONMISES_PROJECTION_H

/*!\file plasticity_2d_vonmises_projection.h
 * \brief Von Mises (J2) plasticity in cone coordinates: yield domain, projection, error and
 * exact solution of the one-element problem.
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

/** Exact solution of the one-element Von Mises problem
 *    sigma in C, -(M sigma + q) in N_C(sigma)
 *  for a 3x3 matrix M (dense, column-major, M_NN > 0 and positive definite symmetric part of
 *  the Schur complement, not necessarily symmetric).
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
int plasticity_2d_vonMises_local_solve(const double M[9], const double q[3], double sigma_y,
                                       double stress[3]);

#if defined(__cplusplus)
}
#endif

#endif /* PLASTICITY_2D_VONMISES_PROJECTION_H */
