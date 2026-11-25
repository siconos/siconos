/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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
#ifndef PLASTICITY_2D_STRESS_COORDINATES_H
#define PLASTICITY_2D_STRESS_COORDINATES_H

/*!\file plasticity_2d_stress_coordinates.h
  Change of variables between Cartesian stresses and the cone coordinates used by the
  Plasticity solvers (see PlasticityStressCoordinates and
  docs/devel/PLASTICITY.md). Cone coordinates are hydrostatic-deviatoric coordinates:
  sigma_N = -(sxx + syy)/sqrt2 is the hydrostatic part (compression positive),
  (sigma_T1, sigma_T2) = ((sxx - syy)/sqrt2, sqrt2 sxy) the deviator.

  On each element, sigma_c = P s and w_c = P^{-T} w_s, with (Mandel notation, P orthogonal)

      P = [ -1/sqrt2  -1/sqrt2  0 ]
          [  1/sqrt2  -1/sqrt2  0 ]
          [  0         0        1 ]

  and P_Voigt = P diag(1, 1, sqrt2). The variational inequality is invariant under this
  change of variables: s solves the Cartesian problem (M, q) iff P s solves the cone
  problem (P^{-T} M P^{-1}, P^{-T} q), with the same yield sets expressed in cone
  coordinates (Drucker-Prager cone, dilatancy and cohesion; Von Mises cylinder).
*/

#include "PlasticityProblem.h"
#include "SiconosConfig.h"  // for BUILD_AS_CPP

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/** sigma_c = P s on nc elements (in place allowed: s == sigma_c).
    \param coordinates coordinates of s (CONE: copy)
    \param nc number of elements
    \param s stresses in the given coordinates (size 3 nc)
    \param[out] sigma_c stresses in cone coordinates (size 3 nc)
*/
void plasticity_2d_stress_to_cone(PlasticityStressCoordinates coordinates, int nc,
                                  const double* s, double* sigma_c);

/** s = P^{-1} sigma_c on nc elements (in place allowed). */
void plasticity_2d_stress_from_cone(PlasticityStressCoordinates coordinates, int nc,
                                    const double* sigma_c, double* s);

/** w_c = P^{-T} w_s on nc elements (strain rate, dual of the stress; in place allowed). */
void plasticity_2d_strainrate_to_cone(PlasticityStressCoordinates coordinates, int nc,
                                      const double* w_s, double* w_c);

/** w_s = P^T w_c on nc elements (in place allowed). */
void plasticity_2d_strainrate_from_cone(PlasticityStressCoordinates coordinates, int nc,
                                        const double* w_c, double* w_s);

/** M_c = P^{-T} M P^{-1}, P block diagonal (one 3x3 block per element).
    Dense M: dense result. Sparse-block M whose block sizes are multiples of 3: sparse-block
    result with the same structure. Otherwise: NM_SPARSE result.
    \param coordinates coordinates of M (CONE: copy)
    \param M the matrix (3 nc x 3 nc)
    \return a new matrix (to be freed with NM_free)
*/
NumericsMatrix* plasticity_2d_matrix_to_cone(PlasticityStressCoordinates coordinates,
                                             NumericsMatrix* M);

/** The problem in cone coordinates (deep copy, with M_c = P^{-T} M P^{-1},
    q_c = P^{-T} q and stress_coordinates = CONE).
    \param problem a problem in any coordinates
    \return a new problem (to be freed with plasticityProblem_free)
*/
PlasticityProblem* plasticity_2d_problem_to_cone(PlasticityProblem* problem);

/** Change of variables of the problem solved by plasticity_2d_driver: on every element,
    s = T y, s the stress in the coordinates of the problem and y the unknown of the solved
    problem, with T = P^{-1} D:
    - P the Cartesian -> cone map (identity in cone coordinates);
    - D = diag(sqrt(2), sqrt(2/3), sqrt(2/3)) for the plane stress Von Mises model, which maps
      the ellipsoid sigma_N^2/2 + 3/2 ||sigma_T||^2 <= sigma_y^2 to the ball
      ||y|| <= sigma_y (identity otherwise).
    The solved problem is M_y = T^T M T, q_y = T^T q, w_y = T^T w.
    \param problem the problem
    \param[out] T the 3x3 map (row-major)
    \param[out] Tinv its inverse (row-major)
    \return 0 if T = I (nothing to do), 1 otherwise
*/
int plasticity_2d_change_of_variables(const PlasticityProblem* problem, double T[9],
                                      double Tinv[9]);

/** y_e = A x_e (or A^T x_e if transpose) on nc elements, A 3x3 row-major (in place allowed).
 */
void plasticity_2d_block_map(const double A[9], int transpose, int nc, const double* x,
                             double* y);

/** T^T M T, T block diagonal with the 3x3 block T (row-major); same storage rules as
    plasticity_2d_matrix_to_cone.
    \return a new matrix (to be freed with NM_free) */
NumericsMatrix* plasticity_2d_matrix_change_of_variables(const double T[9], NumericsMatrix* M);

/** The problem solved by plasticity_2d_driver (deep copy): cone coordinates and, for the
    plane stress Von Mises model, the scaled variables y (model PLASTICITY_VONMISES_BALL),
    M_y = T^T M T, q_y = T^T q (see plasticity_2d_change_of_variables).
    \return a new problem (to be freed with plasticityProblem_free) */
PlasticityProblem* plasticity_2d_problem_change_of_variables(PlasticityProblem* problem);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
