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

#ifndef JORDAN_ALGEBRA_H
#define JORDAN_ALGEBRA_H

/*!\file JordanAlgebra.h
  \brief Functions of the Jordan algebra.
*/

#include "NumericsMatrix.h"  // IWYU pragma: keep - For NM_Types

typedef long double float_type;
/* typedef double float_type; */

/** \return the Arrow representation matrix from vector.
 * \param[in] vec pointer to the vector data.
 * \param[in] vecSize the length of the vector.
 * \param[in] varsCount the count of variables (subvectors) in vec.
 */
NumericsMatrix* Arrow_repr(const double* const vec, size_t vecSize, size_t varsCount);

/** Replace the content of an existing Arrow representation matrix from vector.
 * This is an alternative of Arrow_repr, without creation of matrix
 * (not allocation of memory)
 */
int Arrow_repr_replace(NumericsMatrix* Arw_mat, const double* const vec, const size_t vecSize,
                       const size_t varsCount);

/**
 * Returns reflection matrix.

            |1  0 ...  0|
        R = |0 -1 ...  0|
            | ........ 0|
            |0  0 ... -1|

 * \param size is the size of rectangular metrix
 * \return reflection matrix.
 */
NumericsMatrix* Reflect_mat(const size_t size, NM_types type);

/**
 * Returns quadratic representation of the vector x by formula 2*xx^T - det(x)R
 * \param vec pointer to the vector data.
 * \param vecSize the length of the vector.
 * \param varsCount the count of variables (subvectors) in vec.
 * \return a pointer to a NumericsMatrix
 */
NumericsMatrix* Quad_repr(const double* const vec, const size_t vecSize,
                          const size_t varsCount);

/** Nesterov-Todd vector for the pair (vec1, vec2).
 * \param vec1 first input vector
 * \param vec2 second input vector
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors)
 * \param out result NT vector
 */
void NesterovToddVector(const double* const vec1, const double* const vec2,
                        const size_t vecSize, const size_t varsCount, double* out);

/** Create the Arrow representation matrix from vector.
 * \param vecSize the length of the vector.
 * \param varsCount the count of variables (subvectors) in vec.
 * \return a pointer to the identity element of the Jordan algebra
 */
double* JA_iden(const size_t vecSize, const size_t varsCount);

/** Jordan product of two vectors.
 * \param vec1 is a vector;
 * \param vec2 is a vector;
 * \param vecSize is the length of the vectors;
 * \param varsCount is the count of variables (subvectors) in x and y.
 * \param out is the result vector of the Jordan product.
 */
void JA_prod(const double* const vec1, const double* const vec2, const size_t vecSize,
             const size_t varsCount, double* out);

/** Returns the eigenvalues of each element in the vector.
 * \param vec is the pointer to the vector data.
 * \param vecSize is the length of the vector.
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the result vector of eigenvalues (out\in\mathbb{R}^{2*varsCount})..
 */
void JA_eigenvals(const double* const vec, const size_t vecSize, const size_t varsCount,
                  double* out);

/** Returns the eigenvectors of each element in the vector.
 * \param vec is the pointer to the vector data.
 * \param vecSize is the length of the vector.
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the result matrix of eigenvactors (out\in\mathbb{R}^{vecSize\times
 * 2*varsCount}).
 */
void JA_eigenvecs(const double* const vec, const size_t vecSize, const size_t varsCount,
                  double** out);

/** Compute element by element square root
 * \param vec is the vector
 * \param vecSize is the size of the vector vec
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the sqrt vector
 */
void JA_sqrt(const double* const vec, const size_t vecSize, const size_t varsCount,
             double* out);

/** Compute element by element inverse of square root
 * \param vec is the vector
 * \param vecSize is the size of the vector vec
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the inverse of sqrt vector
 */
void JA_sqrt_inv(const double* const vec, const size_t vecSize, const size_t varsCount,
                 double* out);

/** Compute element by element the square of the vector
 * \param vec is the vector
 * \param vecSize is the size of the vector vec
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the square of the vector
 */
void JA_power2(const double* const vec, const size_t vecSize, const size_t varsCount,
               double* out);

/** Compute element by element the inverse of the vector
 * \param vec is the vector
 * \param vecSize is the size of the vector vec
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the inverse of the vector
 */
void JA_inv(const double* const vec, const size_t vecSize, const size_t varsCount,
            double* out);

/** Compute element by element square root
 * \param vec is the vector
 * \param vecSize is the size of the vector vec
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the vector of determinants
 */
void JA_det(const double* const vec, const size_t vecSize, const size_t varsCount,
            double* out);

/** Euclidean 2-norm of a vector, computed in long double precision.
 * \param n size of the vector
 * \param x input vector
 * \return long-double value ||x||_2
 */
float_type dnrm2l(const size_t n, const double* x);

/** Squared Euclidean 2-norm of a vector, computed in long double precision.
 * \param n size of the vector
 * \param x input vector
 * \return long-double value ||x||_2^2
 */
float_type dnrm2sqrl(const size_t n, const double* x);

/** Quadratic representation applied to a vector: out = Q(x)^{1/2} y.
 * \param x input vector defining the quadratic representation
 * \param y input vector being mapped
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors) in x and y
 * \param out result vector Q(x)^{1/2} y
 */
void Qx05y(const double* const x, const double* const y, const size_t vecSize,
           const size_t varsCount, double* out);

/** Inverse quadratic representation applied to a vector: out = Q(x)^{-1/2} y.
 * \param x input vector defining the quadratic representation
 * \param y input vector being mapped
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors) in x and y
 * \param out result vector Q(x)^{-1/2} y
 */
void Qx50y(const double* const x, const double* const y, const size_t vecSize,
           const size_t varsCount, double* out);

/** Jordan algebra inverse: out = x^{-1}.
 * \param x input vector
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors) in x
 * \param out result vector x^{-1}
 */
void Jinv(const double* const x, const size_t vecSize, const size_t varsCount, double* out);

/** Jordan algebra square root: out = x^{1/2}.
 * \param x input vector
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors) in x
 * \param out result vector x^{1/2}
 */
void Jsqrt(const double* const x, const size_t vecSize, const size_t varsCount, double* out);

/** Jordan algebra inverse square root: out = x^{-1/2}.
 * \param x input vector
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors) in x
 * \param out result vector x^{-1/2}
 */
void Jsqrtinv(const double* const x, const size_t vecSize, const size_t varsCount,
              double* out);

/** Nesterov-Todd vector p for the pair (x, y).
 * \param T selector for the formulation variant
 * \param x first input vector
 * \param y second input vector
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors) in x and y
 * \param p result NT vector
 */
void Nesterov_Todd_vector(short T, const double* const x, const double* const y,
                          const size_t vecSize, const size_t varsCount, double* p);

/** Nesterov-Todd vector p for the pair (x, y), alternative formula.
 * \param x first input vector
 * \param y second input vector
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors) in x and y
 * \param p result NT vector
 */
void Nesterov_Todd_vector_b(const double* const x, const double* const y, const size_t vecSize,
                            const size_t varsCount, double* p);

/** Quadratic representation applied to a vector: z = Q(x) y.
 * Uses the formula 2 (x^T y) x - det(x) R y.
 * \param x input vector defining the quadratic representation
 * \param y input vector being mapped
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors) in x and y
 * \param z result vector Q(x) y
 */
void Qxy(const double* const x, const double* const y, const size_t vecSize,
         const size_t varsCount, double* z);

/** Scaling by the Nesterov-Todd vector: out = Q(p) z with p = NT(x, y).
 * \param x first input vector defining the NT scaling point
 * \param y second input vector defining the NT scaling point
 * \param z input vector being scaled
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors)
 * \param out result vector Q(p) z
 */
void QNTpz(const double* const x, const double* const y, const double* const z,
           const size_t vecSize, const size_t varsCount, double* out);

/** Inverse NT scaling: out = Q(p)^{-1} z with p = NT(x, y).
 * \param x first input vector defining the NT scaling point
 * \param y second input vector defining the NT scaling point
 * \param z input vector being scaled
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors)
 * \param out result vector Q(p)^{-1} z
 */
void QNTpinvz(const double* const x, const double* const y, const double* const z,
              const size_t vecSize, const size_t varsCount, double* out);

/** Squared inverse NT scaling: out = Q(p)^{-2} z with p = NT(x, y).
 * \param x first input vector defining the NT scaling point
 * \param y second input vector defining the NT scaling point
 * \param z input vector being scaled
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors)
 * \param out result vector Q(p)^{-2} z
 */
void QNTpinv2z(const double* const x, const double* const y, const double* const z,
               const size_t vecSize, const size_t varsCount, double* out);

/** Jordan product of the inverse of x with y: out = x^{-1} \circ y.
 * Uses the formula x^{-1} = R x / det(x), where R is the reflection matrix.
 * \param x input vector to be inverted
 * \param y input vector for the Jordan product
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors)
 * \param out result vector x^{-1} \circ y
 */
void Jxinvprody(const double* const x, const double* const y, const size_t vecSize,
                const size_t varsCount, double* out);

/** Quadratic representation matrix Q(vec) of a vector.
 * \param vec input vector
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors) in vec
 * \return newly allocated matrix Q(vec)
 */
NumericsMatrix* QRmat(const double* const vec, const size_t vecSize, const size_t varsCount);

/** Square root of the determinant of a Jordan-algebra element.
 * \param x input vector
 * \param dimension dimension of each Jordan-algebra block
 * \return long-double value sqrt(det(x))
 */
float_type ld_gammal(const double* const x, const size_t dimension);

/** Nesterov-Todd scaling matrix W(x, z).
 * Computed following the procedure in Tütüncü, Toh and Todd,
 * "Solving semidefinite-linear programs using SDPT3", Math. Prog. 2003, pp. 195-196.
 * \param x first input vector
 * \param z second input vector
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors)
 * \return newly allocated NT scaling matrix
 */
NumericsMatrix* NTmat(const double* const x, const double* const z, const size_t vecSize,
                      const size_t varsCount);

/** Inverse of the Nesterov-Todd scaling matrix W(x, z)^{-1}.
 * \param x first input vector
 * \param z second input vector
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors)
 * \return newly allocated inverse NT scaling matrix
 */
NumericsMatrix* NTmatinv(const double* const x, const double* const z, const size_t vecSize,
                         const size_t varsCount);

/** Square of the Nesterov-Todd scaling matrix W(x, z)^2.
 * \param x first input vector
 * \param z second input vector
 * \param vecSize length of the vector data
 * \param varsCount number of Jordan-algebra variables (subvectors)
 * \return newly allocated squared NT scaling matrix
 */
NumericsMatrix* NTmatsqr(const double* const x, const double* const z, const size_t vecSize,
                         const size_t varsCount);

#endif
