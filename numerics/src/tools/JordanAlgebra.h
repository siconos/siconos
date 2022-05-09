/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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


#include "NumericsMatrix.h"

typedef long double float_type;
/* typedef double float_type; */

/** Create the Arrow representation matrix from vector.
 * \param vec pointer to the vector data.
 * \param vecSize the length of the vector.
 * \param varsCount the count of variables (subvectors) in vec.
 * \return a pointer to a NumericsMatrix
 */
NumericsMatrix* Arrow_repr(const double* const vec, const unsigned int vecSize, const size_t varsCount);

/**
 * Returns reflection matrix.

            |1  0 ...  0|
        R = |0 -1 ...  0|
            | ........ 0|
            |0  0 ... -1|

 * \param size is the size of rectangular metrix
 * \return reflection matrix.
 */
NumericsMatrix* Reflect_mat(const unsigned int size, NM_types type);

/**
 * Returns quadratic representation of the vector x by formula 2*xx^T - det(x)R
 * \param vec pointer to the vector data.
 * \param vecSize the length of the vector.
 * \param varsCount the count of variables (subvectors) in vec.
 * \return a pointer to a NumericsMatrix
 */
NumericsMatrix* Quad_repr(const double* const vec, const unsigned int vecSize, const size_t varsCount);


void NesterovToddVector(const double* const vec1, const double* const vec2,
                           const unsigned int vecSize, const size_t varsCount, double * out);

/** Create the Arrow representation matrix from vector.
 * \param vecSize the length of the vector.
 * \param varsCount the count of variables (subvectors) in vec.
 * \return a pointer to the identity element of the Jordan algebra
 */
double * JA_iden(const unsigned int vecSize, const size_t varsCount);

/** Jordan product of two vectors.
 * \param vec1 is a vector;
 * \param vec2 is a vector;
 * \param vecSize is the length of the vectors;
 * \param varsCount is the count of variables (subvectors) in x and y.
 * \param out is the result vector of the Jordan product.
 */
void JA_prod(const double * const vec1, const double * const vec2,
             const unsigned int vecSize, const int varsCount, double * out);


/** Returns the eigenvalues of each element in the vector.
 * \param vec is the pointer to the vector data.
 * \param vecSize is the length of the vector.
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the result vector of eigenvalues (out\in\mathbb{R}^{2*varsCount})..
 */
void JA_eigenvals(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out);


/** Returns the eigenvectors of each element in the vector.
 * \param vec is the pointer to the vector data.
 * \param vecSize is the length of the vector.
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the result matrix of eigenvactors (out\in\mathbb{R}^{vecSize\times 2*varsCount}).
 */
void JA_eigenvecs(const double * const vec, const unsigned int vecSize, const size_t varsCount, double ** out);


/** Compute element by element square root
 * \param vec is the vector
 * \param vecSize is the size of the vector vec
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the sqrt vector
 */
void JA_sqrt(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out);


/** Compute element by element inverse of square root
 * \param vec is the vector
 * \param vecSize is the size of the vector vec
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the inverse of sqrt vector
 */
void JA_sqrt_inv(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out);


/** Compute element by element the square of the vector
 * \param vec is the vector
 * \param vecSize is the size of the vector vec
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the square of the vector
 */
void JA_power2(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out);


/** Compute element by element the inverse of the vector
 * \param vec is the vector
 * \param vecSize is the size of the vector vec
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the inverse of the vector
 */
void JA_inv(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out);


/** Compute element by element square root
 * \param vec is the vector
 * \param vecSize is the size of the vector vec
 * \param varsCount is the count of variables (subvectors) in vec.
 * \param out is the vector of determinants
 */
void JA_det(const double * const vec, const unsigned int vecSize, const size_t varsCount, double * out);

/* PA: Return the product Q_sqrt(x)*y */
void Qx05y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out);


/* PA: Return the product Q_inv_sqrt(x)*y */
void Qx50y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out);


/* PA: Jordan algebra, returns inv(x) */
void Jinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out);


/* PA: Return J_sqrt(x) */
void Jsqrt(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out);


/* PA: Return J_sqrtinv(x) */
void Jsqrtinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out);


/* PA: Return the Nesterov-Todd vector */
void Nesterov_Todd_vector(short T, const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * p);


/* PA: Return the Nesterov-Todd vector by means of the second formula */
void Nesterov_Todd_vector_b(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * p);


/* Computation of Qx*y by means of the formula 2*(x'*y)*x - det(x)*R*y */
void Qxy(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * z);


/* Returns the product Q_{p}*z where p is the NT vector related to the pair (x,y) */
void QNTpz(const double * const x, const double * const y,const double * const z, const unsigned int vecSize, const size_t varsCount, double * out);


/* Returns the product Q_{p^{-1}}*z where p is the NT vector related to the pair (x,y) */
void QNTpinvz(const double * const x, const double * const y,const double * const z, const unsigned int vecSize, const size_t varsCount, double * out);


/* returns the Jordan product x^{-1} o y by using the formula x^{-1} = R*x/det(x), where R is the reflection matrix */
void Jxinvprody(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out);


/* returns the quadratic representation of a vector vec */
NumericsMatrix* QRmat(const double* const vec, const unsigned int vecSize, const size_t varsCount);


/* Returns a long double as the square root of determinant of a vector related to the Jordan product */
float_type ld_gammal(const double * const x, const size_t dimension);

float_type dnrm2l(const unsigned int n, const double * x);

/*
   Returns the NT matrix by performing computation as in
   Solving semidefinite-linear programs using SDPT3
   by Tutuncu, Toh and Todd, Math.Prog 2003, pp. 195-196
*/
NumericsMatrix* NTmat(const double* const x, const double* const z, const unsigned int vecSize, const size_t varsCount);


/*
   Returns the inverse of NT matrix by performing computation as in
   Solving semidefinite-linear programs using SDPT3
   by Tutuncu, Toh and Todd, Math.Prog 2003, pp. 195-196
*/
NumericsMatrix* NTmatinv(const double* const x, const double* const z, const unsigned int vecSize, const size_t varsCount);


/*
   Returns the square of NT matrix by performing computation as in
   Solving semidefinite-linear programs using SDPT3
   by Tutuncu, Toh and Todd, Math.Prog 2003, pp. 195-196
*/
NumericsMatrix* NTmatsqr(const double* const x, const double* const z, const unsigned int vecSize, const size_t varsCount);


#endif
