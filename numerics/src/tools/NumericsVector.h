/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#ifndef NumericsVector_H
#define NumericsVector_H

/*!\file NumericsVector.h
  \brief Structure definition and functions related to vector storage in Numerics
*/
#ifndef __cplusplus
#include <stdbool.h>  // for bool
#endif
#include <stdio.h>    // for FILE
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep
#include "NumericsMatrix.h"
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Screen display of the vector content stored as a double * array
      \param m the vector to be displayed
      \param n the size of the vector
   */
  void NV_display(const double * const m, int n);

  void NV_copy(const double * const vec, unsigned int vecSize, double * out);

  void NV_write_in_file_python(double * m,  int nRow, FILE* file);

  /** Test if two vectors are equal up to a given tolerance
      \param x the vector to be tested
      \param y the vector to be tested
      \param n the size of the vector
      \param tol the tolerance
      \return 1 is equal
   */
  bool NV_equal(double * x, double * y, int nRow, double tol);

  /** Insert the vector y into the vector x starting from position i
      \param x the vector
      \param size_x the size of the vector x
      \param y the vector to be inserted
      \param size_y the size of the vector y
      \param i the position of insertion
   */
  void NV_insert(double * x, const unsigned int xSize,
                 const double * const y, const unsigned int ySize,
                 unsigned int i);

  /** Put all elements of vector to the square (element by element)
      \param vec is the vector
      \param vecSize the size of the vector vec
      \param out is the result power2 vector
   */
  void NV_power2(const double * const vec, const unsigned int vecSize, double * out);

  /** Sum all the elements in vector
      \param vec is the vector
      \param vecSize the size of the vector vec
      \return sum of all the elements in vector
   */
  double NV_reduce(const double * const vec, const unsigned int vecSize);

  /** Element by element product of two vectors
      \param vec1 is the vector
      \param vec2 is the vector
      \param vecSize the size of the vector vec
      \param out is the result product vector
   */
  void NV_prod(const double * const vec1, const double * const vec2, const unsigned int vecSize, double * out);

  /** Element by element division of two vectors
      \param x is the vector
      \param y is the vector
      \param vecSize the size of the vector vec
      \return product vector
   */
  double* NV_div(const double * const x, const double * const y, const unsigned int vecSize);

  /** Find a minimum value of vertor
      \param vec is the vector
      \param vecSize the size of the vector vec
      \return a minimum value of vertor
   */
  double NV_min(const double * const vec, const unsigned int vecSize);

  /** Find a maximum value of vertor
      \param vec is the vector
      \param vecSize the size of the vector vec
      \return a minimum value of vertor
   */
  double NV_max(const double * const vec, const unsigned int vecSize);

  /** Compute abs vector
      \param vec is the vector
      \param vecSize the size of the vector vec
      \return elemet by element abs vector
   */
  double * NV_abs(const double * const vec, const unsigned int vecSize);

  /** Compute element by element by element sum
      \param vec1 is the vector
      \param vec2 is the vector
      \param vecSize the size of the vector vec
      \param out is the sum vector
   */
  void NV_add(const double * const vec1, const double * const vec2, const unsigned int vecSize, double * out);

  /** Compute y = alpha * x + beta
      \param x is the vector
      \param vecSize the size of the vector vec
      \param alpha is a scalar
      \param beta is a scalar
      \param out is y = alpha * x + beta
   */

  void NV_const_add(const double * const vec, const unsigned int vecSize, const double alpha, const double beta, double * out);

  /** Compute element by element by element subtraction
      \param vec1 is the vector
      \param vec2 is the vector
      \param vecSize the size of the vector vec
      \return subtract vector
   */

  void NV_sub(const double * const vec1, const double * const vec2, const unsigned int vecSize, double * out);

  /** Find a L-inf norm of vertor ( max(abs(vec)) )
      \param vec is the vector
      \param vecSize the size of the vector vec
      \return a minimum value of vertor
   */
  double NV_norm_inf(const double * const vec, const unsigned int vecSize);

  /** Find a L-2 norm of vertor ( sqrt(sum(vec^2)) )
      \param vec is the vector
      \param vecSize the size of the vector vec
      \return a minimum value of vertor
   */
  double NV_norm_2(const double * const vec, const unsigned int vecSize);

  /** Compute element by element square root
      \param vec is the vector
      \param vecSize is the size of the vector vec
      \param out is the sqrt vector
   */
  void NV_sqrt(const double * const vec, const unsigned int vecSize, double * out);

  /**
   * Compute scalar product x x^T which is a matrix of rank 1.
   * \param vec1 is the vector
   * \param vec2 is the vector
   * \param vecSize is the size of the vector vec
   * \param out is the resut matrix of rank 1.
   */
  void NV_dott(const double * const vec1, const double * const vec2, const unsigned int vecSize, NumericsMatrix* out);
  
  int NV_isnan(const double * const vec,  const unsigned int vecSize );
  
  
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
