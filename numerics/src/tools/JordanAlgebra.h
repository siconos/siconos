/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/** Create the Arrow representation matrix from vector.
 * \param vecSize the length of the vector.
 * \param vec pointer to the vector data.
 * \param varsCount the count of variables (subvectors) in vec.
 * \return a pointer to a NumericsMatrix
 */
RawNumericsMatrix* Arrow_repr(const unsigned int vecSize, const double* const vec, const size_t varsCount);


#endif
