/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
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

#ifndef SparseMatrix_internal_H
#define SparseMatrix_internal_H

/*!\file SparseMatrix_internal.h
 * \brief non-public functions and data structures for CSparseMatrix
 *
 * Include this file from .c/.cpp files if there is a need to access
 * members of CSparseMatrix (treat it as non-opaque type) or call
 * CXSparse functions.
 *
 * To use the correct cs_* functions, privately define CS_LONG
 * according to SICONOS_INT64. User code must do the same, but we do
 * not impose on user code independently using CXSparse by defining
 * CS_LONG in our external headers. Nonetheless CSparseMatrix is
 * typedef'd according to SICONOS_INT64 in SparseMatrix.h, therefore
 * it should throw a type error if there is a mismatch.
 */

#include "SiconosConfig.h"
#ifdef SICONOS_INT64
#define CS_LONG
#endif

/* Siconos does not need "complex" part of CXSparse, so avoid
 * compilation C++-related problems with this flag (complex_t vs
 * std::complex). */
#define NCOMPLEX

#include "cs.h"

#include "SparseMatrix.h"

#endif // SparseMatrix_internal_H
