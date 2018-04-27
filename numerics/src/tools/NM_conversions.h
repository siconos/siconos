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

#ifndef NM_conversions_H
#define NM_conversions_H

/*!\file NM_conversions.h
  \brief Conversion related functions for the various matrix storages in Numerics
*/

#include "CSparseMatrix.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Convert from csc to triplet (aka coo)
   * \param csc the matrix to convert
   * \return the matrix in triplet format
   */
  CSparseMatrix* NM_csc_to_triplet(CSparseMatrix* csc);

  /** Convert from triplet (aka coo) to csr
   * \param triplet the matrix to convert
   * \return the matrix in csr format
   */
  CSparseMatrix* NM_triplet_to_csr(CSparseMatrix* triplet);

  /** Convert from csr to triplet (aka coo)
   * \param csr the matrix to convert
   * \return the matrix in triplet format
   */
  CSparseMatrix* NM_csr_to_triplet(CSparseMatrix* csr);

  /** Convert from csc to csr
   * \param csc the matrix to convert
   * \return the matrix in csr format
   */
  CSparseMatrix* NM_csc_to_csr(CSparseMatrix* csc);

  /** Convert from csr to csc
   * \param csr the matrix to convert
   * \return the matrix in csc format
   */
  CSparseMatrix* NM_csr_to_csc(CSparseMatrix* csr);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
