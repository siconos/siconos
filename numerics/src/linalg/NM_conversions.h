// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef NM_conversions_H
#define NM_conversions_H

/*!\file NM_conversions.h
  \brief Conversion related functions for the various matrix storages in Numerics
*/
#include "CSparseMatrix.h"  // for CSparseMatrix

#if defined(__cplusplus)
extern "C" {
#endif

/** Convert from csc to triplet (aka coo)
 * \param csc the matrix to convert
 * \return the matrix in triplet format
 */
CSparseMatrix* NM_csc_to_triplet(CSparseMatrix* csc);

/** Convert from csc to half triplet for symmetric matrices
 * \param csc the matrix to convert
 * \return the matrix in triplet format
 */
CSparseMatrix* NM_csc_to_half_triplet(CSparseMatrix* csc);

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

#if defined(__cplusplus)
}
#endif

#endif
