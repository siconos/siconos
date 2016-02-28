/* Siconos-Numerics, Copyright INRIA 2005-2016
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#ifndef NM_conversions_H
#define NM_conversions_H

/*!\file NM_conversions.h
  \brief Conversion related functions for the various matrix storages in Numerics
*/

#include "SiconosConfig.h"
#include "SparseMatrix.h"

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
