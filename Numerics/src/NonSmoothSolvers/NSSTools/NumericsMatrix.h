/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

#ifndef NumericsMatrix_H
#define NumericsMatrix_H

/*!\file NumericsMatrix.h
  \brief Structure definition and functions related to matrix storage in Numerics
  \author Franck Perignon
*/

#include "SparseBlockMatrix.h"

/** Structure used to handle with matrix in Numerics (interface to double*, SparseBlockStructuredMatrix and so on) \n
    Warning: one and only one storage is allowed and thus only one of the pointers below can be different from NULL
    \param storageType, int that identifies the type of storage (0: double*, 1:SparseBlockStructuredMatrix)
    \param size0, number of rows
    \param size1, number of columns
    \param double*, matrix saved as a double * (if storage=0, else equal to NULL)
    \param SparseBlockStructuredMatrix* (if storageType =1, else equal to NULL)
*/
typedef struct
{
  int storageType;
  int size0;
  int size1;
  double* matrix0;
  SparseBlockStructuredMatrix* matrix1;
} NumericsMatrix;

#ifdef __cplusplus
extern "C" {
#endif
  /** Matrix - vector product y = A*x or y += Ax
      \param[in] sizeX, dim of the vector x
      \param[in] sizeY, dim of the vector y
      \param[in] A, the matrix to be multiplied
      \param[in] x, the vector to be multiplied
      \param[in-out] y, the resulting vector
      \param[in] init, = 0 for y += Ax, =1 for y = Ax
  */
  void prod(int sizeX, int sizeY, const NumericsMatrix* const A, const double* const x, double* y, int init);

  /** Row of a Matrix - vector product y = rowA*x or y += rowA*x, rowA being a submatrix of A (sizeY rows and sizeX columns)
      \param[in] sizeX, dim of the vector x
      \param[in] sizeY, dim of the vector y
      \param[in] currentRowNumber, position of the first row of rowA in A (warning: real row if A is a double*, block-row if A is a SparseBlockStructuredMatrix)
      \param[in] A, the matrix to be multiplied
      \param[in] x, the vector to be multiplied
      \param[in-out] y, the resulting vector
      \param[in] init, = 0 for y += Ax, =1 for y = Ax
  */
  void subRowProd(int sizeX, int sizeY, int currentRowNumber, const NumericsMatrix* const A, const double* const x, double* y, int init);

#ifdef __cplusplus
}
#endif

#endif
