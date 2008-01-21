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

#ifndef SparseBlockMatrix_H
#define SparseBlockMatrix_H

/*!\file SparseBlockMatrix.h
  \brief Structure definition and functions related to SparseBlockStructuredMatrix
  \author Pascal Denoyelle and Franck Perignon
*/

/** Structure to store sparse block matrices with square diagonal blocks
    \param nbblocks         : the total number of non null blocks
    \param **block          : *block contains the double values of one block in Fortran storage (column by column)
    **block is the list of non null blocks
    \param size             : the number of blocks along a row (or column)
    \param *blocksize       : the list of the sizes of diagonal (square) blocks
    \param *RowIndex        : the list of *block row indices (first row = 0)
    \param *ColumnIndex     : the list of *block column indices (first column = 0)
*/
typedef struct
{
  int nbblocks;
  double **block;
  int size;
  int *blocksize;
  int *RowIndex;
  int *ColumnIndex;
} SparseBlockStructuredMatrix;

#ifdef __cplusplus
extern "C" {
#endif

  /** SparseMatrix - vector product y = A*x or y += Ax
      \param[in] size, dim of the vectors x and y
      \param[in] A, the matrix to be multiplied
      \param[in] x, the vector to be multiplied
      \param[in-out] y, the resulting vector
      \param[in] init, = 0 for y += Ax, =1 for y = Ax
  */
  void prod(int size, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init);

  /** Row of a SparseMatrix - vector product y = rowA*x or y += rowA*x, rowA being a row of blocks of A
      \param[in] sizeX, dim of the vector x
      \param[in] sizeY, dim of the vector y
      \param[in] currentRowNumber, number of the required row of blocks
      \param[in] A, the matrix to be multiplied
      \param[in] x, the vector to be multiplied
      \param[in-out] y, the resulting vector
      \param[in] init, = 0 for y += Ax, =1 for y = Ax
  */
  void subRowProd(int sizeX, int sizeY, int currentRowNumber, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init);

  /** Destructor for SparseBlockStructuredMatrix objects
      \param SparseBlockStructuredMatrix, the matrix to be destroyed.
   */
  void freeSpBlMat(SparseBlockStructuredMatrix *);

#ifdef __cplusplus
}
#endif

#endif /* NSSPACK_H */
