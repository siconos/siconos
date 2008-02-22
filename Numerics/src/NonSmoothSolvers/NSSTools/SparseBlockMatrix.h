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
    \param *blocksize       : the list of the sum of dim of diagonal blocks of M: blocksize[i] = blocksize[i-1] + ni,\n
    ni being the size of the diagonal block at row(block) i
    \param *RowIndex        : the list of *block row indices (first row = 0)
    \param *ColumnIndex     : the list of *block column indices (first column = 0)
    Related functions: prodSBM(), subRowProdSBM(), freeSBM(), printSBM, getDiagonalBlockPos()
*/
typedef struct
{
  int nbblocks;
  double **block;
  int size;
  int *blocksize;
  int *RowIndex;
  int *ColumnIndex;
  double* vec;
} SparseBlockStructuredMatrix;

typedef struct
{
  int nbbldiag;
  int **indic;
  int **indicop;
  double **submatlcp;
  double **submatlcpop;
  int **ipiv;
  int *sizesublcp;
  int *sizesublcpop;
  double **subq;
  double **bufz;
  double **newz;
  double **workspace;
} SparseBlockStructuredMatrixPred;

#ifdef __cplusplus
extern "C" {
#endif

  /** SparseMatrix - vector product y = alpha*A*x + beta*y
      \param[in] size, dim of the vectors x and y
      \param[in] alpha coefficient
      \param[in] A, the matrix to be multiplied
      \param[in] x, the vector to be multiplied
      \param[in] beta coefficient
      \param[in-out] y, the resulting vector
  */
  void prodSBM(int size, double alpha, const SparseBlockStructuredMatrix* const A, const double* const x, double beta, double* y);

  /** Row of a SparseMatrix - vector product y = rowA*x or y += rowA*x, rowA being a row of blocks of A
      \param[in] sizeX, dim of the vector x
      \param[in] sizeY, dim of the vector y
      \param[in] currentRowNumber, number of the required row of blocks
      \param[in] A, the matrix to be multiplied
      \param[in] x, the vector to be multiplied
      \param[in-out] y, the resulting vector
      \param[in] init, = 0 for y += Ax, =1 for y = Ax
  */
  void subRowProdSBM(int sizeX, int sizeY, int currentRowNumber, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init);

  /** Row of a SparseMatrix - vector product y = rowA*x or y += rowA*x, rowA being a row of blocks of A
      \param[in] sizeX, dim of the vector x
      \param[in] sizeY, dim of the vector y
      \param[in] currentRowNumber, number of the required row of blocks
      \param[in] A, the matrix to be multiplied
      \param[in] x, the vector to be multiplied
      \param[in-out] y, the resulting vector
      \param[in] init, = 0 for y += Ax, =1 for y = Ax
  */
  void rowProdNoDiagSBM(int sizeX, int sizeY, int currentRowNumber, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init);

  /** Destructor for SparseBlockStructuredMatrix objects
      \param SparseBlockStructuredMatrix, the matrix to be destroyed.
   */
  void freeSBM(SparseBlockStructuredMatrix *);

  /** Screen display of the matrix content
      \param M the matrix to be displayed
   */
  void printSBM(const SparseBlockStructuredMatrix* const M);

  /** Destructor for SparseBlockStructuredMatrixPred objects
      \param SparseBlockStructuredMatrix, the matrix to be destroyed.
   */
  void freeSpBlMatPred(SparseBlockStructuredMatrixPred *blmatpred);

  /** Find index position in blocks of the diagonal block of row num
      \param M the SparseBlockStructuredMatrix matrix
      \param num the row of the required block
      \return pos the position of the block
  */
  int getDiagonalBlockPos(const SparseBlockStructuredMatrix* const M, int num);


#ifdef __cplusplus
}
#endif

#endif /* NSSPACK_H */

