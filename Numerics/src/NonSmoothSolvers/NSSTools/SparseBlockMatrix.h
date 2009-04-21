/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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

#include <stddef.h>

/*!\file SparseBlockMatrix.h
  \brief Structure definition and functions related to SparseBlockStructuredMatrix
  \author Pascal Denoyelle and Franck Perignon
*/

/** Structure to store sparse block matrices with square diagonal
    blocks
    \param nbblocks         : the total number of non null blocks
    \param **block : *block contains the double values of one block in
                      Fortran storage (column by column) **block is
    the list of non null blocks

    \param blocknumber0           : the first dimension of the block matrix (number of block rows)
    \param blocknumber1           : the second dimension of the block matrix (number of block columns)
    \param *blocksize0       : the list of sums of the number of rows of the first column of blocks of M: blocksize[i] = blocksize[i-1] + ni,\n
    ni being the number of rows of the  block at  row i
    *blocksize0       : the list of sums of the number of columns of the first row of blocks of M: blocksize[i] = blocksize[i-1] + ni,\n
    ni being the number of columns of the block at  column i
    \param filled1 : number of non empty lines + 1
    \param filled2 : number of non null blocks
    \param index1_data : index1_data is of size filled1 = number of non empty lines + 1. A block with number blockNumber inside a row numbered rowNumber verify : index1_data[rowNumber]<= blockNumber <index1_data[rowNumber+1]
    \param index2_data : index2_data is of size filled2  index2_data[blockNumber] -> columnNumber.


    Related functions: prodSBM(), subRowProdSBM(), freeSBM(), printSBM, getDiagonalBlockPos()
 * If we consider the matrix M and the right-hand-side q defined as
 *
 * \f$
 * M=\left[\begin{array}{cccc|cc|cc}
 *          1 & 2 & 0 & 4   & 3 &-1   & 0 & 0\\
 *          2 & 1 & 0 & 0   & 4 & 1   & 0 & 0\\
 *          0 & 0 & 1 &-1   & 0 & 0   & 0 & 0\\
 *          5 & 0 &-1 & 6   & 0 & 6   & 0 & 0\\
 *          \hline
 *          0 & 0 & 0 & 0   & 1 & 0   & 0 & 5\\
 *          0 & 0 & 0 & 0   & 0 & 2   & 0 & 2\\
 *          \hline
 *          0 & 0 & 2 & 1   & 0 & 0   & 2 & 2\\
 *          0 & 0 & 2 & 2   & 0 & 0   & -1 & 2\\
 *        \end{array}\right] \quad, q=\left[\begin{array}{c}-1\\-1\\0\\-1\\\hline 1\\0\\\hline -1\\2\end{array}\right].
 * \f$
 *
 * then
 * - the number of non null blocks is 6 (nbblocks=6)
 * - the number of rows of blocks is 3 (blocknumber0 =3) and the number of columns of blocks is 3 (blocknumber1 =3)
 * - the vector blocksize0  is equal to {4,6,8} and the vector blocksize1  is equal to {4,6,8}
 * - the integer filled1 is equal to 4
 * - the integer filled2 is equal to 6
 * - the vector index1_data is equal to {0,2,4,6}
 * - the vector index2_data is equal to {0,1,1,2,0,2}
 * - the block contains all non null block matrices stored in Fortran order (column by column) as\n
 *   block[0] = {1,2,0,5,2,1,0,0,0,0,1,-1,4,0,-1,6}\n
 *   block[1] = {3,4,0,0,-1,1,0,6}\n
 *   ...\n
 *   block[5] = {2,-1,2,2}
*/


typedef struct
{
  int nbblocks;
  double **block;
  int blocknumber0;
  int blocknumber1;
  int *blocksize0;
  int *blocksize1;
  size_t filled1;
  size_t filled2;
  size_t *index1_data;
  size_t *index2_data;

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

