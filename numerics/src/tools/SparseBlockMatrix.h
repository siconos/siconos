/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#ifndef SparseBlockMatrix_H
#define SparseBlockMatrix_H

#include <stddef.h>
#include <stdio.h>

#include "SiconosConfig.h"
#include "SiconosCompat.h"

#include "SparseMatrix.h"

/*!\file SparseBlockMatrix.h
  \brief Structure definition and functions related to
  SparseBlockStructuredMatrix

  \author Pascal Denoyelle and Franck Perignon and Co
*/

/** Structure to store sparse block matrices with square diagonal
    blocks.

    Note: the sparse format is the same as the one used by Boost C++
    library to store compressed sparse row matrices. The same member
    names have been adopted in order to simplify usage from Siconos
    Kernel : filled1, filled2, index1_data, index2_data.
    Reference :
    http://ublas.sourceforge.net/refdoc/classboost_1_1numeric_1_1ublas_1_1compressed__matrix.html

    \param nbblocks         the total number of non null blocks
    \param **block : *block contains the double values of one block in
                      Fortran storage (column by column) **block is
    the list of non null blocks
    \param blocknumber0 the first dimension of the block matrix
    (number of block rows)
    \param blocknumber1 the second dimension of the block matrix
    (number of block columns)
    \param *blocksize0 the list of sums of the number of rows of the
    first column of blocks of M: blocksize0[i] = blocksize0[i-1] +
    ni,\n ni being the number of rows of the block at row i
    *blocksize1 the list of sums of the number of columns of the
    first row of blocks of M: blocksize1[i] = blocksize1[i-1] + ni,\n
    ni being the number of columns of the block at column i
    \param filled1 index of the last non empty line + 1
    \param filled2 number of non null blocks
    \param index1_data index1_data is of size equal to number of non
    empty lines + 1. A block with number blockNumber inside a row
    numbered rowNumber verify index1_data[rowNumber]<= blockNumber
    <index1_data[rowNumber+1]`

    \param index2_data index2_data is of size filled2
    index2_data[blockNumber] -> columnNumber.


    Related functions: prodSBM(), subRowProdSBM(), freeSBM(),
    printSBM, getDiagonalBlockPos()
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
 * - the number of rows of blocks is 3 (blocknumber0 =3) and the
     number of columns of blocks is 3 (blocknumber1 =3)
 * - the vector blocksize0 is equal to {4,6,8} and the vector
     blocksize1 is equal to {4,6,8}
 * - the integer filled1 is equal to 4
 * - the integer filled2 is equal to 6
 * - the vector index1_data is equal to {0,2,4,6}
 * - the vector index2_data is equal to {0,1,1,2,0,2}
 * - the block contains all non null block matrices stored in Fortran
     order (column by column) as\n
 *   block[0] = {1,2,0,5,2,1,0,0,0,0,1,-1,4,0,-1,6}\n
 *   block[1] = {3,4,0,0,-1,1,0,6}\n
 *   ...\n
 *   block[5] = {2,-1,2,2}
*/


typedef struct
{
  unsigned int nbblocks;
  double **block;
  unsigned int blocknumber0;
  unsigned int blocknumber1;
  unsigned int *blocksize0;
  unsigned int *blocksize1;
  size_t filled1;
  size_t filled2;
  size_t *index1_data;
  size_t *index2_data;

} SparseBlockStructuredMatrix;

typedef struct
{
  /** number of blocks */
  unsigned int nbblocks;

  /** number of rows */
  unsigned int blocknumber0;

  /** number of columns */
  unsigned int blocknumber1;

  /** block pointers */
  double **block;

  /** cumulative number of rows in blocks */
  unsigned int *blocksize0;

  /** cumulative number of columns in blocks */

  unsigned int *blocksize1;

  /** row indices */
  unsigned int *row;

  /** column indices */
  unsigned int *column;
} SparseBlockCoordinateMatrix;

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

#define NUMERICS_SBM_FREE_BLOCK 4
#define NUMERICS_SBM_FREE_SBM 8

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Creation of an empty Sparse Block Matrix.
   * \return a pointer on allocated and initialized space
   */
  SparseBlockStructuredMatrix* newSBM(void);

  /** SparseMatrix - vector product y = alpha*A*x + beta*y
      \param[in] sizeX dim of the vectors x
      \param[in] sizeY dim of the vectors y
      \param[in] alpha coefficient
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in] beta coefficient
      \param[in,out] y the resulting vector
  */
  void prodSBM(unsigned int sizeX, unsigned int sizeY,
               double alpha, const SparseBlockStructuredMatrix* const A,
               const double* const x, double beta, double* y);

  /** SparseMatrix - vector product y = A*x + y for block of size 3x3
      \param[in] sizeX dim of the vectors x
      \param[in] sizeY dim of the vectors y
      \param[in] alpha coefficient
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in] beta coefficient
      \param[in,out] y the resulting vector
  */
  void prodSBM3x3(unsigned int sizeX, unsigned int sizeY,
                  const SparseBlockStructuredMatrix* const A,
                  double* const x, double* y);

  /** SparseMatrix - SparseMatrix product C = alpha*A*B + beta*C

     \param[in] alpha coefficient
     \param[in] A the matrix to be multiplied
     \param[in] B the matrix to be multiplied
    \param[in] beta coefficient
     \param[in,out] C the resulting matrix
  */
  void prodSBMSBM(double alpha, const SparseBlockStructuredMatrix* const A,
                  const SparseBlockStructuredMatrix* const B,  double beta, SparseBlockStructuredMatrix*  C);

  /** Allocating Memory and initialization for  SparseMatrix - SparseMatrix product C = alpha*A*B + beta*C
    \param[in] A the matrix to be multiplied
    \param[in] B the matrix to be multiplied
    \param[in,out] C the resulting matrix
  */
  void allocateMemoryForProdSBMSBM(const SparseBlockStructuredMatrix* const A, const SparseBlockStructuredMatrix* const B, SparseBlockStructuredMatrix*  C);

  /** Row of a SparseMatrix - vector product y = rowA*x or y += rowA*x, rowA being a row of blocks of A
      \param[in] sizeX dim of the vector x
      \param[in] sizeY dim of the vector y
      \param[in] currentRowNumber number of the required row of blocks
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in,out] y the resulting vector
      \param[in] init = 0 for y += Ax, =1 for y = Ax
  */
  void subRowProdSBM(unsigned int sizeX, unsigned int sizeY, unsigned int currentRowNumber, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init);

  /** Row of a SparseMatrix - vector product y = rowA*x or y += rowA*x, rowA being a row of blocks of A
      \param[in] sizeX dim of the vector x
      \param[in] sizeY dim of the vector y
      \param[in] currentRowNumber number of the required row of blocks
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in,out] y the resulting vector
      \param[in] init = 0 for y += Ax, =1 for y = Ax
  */
  void rowProdNoDiagSBM(unsigned int sizeX, unsigned int sizeY, unsigned int currentRowNumber, const SparseBlockStructuredMatrix* const A, const double* const x, double* y, int init);

  /** Row of a SparseMatrix - vector product y = rowA*x or y += rowA*x, rowA being a row of blocks of A of size 3x3
      \param[in] sizeX dim of the vector x
      \param[in] sizeY dim of the vector y
      \param[in] currentRowNumber number of the required row of blocks
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in,out] y the resulting vector
      \param[in] init = 0 for y += Ax
  */
  void rowProdNoDiagSBM3x3(unsigned int sizeX, unsigned int sizeY, unsigned int currentRowNumber, const SparseBlockStructuredMatrix* const A, double* const x, double* y);
  void rowProdNoDiagSBM3x3_index_block(unsigned int sizeX, unsigned int sizeY, unsigned int currentRowNumber,
                                     const SparseBlockStructuredMatrix* const A, double* const x, double* y,
                                       unsigned int * index_block, int index_block_size);

  /** Destructor for SparseBlockStructuredMatrix objects
      \param blmat SparseBlockStructuredMatrix the matrix to be destroyed.
   */
  void freeSBM(SparseBlockStructuredMatrix * blmat);

  /** Screen display of the matrix content
      \param m the matrix to be displayed
   */
  void printSBM(const SparseBlockStructuredMatrix* const m);

  /** print in file  of the matrix content
  \param m the matrix to be displayed
  \param file the corresponding  file
  */
  void printInFileSBM(const SparseBlockStructuredMatrix* const m, FILE* file);

  /** read in file  of the matrix content without performing memory allocation
  \param M the matrix to be displayed
  \param file the corresponding name of the file
  */
  void readInFileSBM(SparseBlockStructuredMatrix* const M, FILE *file);

  /** Create from file a SparseBlockStructuredMatrix with  memory allocation
      \param outSBM the matrix to be displayed
      \param file the corresponding name of the file
   */
  void newFromFileSBM(SparseBlockStructuredMatrix* const outSBM, FILE *file);

  /** print in file  of the matrix content in Scilab format for each block
      \param M the matrix to be displayed
      \param file the corresponding  file
  */
  void printInFileSBMForScilab(const SparseBlockStructuredMatrix* const M, FILE* file);


  /** print in file  of the matrix content
   \param M the matrix to be displayed
   \param filename the corresponding file
     */
  void printInFileNameSBM(const SparseBlockStructuredMatrix* const M, const char *filename);

  /** read in file  of the matrix content
  \param M the matrix to be displayed
  \param filename the corresponding name of the file
  */
  void readInFileNameSBM(SparseBlockStructuredMatrix* const M, const char *filename);

  /** Destructor for SparseBlockStructuredMatrixPred objects
   *   \param blmatpred SparseBlockStructuredMatrix, the matrix to be destroyed.
   */
  void freeSpBlMatPred(SparseBlockStructuredMatrixPred *blmatpred);

  /** Find index position in blocks of the diagonal block of row num
      \param M the SparseBlockStructuredMatrix matrix
      \param num the row of the required block
      \return pos the position of the block
  */
  unsigned int getDiagonalBlockPos(const SparseBlockStructuredMatrix* const M, unsigned int num);

  /** get the element of row i and column j of the matrix M
     \param M the SparseBlockStructuredMatrix matrix
     \param row the row index
     \param col the column index
     \return the value
  */
  double getValueSBM(const SparseBlockStructuredMatrix* const M, unsigned int row, unsigned int col);

  /** Copy of a SBM  A into B
    \param[in] A the SparseBlockStructuredMatrix matrix to be copied
    \param[out]  B the SparseBlockStructuredMatrix matrix copy of A
    \param[in] copyBlock if copyBlock then the content of block are copied, else only the pointers are copied.
    \return 0 if ok
  */
  int copySBM(const SparseBlockStructuredMatrix* const A, SparseBlockStructuredMatrix*  B, unsigned int copyBlock);


  /** Transpose  by copy of a SBM  A into B
    \param[in] A the SparseBlockStructuredMatrix matrix to be copied
    \param[out]  B the SparseBlockStructuredMatrix matrix copy of transpose A
    \return 0 if ok
  */
  int transposeSBM(const SparseBlockStructuredMatrix* const A, SparseBlockStructuredMatrix*  B);

  /** Inverse (in place) a square diagonal block matrix
  \param[in,out] M the SparseBlockStructuredMatrix matrix to be inversed
  \return 0 ik ok
  */
  int inverseDiagSBM(const SparseBlockStructuredMatrix*  M);

  /** allocate a SparseBlockCoordinateMatrix from a list of 3x3
   * blocks
   * \param[in] m the number of rows
   * \param[in] n the number of colums
   * \param[in] nbblocks the number of blocks
   * \param[in] row a pointer to row of each block
   * \param[in] column a pointer to column of each block
   * \param[in] block a pointer to each block
   * \return a pointer to a SparseBlockCoordinateMatrix structure
   */
  SparseBlockCoordinateMatrix* newSparseBlockCoordinateMatrix3x3fortran(unsigned int m, unsigned int n,
      unsigned int nbblocks,
      unsigned int *row, unsigned int *column, double *block);


  /** free allocated memory in newSparseBlockCoordinateMatrix functions
   * \param[in] MC matrix pointer */
  void freeSparseBlockCoordinateMatrix3x3fortran(SparseBlockCoordinateMatrix *MC);

  /** copy a SparseBlockCoordinateMatrix to a SparseBlockStructuredMatrix
   * \param[in] MC the SparseBlockCoordinateMatrix matrix
   * \return a pointer to a SparseBlockCoordinateMatrix structure
   */
  SparseBlockStructuredMatrix* SBCMToSBM(SparseBlockCoordinateMatrix* MC);

  /** free a SparseBlockStructuredMatrix created with SBCMToSBM
   * \param[in,out] M a SparseBlockStructuredMatrix to free*/
  void freeSBMFromSBCM(SparseBlockStructuredMatrix* M);

  /** Copy a Sparse Matrix into a SBM, with fixed blocksize
      \param[in] blocksize the blocksize
      \param[in] sparseMat pointer on the Sparse Matrix
      \param[in,out] outSBM pointer on an empty SparseBlockStructuredMatrix
      \return 0 in ok
  */
  int sparseToSBM(int blocksize, const CSparseMatrix* const sparseMat, SparseBlockStructuredMatrix* outSBM);

  /** Copy a SBM into a Dense Matrix
  \param[in] A the SparseBlockStructuredMatrix matrix
  \param[in] denseMat pointer on the filled dense Matrix
  */
  void SBMtoDense(const SparseBlockStructuredMatrix* const A, double *denseMat);

  /** Copy a SBM into a Sparse (CSR) Matrix
  \param[in] A the SparseBlockStructuredMatrix matrix
  \param[in] outSparseMat pointer on the filled sparse Matrix
  \return 0 if ok
  */
  int SBMtoSparse(const SparseBlockStructuredMatrix* const A, CSparseMatrix *outSparseMat);

  /** initMemory of a Sparse (CSR) Matrix form a SBM matrix
  \param[in] A the SparseBlockStructuredMatrix matrix
  \param[in] sparseMat pointer on the initialized sparse Matrix
  \return 0 if ok
  */
  int SBMtoSparseInitMemory(const SparseBlockStructuredMatrix* const A, CSparseMatrix *sparseMat);

  /**Copy a block row of the SBM into a Dense Matrix
  \param[in] A the SparseBlockStructuredMatrix matrix to be inversed.
  \param[in] row the block row index copied.
  \param[in] denseMat pointer on the filled dense Matrix.
  \param[in] rowPos line pos in the dense matrix.
  \param[in] rowNb total number of line of the dense matrix.
  The number of line copied is contained in M.
  \
   */
  void SBMRowToDense(const SparseBlockStructuredMatrix* const A, int row, double *denseMat, int rowPos, int rowNb);


  /** To free a SBM matrix (for example allocated by newFromFile).
   * \param[in] A the SparseBlockStructuredMatrix that mus be de-allocated.
   * \param[in] level use NUMERICS_SBM_FREE_BLOCK | NUMERICS_SBM_FREE_SBM
   */
  void SBMfree(SparseBlockStructuredMatrix* A, unsigned int level);

  /*
   * \param [in] rowIndex: permutation: the row numC of C is the row rowIndex[numC] of A.
   * \param [in] A The source SBM.
   * \param [out] C The target SBM. It assumes the structure SBM has been allocated.
   * The memory allocation for its menber is done inside.
   * NB : The blocks are not copied.
   */
  void RowPermutationSBM(unsigned int *rowIndex, SparseBlockStructuredMatrix* A, SparseBlockStructuredMatrix*  C);

  /*
  * \param [in] colIndex: permutation: the col numC of C is the col colIndex[numC] of A.
  * \param [in] A The source SBM.
  * \param [out] C The target SBM. It assumes the structure SBM has been allocated.
  * The memory allocation for its menber is done inside.
  * NB : The blocks are not copied.
  */
  void ColPermutationSBM(unsigned int *colIndex, SparseBlockStructuredMatrix* A, SparseBlockStructuredMatrix*  C);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif /* NSSPACK_H */

