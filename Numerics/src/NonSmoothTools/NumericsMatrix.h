/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#ifndef NumericsMatrix_H
#define NumericsMatrix_H

/*! \page NumericsMatrixPage Matrix Storage in Numerics

\section NumericsMatrixDef What is a NumericsMatrix?

To store matrix objects, Numerics functions use a structure named NumericsMatrix.\n
This objects handles:
 - a number that identify the type of storage (0: double*, 1: sparse block, 2: sparse triplet, 3: sparse compressed columns, 4: sparse compressed transpose )
 - the dimensions of the matrix (rows: size0, columns: size1)
 - a list of pointers among which only one is non NULL and represents the matrix itself.

At the time, the following storage are available: \n
- "classical" column-major storage in a double* (field named matrix0)
- sparse block storage, in a structure of type SparseBlockStructuredMatrix (warning: only for square matrices!!) (field named matrix1)
- triplet storage, to be used with cs_sparse  (field named matrix2)
- compressed columns storage, to be used with cs_sparse  (field named matrix3)
- compressed transposed storage  (field named matrix4)


As an example, consider the following matrix A of size 8X8:\n
\f{equation*}
  A=\left[\begin{array}{cccc|cc|cc}
           1 & 2 & 0 & 4   & 3 &-1   & 0 & 0\\
           2 & 1 & 0 & 0   & 4 & 1   & 0 & 0\\
           0 & 0 & 1 &-1   & 0 & 0   & 0 & 0\\
           5 & 0 &-1 & 6   & 0 & 6   & 0 & 0\\
           \hline
           0 & 0 & 0 & 0   & 1 & 0   & 0 & 5\\
           0 & 0 & 0 & 0   & 0 & 2   & 0 & 2\\
           \hline
           0 & 0 & 2 & 1   & 0 & 0   & 2 & 2\\
           0 & 0 & 2 & 2   & 0 & 0   & -1& 2\\
 \end{array}\right]
\f}

How can we store this matrix ?\n
The first classical storage results in: \n
 - M.storageType = 0
 - M.size0 = 8, M.size1 = 8
 - M.matrix0 = [1 2 0 5 0 0 0 0 2 1 0 0 ...]\n
matrix0 being a double* of size 64.

For the second way of storage, SparseBlockStructuredMatrix we have:
 - M.storageType = 1
 - M.size0 = 8, M.size1 = 8
 - M.matrix1 a SparseBlockStructuredMatrix in which we save:
  - the number of non null blocks, 6 (matrix1->nbblocks) and the number of diagonal blocks, 3 (matrix1->size).
  - the vector matrix1->blocksize which collects the sum of diagonal blocks sizes until the present one, is equal to [4,6,8], \n
  blocksize[i] = blocksize[i-1] + ni,  ni being the size of the diagonal block at row(block) i. \n
  Note that the last element of blocksize corresponds to the real size of the matrix.
  - the list of positions of non null blocks in vectors matrix1->ColumnIndex and matrix1->RowIndex, equal to [0,1,1,2,0,2] and [0,0,1,1,2,2]
  - the list of non null blocks, in matrix1->block, stored in Fortran order (column-major) as\n
    matrix1->block[0] = [1,2,0,5,2,1,0,0,0,0,1,-1,4,0,-1,6]\n
    matrix1->block[1] = [3,4,0,0,-1,1,0,6]\n
    ...\n
    matrix1->block[5] = [2,-1,2,2]

\section NumericsMatrixTools Functions on NumericsMatrix

\subsection NMAlloc Create and delete NumericsMatrix

NumericsMatrix of any kind can be created using either createNumericsMatrix(), which allocates memory or createNumericsMatrixFromData() or fillNumericsMatrix(), which just fills a structure.
These last two functions accept a <i>data</i> parameter, which if non-NULL contains the matrix data.
The function freeNumericsMatrix() is used to clean the fields properly.

\subsection NM_LA Linear Algebra

The following linear algebra operation are supported:

  - product matrix - vector: prodNumericsMatrix()
  - product matrix - matrix: prodNumericsMatrixNumericsMatrix()
  - partial product matrix - vector: subRowProd()

\subsection NM_IO Input / Output

  - display(): display a NumericsMatrix
  - displayRowbyRow(): display a NumericsMatrix row by row
  - printInFileName(), printInFile(): save to filesystem
  - readInFileName(), readInFile(): fill a NumericsMatrix from a file
  - newFromFile(): create new NumericsMatrix from a file

*/




/*!\file NumericsMatrix.h
  \brief Structure definition and functions related to matrix storage in Numerics
  \author Franck Perignon
*/
#include <stdio.h>
#include "NumericsConfig.h"
#include "SparseMatrix.h"
#include "SparseBlockMatrix.h"

/** \struct NumericsMatrix NumericsMatrix.h
 * Structure used to handle with matrix in Numerics (interface to double*, SparseBlockStructuredMatrix and so on) \n
    Warning: one and only one storage is allowed and thus only one of the pointers below can be different from NULL
    Related functions: prod(), subRowProd(), freeNumericsMatrix(), display()
*/
typedef struct
{
  int storageType; /**< the type of storage (0: dense (double*), 1: SparseBlockStructuredMatrix, 2: triplet, 3: compressed col, 4: compressed transpose) */
  int size0; /**< number of rows */
  int size1; /**< number of columns */
  double* matrix0; /**< dense storage */
  SparseBlockStructuredMatrix* matrix1; /**< sparse block storage */
  SparseMatrix* matrix2; /**< triplet storage */
  SparseMatrix* matrix3; /**< compressed column storage */
  SparseMatrix* matrix4; /**< compressed transpose storage */
} NumericsMatrix;

/** Possible types of matrices for NumericsMatrix */
enum NumericsMatrix_types {
  NUMERICS_MATRIX_DENSE,                   /**< dense format */
  NUMERICS_MATRIX_SPARSE_BLOCK,            /**< sparse block format */
  NUMERICS_MATRIX_SPARSE_TRIPLET,          /**< triplet format */
  NUMERICS_MATRIX_SPARSE_COMPRESSED_COL,   /**< compressed columns format */
  NUMERICS_MATRIX_SPARSE_COMPRESSED_TRANS  /**< compressed transpose format */
};

#include "stdio.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** Matrix - vector product y = alpha*A*x + beta*y
      \param[in] sizeX dim of the vector x
      \param[in] sizeY dim of the vector y
      \param[in] alpha coefficient
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in] beta coefficient
      \param[in,out] y the resulting vector
  */
  void prodNumericsMatrix(int sizeX, int sizeY, double alpha, const NumericsMatrix* const A, const double* const x, double beta, double* y);

  /** Matrix - Matrix product C = alpha*A*B + beta*B
      \param[in] alpha coefficient
      \param[in] A the matrix to be multiplied
      \param[in] B the matrix to be multiplied
      \param[in] beta coefficient
      \param[in,out] C the resulting matrix
  */
  void prodNumericsMatrixNumericsMatrix(double alpha, const NumericsMatrix* const A, const NumericsMatrix* const B, double beta, NumericsMatrix* C);



  /** Row of a Matrix - vector product y = rowA*x or y += rowA*x, rowA being a submatrix of A (sizeY rows and sizeX columns)
      \param[in] sizeX dim of the vector x
      \param[in] sizeY dim of the vector y
      \param[in] currentRowNumber position of the first row of rowA in A (warning: real row if A is a double*, block-row if A is a SparseBlockStructuredMatrix)
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in,out] y the resulting vector
      \param[in] init = 0 for y += Ax, =1 for y = Ax
  */
  void subRowProd(int sizeX, int sizeY, int currentRowNumber, const NumericsMatrix* const A, const double* const x, double* y, int init);

  /** Row of a Matrix - vector product y = rowA*x or y += rowA*x, rowA being a submatrix of A (sizeY rows and sizeX columns)
      \param[in] sizeX dim of the vector x
      \param[in] sizeY dim of the vector y
      \param[in] currentRowNumber position of the first row of rowA in A (warning: real row if A is a double*, block-row if A is a SparseBlockStructuredMatrix)
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in,out] y the resulting vector
      \param[in] init = 0 for y += Ax, =1 for y = Ax
  */
  void rowProdNoDiag(int sizeX, int sizeY, int currentRowNumber, const NumericsMatrix* const A, const double* const x, double* y, int init);

  /** Free memory for a NumericsMatrix. Warning: call this function only if you are sure that
      memory has been allocated for the structure in Numerics. This function is assumed that the memory is "owned" by this structure.
      \param m the matrix to be deleted.
   */
  void freeNumericsMatrix(NumericsMatrix* m);

  /** Screen display of the matrix content stored as a double * array in Fortran style
      \param m the matrix to be displayed
      \param nRow the number of rows
      \param nCol the number of columns
      \param lDim the leading dimension of M
   */
  void displayMat(double * m, int nRow, int nCol, int lDim);


  /** Screen display of the matrix content
      \param M the matrix to be displayed
   */
  void display(const NumericsMatrix* const M);

  /** PrintInFile  of the matrix content
     \param M the matrix to be printed
     \param filename the corresponding name of the file
  */
  void printInFileName(const NumericsMatrix* const M, const char *filename);

  /** Read in file  of the matrix content
     \param M the matrix to be read
     \param filename the corresponding name of the file
  */
  void readInFileName(NumericsMatrix* const M, const char *filename);

  /** PrintInFile  of the matrix content
     \param M the matrix to be printed
     \param file filename the corresponding file
  */

  void printInFile(const NumericsMatrix* const M, FILE* file);

  /** Read in file  of the matrix content without performing memory allocation
     \param M the matrix to be read
     \param file the corresponding  file
  */
  void readInFile(NumericsMatrix* const M, FILE *file);

  /** Create from file a NumericsMatrix with  memory allocation
     \param M the matrix to be read
     \param file the corresponding  file
  */
  void newFromFile(NumericsMatrix* const M, FILE *file);

  /** PrintInFileForScilab  of the matrix content
   \param M the matrix to be printed
   \param file the corresponding file
  */
  void printInFileForScilab(const NumericsMatrix* const M, FILE* file);

  /** Read in file for scilab  of the matrix content
     \param M the matrix to be read
     \param file the corresponding  file
  */
  void readInFileForScilab(NumericsMatrix* const M, FILE *file);

  /** Screen display raw by raw of the matrix content
      \param m the matrix to be displayed
  */
  void displayRowbyRow(const NumericsMatrix* const m);

  /** get the diagonal block of a NumericsMatrix. No allocation is done.
   * \param[in] m a NumericsMatrix
   * \param[in] numBlockRow the number of the block Row. Useful only in sparse case
   * \param[in] numRow the starting row. Useful only in dense case.
   * \param[in] size of the diag block.Useful only in dense case.
   * \param[out] Bout the target. In the dense case (*Bout) must be allocated by caller.
   *   In case of sparse case **Bout contains the resulting block (from the SBM).
   */
  void getDiagonalBlock(NumericsMatrix* m,
                        int numBlockRow,
                        int numRow,
                        int size,
                        double **Bout);


  /** create a NumericsMatrix and allocate the memory according to the matrix type
   * \param storageType the type of storage
   * \param size0 number of rows
   * \param size1 number of columns
   * \return a pointer to a NumericsMatrix
   */
  NumericsMatrix* createNumericsMatrix(int storageType, int size0, int size1);

  /** create a NumericsMatrix and possibly set the data
   * \param storageType the type of storage
   * \param size0 number of rows
   * \param size1 number of columns
   * \param data pointer to the matrix data. If NULL, all matrixX fields are
   * set to NULL
   * \return a pointer to a NumericsMatrix
   */
  NumericsMatrix* createNumericsMatrixFromData(int storageType, int size0, int size1, void* data);

  /** fill an existing NumericsMatrix struct
   * \param[in,out] M the struct to fill
   * \param storageType the type of storage
   * \param size0 number of rows
   * \param size1 number of columns
   * \param data pointer to the matrix data. If NULL, all matrixX fields are
   * set to NULL
   */
  void fillNumericsMatrix(NumericsMatrix* M, int storageType, int size0, int size1, void* data);

  /** new NumericsMatrix with sparse storage from minimal set of data
   * \param[in] size0 number of rows
   * \param[in] size1 number of columns
   * \param[in] m1 the SparseBlockStructuredMatrix
   * \return  a pointer to a NumericsMatrix
   */
  NumericsMatrix* newSparseNumericsMatrix(int size0, int size1, SparseBlockStructuredMatrix* m1);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
