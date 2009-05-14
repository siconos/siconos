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

#ifndef NumericsMatrix_H
#define NumericsMatrix_H

/*! \page NumericsMatrixPage Matrix Storage in Numerics

\section NumericsMatrixDef What is a NumericsMatrix?

To store matrix objects, Numerics functions use a structure named NumericsMatrix.\n
This objects handles:
 - a number that identify the type of storage (0: double*, 1: sparse block)
 - the dimensions of the matrix (rows: size0, columns: size1)
 - a list of pointers among which only one is non null and represents the matrix itself.

At the time, two different kind of storage are available: \n
- "classical" column-major storage in a double* (field named matrix0)
- sparse block storage, in a structure of type SparseBlockStructuredMatrix (warning: only for square matrices!!) (field named matrix1)

As an example, consider a matrix A = [aij] of size 8X8:\n
\f$
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
\f$

Let us consider a NumericsMatrix named M used to save the matrix above.\n
The first classical storage results in: \n
 - M.storageType = 0\n
 - M.size0 = 8, M.size1 = 8
 - M.matrix0 = [1 2 0 5 0 0 0 0 2 1 0 0 ...]\n
matrix0 being a double* of size 64.
 - M.matrix1 = NULL

For the second way of storage, SparseBlockStructuredMatrix we have:
 - M.storageType = 1\n
 - M.size0 = 8, M.size1 = 8
 - M.matrix0 = NULL
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

\section NumericsMatrixTools Tools for NumericsMatrix

- product with a vector: function prod()
- product of a sub-part of the matrix with a vector: subRowProd().\n


\section NumericsMatrixBuild How to destroy NumericsMatrix?

In the platform, it's up to the Kernel functions to fill NumericsMatrix. Indeed they are only used in the OneStepNSProblem functions, as
argument of the Numerics driver to solver. In that case NumericsMatrix contains only links to objects created in the OneStepNSProblem. Thus there is no
need to allocate or destroy memory on Numerics' side. \n
However, if required, a function freeNumericsMatrix() is available.


*/




/*!\file NumericsMatrix.h
  \brief Structure definition and functions related to matrix storage in Numerics
  \author Franck Perignon
*/
#include <stdio.h>
#include "SparseBlockMatrix.h"
/** Structure used to handle with matrix in Numerics (interface to double*, SparseBlockStructuredMatrix and so on) \n
    Warning: one and only one storage is allowed and thus only one of the pointers below can be different from NULL
    \param storageType, int that identifies the type of storage (0: double*, 1:SparseBlockStructuredMatrix)
    \param size0, number of rows
    \param size1, number of columns
    \param double*, matrix saved as a double * (if storage=0, else equal to NULL)
    \param SparseBlockStructuredMatrix* (if storageType =1 ,else equal to NULL)\n
    Related functions: prod(), subRowProd(), freeNumericsMatrix(), display()
*/
typedef struct
{
  int storageType;
  int size0;
  int size1;
  double* matrix0;
  SparseBlockStructuredMatrix* matrix1;
} NumericsMatrix;

#include "stdio.h"

#ifdef __cplusplus
extern "C" {
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
      memory has been allocate for the structure in Numerics.
      \param m the matrix to be deleted.
   */
  void freeNumericsMatrix(NumericsMatrix* m);

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
     \param filename the corresponding file
  */

  void printInFile(const NumericsMatrix* const M, FILE* file);

  /** Read in file  of the matrix content
     \param M the matrix to be read
     \param filename the corresponding  file
  */
  void readInFile(NumericsMatrix* const M, FILE *file);

  /** PrintInFileForScilab  of the matrix content
   \param M the matrix to be printed
   \param filename the corresponding file
  */

  void printInFileForScilab(const NumericsMatrix* const M, FILE* file);

  /** Read in file for scilab  of the matrix content
     \param M the matrix to be read
     \param filename the corresponding  file
  */
  void readInFileForScilab(NumericsMatrix* const M, FILE *file);

  /** Screen display raw by raw of the matrix content
      \param M the matrix to be displayed
  */
  void displayRawbyRaw(const NumericsMatrix* const m);

#ifdef __cplusplus
}
#endif

#endif
