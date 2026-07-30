/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
 * limitations under the License.make
 */

#ifndef SparseBlockMatrix_H
#define SparseBlockMatrix_H

#include <stdio.h>  // for size_t, FILE

#include "CSparseMatrix.h"  // for CSparseMatrix
#include "NumericsDataVersion.h"
#include "NumericsFwd.h"  // for SparseBlockStructuredMatrix, SparseBlockC...

/*!\file SparseBlockMatrix.h
  Structure definition and functions related to
  SparseBlockStructuredMatrix

*/

/* Note: the SparseBlockStructuredMatrix format is the same as the one used by Boost C++
   library to store compressed sparse row matrices. The same member
   names have been adopted in order to simplify usage from Siconos
   Kernel : filled1, filled2, index1_data, index2_data.
   Reference :
   http://ublas.sourceforge.net/refdoc/classboost_1_1numeric_1_1ublas_1_1compressed__matrix.html

*/

/** @brief Structure to store sparse block matrices with square diagonal blocks.

  Equivalent to CSR but in block format.

  If we consider the matrix M and the right-hand-side q defined as

   \f[
   M=\left[\begin{array}{cccc|cc|cc}
   1 & 2 & 0 & 4   & 3 &-1   & 0 & 0  \\
   2 & 1 & 0 & 0   & 4 & 1   & 0 & 0  \\
   0 & 0 & 1 &-1   & 0 & 0   & 0 & 0  \\
   5 & 0 &-1 & 6   & 0 & 6   & 0 & 0  \\
   \hline
   0 & 0 & 0 & 0   & 1 & 0   & 0 & 5 \\
   0 & 0 & 0 & 0   & 0 & 2   & 0 & 2 \\
   \hline
   0 & 0 & 2 & 1   & 0 & 0   & 2 & 2 \\
   0 & 0 & 2 & 2   & 0 & 0   & -1 & 2 \\
   \end{array}\right] \quad, q=\left[\begin{array}{c}-1 \\ -1 \\ 0 \\ -1 \\ \hline 1 \\ 0
   \\ \hline -1 \\ 2 \end{array}\right]. \f]

   then

   - the number of non null blocks is 6 (nbblocks=6)
   - the number of rows of blocks is 3 (blocknumber0 =3) and the
   number of columns of blocks is 3 (blocknumber1 =3)
   - the vector blocksize0 is equal to {4,6,8} and the vector
   blocksize1 is equal to {4,6,8}
   - the integer filled1 is equal to 4
   - the integer filled2 is equal to 6
   - the vector index1_data is equal to {0,2,4,6}
   - the vector index2_data is equal to {0,1,1,2,0,2}
   - the block contains all non null block matrices stored in Fortran
   order (column by column) as
   block[0] = {1,2,0,5,2,1,0,0,0,0,1,-1,4,0,-1,6}
   block[1] = {3,4,0,0,-1,1,0,6}
   ...
   block[5] = {2,-1,2,2}
*/
struct SparseBlockStructuredMatrix {
  /** The number of non null blocks */
  size_t nbblocks;

  /**  List of non-null blocks. block[i] = the double values of block i in Fortran storage */
  double** block;

  /** the number of rows of blocks */
  size_t blocknumber0;

  /** the number of columns of blocks */
  size_t blocknumber1;

  /** the vector of cumulated row sizes of blocks. blocksize0[i] = blocksize0[i-1] +
   ni, ni being the number of rows of blocks at block-row i*/
  size_t* blocksize0;

  /** the vector of cumulated column sizes of blocks. blocksize1[i] = blocksize1[i-1] + ni,
   ni being the number of columns of blocks at block-col i */
  size_t* blocksize1;

  /** the index of the last non empty line (of blocks) + 1 */
  size_t filled1;

  /** the size of index2_data that corresponds to the number of non null blocks*/
  size_t filled2;

  /** vector of size "number of lines of blocks" + 1. index1_data[i] = position in
   * index2_data/block of the first non-null block of the line i */
  size_t* index1_data;

  /** vector of size number of blocks. index2_data[j] = column index of block[j] */
  size_t* index2_data;

  /** the indices of the diagonal blocks */
  size_t* diagonal_blocks;

  /** version of storage */
  NumericsDataVersion version;
};

/** @brief Structure to store sparse block matrices with square diagonal blocks.

   Coordinate format.

 */
struct SparseBlockCoordinateMatrix {
  /** number of blocks */
  size_t nbblocks;

  /** number of rows */
  size_t blocknumber0;

  /** number of columns */
  size_t blocknumber1;

  /** block pointers */
  double** block;

  /** cumulative number of rows in blocks */
  size_t* blocksize0;

  /** cumulative number of columns in blocks */
  size_t* blocksize1;

  /** row indices */
  size_t* row;

  /** column indices */
  size_t* column;
};

struct SparseBlockStructuredMatrixPred {
  int nbbldiag;
  int** indic;
  int** indicop;
  double** submatlcp;
  double** submatlcpop;
  int** ipiv;
  int* sizesublcp;
  int* sizesublcpop;
  double** subq;
  double** bufz;
  double** newz;
  double** workspace;
};

struct SBM_index_by_column {
  /* the index of the last non empty line + 1 */
  size_t filled3;
  /* the size of index2_data that corresponds of the number of non null blocks*/
  size_t filled4;

  size_t* index3_data;
  size_t* index4_data;
  size_t* blockMap;
};

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * @brief An enum to define the behavior of SBM_clear and SBM_free
 *
 * When calling SBM_clear everything in the structure is cleared and freed
 * with  possible exceptions for
 *  - block[i] attribute, for which the behavior must be set
 *  - blocksizes1 and 2
 *  - the structure itself (SBM_free only)
 *
 * The current enum is used to set this behavior.
 *
 *  no: not freed. yes: freed.
 * | Flag                          | block[i] | blocksize | sbm | block[] | index1/2 |
 * |-------------------------------|----------|-----------|-----|---------|----------|
 * | SBM_FREE_NONE                 |    no    |     no    |  no |    no   |    no    |
 * | SBM_FREE_BLOCKS               |    yes   |     no    |  no |    no   |    no    |
 * | SBM_FREE_SBM                  |    no    |     no    |  yes|    no   |    no    |
 * | SBM_FREE_SIZES                |    no    |     yes   |  no |    no   |    no    |
 * | SBM_FREE_BLOCK_ARRAY          |    no    |     no    |  no |    yes  |    no    |
 * | SBM_FREE_INDEX                |    no    |     no    |  no |    no   |    yes   |
 * | SBM_FREE_KEEP_BLOCKS          |    no    |     yes   |  yes|    yes  |    yes   |
 * | SBM_FREE_KEEP_BLOCK_ARRAY     |    no    |     yes   |  yes|    no   |    yes   |
 * | SBM_FREE_KEEP_SIZES           |    yes   |     no    |  yes|    yes  |    yes   |
 * | SBM_FREE_KEEP_BLOCK_AND_SIZES |   no     |     no    |  yes|    no   |    yes   |
 * | SBM_FREE_KEEP_BLOCKS_AND_SIZES|   no     |     no    |  yes|    yes  |    yes   |
 * | SBM_FREE_ALL                  |    yes   |     yes   |  yes|    yes  |    yes   |
 *
 *  In all cases:
 *   - all pointers will end up to null
 *   - sbm->_blocks is freed.
 *
 *  Warning FP: combinations (SBM_FREE_KEEP_*) must NOT be combined
 *  with each other or with primitive flags using |.
 *  Only primitive flags (SBM_FREE_BLOCKS ...) can be combined with | else
 *  unexpected results may occur.
 */
typedef enum {
  SBM_FREE_NONE =
      0, /* Do not free sbm.block[i] nor sbm.block, nor sbm, nor blocksizes, nor index*/
  SBM_FREE_BLOCKS = 1,      /**< Free sbm->block[i] */
  SBM_FREE_SBM = 2,         /**< Free the struct sbm itself */
  SBM_FREE_SIZES = 4,       /**< Free blocksize0 and blocksize1 */
  SBM_FREE_BLOCK_ARRAY = 8, /**< Free sbm->block */
  SBM_FREE_INDEX = 16,      /**< Free index1_data and index2_data */
  /* Combinations to ease call of this enum ... */
  SBM_FREE_KEEP_BLOCKS = SBM_FREE_SIZES | SBM_FREE_SBM | SBM_FREE_BLOCK_ARRAY | SBM_FREE_INDEX,
  /**< Free all except block[i] */
  SBM_FREE_KEEP_BLOCK_ARRAY = SBM_FREE_SIZES | SBM_FREE_SBM | SBM_FREE_INDEX,
  /**< Free sizes, sbm, index; keep block[] and block[i] */
  SBM_FREE_KEEP_SIZES = SBM_FREE_BLOCKS | SBM_FREE_SBM | SBM_FREE_BLOCK_ARRAY | SBM_FREE_INDEX,
  /**< Free all except blocksize */
  SBM_FREE_KEEP_BLOCK_AND_SIZES = SBM_FREE_SBM | SBM_FREE_INDEX,
  /**< Free sbm and index; keep block[], block[i] and sizes */
  SBM_FREE_KEEP_BLOCKS_AND_SIZES = SBM_FREE_SBM | SBM_FREE_BLOCK_ARRAY | SBM_FREE_INDEX,
  /**< Free sbm, index and blocl; keep block[i] and sizes */
  SBM_FREE_ALL =
      SBM_FREE_BLOCKS | SBM_FREE_SBM | SBM_FREE_SIZES | SBM_FREE_BLOCK_ARRAY | SBM_FREE_INDEX
  /**< Free everything */
} SBMFreeLevel;

/** Creation of an empty Sparse Block Matrix (all fields to "0")
 *
 *  \return a pointer on allocated and initialized space
 */
SparseBlockStructuredMatrix* SBM_new(void);

/** set all matrix fields to "0"
 *
 *  \param sbm a matrix (must not be nullptr)
 */
void SBM_null(SparseBlockStructuredMatrix* sbm);

/** @brief Free memory allocated for the attribute 'blocks' of a SBM matrix
 *
 *  Warning : the behavior is controlled by the level parameter
 *  @see SBMFreeLevel
 *
 *  @param sbm the matrix to clear
 *  @param level controls what is freed, see comment above.
 */
void SBM_clear_block(SparseBlockStructuredMatrix* sbm, SBMFreeLevel level);

/** @brief Free memory allocated for the SBM matrix
 *
 * Warning : the behavior for the block attribute is controlled by the level parameter
 *  @see SBMFreeLevel
 *
 * @param sbm   the matrix to clear
 * @param level controls what is freed, see comment above.
 */
void SBM_clear(SparseBlockStructuredMatrix* blmat, SBMFreeLevel level);

/** @brief lear and free a  SparseBlockStructuredMatrix
 *
 * Calls SBM_clear then frees sbm itself if level & SBM_FREE_SBM.
 *
 * Warning : the behavior for the block attribute is controlled by the level parameter
 *
 *  @see SBMFreeLevel
 *
 *  @param sbm the SparseBlockStructuredMatrix that must be freed
 *  @param level controls what is freed, see SBMFreeLevel
 *  @return NULL
 */
SparseBlockStructuredMatrix* SBM_free(SparseBlockStructuredMatrix* sbm, SBMFreeLevel level);

/** @brief number of non-null elements in the matrix
 *
 *  @param A input matrix
 */
size_t SBM_nnz(SparseBlockStructuredMatrix* A);

/** SparseMatrix - vector product y = alpha*A*x + beta*y

   \param[in] sizeX dim of the vectors x
   \param[in] sizeY dim of the vectors y
   \param[in] alpha coefficient
   \param[in] A the matrix to be multiplied
   \param[in] x the vector to be multiplied
   \param[in] beta coefficient
   \param[in,out] y the resulting vector
*/
void SBM_gemv(size_t sizeX, size_t sizeY, double alpha,
              const SparseBlockStructuredMatrix* const A, const double* x, double beta,
              double* y);

/**
   SparseMatrix - vector product y = A*x + y for block of size 3x3

   \param[in] sizeX dim of the vectors x
   \param[in] sizeY dim of the vectors y
   \param[in] A the matrix to be multiplied
   \param[in] x the vector to be multiplied
   \param[in,out] y the resulting vector
*/
void SBM_gemv_3x3(size_t sizeX, size_t sizeY, const SparseBlockStructuredMatrix* const A,
                  double* const x, double* y);

/**
   SparseBlockStructuredMatrix - SparseBlockStructuredMatrix product C = alpha*A*B + beta*C
   The routine has to be used with precaution. The allocation of C is not done
   since we want to add beta*C. We assume that the structure and the allocation
   of the matrix C are right. Especially:


   - the blocks C(i,j) must exists
   - the sizes of blocks must be consistent
   - no extra block must be present in C

   \param[in] alpha coefficient
   \param[in] A the matrix to be multiplied
   \param[in] B the matrix to be multiplied
   \param[in] beta coefficient
   \param[in,out] C the resulting matrix
 */
void SBM_gemm_without_allocation(double alpha, const SparseBlockStructuredMatrix* const A,
                                 const SparseBlockStructuredMatrix* const B, double beta,
                                 SparseBlockStructuredMatrix* C);

/**
   SparseBlockStructuredMatrix - SparseBlockStructuredMatrix multiplication C = A * B

   Memory allocation for C.

   \param[in] A the matrix to be multiplied
   \param[in] B the matrix to be multiplied
   \return C the resulting matrix
*/
SparseBlockStructuredMatrix* SBM_multiply(const SparseBlockStructuredMatrix* const A,
                                          const SparseBlockStructuredMatrix* const B);

/** Perform the allocation of a zero matrix that is compatible with multiplication
 */
SparseBlockStructuredMatrix* SBM_zero_matrix_for_multiply(
    const SparseBlockStructuredMatrix* const A, const SparseBlockStructuredMatrix* const B);

/**
   SparseBlockStructuredMatrix - SparseBlockStructuredMatrix addition C = alpha*A + beta*B

   Memory allocation for C.

   \param[in] A the matrix to be added
   \param[in] B the matrix to be added
   \param[in] alpha coefficient
   \param[in] beta coefficient
   \return C the resulting matrix
*/
SparseBlockStructuredMatrix* SBM_add(const SparseBlockStructuredMatrix* A,
                                     const SparseBlockStructuredMatrix* B, double alpha,
                                     double beta);

/**
   SparseBlockStructuredMatrix - SparseBlockStructuredMatrix addition C = alpha*A + beta*B +
   gamma*C without allocation. We assume that C has the correct structure

   \param[in] A the matrix to be added
   \param[in] B the matrix to be added
   \param[in] alpha coefficient
   \param[in] beta coefficient
   \param[in] gamma coefficient
   \param[in,out] C the resulting matrix
*/
void SBM_add_without_allocation(const SparseBlockStructuredMatrix* A,
                                const SparseBlockStructuredMatrix* B, double alpha,
                                double beta, SparseBlockStructuredMatrix* C, double gamma);

/**
   In-place multiply a matrix with a double alpha*A --> A

   \param alpha the  coefficient
   \param A the   matrix
 */
void SBM_scal(double alpha, SparseBlockStructuredMatrix* A);

/**
   Row of a SparseMatrix - vector product y = rowA*x or y += rowA*x, rowA being a row of blocks
   of A

   \param[in] sizeX dim of the vector x
   \param[in] sizeY dim of the vector y
   \param[in] currentRowNumber number of the required row of blocks
   \param[in] A the matrix to be multiplied
   \param[in] x the vector to be multiplied
   \param[in,out] y the resulting vector
   \param[in] init = 0 for y += Ax, =1 for y = Ax
*/
void SBM_row_prod(size_t sizeX, size_t sizeY, size_t currentRowNumber,
                  const SparseBlockStructuredMatrix* const A, const double* const x, double* y,
                  int init);

/**
   Row of a SparseMatrix - vector product y = rowA*x or y += rowA*x, rowA being a row of blocks
   of A

   \param[in] sizeX dim of the vector x
   \param[in] sizeY dim of the vector y
   \param[in] currentRowNumber number of the required row of blocks
   \param[in] A the matrix to be multiplied
   \param[in] x the vector to be multiplied
   \param[in,out] y the resulting vector
   \param[in] init = 0 for y += Ax, =1 for y = Ax
*/
void SBM_row_prod_no_diag(size_t sizeX, size_t sizeY, size_t currentRowNumber,
                          const SparseBlockStructuredMatrix* const A, const double* const x,
                          double* y, int init);

/**
   Row of a SparseMatrix - vector product y = rowA*x or y += rowA*x, rowA being a row of blocks
   of A of size 3x3

   \param[in] sizeX dim of the vector x
   \param[in] sizeY dim of the vector y
   \param[in] currentRowNumber number of the required row of blocks
   \param[in] A the matrix to be multiplied
   \param[in] x the vector to be multiplied
   \param[in,out] y the resulting vector
*/
void SBM_row_prod_no_diag_3x3(size_t sizeX, size_t sizeY, size_t currentRowNumber,
                              const SparseBlockStructuredMatrix* const A, double* const x,
                              double* y);
void SBM_row_prod_no_diag_2x2(size_t sizeX, size_t sizeY, size_t currentRowNumber,
                              const SparseBlockStructuredMatrix* const A, double* const x,
                              double* y);
void SBM_row_prod_no_diag_1x1(size_t sizeX, size_t sizeY, size_t currentRowNumber,
                              const SparseBlockStructuredMatrix* const A, double* const x,
                              double* y);

void SBM_extract_component_3x3(const SparseBlockStructuredMatrix* const A,
                               SparseBlockStructuredMatrix* B, size_t* row_components,
                               size_t row_components_size, size_t* col_components,
                               size_t col_components_size);

/**
    Screen display of the matrix content


    \param m the matrix to be displayed
 */
void SBM_print(const SparseBlockStructuredMatrix* const m);

/**
   print in file  of the matrix content

   \param m the matrix to be displayed
   \param file the corresponding  file
*/
void SBM_write_in_file(const SparseBlockStructuredMatrix* const m, FILE* file);

/**
   read in file  of the matrix content without performing memory allocation

   \param M the matrix to be displayed
   \param file the corresponding name of the file
*/
void SBM_read_from_file(SparseBlockStructuredMatrix* const M, FILE* file);

/**
   Create from file a SparseBlockStructuredMatrix with  memory allocation


   \param file the corresponding name of the file
   \return the matrix to be displayed
 */
SparseBlockStructuredMatrix* SBM_new_from_file(FILE* file);

/**
   print in file  of the matrix content in Scilab format for each block

   \param M the matrix to be displayed
   \param file the corresponding  file
*/
void SBM_write_in_fileForScilab(const SparseBlockStructuredMatrix* const M, FILE* file);

/**
   print in file  of the matrix content

   \param M the matrix to be displayed
   \param filename the corresponding file
*/
void SBM_write_in_filename(const SparseBlockStructuredMatrix* const M, const char* filename);

/**
    read in file  of the matrix content

    \param M the matrix to be displayed
    \param filename the corresponding name of the file
*/
void SBM_read_from_filename(SparseBlockStructuredMatrix* const M, const char* filename);

/**
   Destructor for SparseBlockStructuredMatrixPred objects

   \param blmatpred SparseBlockStructuredMatrix, the matrix to be destroyed.
 */
void SBM_clear_pred(SparseBlockStructuredMatrixPred* blmatpred);

/**
    Compute the indices of blocks of the diagonal block


    \param M the SparseBlockStructuredMatrix matrix
    \return the indices for all the rows
*/
size_t* SBM_compute_diagonal_block_indices(SparseBlockStructuredMatrix* const M);

/**
    Find index of the diagonal block in a row

    \param M the SparseBlockStructuredMatrix matrix
    \param row the row of the required block
    \return pos the position of the block
*/
size_t SBM_diagonal_block_index(SparseBlockStructuredMatrix* const M, size_t row);

/**
    insert an entry into a SparseBlockStructuredMatrix.
    This method is expensive in terms of memory management. For a lot of entries, use
    an alternative

    \param M the SparseBlockStructuredMatrix
    \param i row index
    \param j column index
    \param val the value to be inserted.
    \return true on success, else false
 */
bool SBM_entry(SparseBlockStructuredMatrix* M, size_t row, size_t col, double val);

/**
   get the element of row i and column j of the matrix M

   \param M the SparseBlockStructuredMatrix matrix
   \param row the row index
   \param col the column index
   \return the value
*/
double SBM_get_value(const SparseBlockStructuredMatrix* const M, size_t row, size_t col);

/**
   @brief Full copy of a Sparse Block matrix into another

   @param[in] A the matrix to be copied
   @param[out] B the matrix copy of A
*/
void SBM_copy(const SparseBlockStructuredMatrix* const A, SparseBlockStructuredMatrix* B);

/**
 * @brief Copy of a SBM into another BUT B block data are pointer links to A block data
 *   --> B does not own the blocks
 * Be careful with SBM_free on B matrix.
 * All other attributes of B are reallocated and set.
 * @param[in] A the matrix to be copied
 * @param[out] B the matrix copy of A
 */
void SBM_copy_and_link(const SparseBlockStructuredMatrix* const A,
                       SparseBlockStructuredMatrix* B);

/**
Transpose  by copy of a SBM  A into B

\param[in] A the SparseBlockStructuredMatrix matrix to be copied
\param[out]  B the SparseBlockStructuredMatrix matrix copy of transpose A
\return 0 if ok
*/
int SBM_transpose(const SparseBlockStructuredMatrix* const A, SparseBlockStructuredMatrix* B);

/**
   Inverse (in place) a square diagonal block matrix

   \param[in,out] M the SparseBlockStructuredMatrix matrix to be inversed
   \param ipiv worksapce for storign ipiv
   \return 0 ik ok
*/
int SBM_inverse_diagonal_block_matrix_in_place(const SparseBlockStructuredMatrix* M,
                                               int* ipiv);

/**
   Copy a SBM into a Dense Matrix

   \param[in] A the SparseBlockStructuredMatrix matrix
   \param[in] denseMat pointer on the filled dense Matrix. Must be allocated before the call!
*/
void SBM_to_dense(const SparseBlockStructuredMatrix* const A, double* denseMat);

/**
    Copy a SBM into a Sparse (CSR) Matrix

    Warning: the output matrix must be pre-allocated (use SBM_to_sparse_init_memory)

    \param[in] A the SparseBlockStructuredMatrix matrix
    \param[in,out] outSparseMat pointer on the filled sparse Matrix
    \return 0 if ok
*/
int SBM_to_sparse(const SparseBlockStructuredMatrix* const A, CSparseMatrix* outSparseMat);

/**
   initMemory of a Sparse (CSR) Matrix from a SBM matrix

   \param[in] A the SparseBlockStructuredMatrix matrix
   \param[in] sparseMat pointer on the initialized sparse Matrix
   \return 0 if ok
*/
int SBM_to_sparse_init_memory(const SparseBlockStructuredMatrix* const A,
                              CSparseMatrix* sparseMat);

/**
   Copy a block row of the SBM into a Dense Matrix

   \param[in] A the SparseBlockStructuredMatrix matrix to be inversed.
   \param[in] row the block row index copied.
   \param[in] denseMat pointer on the filled dense Matrix.
   \param[in] rowPos line pos in the dense matrix.
   \param[in] rowNb total number of line of the dense matrix.
   The number of line copied is contained in M.

 */
void SBM_row_to_dense(const SparseBlockStructuredMatrix* const A, size_t row, double* denseMat,
                      size_t rowPos, size_t rowNb);

/**

   \param [in] rowIndex: permutation: the row numC of C is the row rowIndex[numC] of A.
   \param [in] A The source SBM.
   \param [out] C The target SBM. It assumes the structure SBM has been allocated.
   The memory allocation for its menber is done inside.
   Warning: in the end, A and C share their blocks (pointer links). C must be freed
   using SBM_free(..., SBM_FREE_KEEP_BLOCKS).
*/
void SBM_row_permutation(const size_t* const rowIndex,
                         const SparseBlockStructuredMatrix* const A,
                         SparseBlockStructuredMatrix* C);

/**

   \param [in] colIndex: permutation: the col numC of C is the col colIndex[numC] of A.
   \param [in] A The source SBM.
   \param [out] C The target SBM. It assumes the structure SBM has been allocated.
   The memory allocation for its menber is done inside.
   NB : The blocks are not copied.
*/
void SBM_column_permutation(size_t* colIndex, SparseBlockStructuredMatrix* A,
                            SparseBlockStructuredMatrix* C);
/** set all matrix fields to "0"
 *
 *  \param sbm a matrix (must not be nullptr)
 */
void SBCM_null(SparseBlockCoordinateMatrix* MC);

/** @brief Creation of an empty Sparse Block Coordinate Matrix (all fields to "0")
 *
 *  @return a pointer on allocated and initialized space
 */
SparseBlockCoordinateMatrix* SBCM_new(void);

/** @brief Free memory allocated for a Sparse Block Coordinate matrix
 *
 * Warning : the behavior for the block attribute is controlled by the level parameter
 *  @see SBMFreeLevel
 *
 * @param sbm   the matrix to clear
 * @param level controls what is freed, see comment above.
 */
void SBCM_clear(SparseBlockCoordinateMatrix* blmat, SBMFreeLevel level);

/** @brief Free memory allocated for the attribute 'blocks' of a Sparse Block Coordinate
 * matrix
 *
 *  Warning : the behavior is controlled by the level parameter
 *  @see SBMFreeLevel
 *
 *  @param sbm the matrix to clear
 *  @param level controls what is freed, see comment above.
 */
void SBCM_clear_block(SparseBlockCoordinateMatrix* sbm, SBMFreeLevel level);

/** @brief lear and free a  SparseBlockCoordinateMatrix
 *
 * Calls SBCM_clear then frees sbm itself if level & SBM_FREE_SBM.
 *
 * Warning : the behavior for the block attribute is controlled by the level parameter
 *  @see SBMFreeLevel
 *
 *  @param sbm the SparseBlockCoordinateMatrix that must be freed
 *  @param level controls what is freed, see SBMFreeLevel
 *  @return NULL
 */
SparseBlockCoordinateMatrix* SBCM_free(SparseBlockCoordinateMatrix* sbm, SBMFreeLevel level);

/**
   creates a SparseBlockCoordinateMatrix from a list of 3x3
   blocks

   \param[in] m the number of rows
   \param[in] n the number of colums
   \param[in] nbblocks the number of blocks
   \param[in] row a pointer to row of each block
   \param[in] column a pointer to column of each block
   \param[in] block a pointer to each block
   \return a pointer to a SparseBlockCoordinateMatrix structure
*/
SparseBlockCoordinateMatrix* SBCM_new_3x3(size_t m, size_t n, size_t nbblocks, size_t* row,
                                          size_t* column, double* block);

/**
   copy a SparseBlockCoordinateMatrix to a SparseBlockStructuredMatrix

   \param[in] MC the SparseBlockCoordinateMatrix matrix
   \return a pointer to a SparseBlockCoordinateMatrix structure
 */
SparseBlockStructuredMatrix* SBCM_to_SBM(SparseBlockCoordinateMatrix* MC);

/**
   Copy a Sparse Matrix into a SBM, with fixed blocksize

   \param[in] blocksize the blocksize
   \param[in] sparseMat pointer on the Sparse Matrix
   \param[in,out] outSBM pointer on an empty SparseBlockStructuredMatrix
   \return 0 in ok
*/
int SBM_from_csparse(size_t blocksize, const CSparseMatrix* const sparseMat,
                     SparseBlockStructuredMatrix* outSBM);

/**
Copy a Sparse Matrix into a SBM, with fixed blocksize, without big mallocs.

   \param[in] blocksize the blocksize
   \param[in] sparseMat pointer on the Sparse Matrix
   \param[in,out] outSBM pointer on an empty SparseBlockStructuredMatrix
   \return 0 in ok
*/
int SBM_from_csparse_2(size_t blocksize, CSparseMatrix* sparseMat,
                       SparseBlockStructuredMatrix* A);

/**
Copy a dense matrix into a SBM, with fixed blocksize.

   \param[in] blocksize the blocksize
   \param[in] n number of columns
   \param[in] m number of rows
   \param[in] sparseMat pointer on the dense matrix
   \param[in,out] outSBM pointer on an empty SparseBlockStructuredMatrix
   \return 0 in ok
*/
int SBM_from_dense(size_t blocksize, size_t n, size_t m, const double* const denseMat,
                   SparseBlockStructuredMatrix* A);

/** Same as SBM_row_permutation, but copies blocks */
void SBM_row_permutation_copy(size_t* rowIndex, SparseBlockStructuredMatrix* A,
                              SparseBlockStructuredMatrix* C);

/** Same as SBM_row_prod_no_diag_2x2 but allows to choose which column is ignored.
 * Useful when permutation.
 */
void SBM_row_prod_no_diag_2x2_permut(size_t sizeX, size_t sizeY, size_t currentRowNumber,
                                     size_t ignoredCol,
                                     const SparseBlockStructuredMatrix* const A,
                                     double* const x, double* y);

#if defined(__cplusplus)
}
#endif

#endif /* NSSPACK_H */
