/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#ifndef NumericsMatrix_H
#define NumericsMatrix_H

/*!\file NumericsMatrix.h
  \brief Structure definition and functions related to matrix storage in Numerics
*/

#include <assert.h>         // for assert
#include <stdio.h>          // for size_t, FILE, NULL
#include <stdlib.h>         // for malloc
#include "CSparseMatrix.h"  // for CS_INT, CSparseMatrix
#include "NumericsFwd.h"    // for NumericsMatrix, NumericsSparseMatrix, Spa...
#include "NumericsDataVersion.h" // Versioning
#include "NumericsSparseMatrix.h" // for NSM_linear_solver typedef
#include "SiconosConfig.h" // for BUILD_AS_CPP, SICONOS_HAS_MP // IWYU pragma: keep
#include "NM_MPI.h"
#ifndef __cplusplus
#include <stdbool.h>        // for bool
#endif

#ifdef WITH_OPENSSL
#include <openssl/sha.h>
#endif

/** \struct NumericsMatrixInternalData NumericsMatrix.h
 * Structure for simple workspaces
 */
typedef struct
{
  size_t iWorkSize; /**< size of iWork */
  void *iWork; /**< integer workspace */
  size_t sizeof_elt ; /**< sizeof_elt of an element in bytes (result of sizeof for instance)*/
  size_t dWorkSize; /**< size of dWork */
  double *dWork; /**< double workspace */
  bool isLUfactorized; /**<  true if the matrix has already been LU-factorized */
  bool isCholeskyfactorized; /**<  true if the matrix has already been Cholesky factorized */
  bool isLDLTfactorized; /**<  true if the matrix has already been LDLT factorized */
  bool isInversed; /**<  true if the matrix contains its inverse (in place inversion) */
#ifdef SICONOS_HAS_MPI
  MPI_Comm mpi_comm; /**< optional mpi communicator */
#endif
#ifdef WITH_OPENSSL
  unsigned int values_sha1_count; /**< counter for sha1 */
  unsigned char values_sha1[SHA_DIGEST_LENGTH]; /**< sha1 hash of
                                                 * values. Matrices of
                                                 * differents sizes may have
                                                 * the same hash */
#endif
} NumericsMatrixInternalData;

/*! Available types of storage for NumericsMatrix */
typedef enum NumericsMatrix_types {
  NM_DENSE,        /**< dense format */
  NM_SPARSE_BLOCK, /**< sparse block format */
  NM_SPARSE,          /**< compressed column format */
  NM_UNKNOWN, /** unset. Used in NM_null */
} NM_types;

/** \struct NumericsMatrix NumericsMatrix.h
    Interface to different type of matrices in numerics component.

    See NM_* functions for linear algebra operations on dense, sparse block and sparse storage.
*/
struct NumericsMatrix
{
  NM_types storageType; /**< the type of storage:
                      0: dense (double*),
                      1: SparseBlockStructuredMatrix,
                      2: classical sparse (csc, csr or triplet) via CSparse (from T. Davis)*/
  int size0; /**< number of rows */
  int size1; /**< number of columns */
  double* matrix0; /**< dense storage */
  SparseBlockStructuredMatrix* matrix1; /**< sparse block storage */
  NumericsSparseMatrix* matrix2; /**< csc, csr or triplet storage */

  NumericsMatrixInternalData* internalData; /**< internal storage, used for workspace among other things */

  NumericsDataVersion version; /*< version of dense storage */

  NumericsMatrix* destructible; /**< pointer on the destructible
                                 * matrix, by default points toward
                                 * the matrix itself */
};

typedef struct
{
  int size0;
  int size1;
  NumericsMatrix* D1;
  NumericsMatrix* D2;
  NumericsMatrix* A;
} BalancingMatrices;

/*! RawNumericsMatrix is used without conversion in python */
typedef NumericsMatrix RawNumericsMatrix;


typedef enum {
  NM_NONE,          /**< keep nothing */
  NM_KEEP_FACTORS,  /**< keep all the factorization data (useful to reuse the factorization) */
  NM_PRESERVE       /**< keep the matrix as-is (useful for the dense case) */
} NM_gesv_opts;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /**************************************************/
  /** Constructors and destructors   ****************/
  /**************************************************/

  /** Creation of an empty NumericsMatrix.
   * \return a pointer to allocated space
   */
  RawNumericsMatrix* NM_new(void);
  RawNumericsMatrix* NM_eye(int size);
  RawNumericsMatrix* NM_scalar(int size, double s);

  /** create a NumericsMatrix and allocate the memory according to the matrix type
   * \param storageType the type of storage
   * \param size0 number of rows
   * \param size1 number of columns
   * \return a pointer to a NumericsMatrix
   */
  RawNumericsMatrix* NM_create(NM_types storageType, int size0, int size1);

  /** create a NumericsMatrix and possibly set the data
   * \param storageType the type of storage
   * \param size0 number of rows
   * \param size1 number of columns
   * \param data pointer to the matrix data. If NULL, all matrixX fields are
   * set to NULL
   * \return a pointer to a NumericsMatrix
   */
  RawNumericsMatrix* NM_create_from_data(int storageType, int size0, int size1, void* data);
  RawNumericsMatrix* NM_create_from_filename(const char *filename);
  RawNumericsMatrix* NM_create_from_file(FILE *file);

  /** Copy NumericsMatrix version.
   * \param[in] A a NumericsMatrix
   * \param[in,out] B a NumericsMatrix
   */
  void NM_version_copy(const NumericsMatrix* const A, NumericsMatrix* B);

  /** Copy a NumericsMatrix inside another NumericsMatrix (deep).
   *  Reallocations are performed if B cannot hold a copy of A
   * \param[in] A a NumericsMatrix
   * \param[in,out] B a NumericsMatrix
   */
  void NM_copy(const NumericsMatrix* const A, NumericsMatrix* B);

  /** Copy a NumericsMatrix to s sparse one.
   *  Allocation or reallocation are performed on B
   *  \warning It is assumed that B has been properly initialized: its storageType must
   *  be set to NM_SPARSE.
   * \param[in] A a NumericsMatrix
   * \param[in,out] B a NumericsMatrix
   * \param threshold if the original matrix is dense, a threshold can be applied
   * on the absolute value of the entries
   */
  void NM_copy_to_sparse(const NumericsMatrix* const A, NumericsMatrix* B, double threshold);

  /** create a NumericsMatrix similar to the another one. The structure is the same
   * \param mat the model matrix
   * \return a pointer to a NumericsMatrix
   */
  RawNumericsMatrix* NM_duplicate(NumericsMatrix* mat);


  /** Creation, if needed, of sparse matrix storage.
   * \param[in,out] A a NumericsMatrix
   * \return a pointer on the sparse matrix storage
   */
  NumericsSparseMatrix* numericsSparseMatrix(NumericsMatrix* A);

  /** Creation, if needed, of triplet storage from sparse block storage.
   * \param[in,out] A a NumericsMatrix initialized with sparsed block storage.
   * \return the triplet sparse Matrix created in A.
   */
  CSparseMatrix* NM_triplet(NumericsMatrix* A);

   /** Creation, if needed, of half triplet storage from sparse block storage.
   * \param[in,out] A a NumericsMatrix initialized with sparsed block storage.
   * \return the triplet sparse Matrix created in A.
   */
  CSparseMatrix* NM_half_triplet(NumericsMatrix* A);

  /** Creation, if needed, of compress column storage of a NumericsMatrix.
   * \param[in,out] A a NumericsMatrix with sparse block storage initialized
   * \return the compressed column CSparseMatrix created in A.
   */
  CSparseMatrix* NM_csc(NumericsMatrix *A);

  /** Creation, if needed, of the transposed compress column storage
   * from compress column storage.
   * \param[in,out] A a NumericsMatrix with sparse block storage.
   * \return the transposed compressed column matrix created in A.
   */
  CSparseMatrix* NM_csc_trans(NumericsMatrix* A);

  /** Creation, if needed, of compress row storage of a NumericsMatrix
   * \warning This rely on the MKL
   * \param[in,out] A a NumericsMatrix with sparse block storage initialized
   * \return the compressed row CSparseMatrix created in A.
   */
  CSparseMatrix* NM_csr(NumericsMatrix *A);

  /** fill an existing NumericsMatrix struct
   * \param[in,out] M the struct to fill
   * \param storageType the type of storage
   * \param size0 number of rows
   * \param size1 number of columns
   * \param data pointer to the matrix data. If NULL, all matrixX fields are
   * set to NULL
   */
  void NM_fill(NumericsMatrix* M, NM_types storageType, int size0, int size1, void* data);

  /** new NumericsMatrix with sparse storage from minimal set of data
   * \param[in] size0 number of rows
   * \param[in] size1 number of columns
   * \param[in] m1 the SparseBlockStructuredMatrix
   * \return  a pointer to a NumericsMatrix
   */
  RawNumericsMatrix* NM_new_SBM(int size0, int size1, SparseBlockStructuredMatrix* m1);

  /** new NumericsMatrix equal to the transpose of a given matrix
   * \param[in] A
   * \return  a pointer to a NumericsMatrix
   */
  RawNumericsMatrix* NM_transpose(NumericsMatrix * A);

 /** set NumericsMatrix fields to NULL
   * \param A a matrix
   */
  void NM_null(NumericsMatrix* A);

  /** Check if a matrix is destructible.
   * \param[in] A the NumericsMatrix
   * \return true if the matrix is destructible */
  bool NM_destructible(const NumericsMatrix* A);

  /** Preservation of a matrix before in-place transformations such as
   * factorizations.
   * \param[in] A the NumericsMatrix
   * \return a pointer on the preserved Matrix;
   */
  RawNumericsMatrix* NM_preserve(NumericsMatrix* A);

  /** Set the matrix as destructible, clear the preserved data.
   * \param[in] A the NumericsMatrix
   * \return a pointer on the Matrix;
   */
  RawNumericsMatrix* NM_unpreserve(NumericsMatrix* A);

  /** Check for a previous LU factorization.
   * \param[in] A the NumericsMatrix
   * \return true if the matrix has been LU factorized.
   */
  bool NM_LU_factorized(const NumericsMatrix* const A);

  /** Check for a previous Cholesky factorization.
   * \param[in] A the NumericsMatrix
   * \return true if the matrix has been Cholesky factorized.
   */
  bool NM_Cholesky_factorized(const NumericsMatrix* const A);

  /** Check for a previous LDLT factorization.
   * \param[in] A the NumericsMatrix
   * \return true if the matrix has been Cholesky factorized.
   */
  bool NM_LDLT_factorized(const NumericsMatrix* const A);

  /** Set the factorization flag.
   * \param[in] A the NumericsMatrix,
   * \param[in] flag a boolean.
   */
  void NM_set_LU_factorized(NumericsMatrix* A, bool flag);
  void NM_set_Cholesky_factorized(NumericsMatrix* A, bool flag);
  void NM_set_LDLT_factorized(NumericsMatrix* A, bool flag);

  /** update the size of the matrix based on the matrix data
   * \param[in,out] A the matrix which size is updated*/
  void NM_update_size(NumericsMatrix* A);

  /** Allocate a csc matrix in A
   * \param A the matrix
   * \param nzmax number of non-zero elements
   */
  void NM_csc_alloc(NumericsMatrix* A, CS_INT nzmax);

  /** Allocate a csc matrix in A and set the vector of
   * column pointers to 0 such that the matrix is empty.
   * \param A the matrix
   * \param nzmax number of non-zero elements
   */
  void NM_csc_empty_alloc(NumericsMatrix* A, CS_INT nzmax);

  /** Allocate a triplet matrix in A
   * \param A the matrix
   * \param nzmax maximum number of non-zero elements
   */
  void NM_triplet_alloc(NumericsMatrix* A, CS_INT nzmax);


  /** Free memory for a NumericsMatrix. Warning: call this function only if you are sure that
      memory has been allocated for the structure in Numerics. This function is assumed that the memory is "owned" by this structure.
      Note that this function does not free m.
      \param m the matrix to be deleted.
   */

  void NM_clear(NumericsMatrix* m);

  NumericsMatrix *  NM_free(NumericsMatrix* m);

  /** Free memory for a NumericsMatrix except the dense matrix that is assumed not to be owned.
      \param m the matrix to be cleared.
   */
  void NM_clear_not_dense(NumericsMatrix* m);
  NumericsMatrix *  NM_free_not_dense(NumericsMatrix* m);
  /** Free memory for a NumericsMatrix except the SBM matrix that is assumed not to be owned.
      Note that this function does not free m.
      \param m the matrix to be cleared.
   */
  void NM_clear_not_SBM(NumericsMatrix* m);
  NumericsMatrix * NM_free_not_SBM(NumericsMatrix* m);
  



 
  /** Free memory for a NumericsMatrix except for a given storage. Warning: call this function only if you are sure that
      memory has been allocated for the structure in Numerics. This function is assumed that the memory is "owned" by this structure.
      Note that this function does not free m.
      \param m the matrix to be deleted.
      \param storageType to be kept.
   */
  void NM_clear_other_storages(NumericsMatrix* M, NM_types storageType);

  /**************************************************/
  /** setters and getters               *************/
  /**************************************************/

  /** insert an non zero entry into a NumericsMatrix.
   * for storageType = NM_SPARSE, a conversion to triplet is done for performing the entry in the
   * matrix. This method is expensive in terms of memory management. For a lot of entries, use
   * preferably a triplet matrix.
   * \param M the NumericsMatrix
   * \param i row index
   * \param j column index
   * \param val the value to be inserted.
   * \param threshold a threshold to filter the small value in magnitude (useful for dense to sparse conversion)
   */
  void NM_zentry(NumericsMatrix* M, int i, int j, double val, double threshold);

  /** insert an entry into a NumericsMatrix.
   * for storageType = NM_SPARSE, a conversion to triplet is done for performing the entry in the
   * matrix. This method is expensive in terms of memory management. For a lot of entries, use
   * preferably a triplet matrix.
   * \param M the NumericsMatrix
   * \param i row index
   * \param j column index
   * \param val the value to be inserted.
   */
  void NM_entry(NumericsMatrix* M, int i, int j, double val);

  /** get the value of a NumericsMatrix.
   * \param M the NumericsMatrix
   * \param i row index
   * \param j column index
   * \return  the value to be inserted.
   */
  double NM_get_value(const NumericsMatrix* const M, int i, int j);

  /** compare to NumericsMatrix up to machine accuracy (DBL_EPSILON)
   * \param A the NumericsMatrix
   * \param B the NumericsMatrix
   */
  bool NM_equal(NumericsMatrix* A, NumericsMatrix* B);

  /** compare to NumericsMatrix up to a given tolerance
   * \param A the NumericsMatrix
   * \param B the NumericsMatrix
   * \param tol the tolerance
   */
  bool NM_compare(NumericsMatrix* A, NumericsMatrix* B, double tol);

  /** return the number of non-zero element. For a dense matrix, it is the
   * product of the dimensions (e.g. an upper bound). For a sparse matrix, it is the true number
   * \param M the matrix
   * \return the number (or an upper bound) of non-zero elements in the matrix
   */
  size_t NM_nnz(const NumericsMatrix* M);

  /** get the (square) diagonal block of a NumericsMatrix. No allocation is done.
   * \param[in] M a NumericsMatrix
   * \param[in] block_row_nb the number of the block Row. Useful only in sparse case
   * \param[in] start_row the starting row. Useful only in dense case.
   * \param[in] size of the diag block. Only useful in dense case.
   * \param[out] Block the target. In the dense and sparse case (*Block) must be allocated by caller.
   *   In case of SBM case **Bout contains the resulting block (from the SBM).
   */
  void NM_extract_diag_block(NumericsMatrix* M, int block_row_nb, size_t start_row,
                             int size, double **Block);

  /** get a 3x3 diagonal block of a NumericsMatrix. No allocation is done.
   * \param[in] M a NumericsMatrix
   * \param[in] block_row_nb the number of the block row
   * \param[out] Block the target. In the dense and sparse case (*Block) must be allocated by caller.
   *   In case of SBM case **Bout contains the resulting block (from the SBM).
   */

  void NM_extract_diag_block3(NumericsMatrix* M, int block_row_nb, double **Block);
  
  /** get a 2x2 diagonal block of a NumericsMatrix. No allocation is done.
   * \param[in] M a NumericsMatrix
   * \param[in] block_row_nb the number of the block row
   * \param[out] Block the target. In the dense and sparse case (*Block) must be allocated by caller.
   *   In case of SBM case **Bout contains the resulting block (from the SBM).
   */

  void NM_extract_diag_block2(NumericsMatrix* M, int block_row_nb, double **Block);

  /** get a 5x5 diagonal block of a NumericsMatrix. No allocation is done.
   * \param[in] M a NumericsMatrix
   * \param[in] block_row_nb the number of the block row
   * \param[out] Block the target. In the dense and sparse case (*Block) must be allocated by caller.
   *   In case of SBM case **Bout contains the resulting block (from the SBM).
   */
  void NM_extract_diag_block5(NumericsMatrix* M, int block_row_nb, double **Block);

  /** get a 3x3 diagonal block of a NumericsMatrix. No allocation is done.
   * \param[in] M a NumericsMatrix
   * \param[in] block_row_nb the number of the block row
   * \param[out] Block the target.
   * In all cases (dense, sbm, and sparse) (*Block) must be allocated by caller.
   *  A copy is always performed
   */
  void NM_copy_diag_block3(NumericsMatrix* M, int block_row_nb, double **Block);


  /** Set the submatrix B into the matrix A on the position defined in
   *  (start_i, start_j) position.
   * \param[in] A a pointer to NumerixMatrix
   * \param[in] B a pointer toNumericsMatrix
   * \param[in] start_i a start row index
   * \param[in] start_j a start column index
   */
  void NM_insert(NumericsMatrix* A, const NumericsMatrix* const B,
                 const unsigned int start_i, const unsigned int start_j);

  /**************************************************/
  /** Matrix - vector product           *************/
  /**************************************************/

  /** Matrix - vector product y = A*x + y
      \param[in] sizeX dim of the vector x
      \param[in] sizeY dim of the vector y
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in,out] y the resulting vector
  */
  void NM_prod_mv_3x3(int sizeX, int sizeY,  NumericsMatrix* A,
                             double* const x, double* y);

  /** Row of a Matrix - vector product y = rowA*x or y += rowA*x, rowA being a submatrix of A (sizeY rows and sizeX columns)
      \param[in] sizeX dim of the vector x
      \param[in] sizeY dim of the vector y
      \param[in] currentRowNumber position of the first row of rowA in A (warning: real row if A is a double*, block-row if A is a SparseBlockStructuredMatrix)
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in,out] y the resulting vector
      \param[in] init = 0 for y += Ax, =1 for y = Ax
  */
  void NM_row_prod(int sizeX, int sizeY, int currentRowNumber, const NumericsMatrix* const A, const double* const x, double* y, int init);

  /** Row of a Matrix - vector product y = rowA*x or y += rowA*x, rowA being a submatrix of A (sizeY rows and sizeX columns)
      \param[in] sizeX dim of the vector x
      \param[in] sizeY dim of the vector y
      \param[in] block_start block number (only used for SBM)
      \param[in] row_start position of the first row of A (unused if A is SBM)
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in,out] y the resulting vector
      \param[in] xsave storage for saving the part of x set to 0
      \param[in] init if True y = Ax, else y += Ax
  */
  void NM_row_prod_no_diag(size_t sizeX, size_t sizeY, int block_start, size_t row_start, NumericsMatrix* A, double* x, double* y, double* xsave, bool init);

  /** Row of a Matrix - vector product y = rowA*x or y += rowA*x, rowA being a submatrix of A (3 rows and sizeX columns)
      \param[in] sizeX dim of the vector x
      \param[in] block_start block number (only used for SBM)
      \param[in] row_start position of the first row of A (unused if A is SBM)
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in,out] y the resulting vector
      \param[in] init if True y = Ax, else y += Ax
  */
  void NM_row_prod_no_diag3(size_t sizeX, int block_start, size_t row_start, NumericsMatrix* A, double* x, double* y, bool init);
  
  /** Row of a Matrix - vector product y = rowA*x or y += rowA*x, rowA being a submatrix of A (2 rows and sizeX columns)
      \param[in] sizeX dim of the vector x
      \param[in] block_start block number (only used for SBM)
      \param[in] row_start position of the first row of A (unused if A is SBM)
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in,out] y the resulting vector
      \param[in] init if True y = Ax, else y += Ax
  */
  void NM_row_prod_no_diag2(size_t sizeX, int block_start, size_t row_start, NumericsMatrix* A, double* x, double* y, bool init);


  void NM_row_prod_no_diag1x1(size_t sizeX, int block_start, size_t row_start, NumericsMatrix* A, double* x, double* y, bool init);

  /** Matrix vector multiplication : y = alpha A x + beta y
   * \param[in] alpha scalar
   * \param[in] A a NumericsMatrix
   * \param[in] x pointer on a dense vector of size A->size1
   * \param[in] beta scalar
   * \param[in,out] y pointer on a dense vector of size A->size1
   */
  void NM_gemv(const double alpha, NumericsMatrix* A, const double *x,
               const double beta,
               double *y);

  /** Matrix matrix multiplication : C = alpha A B + beta C
   * \param[in] alpha scalar
   * \param[in] A a NumericsMatrix
   * \param[in] B a NumericsMatrix
   * \param[in] beta scalar
   * \param[in,out] C a NumericsMatrix
   */
  void NM_gemm(const double alpha, NumericsMatrix* A, NumericsMatrix* B,
               const double beta, NumericsMatrix *C);

   /** Matrix matrix multiplication : C = A B
   * \param[in] A a NumericsMatrix
   * \param[in] B a NumericsMatrix
   * \param[in,out] C a NumericsMatrix
   */
  RawNumericsMatrix * NM_multiply(NumericsMatrix* A, NumericsMatrix* B);

  /** Transposed matrix multiplication : y += alpha transpose(A) x + y
   * \param[in] alpha scalar
   * \param[in] A a NumericsMatrix
   * \param[in] x pointer on a dense vector of size A->size1
   * \param[in] beta scalar
   * \param[in,out] y pointer on a dense vector of size A->size1
   */
  void NM_tgemv(const double alpha, NumericsMatrix* A, const double *x,
                const double beta,
                double *y);

  /**************************************************/
  /** matrix conversion display *********************/
  /**************************************************/


  /**************************************************/
  /** matrix and vector display *********************/
  /**************************************************/

  void NM_dense_to_sparse(const NumericsMatrix* const A, NumericsMatrix* B, double threshold);

  /** Copy a NumericsMatrix into another with dense storage.
      \param A source matrix (any kind of storage)
      \param B targeted matrix, must be dense with the same dimension as A
  */
  int NM_to_dense(const NumericsMatrix* const A, NumericsMatrix* B);

  /** Screen display of the matrix content stored as a double * array in Fortran style
      \param m the matrix to be displayed
      \param nRow the number of rows
      \param nCol the number of columns
      \param lDim the leading dimension of M
   */
  void NM_dense_display_matlab(double * m, int nRow, int nCol, int lDim);

  /** Screen display of the matrix content stored as a double * array in Fortran style
      \param m the matrix to be displayed
      \param nRow the number of rows
      \param nCol the number of columns
      \param lDim the leading dimension of M
   */
  void NM_dense_display(double * m, int nRow, int nCol, int lDim);

  /** Screen display of the vector content stored as a double * array
      \param m the vector to be displayed
      \param nRow the number of rows
   */
  void NM_vector_display(double * m, int nRow);


  /** Screen display of the matrix content
      \param M the matrix to be displayed
   */
  void NM_display(const NumericsMatrix* const M);

  /** Screen display of the matrix storage
      \param M the matrix to be displayed
   */
  void NM_display_storageType(const NumericsMatrix* const M);


  /** Screen display raw by raw of the matrix content
      \param m the matrix to be displayed
  */
  void NM_display_row_by_row(const NumericsMatrix* const m);

  /**************************************************/
  /** matrix I/O                *********************/
  /**************************************************/

  /** PrintInFile  of the matrix content
     \param M the matrix to be printed
     \param filename the corresponding name of the file
  */
  void NM_write_in_filename(const NumericsMatrix* const M, const char *filename);

  /** Read in file  of the matrix content
     \param M the matrix to be read
     \param filename the corresponding name of the file
  */
  void NM_read_in_filename(NumericsMatrix* const M, const char *filename);

  /** PrintInFile  of the matrix content
     \param M the matrix to be printed
     \param file filename the corresponding file
  */

  void NM_write_in_file(const NumericsMatrix* const M, FILE* file);

  /** Read in file  of the matrix content without performing memory allocation
     \param M the matrix to be read
     \param file the corresponding  file
  */
  void NM_read_in_file(NumericsMatrix* const M, FILE *file);

  /** Create from file a NumericsMatrix with  memory allocation
     \param file the corresponding  file
     \return 0 if the matrix
  */
  RawNumericsMatrix*  NM_new_from_file(FILE *file);
  RawNumericsMatrix*  NM_new_from_filename(const char * filename);

  /**  NM_write_in_file_scilab of the matrix content
   \param M the matrix to be printed
   \param file the corresponding file
  */
  void NM_write_in_file_scilab(const NumericsMatrix* const M, FILE* file);

 /**   NM_write_in_file_python of the matrix content
   \param M the matrix to be printed
   \param file the corresponding file
  */
  void NM_write_in_file_python(const NumericsMatrix* const M, FILE* file);

  /** Read in file for scilab  of the matrix content
     \param M the matrix to be read
     \param file the corresponding  file
  */
  void NM_read_in_file_scilab(NumericsMatrix* const M, FILE *file);

  /** Clear dense storage, if it is existent.
   * \param[in,out] A a Numericsmatrix
   */
  void NM_clearDense(NumericsMatrix* A);

  /** Clear sparse block storage, if it is existent.
   * \param[in,out] A a Numericsmatrix
   */
  void NM_clearSparseBlock(NumericsMatrix* A);

  /** Clear sparse data, if it is existent.
      The linear solver parameters are also cleared.
   * \param[in,out] A a Numericsmatrix
   */
  void NM_clearSparse(NumericsMatrix* A);

  /** Clear triplet storage, if it is existent.
   * \param[in,out] A a Numericsmatrix
   */
  void NM_clearTriplet(NumericsMatrix* A);

  /** Clear half triplet storage, if it is existent.
   * \param[in,out] A a Numericsmatrix
   */
  void NM_clearHalfTriplet(NumericsMatrix* A);

  /** Clear compressed column storage, if it is existent.
   * \param[in,out] A a Numericsmatrix
   */
  void NM_clearCSC(NumericsMatrix* A);

  /** Clear transposed compressed column storage, if it is existent.
    * \param[in,out] A a Numericsmatrix
    */
  void NM_clearCSCTranspose(NumericsMatrix* A);

  /** Clear compressed row storage, if it is existent.
   * \param[in,out] A a Numericsmatrix
   */
  void NM_clearCSR(NumericsMatrix* A);

  /** Clear triplet, csc, csc transposed storage, if they are existent.
    * Linear solver parameters are preserved.
    * \param[in,out] A a Numericsmatrix
    */
  void NM_clearSparseStorage(NumericsMatrix *A);



  /** XXXXXX: to be rewritten Direct computation of the solution of a real system of linear
   * equations: A x = b. The factorized matrix A is kept for future solve.
   * If A is already factorized, the solve the linear system from it
   * \warning this is not enable for all the solvers, your mileage may vary
   * \param[in,out] A a NumericsMatrix. On a dense factorisation
   * A.iWork is initialized.
   * \param[in,out] b pointer on a dense vector of size A->size1
   * \param keep if set to NM_KEEP_FACTORS, keep all the info related to the factorization to
   * allow for future solves. If A is already factorized, just solve the linear
   * system. If set to NM_PRESERVE, preserve the original matrix (just used in
   * the dense case). if NM_NONE, discard everything.
   * \return 0 if successful, else the error is specific to the backend solver
   * used
   */

  /* LU factorization of the matrix. If the matrix has already been
   * factorized (i.e if NM_LU_factorized(A) return true), nothing is
   * done. To force a new factorization one has to set factorization
   * flag to false : NM_set_LU_factorized(A, false) before the call to
   * NM_LU_factorize.
   * If the matrix is preserved, that means that a call to
   * NM_preserve(A) has been done before the call to NM_LU_factorize,
   * it is not destroyed, but the factorized part remains accessible for
   * subsequent calls to NM_LU_solve.
   * If the matrix is not preserved, then it is replaced by the
   * factorized part.
   * \param[in] A the NumericsMatrix
   * \return an int, 0 means the matrix has been factorized. */
  int NM_LU_factorize(NumericsMatrix* A);
  int NM_Cholesky_factorize(NumericsMatrix* A);
  int NM_LDLT_factorize(NumericsMatrix* A);

  /* Solve linear system with multiple right hand size. A call to
   * NM_LU_factorize is done at the beginning.

   * \param[in] A the NumericsMatrix. A is not destroyed if it has
   * been preserved by a call to NM_preserve(A).

   * \param[in,out] b the right hand size which is a pointer on a
   * matrix of double. It is replaced by the solutions

   * \param[in] nrhs the number of right hand side.
   * \return 0 if the solve succeeded.
   */
  int NM_LU_solve(NumericsMatrix* A,  double *b, unsigned int nrhs);
  int NM_LU_solve_matrix_rhs(NumericsMatrix* Ao, NumericsMatrix* B);
  int NM_Cholesky_solve(NumericsMatrix* A,  double *b, unsigned int nrhs);
  int NM_Cholesky_solve_matrix_rhs(NumericsMatrix* Ao, NumericsMatrix* B);
  int NM_LDLT_solve(NumericsMatrix* A,  double *b, unsigned int nrhs);
  int NM_LDLT_refine(NumericsMatrix* Ao, double *x , double *b, unsigned int nrhs, double tol, int maxitref, int job );

  int NM_gesv_expert(NumericsMatrix* A, double *b, unsigned keep);
  int NM_posv_expert(NumericsMatrix* A, double *b, unsigned keep);

  int NM_gesv_expert_multiple_rhs(NumericsMatrix* A, double *b, unsigned int n_rhs, unsigned keep);



  /**  Computation of the inverse of a NumericsMatrix A usinf NM_gesv_expert
   * \param[in,out] A a NumericsMatrix.
   * \return the matrix inverse.
   */
  NumericsMatrix* NM_LU_inv(NumericsMatrix* A);

  int NM_inverse_diagonal_block_matrix_in_place(NumericsMatrix* A);

  /** Direct computation of the solution of a real system of linear
   * equations: A x = b.
   * \param[in,out] A a NumericsMatrix. On a dense factorisation
   * A.iWork is initialized.
   * \param[in,out] b pointer on a dense vector of size A->size1
   * \param preserve preserve the original matrix data. Only useful in the
   * dense case, where the LU factorization is done in-place.
   * \return 0 if successful, else the error is specific to the backend solver
   * used
   */
  static inline int NM_gesv(NumericsMatrix* A, double *b, bool preserve)
  {
    return NM_gesv_expert(A, b, preserve ? NM_PRESERVE : NM_NONE);
  }

  /**  Computation of the inverse of a NumericsMatrix A usinf NM_gesv_expert
   * \param[in,out] A a NumericsMatrix.
   * \return the matrix inverse.
   */
  NumericsMatrix* NM_gesv_inv(NumericsMatrix* A);

  /** Set the linear solver
   * \param A the matrix
   * \param solver_id id of the solver
   */
  void NM_setSparseSolver(NumericsMatrix* A, NSM_linear_solver solver_id);

  /** Get Matrix internal data with initialization if needed.
   * \param[in,out] A a NumericsMatrix.
   * \return a pointer on internal data.
   */
  NumericsMatrixInternalData* NM_internalData(NumericsMatrix* A);

  /** Allocate the internalData structure (but not its content!)
   * \param M the matrix to modify
   */
  void NM_internalData_new(NumericsMatrix* M);

  /** Copy the internalData structure
   * \param M the matrix to modify
   */
  void NM_internalData_copy(const NumericsMatrix* const A, NumericsMatrix* B );

  /** Integer work vector initialization, if needed.
   * \param[in,out] A pointer on a NumericsMatrix.
   * \param[in] size number of element to allocate
   * \param[in] sizeof_elt of an element in bytes (result of sizeof for instance)
   * \return pointer on A->iWork allocated space of with the right size
   */
  void* NM_iWork(NumericsMatrix *A, size_t size, size_t sizeof_elt);

  /** Double workspace initialization, if needed.
   * \param[in,out] A pointer on a NumericsMatrix.
   * \param[in] size the size of needed space.
   * \return pointer on A->dWork allocated space of with the right size
   */
  double* NM_dWork(NumericsMatrix *A, int size);

  /** Add a constant term to the diagonal elements, when the block of the SBM
   * are 3x3
   * \param M the matrix
   * \param alpha the term to add
   */
  void NM_add_to_diag3(NumericsMatrix* M, double alpha);

  /** Add a constant term to the diagonal elements, when the block of the SBM
   * are 5x5
   * \param M the matrix
   * \param alpha the term to add
   */
  void NM_add_to_diag5(NumericsMatrix* M, double alpha);

  /** Add two matrices with coefficients C = alpha*A + beta*B
   * \param alpha the first coefficient
   * \param A the first  matrix
   * \param beta the second coefficient
   * \param B the second  matrix
   * \return C a new NumericsMatrix
   */
  RawNumericsMatrix *  NM_add(double alpha, NumericsMatrix* A, double beta, NumericsMatrix* B);

  /** Multiply a matrix with a double alpha*A --> A
   * \param alpha the  coefficient
   * \param A the   matrix
   */
  void  NM_scal(double alpha, NumericsMatrix* A);

  /** assert that a NumericsMatrix has the right structure given its type
   * \param type expected type
   * \param M the matrix to check
   */
    static inline void NM_assert(NM_types type, NumericsMatrix* M)
  {
#ifndef NDEBUG
    assert(M && "NM_assert :: the matrix is NULL");
    assert(M->storageType == type && "NM_assert :: the matrix has the wrong type");
    switch(type)
    {
      case NM_DENSE:
        assert(M->matrix0);
        break;
      case NM_SPARSE_BLOCK:
        assert(M->matrix1);
        break;
      case NM_SPARSE:
        assert(M->matrix2);
        break;
      default:
        assert(0 && "NM_assert :: unknown storageType");
    }
#endif
  }

  /** Check the matrix (the sparse format for now)
   * \param A the matrix to check
   * \return 0 if the matrix storage is fine, 1 if not*/
  int NM_check(const NumericsMatrix* const A);

 /** Compute the  1-norm of a sparse matrix = max (sum (abs (A))), largest column sum of a matrix (the sparse format for now)
   * \param A the matrix
   * \return the norm*/
  double NM_norm_1(NumericsMatrix* const A);

 /** Compute the  inf-norm of a sparse matrix = max (sum (abs (A^T))), largest row  sum of a matrix (the sparse format for now)
   * \param A the matrix
   * \return the norm*/
  double NM_norm_inf(NumericsMatrix* const A);

  int NM_is_symmetric(NumericsMatrix* A);
  double NM_symmetry_discrepancy(NumericsMatrix* A);


  /** Pass a NumericsMatrix through swig typemaps.
   * This is only useful in python.
   * \param A the matrix
   * \return a NumericsMatrix
   */
  static inline NumericsMatrix* NM_convert(NumericsMatrix* A)
  {
    return A;
  }

  /** Compute the  maximum eigenvalue with the iterated power method
   * \param A the matrix
   * \return the maximum eigenvalue*/
  double NM_iterated_power_method(NumericsMatrix* A, double tol, int itermax);

  /* Compute the maximum values by columns
   *  \param A the matrix
   *  \param max the vector of max that must be preallocated
   *  \return info
   */
  int NM_max_by_columns(NumericsMatrix *A, double * max);

  /* Compute the maximum values by rows
   *  \param A the matrix
   *  \param max the vector of max that must be preallocated
   *  \return info
   */
  int NM_max_by_rows(NumericsMatrix *A, double * max);

  /* Compute the maximum absolute values by columns
   *  \param A the matrix
   *  \param max the vector of max that must be preallocated
   *  \return info
   */
  int NM_max_abs_by_columns(NumericsMatrix *A, double * max);

  /* Compute the maximum absolute values by rows
   *  \param A the matrix
   *  \param max the vector of max that must be preallocated
   *  \return info
   */
  int NM_max_abs_by_rows(NumericsMatrix *A, double * max);

  /* Compute the balancing matrices for a given matrix by iteration
   *  \param A the matrix
   *  \param tol tolerance on the balanced matrix
   *  \param itermax max number of iterations
   *  \param alloated structure for the balancing matrices and the balanced matrix
   * \return 0 if succeed.
   */
  int NM_compute_balancing_matrices(NumericsMatrix* A, double tol, int itermax, BalancingMatrices * B);

  /* Create a Balancing Matrices structure
   *  \param A the matrix  to be balanced
   */
  BalancingMatrices * NM_BalancingMatrices_new(NumericsMatrix* A);


  /* free a Balancing Matrices structure
   */
  BalancingMatrices * NM_BalancingMatrices_free(BalancingMatrices* A);

  /* Reset the version of a NM_types storage.
   *\param M the NumericsMatrix,
   *\param id the NM_types storage
   */
  void NM_reset_version(NumericsMatrix* M, NM_types id);

  /* Reset versions of all storages.
   *\param M the NumericsMatrix
   */
  void NM_reset_versions(NumericsMatrix* M);


#ifdef WITH_OPENSSL
  /* Compute sha1 hash of matrix values. Matrices of differents size and same
   * values have the same hash.
   * \param[in] A the matrix
   * \param[in,out] digest an allocated space of size SHA_DIGEST_LENGTH
   */
  void NM_compute_values_sha1(NumericsMatrix* A, unsigned char * digest);

  /* Get stored sha1 hash of a matrix.
   * \param[in] A the matrix
   * \return a pointer on the current raw sha1 hash
   */
  unsigned char* NM_values_sha1(NumericsMatrix* A);

   /* Set sha1 hash of a matrix. Matrices of differents size and same
   * values have the same hash.
   * \param[in] A the matrix
   */
  void NM_set_values_sha1(NumericsMatrix* A);

  /* Clear sha1 hash of a matrix.
   * \param[in] A the matrix
   */
  void NM_clear_values_sha1(NumericsMatrix* A);

  /* Check if matrix has beend modified after a previous NM_set_values_sha1.
   * \param[in] A the NumericsMatrix
   * \return true if the matrix is the same
   */
  bool NM_check_values_sha1(NumericsMatrix* A);

  /* Compare two matrices with sha1. NM_set_values_sha1 must be called
   * on the two matrices before.
   * \param[in] A a NumericsMatrix
   * \param[in] B a NumericsMatrix
   * \return true if the matrices have the same values
   */
  bool NM_equal_values_sha1(NumericsMatrix* A, NumericsMatrix* B);

#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
