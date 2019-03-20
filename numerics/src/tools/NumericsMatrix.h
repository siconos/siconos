/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>

#include "NumericsFwd.h"
#include "SiconosConfig.h"
#include "CSparseMatrix.h" // for CS_INT
#include "NM_MPI.h"

/** \struct NumericsMatrixInternalData NumericsMatrix.h
 * Structure for simple workspaces
 */
typedef struct
{
  size_t iWorkSize; /**< size of iWork */
  void *iWork; /**< integer workspace */
  size_t sizeof_elt ; ; /**< sizeof_elt of an element in bytes (result of sizeof for instance)*/
  size_t dWorkSize; /**< size of dWork */
  double *dWork; /**< double workspace */
  bool isLUfactorized; /**<  true if the matrix has already been LU-factorized */
  bool isInversed; /**<  true if the matrix containes its inverse (in place inversion) */
#ifdef HAVE_MPI
  MPI_Comm mpi_comm; /**< optional mpi communicator */
#endif
} NumericsMatrixInternalData;

/** \struct NumericsMatrix NumericsMatrix.h
    Interface to different type of matrices in numerics component.

    See NM_* functions for linear algebra operations on dense, sparse block and sparse storage.
*/
struct NumericsMatrix
{
  int storageType; /**< the type of storage:
                      0: dense (double*),
                      1: SparseBlockStructuredMatrix,
                      2: classical sparse (csc, csr or triplet) via CSparse (from T. Davis)*/
  int size0; /**< number of rows */
  int size1; /**< number of columns */
  double* matrix0; /**< dense storage */
  SparseBlockStructuredMatrix* matrix1; /**< sparse block storage */
  NumericsSparseMatrix* matrix2; /**< csc, csr or triplet storage */

  NumericsMatrixInternalData* internalData; /**< internal storage, used for workspace among other things */

};

/*! RawNumericsMatrix is used without conversion in python */
typedef NumericsMatrix RawNumericsMatrix;

/*! Available types of storage for NumericsMatrix */
typedef enum NumericsMatrix_types {
  NM_DENSE,        /**< dense format */
  NM_SPARSE_BLOCK, /**< sparse block format */
  NM_SPARSE,          /**< compressed column format */
} NM_types;

/*! option for gesv factorization */
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
  /** create a NumericsMatrix and allocate the memory according to the matrix type
   * \param storageType the type of storage
   * \param size0 number of rows
   * \param size1 number of columns
   * \return a pointer to a NumericsMatrix
   */
  RawNumericsMatrix* NM_create(int storageType, int size0, int size1);

  /** create a NumericsMatrix and possibly set the data
   * \param storageType the type of storage
   * \param size0 number of rows
   * \param size1 number of columns
   * \param data pointer to the matrix data. If NULL, all matrixX fields are
   * set to NULL
   * \return a pointer to a NumericsMatrix
   */
  RawNumericsMatrix* NM_create_from_data(int storageType, int size0, int size1, void* data);
  RawNumericsMatrix* NM_create_from_filename(char *filename);
  RawNumericsMatrix* NM_create_from_file(FILE *file);


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
   */
  void NM_copy_to_sparse(const NumericsMatrix* const A, NumericsMatrix* B);

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
  void NM_fill(NumericsMatrix* M, int storageType, int size0, int size1, void* data);

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
  static inline void NM_null(NumericsMatrix* A)
  {
    A->matrix0 = NULL;
    A->matrix1 = NULL;
    A->matrix2 = NULL;
    A->internalData = NULL;
  }

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

  /** Allocate a csr matrix in A
   * \param A the matrix
   * \param nzmax number of non-zero elements
   */
  void NM_csr_alloc(NumericsMatrix* A, CS_INT nzmax);

  /** Free memory for a NumericsMatrix. Warning: call this function only if you are sure that
      memory has been allocated for the structure in Numerics. This function is assumed that the memory is "owned" by this structure.
      Note that this function does not free m.
      \param m the matrix to be deleted.
   */
  void NM_free(NumericsMatrix* m);


  /**************************************************/
  /** setters and getters               *************/
  /**************************************************/

  /** insert an non zero entry into a NumericsMatrix.
   * for storageType = NM_SPARSE, a conversion to triplet is done for performing the entry in the
   * matrix. This method is expensice in terms of memory management. For a lot of entries, use
   * preferably a triplet matrix.
   * \param M the NumericsMatrix
   * \param i row index
   * \param j column index
   * \param val the value to be inserted.
   */
  void NM_zentry(NumericsMatrix* M, int i, int j, double val);

  /** get the value of a NumericsMatrix.
   * \param M the NumericsMatrix
   * \param i row index
   * \param j column index
   * \return  the value to be inserted.
   */
  double NM_get_value(NumericsMatrix* M, int i, int j);

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
  /** get a 3x3 diagonal block of a NumericsMatrix. No allocation is done.
   * \param[in] M a NumericsMatrix
   * \param[in] block_row_nb the number of the block row
   * \param[out] Block the target.
   * In all cases (dense, sbm, and sparse) (*Block) must be allocated by caller.
   *  A copy is always performed
   */
  void NM_copy_diag_block3(NumericsMatrix* M, int block_row_nb, double **Block);

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

  void NM_dense_to_sparse(const NumericsMatrix* const A, NumericsMatrix* B);

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
  RawNumericsMatrix*  NM_new_from_filename(char * filename);

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



  /** Direct computation of the solution of a real system of linear
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
  int NM_gesv_expert(NumericsMatrix* A, double *b, unsigned keep);

  int NM_gesv_expert_multiple_rhs(NumericsMatrix* A, double *b, unsigned int n_rhs, unsigned keep);

  /**  Computation of the inverse of a NumericsMatrix A usinf NM_gesv_expert
   * \param[in,out] A a NumericsMatrix.
   * \param[out] Ainv the matrix inverse.
   */
  int NM_inv(NumericsMatrix* A, NumericsMatrix* Ainv);

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


  /** Set the linear solver
   * \param A the matrix
   * \param solver_id the solver
   */
  void NM_setSparseSolver(NumericsMatrix* A, unsigned solver_id);

  /** Get Matrix internal data with initialization if needed.
   * \param[in,out] A a NumericsMatrix.
   * \return a pointer on internal data.
   */
  NumericsMatrixInternalData* NM_internalData(NumericsMatrix* A);

  /** Allocate the internalData structure (but not its content!)
   * \param M the matrix to modify
   */
  static inline void NM_internalData_new(NumericsMatrix* M)
  {
    M->internalData = (NumericsMatrixInternalData *)malloc(sizeof(NumericsMatrixInternalData));
    M->internalData->iWork = NULL;
    M->internalData->iWorkSize = 0;
    M->internalData->dWork = NULL;
    M->internalData->dWorkSize = 0;
    M->internalData->isLUfactorized = 0;
#ifdef HAVE_MPI
    M->internalData->mpi_comm = MPI_COMM_NULL;
#endif
  }
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

  /** Add two matrices with coefficients C = alpha*A + beta*B
   * \param alpha the first coefficient
   * \param A the first  matrix
   * \param beta the second coefficient
   * \param B the second  matrix
   * \return C a new NumericsMatrix
   */
  RawNumericsMatrix *  NM_add(double alpha, NumericsMatrix* A, double beta, NumericsMatrix* B);

  /** assert that a NumericsMatrix has the right structure given its type
   * \param type expected type
   * \param M the matrix to check
   */
  static inline void NM_assert(const int type, NumericsMatrix* M)
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
   * This is only usefull in python.
   * \param A the matrix
   */
  static inline NumericsMatrix* NM_convert(NumericsMatrix* A)
  {
    return A;
  };

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
