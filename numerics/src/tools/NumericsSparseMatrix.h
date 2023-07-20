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

#ifndef NumericsSparseMatrix_H
#define NumericsSparseMatrix_H

/*!\file NumericsSparseMatrix.h
 * Data structures and functions for sparse matrices
 *
 */

#include <stdio.h>          // for size_t, FILE
#include "CSparseMatrix.h"  // for CSparseMatrix, CS_INT
#include "NumericsFwd.h"    // for NumericsSparseMatrix, NSM_linear_solver_p...
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

#include "NumericsDataVersion.h"

/**\struct linalg_data_t NumericsSparseMatrix.h
 * generic data struct for linear algebra operations
 */
typedef struct linalg_data_t
{
  int id;
  void  (*free_fn)(struct linalg_data_t*);
} linalg_data_t;

typedef enum { SN_LINALG_UNKNOWN, SN_LINALG_MKL } linalg_data_id;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif


  /** \enum NSM_linear_solver NumericsSparseMatrix.h
   * id for linear algebra solvers */
  typedef enum { NSM_CSPARSE, NSM_MUMPS, NSM_UMFPACK, NSM_MKL_PARDISO, NSM_SUPERLU, NSM_SUPERLU_MT, NSM_HSL } NSM_linear_solver;

  typedef void (*freeNSLSP)(void* p);

  /** \enum NumericsSparseTypesNZ
   * value of nz for some matrix storage type */
  typedef enum { NSM_CS_CSC = -1, NSM_CS_CSR = -2 } NumericsSparseTypesNZ;

  /** \struct NSM_linear_solver_params NumericsSparseMatrix.h
   * solver-specific parameters*/
  struct NSM_linear_solver_params
  {

    NumericsMatrix * parent_matrix;
    NSM_linear_solver solver;
    NSM_linear_solver LDLT_solver;

    void* linear_solver_data; /**< solver-specific data (or workspace) */
    freeNSLSP solver_free_hook; /**< solver-specific hook to free linear_solver_data  */

    int* iWork; /**< integer work vector array (internal) */
    int iWorkSize; /**< size of integer work vector array */
    double* dWork;
    int dWorkSize;

    linalg_data_t* linalg_data; /**< data for the linear algebra */
  };

  /**\enum NumericsSparseOrigin NumericsSparseMatrix.h
   * matrix storage types */
  typedef enum { NSM_UNKNOWN, NSM_TRIPLET, NSM_CSC, NSM_CSR, NSM_HALF_TRIPLET } NumericsSparseOrigin;

  typedef NumericsSparseOrigin NSM_t;

  /** \struct NumericsSparseMatrix NumericsSparseMatrix.h
   * Sparse matrix representation in Numerics. The supported format are:
   * triplet (aka coordinate, COO), CSC (via CSparse) and CSR if MKL is used */
  struct NumericsSparseMatrix
  {
    CSparseMatrix* triplet;         /**< triplet format, aka coordinate */
    CSparseMatrix* half_triplet;    /**< half triplet format for symmetric matrices */
    CSparseMatrix* csc;             /**< csc matrix */
    CSparseMatrix* trans_csc;       /**< transpose of a csc matrix (used by CSparse) */
    CSparseMatrix* csr;             /**< csr matrix, only supported with mkl */
    CS_INT*        diag_indx;       /**< indices for the diagonal terms.
                                         Very useful for the proximal perturbation */
    NSM_t       origin;          /**< original format of the matrix */
    NSM_linear_solver_params* linearSolverParams;
                                    /**< solver-specific parameters */

    NumericsDataVersion versions[5];
  };


  /** Initialize the fields of a NumericsSparseMatrix
   *
   *  \param A the sparse matrix
   */
  void NSM_null(NumericsSparseMatrix* A);

  /** New and empty NumericsSparseMatrix with correctly initialized fields.
   *
   *  \return a pointer on the allocated space.
   */
  NumericsSparseMatrix* NSM_new(void);

  NumericsSparseMatrix * NSM_triplet_eye(unsigned int size);

  NumericsSparseMatrix * NSM_triplet_scalar(unsigned int size, double s);

  /** Free allocated space for a NumericsSparseMatrix.
   *
   *  \param A a NumericsSparseMatrix
   *  \return NULL on success
   */
  NumericsSparseMatrix* NSM_clear(NumericsSparseMatrix* A);

  /** Copy NumericsSparseMatrix version.
   *
   *  \param A a NumericsSparseMatrix
   *  \param B a NumericsSparseMatrix
   */
  void NSM_version_copy(const NumericsSparseMatrix* const A,
                        NumericsSparseMatrix* B);

  /** Copy a NumericsSparseMatrix.
   *
   *  \param A a NumericsSparseMatrix
   *  \param B a NumericsSparseMatrix
   */
  void NSM_copy(NumericsSparseMatrix* A, NumericsSparseMatrix* B);

   /** Free a workspace related to a LU factorization
   *
   *  \param p the structure to free
   */
  void NSM_clear_p(void *p);

  /** Get the data part of sparse matrix
   *
   *  \param A the sparse matrix
   *  \return a pointer to the data array
   */
  double* NSM_data(NumericsSparseMatrix* A);


  /** Get the LU factors for cs_lusol
   *
   *  \param p the structure holding the data for the solver
   */
  static inline void* NSM_linear_solver_data(NSM_linear_solver_params* p)
  {
    return p->linear_solver_data;
  }
  /** Get the workspace for the sparse solver
   *
   *  \param p the structure holding the data for the solver
   *  \return the (double) workspace
   */
  static inline double* NSM_workspace(NSM_linear_solver_params* p)

  {
    return p->dWork;
  }

  /** get the number of non-zero (nnz) in a sparse matrix
   *
   *  \param A the matrix
   *  \return the number of non-zero elements in the matrix
   */
  size_t NSM_nnz(const CSparseMatrix* const A);

  /** return the set of indices corresponding to the diagonal elements of the
   *  matrix
   *  \warning should be better tested
   *
   *  \param M the matrix
   *  \return the list of indices for the diagonal elements
   */
  CS_INT* NSM_diag_indices(NumericsMatrix* M);


  /** Extract a block from a sparse matrix
   *
   *  \param M matrix
   *  \param blockM dense storage for the block
   *  \param pos_row starting row for the block
   *  \param pos_col starting column for the block
   *  \param block_row_size block width
   *  \param block_col_size block height
   */
  void NSM_extract_block(NumericsMatrix* M, double* blockM, size_t pos_row, size_t pos_col, size_t block_row_size, size_t block_col_size);

  /** Free allocated space for NSM_linear_solver_params.
   *
   *  \param p a NSM_linear_solver_params
   *  \return NULL on success
   */
  NSM_linear_solver_params* NSM_linearSolverParams_free(NSM_linear_solver_params* p);

  /** New and empty NSM_linear_solver_params.
   *
   *  \return a pointer on the allocated space.
   */
  NSM_linear_solver_params* NSM_linearSolverParams_new(void);


  /** Get linear solver parameters with initialization if needed.
   *
   *  \param[in,out] A a NumericsMatrix.
   *  \return a pointer on parameters.
   */
  NSM_linear_solver_params* NSM_linearSolverParams(NumericsMatrix* A);

  /** Check and fix a matrix, if needed
   *
   *  \param A the matrix to check, modified if necessary to have ordered indices
   */
  void NSM_fix_csc(CSparseMatrix* A);

  void NSM_sort_csc(CSparseMatrix* A);

  /** return the origin of a sparse part of a matrix
   *
   *  \param M the matrix
   *  \return -1 if the matrix has no sparse representation, the origin
   * otherwise*/
  unsigned NSM_origin(const NumericsSparseMatrix* M);

  /** return the sparse matrix that has the original label
   *
   *  \param M the matrix
   *  \return the sparse matrix that is at the origin, or NULL if an error occur
   **/
  CSparseMatrix* NSM_get_origin(const NumericsSparseMatrix* M);

  void NSM_write_in_file(const NumericsSparseMatrix* m, FILE* file);

  /** New and empty NumericsSparseMatrix with correctly initialized fields.
   *
   *  \return a pointer on the allocated space.
   */
  NumericsSparseMatrix* NSM_new_from_file(FILE *file);

  int NSM_to_dense(const NumericsSparseMatrix * const A, double * B);

  /** Get current version of a type of csparse matrix.
   *
   *  \param M the NumericsSparseMatrix,
   *  \param type the type of sparse storage from NumericsSparseOrigin
   *  \return a comparable version. */
  version_t NSM_version(const NumericsSparseMatrix* M, NSM_t type);

  /** Get the maximum of versions of csparse matrices.
   *
   *  \param M the NumericsSparseMatrix,
   *  \return a comparable version. */
  version_t NSM_max_version(const NumericsSparseMatrix* M);

  /** Set the version of a NumericsSparseMatrix.
   *
   *  \param M the NumericsSparseMatrix,
   *  \param type the NumericsSparseOrigin of storage,
   *  \param value the new version.
   */
  void NSM_set_version(NumericsSparseMatrix* M, NSM_t type,
                       version_t value);

  /** Reset all versions of a NumericsSparseMatrix.
   *
   *  \param M the NumericsSparseMatrix.
   */
  void NSM_reset_versions(NumericsSparseMatrix *M);

  /** Reset version of a sparse storage.
   *
   *  \param M the NumericsSparseMatrix,
   *  \param type the NumericsSparseOrigin of storage.
   */
  void NSM_reset_version(NumericsSparseMatrix*M, NSM_t type);

  /** Increment the version of a NumericsSparseMatrix.
   *
   *  \param M the NumericsSparseMatrix,
   *  \param type the NumericsSparseOrigin of storage
   */
  void NSM_inc_version(NumericsSparseMatrix* M, NSM_t type);

  /** Get the NumericsSparseOrigin with the latest version.
   *
   *  \param M the NumericsSparseMatrix
   *  \return the NumericsSparseOrigin.
   */
  NSM_t NSM_latest_id(const NumericsSparseMatrix* M);
  
  /** Get most recent CSparseMatrix.
   *
   *  \param M the NumericsSparseMatrix
   *  \return a pointer on a CSparseMatrix.
   */
  CSparseMatrix* NSM_latest(const NumericsSparseMatrix* M);

  /** Sync matrix origin and version
   * \param M the NumericsSparseMatrix
   */
  void NSM_version_sync(NumericsSparseMatrix* M);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
