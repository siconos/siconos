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

#ifndef NumericsSparseMatrix_H
#define NumericsSparseMatrix_H

/*!\file NumericsSparseMatrix.h
 * \brief Data structures and functions for sparse matrices
 *
 * \author Olivier Huber
 */

#include "SiconosConfig.h"
#include "NumericsFwd.h"

/* Private definitions -- must #include csparse.h, not installed, to
 * use internally. */
struct cs_sparse;
typedef struct cs_sparse CSparseMatrix;

typedef void (*freeNSLSP)(void* p);

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


  /** \enum NumericsSparseLinearSolver NumericsSparseMatrix.h
   * id for linear algebra solvers */
  typedef enum { NS_CS_LUSOL, NS_MUMPS, NS_UMFPACK, NS_MKL_PARDISO, NS_SUPERLU, NS_SUPERLU_MT } NumericsSparseLinearSolver;

  /** \enum NumericsSparseTypesNZ
   * value of nz for some matrix storage type */
  typedef enum { NS_CS_CSC = -1, NS_CS_CSR = -2 } NumericsSparseTypesNZ;

  /** \struct NumericsSparseLinearSolverParams NumericsSparseMatrix.h
   * solver-specific parameters*/
  struct NumericsSparseLinearSolverParams
  {
    NumericsSparseLinearSolver solver;

    void* solver_data; /**< solver-specific data (or workspace) */
    freeNSLSP solver_free_hook; /**< solver-specific hook to free solver_data  */

    int* iWork; /**< integer work vector array (internal) */
    int iWorkSize; /**< size of integer work vector array */
    double* dWork;
    int dWorkSize;

    linalg_data_t* linalg_data; /**< data for the linear algebra */
  };

  /**\enum NumericsSparseOrigin NumericsSparseMatrix.h
   * matrix storage types */
  typedef enum { NS_UNKNOWN, NS_TRIPLET, NS_CSC, NS_CSR } NumericsSparseOrigin;

  /** \struct NumericsSparseMatrix NumericsSparseMatrix.h
   * Sparse matrix representation in Numerics. The supported format are:
   * triplet (aka coordinate, COO), CSC (via CSparse) and CSR if MKL is used */
  struct NumericsSparseMatrix
  {
    NumericsSparseLinearSolverParams* linearSolverParams;
                               /**< solver-specific parameters */
    CSparseMatrix* triplet;    /**< triplet format, aka coordinate */
    CSparseMatrix* csc;        /**< csc matrix */
    CSparseMatrix* trans_csc;  /**< transpose of a csc matrix (used by CSparse) */
    CSparseMatrix* csr;        /**< csr matrix, only supported with mkl */
    csi*           diag_indx;  /**< indices for the diagonal terms.
                                    Very useful for the proximal perturbation */
    unsigned       origin;     /**< original format of the matrix */
  };


  /** Initialize the fields of a NumericsSparseMatrix
   * \param A the sparse matrix
   */
  void NM_sparse_null(NumericsSparseMatrix* A);

  /** New and empty NumericsSparseMatrix with correctly initialized fields.
   * \return a pointer on the allocated space.
   */
  NumericsSparseMatrix* newNumericsSparseMatrix(void);


  /** Free allocated space for a NumericsSparseMatrix.
   * \param A a NumericsSparseMatrix
   * \return NULL on success
   */
  NumericsSparseMatrix* freeNumericsSparseMatrix(NumericsSparseMatrix* A);


  /** New and empty NumericsSparseLinearSolverParams.
   * \return a pointer on the allocated space.
   */
  NumericsSparseLinearSolverParams* newNumericsSparseLinearSolverParams(void);

   /** Free a workspace related to a LU factorization
   * \param p the structure to free
   */
  void NM_sparse_free(void *p);

  /** Get the data part of sparse matrix
   * \param A the sparse matrix
   * \return a pointer to the data array
   */
  double* NM_sparse_data(NumericsSparseMatrix* A);

  /** Allocate a CSparse matrix for future copy (as in NM_sparse_copy)
   * \param m the matrix used as model
   * \return an newly allocated matrix
   */
  CSparseMatrix* NM_csparse_alloc_for_copy(const CSparseMatrix* const m);

  /** Get the LU factors for cs_lusol
   * \param p the structure holding the data for the solver
   */
  static inline void* NM_sparse_solver_data(NumericsSparseLinearSolverParams* p)
  {
    return p->solver_data;
  }
  /** Get the workspace for the sparse solver
   * \param p the structure holding the data for the solver
   * \return the (double) workspace
   */
  static inline double* NM_sparse_workspace(NumericsSparseLinearSolverParams* p)

  {
    return p->dWork;
  }

  /** get the number of non-zero (nnz) in a sparse matrix
   * \param A the matrix
   * \return the number of non-zero elements in the matrix
   */
  size_t NM_sparse_nnz(const CSparseMatrix* const A);

  /** Free allocated space for NumericsSparseLinearSolverParams.
   * \param p a NumericsSparseLinearSolverParams
   * \return NULL on success
   */
  NumericsSparseLinearSolverParams* freeNumericsSparseLinearSolverParams(NumericsSparseLinearSolverParams* p);

  /** Check and fix a matrix, if needed
   * \param A the matrix to check, modified if necessary to have ordered indices
   */
  void NM_sparse_fix_csc(CSparseMatrix* A);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
