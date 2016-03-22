#ifndef NumericsSparseMatrix_H
#define NumericsSparseMatrix_H

/*!


 */

#include "SiconosConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif


  typedef enum { NS_CS_LUSOL, NS_MUMPS, NS_UMFPACK, NS_PARDISO } NumericsSparseLinearSolver;

  /** \struct NumericsSparseLinearSolverParams SparseMatrix.h
   * solver-specific parameters*/
  typedef struct
  {
    NumericsSparseLinearSolver solver;

    int* iparam;
    double iSize;
    double* dparam;
    double dSize;

    void* solver_data; /**< solver-specific data (or workspace) */
    freeNSLSP solver_free_hook; /**< solver-specific hook to free solver_data  */

    int* iWork; /**< integer work vector array (internal) */
    int iWorkSize; /**< size of integer work vector array */
    double* dWork;
    int dWorkSize;
  } NumericsSparseLinearSolverParams;

  typedef enum { NS_UNKNOWN, NS_TRIPLET, NS_CSC, NS_CSR } NumericsSparseOrigin;

  /** \struct NumericsSparseMatrix SparseMatrix.h
   * Sparse matrix representation in Numerics. The supported format are:
   * triplet (aka coordinate, COO), CSC (via CSparse) and CSR if MKL is used */
  typedef struct
  {
    NumericsSparseLinearSolverParams* linearSolverParams;
                               /**< solver-specific parameters */
    CSparseMatrix* triplet;    /**< triplet format, aka coordinate */
    CSparseMatrix* csc;        /**< csc matrix */
    CSparseMatrix* trans_csc;  /**< transpose of a csc matrix (used by CSparse) */
    CSparseMatrix* csr;        /**< csr matrix, only supported with mkl */
    unsigned       origin;     /**< original format of the matrix */
  } NumericsSparseMatrix;


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


  /** Get the LU factors for cs_lusol
   * \param p the structure holding the data for the solver
   */

  static inline cs_lu_factors* NM_csparse_lu_factors(NumericsSparseLinearSolverParams* p)
  {
    return (cs_lu_factors*)p->solver_data;
  }
  /** Get the workspace for the sparse solver
   * \param p the structure holding the data for the solver
   */
  static inline double* NM_csparse_workspace(NumericsSparseLinearSolverParams* p)

  {
    return p->dWork;
  }


  /** Free allocated space for NumericsSparseLinearSolverParams.
   * \param p a NumericsSparseLinearSolverParams
   * \return NULL on success
   */
  NumericsSparseLinearSolverParams* freeNumericsSparseLinearSolverParams(NumericsSparseLinearSolverParams* p);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
