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

#ifndef SparseMatrix_H
#define SparseMatrix_H

/*! \page SparseMatrixPage Sparse Matrix Storage in Numerics

Documentation to be done

*/


/*!\file SparseMatrix.h
  \brief Structure definition and functions related to sparse matrix storage in Numerics
*/

#include "NumericsConfig.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
#include "csparse.h"
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

/*  we use csparse from Timothy Davis

    Timothy Davis,
    Direct Methods for Sparse Linear Systems,
    SIAM, 2006,
    ISBN: 0898716136,
    LC: QA188.D386.

   matrix in compressed row/column or triplet form :
{
int nzmax ;   : maximum number of entries 
int m  ;      : number of rows 
int n ;       : number of columns 
int *p ;      : compressed: row (size m+1) or column (size n+1) pointers; triplet: row indices (size nz) 
int *i ;      : compressed: column or row indices, size nzmax; triplet: column indices (size nz) 
double *x ;   :  numerical values, size nzmax 
int nz ;      : # of entries in triplet matrix;
-1 for compressed columns;
-2 for compressed rows 

} */


#define CSparseMatrix struct cs_sparse

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  typedef enum { NS_CS_LUSOL, NS_MUMPS } NumericsSparseLinearSolver;

  typedef struct
  {
    NumericsSparseLinearSolver solver;

    int* iparam;
    double iSize;
    double* dparam;
    double dSize;

#ifdef HAVE_MPI
    int mpi_com_init;
    MPI_Comm mpi_com;
#endif
    void* solver_data;

    int* iWork; /**< integer work vector array (internal) */
    int iWorkSize; /**< size of integer work vector array */
    double* dWork;
    int dWorkSize;
  } NumericsSparseLinearSolverParams;

  typedef struct
  {
    NumericsSparseLinearSolverParams* linearSolverParams;

    CSparseMatrix* triplet;
    CSparseMatrix* csc;
    CSparseMatrix* trans_csc;
  } NumericsSparseMatrix;


  /** Add an entry to a triplet matrix only if the absolute value is
   * greater than DBL_EPSILON.
   * \param T the CSparseMatrix
   * \param i row index
   * \param j column index
   * \param x the value
   * \return integer value : 1 if the absolute value is less than
   * DBL_EPSILON, otherwise the return value of cs_entry.
   */
  csi cs_zentry(CSparseMatrix *T, csi i, csi j, double x);

  /** Create dense matrix from a CSparseMatrix.
   * \param A the CSparseMatrix
   * \return a pointer on A->m * A->n allocated storage
   */
  double* cs_dense(CSparseMatrix *A);

  /** Matrix vector multiplication : y = alpha*A*x+beta*y
   * \param[in] alpha matrix coefficient
   * \param[in] A the CSparseMatrix
   * \param[in] x pointer on a dense vector of size A->n
   * \param[in] beta vector coefficient
   * \param[in, out] y pointer on a dense vector of size A->n
   * \return 0 if A x or y is NULL else 1
   */
  int cs_aaxpy(const double alpha, const cs *A, const double *x,
               const double beta, double *y);

  /** Free space allocated for a SparseMatrix. note : cs_spfree also
   *  free the cs_struct this fails when the struct is allocated on
   *  the stack.
   * \param A a CSparseMatrix
   * \return NULL on success
  */
  CSparseMatrix* cs_spfree_on_stack(CSparseMatrix* A);

  /** New and empty NumericsSparseLinearSolverParams.
   * \return a pointer on the allocated space.
   */
  NumericsSparseLinearSolverParams* newNumericsSparseLinearSolverParams(void);

  /** New and empty NumericsSparseMatrix with correctly initialized fields.
   * \return a pointer on the allocated space.
   */
  NumericsSparseMatrix* newNumericsSparseMatrix(void);

  /** Free allocated space for NumericsSparseLinearSolverParams.
   * \param p a NumericsSparseLinearSolverParams
   * \return NULL on sucess
   */
  NumericsSparseLinearSolverParams* freeNumericsSparseLinearSolverParams(NumericsSparseLinearSolverParams* p);

  /** Free allocated space for a NumericsSparseMatrix.
   * \param A a NumericsSparseMatrix
   * \return NULL on sucess
   */
  NumericsSparseMatrix* freeNumericsSparseMatrix(NumericsSparseMatrix* A);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
