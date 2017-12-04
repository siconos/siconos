/* Siconos, Copyright INRIA 2005-2016.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
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

#include <stdio.h>

/*!\file SparseMatrix.h
  \brief Structure definition and functions related to sparse matrix storage in Numerics
*/

#include "SiconosConfig.h"

/* Compile-time assertion: users of SparseMatrix.h must have CS_LONG
 * set if and only if SICONOS_INT64 is also set. If it is unset, it is
 * set here. */
#ifdef SICONOS_INT64
#ifndef CS_LONG
#define CS_LONG
#endif
#else
#ifdef CS_LONG
#error "CS_LONG (set) does not correspond with SICONOS_INT64 (unset)"
#endif
#endif

/* From cs.h */
#ifndef CS_INT
#ifdef CS_LONG
#define CS_INT long
#else
#define CS_INT int
#endif
#endif

/* Treat CXSparse structs as opaque types.  Users may #include "cs.h"
 * to use them outside Siconos. */
struct cs_dl_sparse;
struct cs_di_sparse;
struct cs_dl_symbolic;
struct cs_di_symbolic;
struct cs_dl_numeric;
struct cs_di_numeric;
typedef struct cs_dl_symbolic cs_dls;
typedef struct cs_di_symbolic cs_dis;
typedef struct cs_dl_numeric cs_dln;
typedef struct cs_di_numeric cs_din;

#ifdef SICONOS_INT64 // SWIG gives syntax error for CS_NAME(_sparse)
typedef struct cs_dl_sparse CSparseMatrix;
#ifndef css
#define css cs_dls
#endif
#ifndef csn
#define csn cs_dln
#endif
#else
typedef struct cs_di_sparse CSparseMatrix;
#ifndef css
#define css cs_dis
#endif
#ifndef csn
#define csn cs_din
#endif
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

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
csi nzmax ;   : maximum number of entries
csi m  ;      : number of rows
csi n ;       : number of columns
csi *p ;      : compressed: row (size m+1) or column (size n+1) pointers; triplet: row indices (size nz)
csi *i ;      : compressed: column or row indices, size nzmax; triplet: column indices (size nz)
double *x ;   :  numerical values, size nzmax
csi nz ;      : # of entries in triplet matrix;
-1 for compressed columns;
-2 for compressed rows

}

csi is either int64_t or int32_t and this is controlled at compile time*/


#define NS_UNKNOWN_ERR(func, orig) \
fprintf(stderr, #func ": unknown origin %d for sparse matrix\n", orig);


#define NS_NROW_CSR(mat) mat->n
#define NS_NCOL_CSR(mat) mat->m

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** \struct cs_lu_factors SparseMatrix.h
   * Information used and produced by CSparse for an LU factorization*/
  typedef struct {
    CS_INT n;       /**< size of linear system */
    css* S;      /**< symbolic analysis */
    csn* N;      /**< numerics factorization */
  } cs_lu_factors;

  typedef void (*freeNSLSP)(void* p);


  /** Add an entry to a triplet matrix only if the absolute value is
   * greater than DBL_EPSILON.
   * \param T the sparse matrix
   * \param i row index
   * \param j column index
   * \param x the value
   * \return integer value : 1 if the absolute value is less than
   * DBL_EPSILON, otherwise the return value of cs_entry.
   */
  CS_INT cs_zentry(CSparseMatrix *T, CS_INT i, CS_INT j, double x);

  /** Check if the given triplet matrix is properly constructed (col and row indices are correct)
   * \param T the sparse matrix to check
   * \return 0 if the matrix is fine, 1 otherwise
   * */
  int cs_check_triplet(CSparseMatrix *T);

  /** Check if the given triplet matrix is properly constructed (col and row indices are correct)
   * \param T the sparse matrix to check
   * \return 0 if the matrix is fine, 1 otherwise
   * */
  int cs_check_csc(CSparseMatrix *T);

  /** Create dense matrix from a CSparseMatrix.
   * \param A the sparse matrix
   * \return a pointer on A->m * A->n allocated storage
   */
  double* cs_dense(CSparseMatrix *A);

  /** Matrix vector multiplication : y = alpha*A*x+beta*y
   * \param[in] alpha matrix coefficient
   * \param[in] A the sparse matrix
   * \param[in] x pointer on a dense vector of size A->n
   * \param[in] beta vector coefficient
   * \param[in, out] y pointer on a dense vector of size A->n
   * \return 0 if A x or y is NULL else 1
   */
  int cs_aaxpy(const double alpha, const CSparseMatrix *A, const double *x,
               const double beta, double *y);

  /** Free space allocated for a SparseMatrix. note : cs_spfree also
   *  free the cs_struct this fails when the struct is allocated on
   *  the stack.
   * \param A the sparse matrix
   * \return NULL on success
  */
  CSparseMatrix* cs_spfree_on_stack(CSparseMatrix* A);







  /** reuse a LU factorization (stored in the cs_lu_A) to solve a linear system Ax = b
   * \param cs_lu_A contains the LU factors of A, permutation information
   * \param x workspace
   * \param[in,out] b on input RHS of the linear system; on output the solution
   * \return 0 if failed, 1 otherwise*/
  CS_INT cs_solve (cs_lu_factors* cs_lu_A, double* x, double *b);

  /** compute a LU factorization of A and store it in a workspace
   * \param order control if ordering is used
   * \param A the sparse matrix
   * \param tol the tolerance
   * \param cs_lu_A the parameter structure that eventually holds the factors
   * \return 1 if the factorization was successful, 1 otherwise
   */
  int cs_lu_factorization(CS_INT order, const CSparseMatrix *A, double tol, cs_lu_factors * cs_lu_A);

  /** Free a workspace related to a LU factorization
   * \param cs_lu_A the structure to free
   */
  void cs_sparse_free(cs_lu_factors* cs_lu_A);

  /** print a matrix to a text file
   * \param A matrix to print
   * \param brief if positive, print only a portion of the matrix
   * \param file file descriptor*/
  int cs_printInFile(const CSparseMatrix *A, int brief, FILE* file);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
