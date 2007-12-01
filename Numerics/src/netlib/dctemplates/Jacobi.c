/* Jacobi.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b2 = 1.;
static integer c__1 = 1;
static doublereal c_b13 = -1.;

/*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*  Purpose
*  =======
*
*  JACOBI solves the linear system Ax = b using the Jacobi iterative
*  method. The matrix splitting should be accomplished before calling
*  this routine. The diagonal elements of the matrix must be passed into
*  this routine in the first column of matrix WORK.
*
*  Relative error measured: norm( X - X_1 ) / norm( X ).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,4).
*          Workspace for residual, direction vector, etc.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter.
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*  BLAS CALLS:   DAXPY, DCOPY, DNRM2
*  ============================================================ */

int jacobi_(n, b, x, work, ldw, iter, resid, matvec, info)

integer *n, *ldw, *iter, *info;
doublereal *b, *x, *work, *resid;
int (*matvec)();
{
  /* System generated locals */
  integer work_dim1, work_offset, i__1;

  /* Local variables */
  static integer temp;
  extern /* Subroutine */ int matsplit_();
  extern doublereal dnrm2_();
  static integer i;
  extern /* Subroutine */ int dcopy_();
  static integer maxit;
  extern /* Subroutine */ int daxpy_();
  static integer x1, mm;
  static doublereal tol;

  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  work_dim1 = *ldw;
  work_offset = work_dim1 + 1;
  work -= work_offset;
  --x;
  --b;

  /* Function Body */
  *info = 0;

  /*     Test the input parameters. */

  if (*n < 0)
  {
    *info = -1;
  }
  else if (*ldw < max(1, *n))
  {
    *info = -2;
  }
  else if (*iter <= 0)
  {
    *info = -3;
  }
  if (*info != 0)
  {
    return 0;
  }

  maxit = *iter;
  tol = *resid;

  /*     Alias workspace columns. */

  mm = 1;
  x1 = 2;
  temp = 3;

  *iter = 0;

  /*     Form matrix splitting inv(M) and N. */

  matsplit_(&c_b2, &b[1], &work[mm * work_dim1 + 1], ldw, "JACOBI", "SPLIT",
            6L, 5L);

L10:

  /*        Perform Jacobi iteration */

  ++(*iter);

  /*        Save the current approximation to X in X1. */

  dcopy_(n, &x[1], &c__1, &work[x1 * work_dim1 + 1], &c__1);

  /*        Apply iteration; result is updated approximation vector x. */

  dcopy_(n, &b[1], &c__1, &work[temp * work_dim1 + 1], &c__1);
  (*matvec)(&c_b2, &x[1], &c_b2, &work[temp * work_dim1 + 1]);
  i__1 = *n;
  for (i = 1; i <= i__1; ++i)
  {
    x[i] = work[i + mm * work_dim1] * work[i + temp * work_dim1];
    /* L15: */
  }

  /*        Compute error and check for acceptable convergence. */

  daxpy_(n, &c_b13, &x[1], &c__1, &work[x1 * work_dim1 + 1], &c__1);
  *resid = dnrm2_(n, &work[x1 * work_dim1 + 1], &c__1) / dnrm2_(n, &x[1], &
           c__1);

  if (*resid <= tol)
  {
    goto L30;
  }
  if (*iter == maxit)
  {
    goto L20;
  }

  goto L10;

L20:

  /*     Iteration fails */

  *info = 1;
  goto L30;

L30:

  /*     Iteration successful. Reconstruct matrix A. */

  matsplit_(&c_b2, &b[1], &work[mm * work_dim1 + 1], ldw, "JACOBI", "RECON\
STRUCT", 6L, 11L);

  return 0;

  /*     End of JACOBI */

}

