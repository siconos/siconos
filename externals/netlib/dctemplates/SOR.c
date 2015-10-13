/* SOR.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b5 = -1.;
static doublereal c_b6 = 1.;

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
*  SOR solves the linear system Ax = b using the Successive
*  Over-Relaxation iterative method.
*  The matrix splitting is formed by copying the strict upper
*  triangular portion of A onto matrix N, stored in WORK. Matrix M
*  is the lower triangular portion of A.
*  On exit, matrix A and right hand side b are reset to their
*  original form.
*
*  Relative error measured: norm( X - X_1 ) / norm( X ).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (input/workspace) DOUBLE PRECISION array, dimension (N*(N+3)).
*          The relaxation parameter, OMEGA, should be input in WORK(1).
*          The amount of workspace can be significantly reduced (to 2*N)
*          by customizing the matrix-vector product and backsolve.
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
*          norm( x - x_1 ) / norm( x ).
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
*  BACKSOLVE (external subroutine)
*          The user must provide a subroutine to perform the
*          linear system solve
*
*               x := M*x,
*
*          where x is a vector and M is a lower triangular matrix.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL BACKSOLVE( N, M, LDM, X )
*
*          The matrix is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -4: Relaxation parameter OMEGA not in interval (0,2).
*
*  BLAS CALLS:   DAXPY, DCOPY, DNRM2
*  ==========================================================
*/

int sor_(n, b, x, work, ldw, iter, resid, matvec, backsolve, info)
integer *n, *ldw, *iter, *info;
doublereal *b, *x, *work, *resid;
int (*matvec)(), (*backsolve)();
{
  /* System generated locals */
  integer work_dim1, work_offset;

  /* Local variables */
  static integer temp;
  extern /* Subroutine */ int matsplit_();
  static doublereal bnrm2;
  extern doublereal dnrm2_();
  static doublereal omega;
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
  else if (work[work_dim1 + 1] <= 0. || work[work_dim1 + 1] >= 2.)
  {
    *info = -4;
  }
  if (*info != 0)
  {
    return 0;
  }

  maxit = *iter;
  tol = *resid;

  /*     Alias workspace columns. */

  x1 = 1;
  temp = 2;
  mm = 3;

  /*     Set relaxation parameter. */

  omega = work[work_dim1 + 1];
  if (omega == 0.)
  {
    omega = 1.;
  }

  /*     Compute initial residual for ( convergence criteria ). */

  dcopy_(n, &b[1], &c__1, &work[x1 * work_dim1 + 1], &c__1);
  if (dnrm2_(n, &x[1], &c__1) != 0.)
  {
    (*matvec)(&c_b5, &x[1], &c_b6, &work[x1 * work_dim1 + 1]);
    if (dnrm2_(n, &work[x1 * work_dim1 + 1], &c__1) < tol)
    {
      goto L30;
    }
  }
  bnrm2 = dnrm2_(n, &b[1], &c__1);
  if (bnrm2 == 0.)
  {
    bnrm2 = 1.;
  }

  /*     Matrix A is set to N. WORK(1:N,1:N) is set to MM. */

  matsplit_(&omega, &b[1], &work[mm * work_dim1 + 1], ldw, "SOR", "SPLIT",
            3L, 5L);

  *iter = 0;

L10:

  /*     Perform SOR iteration */

  ++(*iter);

  /*        Save the current approximation to X in X1, */

  dcopy_(n, &x[1], &c__1, &work[x1 * work_dim1 + 1], &c__1);

  /*        Apply iteration; result is updated approximation vector x */

  dcopy_(n, &b[1], &c__1, &work[temp * work_dim1 + 1], &c__1);
  (*matvec)(&c_b6, &x[1], &c_b6, &work[temp * work_dim1 + 1]);
  dcopy_(n, &work[temp * work_dim1 + 1], &c__1, &x[1], &c__1);
  (*backsolve)(n, &work[mm * work_dim1 + 1], ldw, &x[1]);

  /*        Compute error and check for acceptable convergence. */

  daxpy_(n, &c_b5, &x[1], &c__1, &work[x1 * work_dim1 + 1], &c__1);
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

L30:

  /*     Iteration successful; restore A and B to original form, */
  /*     compute residual norm, and return */

  matsplit_(&omega, &b[1], &work[mm * work_dim1 + 1], ldw, "SOR", "RESTORE",
            3L, 7L);

  return 0;

  /*     End of SOR */

}
