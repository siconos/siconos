/* BiCG.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b5 = -1.;
static doublereal c_b6 = 0.;
static doublereal c_b28 = 1.;

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
*  BiCG solves the linear system Ax = b using the
*  BiConjugate Gradient iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
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
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,6).
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
*  MATVECTRANS  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A'*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A' is the tranpose of a matrix A. Vector x must remain
*          unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVECTRANS( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*         The preconditioner is passed into the routine in a common block
*
*  PSOLVETRANS  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M'*x = b,
*
*          where x and y are vectors, and M' is the tranpose of a
*          matrix M. Vector b must remain unchanged. The solution
*          is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVETRANS( X, B )
*
*         The preconditioner is passed into the routine in a common block
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
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO < BREAKTOL: RHO and RTLD have become
*                                       orthogonal.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  BLAS CALLS:   DAXPY, DCOPY, DDOT, DNRM2
*  ============================================================ */

int bicg_(n, b, x, work, ldw, iter, resid, matvec, matvectrans,
          psolve, psolvetrans, info)
integer *n, *ldw, *iter, *info;
doublereal *b, *x, *work, *resid;
int (*matvec)(), (*matvectrans)(), (*psolve)(), (*psolvetrans)();
{
  /* System generated locals */
  integer work_dim1, work_offset;
  doublereal d__1;

  /* Local variables */
  static doublereal beta;
  extern doublereal ddot_();
  static integer ptld, rtld, qtld;
  extern doublereal getbreak_();
  static integer ztld;
  static doublereal bnrm2;
  extern doublereal dnrm2_();
  static integer p, q, r;
  static doublereal alpha;
  static integer z;
  extern /* Subroutine */ int dcopy_();
  static integer maxit;
  extern /* Subroutine */ int daxpy_();
  static doublereal rhotol, rho, tol, rho1;

  /* Parameter adjustments */
  work_dim1 = *ldw;
  work_offset = work_dim1 + 1;
  work -= work_offset;
  --x;
  --b;

  /*  Executable Statements */

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

  r = 1;
  rtld = 2;
  z = 3;
  ztld = 4;
  p = 5;
  ptld = 6;
  q = 3;
  qtld = 4;

  /*     Set breakdown parameters. */

  rhotol = getbreak_();

  /*     Set initial residual. */

  dcopy_(n, &b[1], &c__1, &work[r * work_dim1 + 1], &c__1);
  if (dnrm2_(n, &x[1], &c__1) != 0.)
  {
    (*matvec)(&c_b5, &x[1], &c_b28, &work[r * work_dim1 + 1]);
    if (dnrm2_(n, &work[r * work_dim1 + 1], &c__1) <= tol)
    {
      goto L30;
    }
  }
  dcopy_(n, &work[r * work_dim1 + 1], &c__1, &work[rtld * work_dim1 + 1], &
         c__1);
  bnrm2 = dnrm2_(n, &b[1], &c__1);
  if (bnrm2 == 0.)
  {
    bnrm2 = 1.;
  }

  *iter = 0;

L10:

  /*     Perform BiConjugate Gradient iteration. */

  ++(*iter);

  /*        Compute direction vectors PK and PTLD. */

  (*psolve)(&work[z * work_dim1 + 1], &work[r * work_dim1 + 1]);
  (*psolvetrans)(&work[ztld * work_dim1 + 1], &work[rtld * work_dim1 + 1]);

  rho = ddot_(n, &work[z * work_dim1 + 1], &c__1, &work[rtld * work_dim1 +
              1], &c__1);
  if (abs(rho) < rhotol)
  {
    goto L25;
  }

  if (*iter > 1)
  {
    beta = rho / rho1;
    daxpy_(n, &beta, &work[p * work_dim1 + 1], &c__1, &work[z * work_dim1
           + 1], &c__1);
    daxpy_(n, &beta, &work[ptld * work_dim1 + 1], &c__1, &work[ztld *
           work_dim1 + 1], &c__1);
    dcopy_(n, &work[z * work_dim1 + 1], &c__1, &work[p * work_dim1 + 1], &
           c__1);
    dcopy_(n, &work[ztld * work_dim1 + 1], &c__1, &work[ptld * work_dim1
           + 1], &c__1);
  }
  else
  {
    dcopy_(n, &work[z * work_dim1 + 1], &c__1, &work[p * work_dim1 + 1], &
           c__1);
    dcopy_(n, &work[ztld * work_dim1 + 1], &c__1, &work[ptld * work_dim1
           + 1], &c__1);
  }

  (*matvec)(&c_b28, &work[p * work_dim1 + 1], &c_b6, &work[q * work_dim1 +
            1]);
  (*matvectrans)(&c_b28, &work[ptld * work_dim1 + 1], &c_b6, &work[qtld *
                 work_dim1 + 1]);
  alpha = rho / ddot_(n, &work[ptld * work_dim1 + 1], &c__1, &work[q *
                      work_dim1 + 1], &c__1);

  /*        Compute current solution vector x. */

  daxpy_(n, &alpha, &work[p * work_dim1 + 1], &c__1, &x[1], &c__1);

  /*        Compute residual vector rk, find norm, */
  /*        then check for tolerance. */

  d__1 = -alpha;
  daxpy_(n, &d__1, &work[q * work_dim1 + 1], &c__1, &work[r * work_dim1 + 1]
         , &c__1);

  *resid = dnrm2_(n, &work[r * work_dim1 + 1], &c__1) / bnrm2;
  if (*resid <= tol)
  {
    goto L30;
  }
  if (*iter == maxit)
  {
    goto L20;
  }

  d__1 = -alpha;
  daxpy_(n, &d__1, &work[qtld * work_dim1 + 1], &c__1, &work[rtld *
         work_dim1 + 1], &c__1);
  rho1 = rho;

  goto L10;

L20:

  /*     Iteration fails. */

  *info = 1;
  return 0;

L25:

  /*     Set breakdown flag. */

  *info = -10;
  return 0;

L30:

  /*     Iteration successful; return. */

  return 0;

  /*     End of BICG */

}
