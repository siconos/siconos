/* BiCGSTAB.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b5 = -1.;
static doublereal c_b6 = 1.;
static doublereal c_b25 = 0.;

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
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with
*  preconditioning.
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
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors R and S shared the same workspace.
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
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must
*          remain unchanged.
*          The solution is over-written on vector b.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common block.
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
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO < BREAKTOL: RHO and RTLD have become
*                                       orthogonal.
*                  -11: OMEGA < BREAKTOL: S and T have become
*                                         orthogonal relative to T'*T.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  BLAS CALLS: DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*  ==============================================================
*/

int bicgstab_(n, b, x, work, ldw, iter, resid, matvec, psolve, info)

integer *n, *ldw, *iter, *info;
doublereal *b, *x, *work, *resid;
int (*matvec)(), (*psolve)();
{
  /* System generated locals */
  integer work_dim1, work_offset;
  doublereal d__1;

  /* Local variables */
  static doublereal beta;
  extern doublereal ddot_();
  static integer phat, shat;
  extern doublereal getbreak_();
  static integer rtld;
  static doublereal omegatol, bnrm2;
  extern doublereal dnrm2_();
  static integer p, r, s, t;
  static doublereal alpha;
  static integer v;
  extern /* Subroutine */ int dscal_();
  static doublereal omega;
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

  /* Executable Statements */
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
  p = 3;
  v = 4;
  t = 5;
  phat = 6;
  shat = 7;
  s = 1;

  /*     Set parameter tolerances. */

  rhotol = getbreak_();
  omegatol = getbreak_();

  /*     Set initial residual. */

  dcopy_(n, &b[1], &c__1, &work[r * work_dim1 + 1], &c__1);
  if (dnrm2_(n, &x[1], &c__1) != 0.)
  {
    (*matvec)(&c_b5, &x[1], &c_b6, &work[r * work_dim1 + 1]);
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

  /*     Perform BiConjugate Gradient Stabilized iteration. */

  ++(*iter);

  rho = ddot_(n, &work[rtld * work_dim1 + 1], &c__1, &work[r * work_dim1 +
              1], &c__1);
  if (abs(rho) < rhotol)
  {
    goto L25;
  }

  /*        Compute vector P. */

  if (*iter > 1)
  {
    beta = rho / rho1 * (alpha / omega);
    d__1 = -omega;
    daxpy_(n, &d__1, &work[v * work_dim1 + 1], &c__1, &work[p * work_dim1
           + 1], &c__1);
    dscal_(n, &beta, &work[p * work_dim1 + 1], &c__1);
    daxpy_(n, &c_b6, &work[r * work_dim1 + 1], &c__1, &work[p * work_dim1
           + 1], &c__1);
  }
  else
  {
    dcopy_(n, &work[r * work_dim1 + 1], &c__1, &work[p * work_dim1 + 1], &
           c__1);
  }

  /*        Compute direction adjusting vector PHAT and scalar ALPHA. */

  (*psolve)(&work[phat * work_dim1 + 1], &work[p * work_dim1 + 1]);
  (*matvec)(&c_b6, &work[phat * work_dim1 + 1], &c_b25, &work[v * work_dim1
            + 1]);
  alpha = rho / ddot_(n, &work[rtld * work_dim1 + 1], &c__1, &work[v *
                      work_dim1 + 1], &c__1);

  /*        Early check for tolerance. */

  d__1 = -alpha;
  daxpy_(n, &d__1, &work[v * work_dim1 + 1], &c__1, &work[r * work_dim1 + 1]
         , &c__1);
  dcopy_(n, &work[r * work_dim1 + 1], &c__1, &work[s * work_dim1 + 1], &
         c__1);
  if (dnrm2_(n, &work[s * work_dim1 + 1], &c__1) <= tol)
  {
    daxpy_(n, &alpha, &work[phat * work_dim1 + 1], &c__1, &x[1], &c__1);
    *resid = dnrm2_(n, &work[s * work_dim1 + 1], &c__1) / bnrm2;
    goto L30;
  }
  else
  {

    /*           Compute stabilizer vector SHAT and scalar OMEGA. */

    (*psolve)(&work[shat * work_dim1 + 1], &work[s * work_dim1 + 1]);
    (*matvec)(&c_b6, &work[shat * work_dim1 + 1], &c_b25, &work[t *
              work_dim1 + 1]);
    omega = ddot_(n, &work[t * work_dim1 + 1], &c__1, &work[s * work_dim1
                  + 1], &c__1) / ddot_(n, &work[t * work_dim1 + 1], &c__1, &
                                       work[t * work_dim1 + 1], &c__1);

    /*           Compute new solution approximation vector X. */

    daxpy_(n, &alpha, &work[phat * work_dim1 + 1], &c__1, &x[1], &c__1);
    daxpy_(n, &omega, &work[shat * work_dim1 + 1], &c__1, &x[1], &c__1);

    /*           Compute residual R, check for tolerance. */

    d__1 = -omega;
    daxpy_(n, &d__1, &work[t * work_dim1 + 1], &c__1, &work[r * work_dim1
           + 1], &c__1);
    *resid = dnrm2_(n, &work[r * work_dim1 + 1], &c__1) / bnrm2;
    if (*resid <= tol)
    {
      goto L30;
    }
    if (*iter == maxit)
    {
      goto L20;
    }
  }

  if (abs(omega) < omegatol)
  {
    goto L25;
  }
  else
  {
    rho1 = rho;
    goto L10;
  }

L20:

  /*     Iteration fails. */

  *info = 1;
  return 0;

L25:

  /*     Set breakdown flag. */

  if (abs(rho) < rhotol)
  {
    *info = -10;
  }
  else if (abs(omega) < omegatol)
  {
    *info = -11;
  }
  return 0;

L30:

  /*     Iteration successful; return. */

  return 0;

  /*     End of BICGSTAB */

}
