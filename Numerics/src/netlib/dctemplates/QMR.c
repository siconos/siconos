/* QMR.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b5 = -1.;
static doublereal c_b6 = 0.;
static doublereal c_b44 = 1.;

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
*  PSOLVEQ  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. As QMR uses left
*          and right preconditioning and the preconditioners are in
*          common, we must specify in the call which to use. Vector b
*          must remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVEQ( X, B, 'LEFT' )
*
*         The preconditioner is passed into the routine in a common block
*
*
*  PSOLVETRANSQ  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M'*x = b,
*
*          where x and y are vectors, and M' is the tranpose of a
*          matrix M. As QMR uses left and right preconditioning and
*          the preconditioners are in common, we must specify in the
*          call which to use. Vector b must remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVETRANSQ( X, B, 'LEFT' )
*
*         The preconditioner is passed into the routine in a common block.
*
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
*                  -10: RHO   < BREAKTOL: RHO and RTLD have become
*                                         orthogonal.
*                  -11: BETA  < BREAKTOL: EPS too small in relation to DELTA.
*                                         Convergence has stalled.
*                  -12: GAMMA < BREAKTOL: THETA too large.
*                                         Convergence has stalled.
*                  -13: DELTA < BREAKTOL: Y and Z have become
*                                         orthogonal.
*                  -14: EPS   < BREAKTOL: Q and PTLD have become
*                                         orthogonal.
*                  -15: XI    < BREAKTOL: Z too small. Convergence has stalled.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  BLAS CALLS:   DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*  ============================================================ */

int qmr_(n, b, x, work, ldw, iter, resid, matvec,
         matvectrans, psolveq, psolvetransq, info)

integer *n, *ldw, *iter, *info;
doublereal *b, *x, *work, *resid;
int (*matvec)(), (*matvectrans)(), (*psolveq)(), (*psolvetransq)();
{
  /* System generated locals */
  integer work_dim1, work_offset;
  doublereal d__1, d__2;

  /* Builtin functions */
  double sqrt();

  /* Local variables */
  static doublereal beta;
  extern doublereal ddot_();
  static integer ptld;
  extern doublereal getbreak_();
  static integer vtld, wtld, ytld, ztld;
  static doublereal gammatol, deltatol, bnrm2;
  extern doublereal dnrm2_();
  static integer d, p, q, r, s;
  static doublereal gamma;
  static integer v, w, y, z;
  static doublereal delta;
  extern /* Subroutine */ int dscal_();
  static doublereal theta;
  extern /* Subroutine */ int dcopy_();
  static integer maxit;
  static doublereal c1;
  extern /* Subroutine */ int daxpy_();
  static doublereal xitol, gamma1, theta1, xi, epstol, rhotol, eta, eps,
         rho, tol, betatol, rho1;

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
  d = 2;
  p = 3;
  ptld = 4;
  q = 5;
  s = 6;
  v = 7;
  vtld = 8;
  w = 9;
  wtld = 9;
  y = 10;
  ytld = 10;
  z = 11;
  ztld = 11;

  /*     Set breakdown tolerances. */

  rhotol = getbreak_();
  betatol = getbreak_();
  gammatol = getbreak_();
  deltatol = getbreak_();
  epstol = getbreak_();
  xitol = getbreak_();

  /*     Set initial residual. */

  dcopy_(n, &b[1], &c__1, &work[r * work_dim1 + 1], &c__1);
  if (dnrm2_(n, &x[1], &c__1) != 0.)
  {
    (*matvec)(&c_b5, &x[1], &c_b44, &work[r * work_dim1 + 1]);
    if (dnrm2_(n, &work[r * work_dim1 + 1], &c__1) < tol)
    {
      goto L30;
    }
  }

  bnrm2 = dnrm2_(n, &b[1], &c__1);
  if (bnrm2 == 0.)
  {
    bnrm2 = 1.;
  }

  dcopy_(n, &work[r * work_dim1 + 1], &c__1, &work[vtld * work_dim1 + 1], &
         c__1);
  (*psolveq)(&work[y * work_dim1 + 1], &work[vtld * work_dim1 + 1], "LEFT",
             4L);
  rho = dnrm2_(n, &work[y * work_dim1 + 1], &c__1);

  dcopy_(n, &work[r * work_dim1 + 1], &c__1, &work[wtld * work_dim1 + 1], &
         c__1);
  (*psolvetransq)(&work[z * work_dim1 + 1], &work[wtld * work_dim1 + 1],
                  "RIGHT", 5L);
  xi = dnrm2_(n, &work[z * work_dim1 + 1], &c__1);

  gamma = 1.;
  eta = -1.;
  theta = 0.;

  *iter = 0;

L10:

  /*     Perform Preconditioned QMR iteration. */

  ++(*iter);

  if (abs(rho) < rhotol || abs(xi) < xitol)
  {
    goto L25;
  }

  dcopy_(n, &work[vtld * work_dim1 + 1], &c__1, &work[v * work_dim1 + 1], &
         c__1);
  d__1 = 1. / rho;
  dscal_(n, &d__1, &work[v * work_dim1 + 1], &c__1);
  d__1 = 1. / rho;
  dscal_(n, &d__1, &work[y * work_dim1 + 1], &c__1);

  dcopy_(n, &work[wtld * work_dim1 + 1], &c__1, &work[w * work_dim1 + 1], &
         c__1);
  d__1 = 1. / xi;
  dscal_(n, &d__1, &work[w * work_dim1 + 1], &c__1);
  d__1 = 1. / xi;
  dscal_(n, &d__1, &work[z * work_dim1 + 1], &c__1);

  delta = ddot_(n, &work[z * work_dim1 + 1], &c__1, &work[y * work_dim1 + 1]
                , &c__1);
  if (abs(delta) < deltatol)
  {
    goto L25;
  }

  (*psolveq)(&work[ytld * work_dim1 + 1], &work[y * work_dim1 + 1], "RIGHT",
             5L);
  (*psolvetransq)(&work[ztld * work_dim1 + 1], &work[z * work_dim1 + 1],
                  "LEFT", 4L);

  if (*iter > 1)
  {
    c1 = -(xi * delta / eps);
    daxpy_(n, &c1, &work[p * work_dim1 + 1], &c__1, &work[ytld *
           work_dim1 + 1], &c__1);
    dcopy_(n, &work[ytld * work_dim1 + 1], &c__1, &work[p * work_dim1 + 1]
           , &c__1);
    d__1 = -(rho * delta / eps);
    daxpy_(n, &d__1, &work[q * work_dim1 + 1], &c__1, &work[ztld *
           work_dim1 + 1], &c__1);
    dcopy_(n, &work[ztld * work_dim1 + 1], &c__1, &work[q * work_dim1 + 1]
           , &c__1);
  }
  else
  {
    dcopy_(n, &work[ytld * work_dim1 + 1], &c__1, &work[p * work_dim1 + 1]
           , &c__1);
    dcopy_(n, &work[ztld * work_dim1 + 1], &c__1, &work[q * work_dim1 + 1]
           , &c__1);
  }

  (*matvec)(&c_b44, &work[p * work_dim1 + 1], &c_b6, &work[ptld * work_dim1
            + 1]);

  eps = ddot_(n, &work[q * work_dim1 + 1], &c__1, &work[ptld * work_dim1 +
              1], &c__1);
  if (abs(eps) < epstol)
  {
    goto L25;
  }

  beta = eps / delta;
  if (abs(beta) < betatol)
  {
    goto L25;
  }

  dcopy_(n, &work[ptld * work_dim1 + 1], &c__1, &work[vtld * work_dim1 + 1],
         &c__1);
  d__1 = -beta;
  daxpy_(n, &d__1, &work[v * work_dim1 + 1], &c__1, &work[vtld * work_dim1
         + 1], &c__1);
  (*psolveq)(&work[y * work_dim1 + 1], &work[vtld * work_dim1 + 1], "LEFT",
             4L);

  rho1 = rho;
  rho = dnrm2_(n, &work[y * work_dim1 + 1], &c__1);

  dcopy_(n, &work[w * work_dim1 + 1], &c__1, &work[wtld * work_dim1 + 1], &
         c__1);
  d__1 = -beta;
  (*matvectrans)(&c_b44, &work[q * work_dim1 + 1], &d__1, &work[wtld *
                 work_dim1 + 1]);
  (*psolvetransq)(&work[z * work_dim1 + 1], &work[wtld * work_dim1 + 1],
                  "RIGHT", 5L);

  xi = dnrm2_(n, &work[z * work_dim1 + 1], &c__1);

  gamma1 = gamma;
  theta1 = theta;

  theta = rho / (gamma1 * abs(beta));
  /* Computing 2nd power */
  d__1 = theta;
  gamma = 1. / sqrt(d__1 * d__1 + 1.);
  if (abs(gamma) < gammatol)
  {
    goto L25;
  }

  /* Computing 2nd power */
  d__1 = gamma;
  /* Computing 2nd power */
  d__2 = gamma1;
  eta = -eta * rho1 * (d__1 * d__1) / (beta * (d__2 * d__2));

  if (*iter > 1)
  {
    /* Computing 2nd power */
    d__2 = theta1 * gamma;
    d__1 = d__2 * d__2;
    dscal_(n, &d__1, &work[d * work_dim1 + 1], &c__1);
    daxpy_(n, &eta, &work[p * work_dim1 + 1], &c__1, &work[d * work_dim1
           + 1], &c__1);
    /* Computing 2nd power */
    d__2 = theta1 * gamma;
    d__1 = d__2 * d__2;
    dscal_(n, &d__1, &work[s * work_dim1 + 1], &c__1);
    daxpy_(n, &eta, &work[ptld * work_dim1 + 1], &c__1, &work[s *
           work_dim1 + 1], &c__1);
  }
  else
  {
    dcopy_(n, &work[p * work_dim1 + 1], &c__1, &work[d * work_dim1 + 1], &
           c__1);
    dscal_(n, &eta, &work[d * work_dim1 + 1], &c__1);
    dcopy_(n, &work[ptld * work_dim1 + 1], &c__1, &work[s * work_dim1 + 1]
           , &c__1);
    dscal_(n, &eta, &work[s * work_dim1 + 1], &c__1);
  }

  /*        Compute current solution vector x. */

  daxpy_(n, &c_b44, &work[d * work_dim1 + 1], &c__1, &x[1], &c__1);

  /*        Compute residual vector rk, find norm, */
  /*        then check for tolerance. */

  daxpy_(n, &c_b5, &work[s * work_dim1 + 1], &c__1, &work[r * work_dim1 + 1]
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

  goto L10;

L20:

  /*     Iteration fails. */

  *info = 1;
  return 0;

L25:

  /*     Method breakdown. */

  if (abs(rho) < rhotol)
  {
    *info = -10;
  }
  else if (abs(beta) < betatol)
  {
    *info = -11;
  }
  else if (abs(gamma) < gammatol)
  {
    *info = -12;
  }
  else if (abs(delta) < deltatol)
  {
    *info = -13;
  }
  else if (abs(eps) < epstol)
  {
    *info = -14;
  }
  else if (abs(xi) < xitol)
  {
    *info = -15;
  }

  return 0;

L30:

  /*     Iteration successful; return. */

  return 0;

  /*     End of QMR */

}
