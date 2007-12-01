/* dblas.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b91 = 1.;

/* Subroutine */
int daxpy_(n, da, dx, incx, dy, incy)
integer *n;
doublereal *da, *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
  /* System generated locals */
  integer i__1;

  /* Local variables */
  static integer i, m, ix, iy, mp1;


  /*     constant times a vector plus a vector. */
  /*     uses unrolled loops for increments equal to one. */
  /*     jack dongarra, linpack, 3/11/78. */


  /* Parameter adjustments */
  --dy;
  --dx;

  /* Function Body */
  if (*n <= 0)
  {
    return 0;
  }
  if (*da == 0.)
  {
    return 0;
  }
  if (*incx == 1 && *incy == 1)
  {
    goto L20;
  }

  /*        code for unequal increments or equal increments */
  /*          not equal to 1 */

  ix = 1;
  iy = 1;
  if (*incx < 0)
  {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0)
  {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i = 1; i <= i__1; ++i)
  {
    dy[iy] += *da * dx[ix];
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;

  /*        code for both increments equal to 1 */


  /*        clean-up loop */

L20:
  m = *n % 4;
  if (m == 0)
  {
    goto L40;
  }
  i__1 = m;
  for (i = 1; i <= i__1; ++i)
  {
    dy[i] += *da * dx[i];
    /* L30: */
  }
  if (*n < 4)
  {
    return 0;
  }
L40:
  mp1 = m + 1;
  i__1 = *n;
  for (i = mp1; i <= i__1; i += 4)
  {
    dy[i] += *da * dx[i];
    dy[i + 1] += *da * dx[i + 1];
    dy[i + 2] += *da * dx[i + 2];
    dy[i + 3] += *da * dx[i + 3];
    /* L50: */
  }
  return 0;
} /* daxpy_ */

/* Subroutine */ int dcopy_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
  /* System generated locals */
  integer i__1;

  /* Local variables */
  static integer i, m, ix, iy, mp1;


  /*     copies a vector, x, to a vector, y. */
  /*     uses unrolled loops for increments equal to one. */
  /*     jack dongarra, linpack, 3/11/78. */


  /* Parameter adjustments */
  --dy;
  --dx;

  /* Function Body */
  if (*n <= 0)
  {
    return 0;
  }
  if (*incx == 1 && *incy == 1)
  {
    goto L20;
  }

  /*        code for unequal increments or equal increments */
  /*          not equal to 1 */

  ix = 1;
  iy = 1;
  if (*incx < 0)
  {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0)
  {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i = 1; i <= i__1; ++i)
  {
    dy[iy] = dx[ix];
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;

  /*        code for both increments equal to 1 */


  /*        clean-up loop */

L20:
  m = *n % 7;
  if (m == 0)
  {
    goto L40;
  }
  i__1 = m;
  for (i = 1; i <= i__1; ++i)
  {
    dy[i] = dx[i];
    /* L30: */
  }
  if (*n < 7)
  {
    return 0;
  }
L40:
  mp1 = m + 1;
  i__1 = *n;
  for (i = mp1; i <= i__1; i += 7)
  {
    dy[i] = dx[i];
    dy[i + 1] = dx[i + 1];
    dy[i + 2] = dx[i + 2];
    dy[i + 3] = dx[i + 3];
    dy[i + 4] = dx[i + 4];
    dy[i + 5] = dx[i + 5];
    dy[i + 6] = dx[i + 6];
    /* L50: */
  }
  return 0;
} /* dcopy_ */

doublereal ddot_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
  /* System generated locals */
  integer i__1;
  doublereal ret_val;

  /* Local variables */
  static integer i, m;
  static doublereal dtemp;
  static integer ix, iy, mp1;


  /*     forms the dot product of two vectors. */
  /*     uses unrolled loops for increments equal to one. */
  /*     jack dongarra, linpack, 3/11/78. */


  /* Parameter adjustments */
  --dy;
  --dx;

  /* Function Body */
  ret_val = 0.;
  dtemp = 0.;
  if (*n <= 0)
  {
    return ret_val;
  }
  if (*incx == 1 && *incy == 1)
  {
    goto L20;
  }

  /*        code for unequal increments or equal increments */
  /*          not equal to 1 */

  ix = 1;
  iy = 1;
  if (*incx < 0)
  {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0)
  {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i = 1; i <= i__1; ++i)
  {
    dtemp += dx[ix] * dy[iy];
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  ret_val = dtemp;
  return ret_val;

  /*        code for both increments equal to 1 */


  /*        clean-up loop */

L20:
  m = *n % 5;
  if (m == 0)
  {
    goto L40;
  }
  i__1 = m;
  for (i = 1; i <= i__1; ++i)
  {
    dtemp += dx[i] * dy[i];
    /* L30: */
  }
  if (*n < 5)
  {
    goto L60;
  }
L40:
  mp1 = m + 1;
  i__1 = *n;
  for (i = mp1; i <= i__1; i += 5)
  {
    dtemp = dtemp + dx[i] * dy[i] + dx[i + 1] * dy[i + 1] + dx[i + 2] *
            dy[i + 2] + dx[i + 3] * dy[i + 3] + dx[i + 4] * dy[i + 4];
    /* L50: */
  }
L60:
  ret_val = dtemp;
  return ret_val;
} /* ddot_ */

/* Subroutine */ int dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb,
                            beta, c, ldc, transa_len, transb_len)
char *transa, *transb;
integer *m, *n, *k;
doublereal *alpha, *a;
integer *lda;
doublereal *b;
integer *ldb;
doublereal *beta, *c;
integer *ldc;
ftnlen transa_len;
ftnlen transb_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
          i__3;

  /* Local variables */
  static integer info;
  static logical nota, notb;
  static doublereal temp;
  static integer i, j, l, ncola;
  extern logical lsame_();
  static integer nrowa, nrowb;
  extern /* Subroutine */ int xerbla_();

  /*     .. Scalar Arguments .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DGEMM  performs one of the matrix-matrix operations */

  /*     C := alpha*op( A )*op( B ) + beta*C, */

  /*  where  op( X ) is one of */

  /*     op( X ) = X   or   op( X ) = X', */

  /*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
  */
  /*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  */

  /*  Parameters */
  /*  ========== */

  /*  TRANSA - CHARACTER*1. */
  /*           On entry, TRANSA specifies the form of op( A ) to be used in
  */
  /*           the matrix multiplication as follows: */

  /*              TRANSA = 'N' or 'n',  op( A ) = A. */

  /*              TRANSA = 'T' or 't',  op( A ) = A'. */

  /*              TRANSA = 'C' or 'c',  op( A ) = A'. */

  /*           Unchanged on exit. */

  /*  TRANSB - CHARACTER*1. */
  /*           On entry, TRANSB specifies the form of op( B ) to be used in
  */
  /*           the matrix multiplication as follows: */

  /*              TRANSB = 'N' or 'n',  op( B ) = B. */

  /*              TRANSB = 'T' or 't',  op( B ) = B'. */

  /*              TRANSB = 'C' or 'c',  op( B ) = B'. */

  /*           Unchanged on exit. */

  /*  M      - INTEGER. */
  /*           On entry,  M  specifies  the number  of rows  of the  matrix
  */
  /*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
  */
  /*           Unchanged on exit. */

  /*  N      - INTEGER. */
  /*           On entry,  N  specifies the number  of columns of the matrix
  */
  /*           op( B ) and the number of columns of the matrix C. N must be
  */
  /*           at least zero. */
  /*           Unchanged on exit. */

  /*  K      - INTEGER. */
  /*           On entry,  K  specifies  the number of columns of the matrix
  */
  /*           op( A ) and the number of rows of the matrix op( B ). K must
  */
  /*           be at least  zero. */
  /*           Unchanged on exit. */

  /*  ALPHA  - DOUBLE PRECISION. */
  /*           On entry, ALPHA specifies the scalar alpha. */
  /*           Unchanged on exit. */

  /*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
  */
  /*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
  /*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
  */
  /*           part of the array  A  must contain the matrix  A,  otherwise
  */
  /*           the leading  k by m  part of the array  A  must contain  the
  */
  /*           matrix A. */
  /*           Unchanged on exit. */

  /*  LDA    - INTEGER. */
  /*           On entry, LDA specifies the first dimension of A as declared
  */
  /*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
  */
  /*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
  */
  /*           least  max( 1, k ). */
  /*           Unchanged on exit. */

  /*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
  */
  /*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
  /*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
  */
  /*           part of the array  B  must contain the matrix  B,  otherwise
  */
  /*           the leading  n by k  part of the array  B  must contain  the
  */
  /*           matrix B. */
  /*           Unchanged on exit. */

  /*  LDB    - INTEGER. */
  /*           On entry, LDB specifies the first dimension of B as declared
  */
  /*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
  */
  /*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
  */
  /*           least  max( 1, n ). */
  /*           Unchanged on exit. */

  /*  BETA   - DOUBLE PRECISION. */
  /*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
  */
  /*           supplied as zero then C need not be set on input. */
  /*           Unchanged on exit. */

  /*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ). */
  /*           Before entry, the leading  m by n  part of the array  C must
  */
  /*           contain the matrix  C,  except when  beta  is zero, in which
  */
  /*           case C need not be set on entry. */
  /*           On exit, the array  C  is overwritten by the  m by n  matrix
  */
  /*           ( alpha*op( A )*op( B ) + beta*C ). */

  /*  LDC    - INTEGER. */
  /*           On entry, LDC specifies the first dimension of C as declared
  */
  /*           in  the  calling  (sub)  program.   LDC  must  be  at  least
  */
  /*           max( 1, m ). */
  /*           Unchanged on exit. */


  /*  Level 3 Blas routine. */

  /*  -- Written on 8-February-1989. */
  /*     Jack Dongarra, Argonne National Laboratory. */
  /*     Iain Duff, AERE Harwell. */
  /*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
  /*     Sven Hammarling, Numerical Algorithms Group Ltd. */


  /*     .. External Functions .. */
  /*     .. External Subroutines .. */
  /*     .. Intrinsic Functions .. */
  /*     .. Local Scalars .. */
  /*     .. Parameters .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
  */
  /*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
  */
  /*     and  columns of  A  and the  number of  rows  of  B  respectively.
  */

  /* Parameter adjustments */
  c_dim1 = *ldc;
  c_offset = c_dim1 + 1;
  c -= c_offset;
  b_dim1 = *ldb;
  b_offset = b_dim1 + 1;
  b -= b_offset;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  nota = lsame_(transa, "N", 1L, 1L);
  notb = lsame_(transb, "N", 1L, 1L);
  if (nota)
  {
    nrowa = *m;
    ncola = *k;
  }
  else
  {
    nrowa = *k;
    ncola = *m;
  }
  if (notb)
  {
    nrowb = *k;
  }
  else
  {
    nrowb = *n;
  }

  /*     Test the input parameters. */

  info = 0;
  if (! nota && ! lsame_(transa, "C", 1L, 1L) && ! lsame_(transa, "T", 1L,
      1L))
  {
    info = 1;
  }
  else if (! notb && ! lsame_(transb, "C", 1L, 1L) && ! lsame_(transb,
           "T", 1L, 1L))
  {
    info = 2;
  }
  else if (*m < 0)
  {
    info = 3;
  }
  else if (*n < 0)
  {
    info = 4;
  }
  else if (*k < 0)
  {
    info = 5;
  }
  else if (*lda < max(1, nrowa))
  {
    info = 8;
  }
  else if (*ldb < max(1, nrowb))
  {
    info = 10;
  }
  else if (*ldc < max(1, *m))
  {
    info = 13;
  }
  if (info != 0)
  {
    xerbla_("DGEMM ", &info, 6L);
    return 0;
  }

  /*     Quick return if possible. */

  if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.)
  {
    return 0;
  }

  /*     And if  alpha.eq.zero. */

  if (*alpha == 0.)
  {
    if (*beta == 0.)
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        i__2 = *m;
        for (i = 1; i <= i__2; ++i)
        {
          c[i + j * c_dim1] = 0.;
          /* L10: */
        }
        /* L20: */
      }
    }
    else
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        i__2 = *m;
        for (i = 1; i <= i__2; ++i)
        {
          c[i + j * c_dim1] = *beta * c[i + j * c_dim1];
          /* L30: */
        }
        /* L40: */
      }
    }
    return 0;
  }

  /*     Start the operations. */

  if (notb)
  {
    if (nota)
    {

      /*           Form  C := alpha*A*B + beta*C. */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        if (*beta == 0.)
        {
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = 0.;
            /* L50: */
          }
        }
        else if (*beta != 1.)
        {
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = *beta * c[i + j * c_dim1];
            /* L60: */
          }
        }
        i__2 = *k;
        for (l = 1; l <= i__2; ++l)
        {
          if (b[l + j * b_dim1] != 0.)
          {
            temp = *alpha * b[l + j * b_dim1];
            i__3 = *m;
            for (i = 1; i <= i__3; ++i)
            {
              c[i + j * c_dim1] += temp * a[i + l * a_dim1];
              /* L70: */
            }
          }
          /* L80: */
        }
        /* L90: */
      }
    }
    else
    {

      /*           Form  C := alpha*A'*B + beta*C */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        i__2 = *m;
        for (i = 1; i <= i__2; ++i)
        {
          temp = 0.;
          i__3 = *k;
          for (l = 1; l <= i__3; ++l)
          {
            temp += a[l + i * a_dim1] * b[l + j * b_dim1];
            /* L100: */
          }
          if (*beta == 0.)
          {
            c[i + j * c_dim1] = *alpha * temp;
          }
          else
          {
            c[i + j * c_dim1] = *alpha * temp + *beta * c[i + j *
                                c_dim1];
          }
          /* L110: */
        }
        /* L120: */
      }
    }
  }
  else
  {
    if (nota)
    {

      /*           Form  C := alpha*A*B' + beta*C */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        if (*beta == 0.)
        {
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = 0.;
            /* L130: */
          }
        }
        else if (*beta != 1.)
        {
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = *beta * c[i + j * c_dim1];
            /* L140: */
          }
        }
        i__2 = *k;
        for (l = 1; l <= i__2; ++l)
        {
          if (b[j + l * b_dim1] != 0.)
          {
            temp = *alpha * b[j + l * b_dim1];
            i__3 = *m;
            for (i = 1; i <= i__3; ++i)
            {
              c[i + j * c_dim1] += temp * a[i + l * a_dim1];
              /* L150: */
            }
          }
          /* L160: */
        }
        /* L170: */
      }
    }
    else
    {

      /*           Form  C := alpha*A'*B' + beta*C */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        i__2 = *m;
        for (i = 1; i <= i__2; ++i)
        {
          temp = 0.;
          i__3 = *k;
          for (l = 1; l <= i__3; ++l)
          {
            temp += a[l + i * a_dim1] * b[j + l * b_dim1];
            /* L180: */
          }
          if (*beta == 0.)
          {
            c[i + j * c_dim1] = *alpha * temp;
          }
          else
          {
            c[i + j * c_dim1] = *alpha * temp + *beta * c[i + j *
                                c_dim1];
          }
          /* L190: */
        }
        /* L200: */
      }
    }
  }

  return 0;

  /*     End of DGEMM . */

} /* dgemm_ */

/* Subroutine */ int dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y,
                            incy, trans_len)
char *trans;
integer *m, *n;
doublereal *alpha, *a;
integer *lda;
doublereal *x;
integer *incx;
doublereal *beta, *y;
integer *incy;
ftnlen trans_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;

  /* Local variables */
  static integer info;
  static doublereal temp;
  static integer lenx, leny, i, j;
  extern logical lsame_();
  static integer ix, iy, jx, jy, kx, ky;
  extern /* Subroutine */ int xerbla_();

  /*     .. Scalar Arguments .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DGEMV  performs one of the matrix-vector operations */

  /*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, */

  /*  where alpha and beta are scalars, x and y are vectors and A is an */
  /*  m by n matrix. */

  /*  Parameters */
  /*  ========== */

  /*  TRANS  - CHARACTER*1. */
  /*           On entry, TRANS specifies the operation to be performed as */
  /*           follows: */

  /*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */

  /*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y. */

  /*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y. */

  /*           Unchanged on exit. */

  /*  M      - INTEGER. */
  /*           On entry, M specifies the number of rows of the matrix A. */
  /*           M must be at least zero. */
  /*           Unchanged on exit. */

  /*  N      - INTEGER. */
  /*           On entry, N specifies the number of columns of the matrix A.
  */
  /*           N must be at least zero. */
  /*           Unchanged on exit. */

  /*  ALPHA  - DOUBLE PRECISION. */
  /*           On entry, ALPHA specifies the scalar alpha. */
  /*           Unchanged on exit. */

  /*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
  /*           Before entry, the leading m by n part of the array A must */
  /*           contain the matrix of coefficients. */
  /*           Unchanged on exit. */

  /*  LDA    - INTEGER. */
  /*           On entry, LDA specifies the first dimension of A as declared
  */
  /*           in the calling (sub) program. LDA must be at least */
  /*           max( 1, m ). */
  /*           Unchanged on exit. */

  /*  X      - DOUBLE PRECISION array of DIMENSION at least */
  /*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
  /*           and at least */
  /*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
  /*           Before entry, the incremented array X must contain the */
  /*           vector x. */
  /*           Unchanged on exit. */

  /*  INCX   - INTEGER. */
  /*           On entry, INCX specifies the increment for the elements of */
  /*           X. INCX must not be zero. */
  /*           Unchanged on exit. */

  /*  BETA   - DOUBLE PRECISION. */
  /*           On entry, BETA specifies the scalar beta. When BETA is */
  /*           supplied as zero then Y need not be set on input. */
  /*           Unchanged on exit. */

  /*  Y      - DOUBLE PRECISION array of DIMENSION at least */
  /*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
  /*           and at least */
  /*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
  /*           Before entry with BETA non-zero, the incremented array Y */
  /*           must contain the vector y. On exit, Y is overwritten by the
  */
  /*           updated vector y. */

  /*  INCY   - INTEGER. */
  /*           On entry, INCY specifies the increment for the elements of */
  /*           Y. INCY must not be zero. */
  /*           Unchanged on exit. */


  /*  Level 2 Blas routine. */

  /*  -- Written on 22-October-1986. */
  /*     Jack Dongarra, Argonne National Lab. */
  /*     Jeremy Du Croz, Nag Central Office. */
  /*     Sven Hammarling, Nag Central Office. */
  /*     Richard Hanson, Sandia National Labs. */


  /*     .. Parameters .. */
  /*     .. Local Scalars .. */
  /*     .. External Functions .. */
  /*     .. External Subroutines .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  --y;
  --x;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  info = 0;
  if (! lsame_(trans, "N", 1L, 1L) && ! lsame_(trans, "T", 1L, 1L) && !
      lsame_(trans, "C", 1L, 1L))
  {
    info = 1;
  }
  else if (*m < 0)
  {
    info = 2;
  }
  else if (*n < 0)
  {
    info = 3;
  }
  else if (*lda < max(1, *m))
  {
    info = 6;
  }
  else if (*incx == 0)
  {
    info = 8;
  }
  else if (*incy == 0)
  {
    info = 11;
  }
  if (info != 0)
  {
    xerbla_("DGEMV ", &info, 6L);
    return 0;
  }

  /*     Quick return if possible. */

  if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.)
  {
    return 0;
  }

  /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  */
  /*     up the start points in  X  and  Y. */

  if (lsame_(trans, "N", 1L, 1L))
  {
    lenx = *n;
    leny = *m;
  }
  else
  {
    lenx = *m;
    leny = *n;
  }
  if (*incx > 0)
  {
    kx = 1;
  }
  else
  {
    kx = 1 - (lenx - 1) * *incx;
  }
  if (*incy > 0)
  {
    ky = 1;
  }
  else
  {
    ky = 1 - (leny - 1) * *incy;
  }

  /*     Start the operations. In this version the elements of A are */
  /*     accessed sequentially with one pass through A. */

  /*     First form  y := beta*y. */

  if (*beta != 1.)
  {
    if (*incy == 1)
    {
      if (*beta == 0.)
      {
        i__1 = leny;
        for (i = 1; i <= i__1; ++i)
        {
          y[i] = 0.;
          /* L10: */
        }
      }
      else
      {
        i__1 = leny;
        for (i = 1; i <= i__1; ++i)
        {
          y[i] = *beta * y[i];
          /* L20: */
        }
      }
    }
    else
    {
      iy = ky;
      if (*beta == 0.)
      {
        i__1 = leny;
        for (i = 1; i <= i__1; ++i)
        {
          y[iy] = 0.;
          iy += *incy;
          /* L30: */
        }
      }
      else
      {
        i__1 = leny;
        for (i = 1; i <= i__1; ++i)
        {
          y[iy] = *beta * y[iy];
          iy += *incy;
          /* L40: */
        }
      }
    }
  }
  if (*alpha == 0.)
  {
    return 0;
  }
  if (lsame_(trans, "N", 1L, 1L))
  {

    /*        Form  y := alpha*A*x + y. */

    jx = kx;
    if (*incy == 1)
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        if (x[jx] != 0.)
        {
          temp = *alpha * x[jx];
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            y[i] += temp * a[i + j * a_dim1];
            /* L50: */
          }
        }
        jx += *incx;
        /* L60: */
      }
    }
    else
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        if (x[jx] != 0.)
        {
          temp = *alpha * x[jx];
          iy = ky;
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            y[iy] += temp * a[i + j * a_dim1];
            iy += *incy;
            /* L70: */
          }
        }
        jx += *incx;
        /* L80: */
      }
    }
  }
  else
  {

    /*        Form  y := alpha*A'*x + y. */

    jy = ky;
    if (*incx == 1)
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        temp = 0.;
        i__2 = *m;
        for (i = 1; i <= i__2; ++i)
        {
          temp += a[i + j * a_dim1] * x[i];
          /* L90: */
        }
        y[jy] += *alpha * temp;
        jy += *incy;
        /* L100: */
      }
    }
    else
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        temp = 0.;
        ix = kx;
        i__2 = *m;
        for (i = 1; i <= i__2; ++i)
        {
          temp += a[i + j * a_dim1] * x[ix];
          ix += *incx;
          /* L110: */
        }
        y[jy] += *alpha * temp;
        jy += *incy;
        /* L120: */
      }
    }
  }

  return 0;

  /*     End of DGEMV . */

} /* dgemv_ */

/* Subroutine */ int dger_(m, n, alpha, x, incx, y, incy, a, lda)
integer *m, *n;
doublereal *alpha, *x;
integer *incx;
doublereal *y;
integer *incy;
doublereal *a;
integer *lda;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;

  /* Local variables */
  static integer info;
  static doublereal temp;
  static integer i, j, ix, jy, kx;
  extern /* Subroutine */ int xerbla_();

  /*     .. Scalar Arguments .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DGER   performs the rank 1 operation */

  /*     A := alpha*x*y' + A, */

  /*  where alpha is a scalar, x is an m element vector, y is an n element
  */
  /*  vector and A is an m by n matrix. */

  /*  Parameters */
  /*  ========== */

  /*  M      - INTEGER. */
  /*           On entry, M specifies the number of rows of the matrix A. */
  /*           M must be at least zero. */
  /*           Unchanged on exit. */

  /*  N      - INTEGER. */
  /*           On entry, N specifies the number of columns of the matrix A.
  */
  /*           N must be at least zero. */
  /*           Unchanged on exit. */

  /*  ALPHA  - DOUBLE PRECISION. */
  /*           On entry, ALPHA specifies the scalar alpha. */
  /*           Unchanged on exit. */

  /*  X      - DOUBLE PRECISION array of dimension at least */
  /*           ( 1 + ( m - 1 )*abs( INCX ) ). */
  /*           Before entry, the incremented array X must contain the m */
  /*           element vector x. */
  /*           Unchanged on exit. */

  /*  INCX   - INTEGER. */
  /*           On entry, INCX specifies the increment for the elements of */
  /*           X. INCX must not be zero. */
  /*           Unchanged on exit. */

  /*  Y      - DOUBLE PRECISION array of dimension at least */
  /*           ( 1 + ( n - 1 )*abs( INCY ) ). */
  /*           Before entry, the incremented array Y must contain the n */
  /*           element vector y. */
  /*           Unchanged on exit. */

  /*  INCY   - INTEGER. */
  /*           On entry, INCY specifies the increment for the elements of */
  /*           Y. INCY must not be zero. */
  /*           Unchanged on exit. */

  /*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
  /*           Before entry, the leading m by n part of the array A must */
  /*           contain the matrix of coefficients. On exit, A is */
  /*           overwritten by the updated matrix. */

  /*  LDA    - INTEGER. */
  /*           On entry, LDA specifies the first dimension of A as declared
  */
  /*           in the calling (sub) program. LDA must be at least */
  /*           max( 1, m ). */
  /*           Unchanged on exit. */


  /*  Level 2 Blas routine. */

  /*  -- Written on 22-October-1986. */
  /*     Jack Dongarra, Argonne National Lab. */
  /*     Jeremy Du Croz, Nag Central Office. */
  /*     Sven Hammarling, Nag Central Office. */
  /*     Richard Hanson, Sandia National Labs. */


  /*     .. Parameters .. */
  /*     .. Local Scalars .. */
  /*     .. External Subroutines .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;
  --y;
  --x;

  /* Function Body */
  info = 0;
  if (*m < 0)
  {
    info = 1;
  }
  else if (*n < 0)
  {
    info = 2;
  }
  else if (*incx == 0)
  {
    info = 5;
  }
  else if (*incy == 0)
  {
    info = 7;
  }
  else if (*lda < max(1, *m))
  {
    info = 9;
  }
  if (info != 0)
  {
    xerbla_("DGER  ", &info, 6L);
    return 0;
  }

  /*     Quick return if possible. */

  if (*m == 0 || *n == 0 || *alpha == 0.)
  {
    return 0;
  }

  /*     Start the operations. In this version the elements of A are */
  /*     accessed sequentially with one pass through A. */

  if (*incy > 0)
  {
    jy = 1;
  }
  else
  {
    jy = 1 - (*n - 1) * *incy;
  }
  if (*incx == 1)
  {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j)
    {
      if (y[jy] != 0.)
      {
        temp = *alpha * y[jy];
        i__2 = *m;
        for (i = 1; i <= i__2; ++i)
        {
          a[i + j * a_dim1] += x[i] * temp;
          /* L10: */
        }
      }
      jy += *incy;
      /* L20: */
    }
  }
  else
  {
    if (*incx > 0)
    {
      kx = 1;
    }
    else
    {
      kx = 1 - (*m - 1) * *incx;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j)
    {
      if (y[jy] != 0.)
      {
        temp = *alpha * y[jy];
        ix = kx;
        i__2 = *m;
        for (i = 1; i <= i__2; ++i)
        {
          a[i + j * a_dim1] += x[ix] * temp;
          ix += *incx;
          /* L30: */
        }
      }
      jy += *incy;
      /* L40: */
    }
  }

  return 0;

  /*     End of DGER  . */

} /* dger_ */

doublereal dnrm2_(n, dx, incx)
integer *n;
doublereal *dx;
integer *incx;
{
  /* Initialized data */

  static doublereal zero = 0.;
  static doublereal one = 1.;
  static doublereal cutlo = 8.232e-11;
  static doublereal cuthi = 1.304e19;

  /* Format strings */
  static char fmt_30[] = "";
  static char fmt_50[] = "";
  static char fmt_70[] = "";
  static char fmt_110[] = "";

  /* System generated locals */
  integer i__1;
  doublereal ret_val, d__1;

  /* Builtin functions */
  double sqrt();

  /* Local variables */
  static doublereal xmax;
  static integer next, i, j, ix;
  static doublereal hitest, sum;

  /* Assigned format variables */
  char *next_fmt;

  /* Parameter adjustments */
  --dx;

  /* Function Body */

  /*     euclidean norm of the n-vector stored in dx() with storage */
  /*     increment incx . */
  /*     if    n .le. 0 return with result = 0. */
  /*     if n .ge. 1 then incx must be .ge. 1 */

  /*           c.l.lawson, 1978 jan 08 */
  /*     modified to correct failure to update ix, 1/25/92. */
  /*     modified 3/93 to return if incx .le. 0. */

  /*     four phase method     using two built-in constants that are */
  /*     hopefully applicable to all machines. */
  /*         cutlo = maximum of  dsqrt(u/eps)  over all known machines. */
  /*         cuthi = minimum of  dsqrt(v)      over all known machines. */
  /*     where */
  /*         eps = smallest no. such that eps + 1. .gt. 1. */
  /*         u   = smallest positive no.   (underflow limit) */
  /*         v   = largest  no.            (overflow  limit) */

  /*     brief outline of algorithm.. */

  /*     phase 1    scans zero components. */
  /*     move to phase 2 when a component is nonzero and .le. cutlo */
  /*     move to phase 3 when a component is .gt. cutlo */
  /*     move to phase 4 when a component is .ge. cuthi/m */
  /*     where m = n for x() real and m = 2*n for complex. */

  /*     values for cutlo and cuthi.. */
  /*     from the environmental parameters listed in the imsl converter */
  /*     document the limiting values are as follows.. */
  /*     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
  */
  /*                   univac and dec at 2**(-103) */
  /*                   thus cutlo = 2**(-51) = 4.44089e-16 */
  /*     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec. */
  /*                   thus cuthi = 2**(63.5) = 1.30438e19 */
  /*     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec. */
  /*                   thus cutlo = 2**(-33.5) = 8.23181d-11 */
  /*     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19 */
  /*     data cutlo, cuthi / 8.232d-11,  1.304d19 / */
  /*     data cutlo, cuthi / 4.441e-16,  1.304e19 / */

  if (*n > 0 && *incx > 0)
  {
    goto L10;
  }
  ret_val = zero;
  goto L300;

L10:
  next = 0;
  next_fmt = fmt_30;
  sum = zero;
  i = 1;
  ix = 1;
  /*                                                 begin main loop */
L20:
  switch ((int)next)
  {
  case 0:
    goto L30;
  case 1:
    goto L50;
  case 2:
    goto L70;
  case 3:
    goto L110;
  }
L30:
  if ((d__1 = dx[i], abs(d__1)) > cutlo)
  {
    goto L85;
  }
  next = 1;
  next_fmt = fmt_50;
  xmax = zero;

  /*                        phase 1.  sum is zero */

L50:
  if (dx[i] == zero)
  {
    goto L200;
  }
  if ((d__1 = dx[i], abs(d__1)) > cutlo)
  {
    goto L85;
  }

  /*                                prepare for phase 2. */
  next = 2;
  next_fmt = fmt_70;
  goto L105;

  /*                                prepare for phase 4. */

L100:
  ix = j;
  next = 3;
  next_fmt = fmt_110;
  sum = sum / dx[i] / dx[i];
L105:
  xmax = (d__1 = dx[i], abs(d__1));
  goto L115;

  /*                   phase 2.  sum is small. */
  /*                             scale to avoid destructive underflow. */

L70:
  if ((d__1 = dx[i], abs(d__1)) > cutlo)
  {
    goto L75;
  }

  /*                     common code for phases 2 and 4. */
  /*                     in phase 4 sum is large.  scale to avoid overflow.
  */

L110:
  if ((d__1 = dx[i], abs(d__1)) <= xmax)
  {
    goto L115;
  }
  /* Computing 2nd power */
  d__1 = xmax / dx[i];
  sum = one + sum * (d__1 * d__1);
  xmax = (d__1 = dx[i], abs(d__1));
  goto L200;

L115:
  /* Computing 2nd power */
  d__1 = dx[i] / xmax;
  sum += d__1 * d__1;
  goto L200;


  /*                  prepare for phase 3. */

L75:
  sum = sum * xmax * xmax;


  /*     for real or d.p. set hitest = cuthi/n */
  /*     for complex      set hitest = cuthi/(2*n) */

L85:
  hitest = cuthi / (real)(*n);

  /*                   phase 3.  sum is mid-range.  no scaling. */

  i__1 = *n;
  for (j = ix; j <= i__1; ++j)
  {
    if ((d__1 = dx[i], abs(d__1)) >= hitest)
    {
      goto L100;
    }
    /* Computing 2nd power */
    d__1 = dx[i];
    sum += d__1 * d__1;
    i += *incx;
    /* L95: */
  }
  ret_val = sqrt(sum);
  goto L300;

L200:
  ++ix;
  i += *incx;
  if (ix <= *n)
  {
    goto L20;
  }

  /*              end of main loop. */

  /*              compute square root and adjust for scaling. */

  ret_val = xmax * sqrt(sum);
L300:
  return ret_val;
} /* dnrm2_ */

/* Subroutine */ int drotg_(da, db, c, s)
doublereal *da, *db, *c, *s;
{
  /* System generated locals */
  doublereal d__1, d__2;

  /* Builtin functions */
  double sqrt(), d_sign();

  /* Local variables */
  static doublereal r, scale, z, roe;


  /*     construct givens plane rotation. */
  /*     jack dongarra, linpack, 3/11/78. */


  roe = *db;
  if (abs(*da) > abs(*db))
  {
    roe = *da;
  }
  scale = abs(*da) + abs(*db);
  if (scale != 0.)
  {
    goto L10;
  }
  *c = 1.;
  *s = 0.;
  r = 0.;
  z = 0.;
  goto L20;
L10:
  /* Computing 2nd power */
  d__1 = *da / scale;
  /* Computing 2nd power */
  d__2 = *db / scale;
  r = scale * sqrt(d__1 * d__1 + d__2 * d__2);
  r = d_sign(&c_b91, &roe) * r;
  *c = *da / r;
  *s = *db / r;
  z = 1.;
  if (abs(*da) > abs(*db))
  {
    z = *s;
  }
  if (abs(*db) >= abs(*da) && *c != 0.)
  {
    z = 1. / *c;
  }
L20:
  *da = r;
  *db = z;
  return 0;
} /* drotg_ */

/* Subroutine */ int drot_(n, dx, incx, dy, incy, c, s)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
doublereal *c, *s;
{
  /* System generated locals */
  integer i__1;

  /* Local variables */
  static integer i;
  static doublereal dtemp;
  static integer ix, iy;


  /*     applies a plane rotation. */
  /*     jack dongarra, linpack, 3/11/78. */


  /* Parameter adjustments */
  --dy;
  --dx;

  /* Function Body */
  if (*n <= 0)
  {
    return 0;
  }
  if (*incx == 1 && *incy == 1)
  {
    goto L20;
  }

  /*       code for unequal increments or equal increments not equal */
  /*         to 1 */

  ix = 1;
  iy = 1;
  if (*incx < 0)
  {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0)
  {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i = 1; i <= i__1; ++i)
  {
    dtemp = *c * dx[ix] + *s * dy[iy];
    dy[iy] = *c * dy[iy] - *s * dx[ix];
    dx[ix] = dtemp;
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;

  /*       code for both increments equal to 1 */

L20:
  i__1 = *n;
  for (i = 1; i <= i__1; ++i)
  {
    dtemp = *c * dx[i] + *s * dy[i];
    dy[i] = *c * dy[i] - *s * dx[i];
    dx[i] = dtemp;
    /* L30: */
  }
  return 0;
} /* drot_ */

/* Subroutine */ int dscal_(n, da, dx, incx)
integer *n;
doublereal *da, *dx;
integer *incx;
{
  /* System generated locals */
  integer i__1, i__2;

  /* Local variables */
  static integer i, m, nincx, mp1;


  /*     scales a vector by a constant. */
  /*     uses unrolled loops for increment equal to one. */
  /*     jack dongarra, linpack, 3/11/78. */
  /*     modified 3/93 to return if incx .le. 0. */


  /* Parameter adjustments */
  --dx;

  /* Function Body */
  if (*n <= 0 || *incx <= 0)
  {
    return 0;
  }
  if (*incx == 1)
  {
    goto L20;
  }

  /*        code for increment not equal to 1 */

  nincx = *n * *incx;
  i__1 = nincx;
  i__2 = *incx;
  for (i = 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2)
  {
    dx[i] = *da * dx[i];
    /* L10: */
  }
  return 0;

  /*        code for increment equal to 1 */


  /*        clean-up loop */

L20:
  m = *n % 5;
  if (m == 0)
  {
    goto L40;
  }
  i__2 = m;
  for (i = 1; i <= i__2; ++i)
  {
    dx[i] = *da * dx[i];
    /* L30: */
  }
  if (*n < 5)
  {
    return 0;
  }
L40:
  mp1 = m + 1;
  i__2 = *n;
  for (i = mp1; i <= i__2; i += 5)
  {
    dx[i] = *da * dx[i];
    dx[i + 1] = *da * dx[i + 1];
    dx[i + 2] = *da * dx[i + 2];
    dx[i + 3] = *da * dx[i + 3];
    dx[i + 4] = *da * dx[i + 4];
    /* L50: */
  }
  return 0;
} /* dscal_ */

/* Subroutine */ int dswap_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
  /* System generated locals */
  integer i__1;

  /* Local variables */
  static integer i, m;
  static doublereal dtemp;
  static integer ix, iy, mp1;


  /*     interchanges two vectors. */
  /*     uses unrolled loops for increments equal one. */
  /*     jack dongarra, linpack, 3/11/78. */


  /* Parameter adjustments */
  --dy;
  --dx;

  /* Function Body */
  if (*n <= 0)
  {
    return 0;
  }
  if (*incx == 1 && *incy == 1)
  {
    goto L20;
  }

  /*       code for unequal increments or equal increments not equal */
  /*         to 1 */

  ix = 1;
  iy = 1;
  if (*incx < 0)
  {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0)
  {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i = 1; i <= i__1; ++i)
  {
    dtemp = dx[ix];
    dx[ix] = dy[iy];
    dy[iy] = dtemp;
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;

  /*       code for both increments equal to 1 */


  /*       clean-up loop */

L20:
  m = *n % 3;
  if (m == 0)
  {
    goto L40;
  }
  i__1 = m;
  for (i = 1; i <= i__1; ++i)
  {
    dtemp = dx[i];
    dx[i] = dy[i];
    dy[i] = dtemp;
    /* L30: */
  }
  if (*n < 3)
  {
    return 0;
  }
L40:
  mp1 = m + 1;
  i__1 = *n;
  for (i = mp1; i <= i__1; i += 3)
  {
    dtemp = dx[i];
    dx[i] = dy[i];
    dy[i] = dtemp;
    dtemp = dx[i + 1];
    dx[i + 1] = dy[i + 1];
    dy[i + 1] = dtemp;
    dtemp = dx[i + 2];
    dx[i + 2] = dy[i + 2];
    dy[i + 2] = dtemp;
    /* L50: */
  }
  return 0;
} /* dswap_ */

/* Subroutine */ int dsymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy,
                            uplo_len)
char *uplo;
integer *n;
doublereal *alpha, *a;
integer *lda;
doublereal *x;
integer *incx;
doublereal *beta, *y;
integer *incy;
ftnlen uplo_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;

  /* Local variables */
  static integer info;
  static doublereal temp1, temp2;
  static integer i, j;
  extern logical lsame_();
  static integer ix, iy, jx, jy, kx, ky;
  extern /* Subroutine */ int xerbla_();

  /*     .. Scalar Arguments .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DSYMV  performs the matrix-vector  operation */

  /*     y := alpha*A*x + beta*y, */

  /*  where alpha and beta are scalars, x and y are n element vectors and */
  /*  A is an n by n symmetric matrix. */

  /*  Parameters */
  /*  ========== */

  /*  UPLO   - CHARACTER*1. */
  /*           On entry, UPLO specifies whether the upper or lower */
  /*           triangular part of the array A is to be referenced as */
  /*           follows: */

  /*              UPLO = 'U' or 'u'   Only the upper triangular part of A */
  /*                                  is to be referenced. */

  /*              UPLO = 'L' or 'l'   Only the lower triangular part of A */
  /*                                  is to be referenced. */

  /*           Unchanged on exit. */

  /*  N      - INTEGER. */
  /*           On entry, N specifies the order of the matrix A. */
  /*           N must be at least zero. */
  /*           Unchanged on exit. */

  /*  ALPHA  - DOUBLE PRECISION. */
  /*           On entry, ALPHA specifies the scalar alpha. */
  /*           Unchanged on exit. */

  /*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
  /*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
  /*           upper triangular part of the array A must contain the upper
  */
  /*           triangular part of the symmetric matrix and the strictly */
  /*           lower triangular part of A is not referenced. */
  /*           Before entry with UPLO = 'L' or 'l', the leading n by n */
  /*           lower triangular part of the array A must contain the lower
  */
  /*           triangular part of the symmetric matrix and the strictly */
  /*           upper triangular part of A is not referenced. */
  /*           Unchanged on exit. */

  /*  LDA    - INTEGER. */
  /*           On entry, LDA specifies the first dimension of A as declared
  */
  /*           in the calling (sub) program. LDA must be at least */
  /*           max( 1, n ). */
  /*           Unchanged on exit. */

  /*  X      - DOUBLE PRECISION array of dimension at least */
  /*           ( 1 + ( n - 1 )*abs( INCX ) ). */
  /*           Before entry, the incremented array X must contain the n */
  /*           element vector x. */
  /*           Unchanged on exit. */

  /*  INCX   - INTEGER. */
  /*           On entry, INCX specifies the increment for the elements of */
  /*           X. INCX must not be zero. */
  /*           Unchanged on exit. */

  /*  BETA   - DOUBLE PRECISION. */
  /*           On entry, BETA specifies the scalar beta. When BETA is */
  /*           supplied as zero then Y need not be set on input. */
  /*           Unchanged on exit. */

  /*  Y      - DOUBLE PRECISION array of dimension at least */
  /*           ( 1 + ( n - 1 )*abs( INCY ) ). */
  /*           Before entry, the incremented array Y must contain the n */
  /*           element vector y. On exit, Y is overwritten by the updated */
  /*           vector y. */

  /*  INCY   - INTEGER. */
  /*           On entry, INCY specifies the increment for the elements of */
  /*           Y. INCY must not be zero. */
  /*           Unchanged on exit. */


  /*  Level 2 Blas routine. */

  /*  -- Written on 22-October-1986. */
  /*     Jack Dongarra, Argonne National Lab. */
  /*     Jeremy Du Croz, Nag Central Office. */
  /*     Sven Hammarling, Nag Central Office. */
  /*     Richard Hanson, Sandia National Labs. */


  /*     .. Parameters .. */
  /*     .. Local Scalars .. */
  /*     .. External Functions .. */
  /*     .. External Subroutines .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  --y;
  --x;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  info = 0;
  if (! lsame_(uplo, "U", 1L, 1L) && ! lsame_(uplo, "L", 1L, 1L))
  {
    info = 1;
  }
  else if (*n < 0)
  {
    info = 2;
  }
  else if (*lda < max(1, *n))
  {
    info = 5;
  }
  else if (*incx == 0)
  {
    info = 7;
  }
  else if (*incy == 0)
  {
    info = 10;
  }
  if (info != 0)
  {
    xerbla_("DSYMV ", &info, 6L);
    return 0;
  }

  /*     Quick return if possible. */

  if (*n == 0 || *alpha == 0. && *beta == 1.)
  {
    return 0;
  }

  /*     Set up the start points in  X  and  Y. */

  if (*incx > 0)
  {
    kx = 1;
  }
  else
  {
    kx = 1 - (*n - 1) * *incx;
  }
  if (*incy > 0)
  {
    ky = 1;
  }
  else
  {
    ky = 1 - (*n - 1) * *incy;
  }

  /*     Start the operations. In this version the elements of A are */
  /*     accessed sequentially with one pass through the triangular part */
  /*     of A. */

  /*     First form  y := beta*y. */

  if (*beta != 1.)
  {
    if (*incy == 1)
    {
      if (*beta == 0.)
      {
        i__1 = *n;
        for (i = 1; i <= i__1; ++i)
        {
          y[i] = 0.;
          /* L10: */
        }
      }
      else
      {
        i__1 = *n;
        for (i = 1; i <= i__1; ++i)
        {
          y[i] = *beta * y[i];
          /* L20: */
        }
      }
    }
    else
    {
      iy = ky;
      if (*beta == 0.)
      {
        i__1 = *n;
        for (i = 1; i <= i__1; ++i)
        {
          y[iy] = 0.;
          iy += *incy;
          /* L30: */
        }
      }
      else
      {
        i__1 = *n;
        for (i = 1; i <= i__1; ++i)
        {
          y[iy] = *beta * y[iy];
          iy += *incy;
          /* L40: */
        }
      }
    }
  }
  if (*alpha == 0.)
  {
    return 0;
  }
  if (lsame_(uplo, "U", 1L, 1L))
  {

    /*        Form  y  when A is stored in upper triangle. */

    if (*incx == 1 && *incy == 1)
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        temp1 = *alpha * x[j];
        temp2 = 0.;
        i__2 = j - 1;
        for (i = 1; i <= i__2; ++i)
        {
          y[i] += temp1 * a[i + j * a_dim1];
          temp2 += a[i + j * a_dim1] * x[i];
          /* L50: */
        }
        y[j] = y[j] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
        /* L60: */
      }
    }
    else
    {
      jx = kx;
      jy = ky;
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        temp1 = *alpha * x[jx];
        temp2 = 0.;
        ix = kx;
        iy = ky;
        i__2 = j - 1;
        for (i = 1; i <= i__2; ++i)
        {
          y[iy] += temp1 * a[i + j * a_dim1];
          temp2 += a[i + j * a_dim1] * x[ix];
          ix += *incx;
          iy += *incy;
          /* L70: */
        }
        y[jy] = y[jy] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
        jx += *incx;
        jy += *incy;
        /* L80: */
      }
    }
  }
  else
  {

    /*        Form  y  when A is stored in lower triangle. */

    if (*incx == 1 && *incy == 1)
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        temp1 = *alpha * x[j];
        temp2 = 0.;
        y[j] += temp1 * a[j + j * a_dim1];
        i__2 = *n;
        for (i = j + 1; i <= i__2; ++i)
        {
          y[i] += temp1 * a[i + j * a_dim1];
          temp2 += a[i + j * a_dim1] * x[i];
          /* L90: */
        }
        y[j] += *alpha * temp2;
        /* L100: */
      }
    }
    else
    {
      jx = kx;
      jy = ky;
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        temp1 = *alpha * x[jx];
        temp2 = 0.;
        y[jy] += temp1 * a[j + j * a_dim1];
        ix = jx;
        iy = jy;
        i__2 = *n;
        for (i = j + 1; i <= i__2; ++i)
        {
          ix += *incx;
          iy += *incy;
          y[iy] += temp1 * a[i + j * a_dim1];
          temp2 += a[i + j * a_dim1] * x[ix];
          /* L110: */
        }
        y[jy] += *alpha * temp2;
        jx += *incx;
        jy += *incy;
        /* L120: */
      }
    }
  }

  return 0;

  /*     End of DSYMV . */

} /* dsymv_ */

/* Subroutine */ int dsyr2_(uplo, n, alpha, x, incx, y, incy, a, lda,
                            uplo_len)
char *uplo;
integer *n;
doublereal *alpha, *x;
integer *incx;
doublereal *y;
integer *incy;
doublereal *a;
integer *lda;
ftnlen uplo_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;

  /* Local variables */
  static integer info;
  static doublereal temp1, temp2;
  static integer i, j;
  extern logical lsame_();
  static integer ix, iy, jx, jy, kx, ky;
  extern /* Subroutine */ int xerbla_();

  /*     .. Scalar Arguments .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DSYR2  performs the symmetric rank 2 operation */

  /*     A := alpha*x*y' + alpha*y*x' + A, */

  /*  where alpha is a scalar, x and y are n element vectors and A is an n
  */
  /*  by n symmetric matrix. */

  /*  Parameters */
  /*  ========== */

  /*  UPLO   - CHARACTER*1. */
  /*           On entry, UPLO specifies whether the upper or lower */
  /*           triangular part of the array A is to be referenced as */
  /*           follows: */

  /*              UPLO = 'U' or 'u'   Only the upper triangular part of A */
  /*                                  is to be referenced. */

  /*              UPLO = 'L' or 'l'   Only the lower triangular part of A */
  /*                                  is to be referenced. */

  /*           Unchanged on exit. */

  /*  N      - INTEGER. */
  /*           On entry, N specifies the order of the matrix A. */
  /*           N must be at least zero. */
  /*           Unchanged on exit. */

  /*  ALPHA  - DOUBLE PRECISION. */
  /*           On entry, ALPHA specifies the scalar alpha. */
  /*           Unchanged on exit. */

  /*  X      - DOUBLE PRECISION array of dimension at least */
  /*           ( 1 + ( n - 1 )*abs( INCX ) ). */
  /*           Before entry, the incremented array X must contain the n */
  /*           element vector x. */
  /*           Unchanged on exit. */

  /*  INCX   - INTEGER. */
  /*           On entry, INCX specifies the increment for the elements of */
  /*           X. INCX must not be zero. */
  /*           Unchanged on exit. */

  /*  Y      - DOUBLE PRECISION array of dimension at least */
  /*           ( 1 + ( n - 1 )*abs( INCY ) ). */
  /*           Before entry, the incremented array Y must contain the n */
  /*           element vector y. */
  /*           Unchanged on exit. */

  /*  INCY   - INTEGER. */
  /*           On entry, INCY specifies the increment for the elements of */
  /*           Y. INCY must not be zero. */
  /*           Unchanged on exit. */

  /*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
  /*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
  /*           upper triangular part of the array A must contain the upper
  */
  /*           triangular part of the symmetric matrix and the strictly */
  /*           lower triangular part of A is not referenced. On exit, the */
  /*           upper triangular part of the array A is overwritten by the */
  /*           upper triangular part of the updated matrix. */
  /*           Before entry with UPLO = 'L' or 'l', the leading n by n */
  /*           lower triangular part of the array A must contain the lower
  */
  /*           triangular part of the symmetric matrix and the strictly */
  /*           upper triangular part of A is not referenced. On exit, the */
  /*           lower triangular part of the array A is overwritten by the */
  /*           lower triangular part of the updated matrix. */

  /*  LDA    - INTEGER. */
  /*           On entry, LDA specifies the first dimension of A as declared
  */
  /*           in the calling (sub) program. LDA must be at least */
  /*           max( 1, n ). */
  /*           Unchanged on exit. */


  /*  Level 2 Blas routine. */

  /*  -- Written on 22-October-1986. */
  /*     Jack Dongarra, Argonne National Lab. */
  /*     Jeremy Du Croz, Nag Central Office. */
  /*     Sven Hammarling, Nag Central Office. */
  /*     Richard Hanson, Sandia National Labs. */


  /*     .. Parameters .. */
  /*     .. Local Scalars .. */
  /*     .. External Functions .. */
  /*     .. External Subroutines .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;
  --y;
  --x;

  /* Function Body */
  info = 0;
  if (! lsame_(uplo, "U", 1L, 1L) && ! lsame_(uplo, "L", 1L, 1L))
  {
    info = 1;
  }
  else if (*n < 0)
  {
    info = 2;
  }
  else if (*incx == 0)
  {
    info = 5;
  }
  else if (*incy == 0)
  {
    info = 7;
  }
  else if (*lda < max(1, *n))
  {
    info = 9;
  }
  if (info != 0)
  {
    xerbla_("DSYR2 ", &info, 6L);
    return 0;
  }

  /*     Quick return if possible. */

  if (*n == 0 || *alpha == 0.)
  {
    return 0;
  }

  /*     Set up the start points in X and Y if the increments are not both
  */
  /*     unity. */

  if (*incx != 1 || *incy != 1)
  {
    if (*incx > 0)
    {
      kx = 1;
    }
    else
    {
      kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0)
    {
      ky = 1;
    }
    else
    {
      ky = 1 - (*n - 1) * *incy;
    }
    jx = kx;
    jy = ky;
  }

  /*     Start the operations. In this version the elements of A are */
  /*     accessed sequentially with one pass through the triangular part */
  /*     of A. */

  if (lsame_(uplo, "U", 1L, 1L))
  {

    /*        Form  A  when A is stored in the upper triangle. */

    if (*incx == 1 && *incy == 1)
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        if (x[j] != 0. || y[j] != 0.)
        {
          temp1 = *alpha * y[j];
          temp2 = *alpha * x[j];
          i__2 = j;
          for (i = 1; i <= i__2; ++i)
          {
            a[i + j * a_dim1] = a[i + j * a_dim1] + x[i] * temp1
                                + y[i] * temp2;
            /* L10: */
          }
        }
        /* L20: */
      }
    }
    else
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        if (x[jx] != 0. || y[jy] != 0.)
        {
          temp1 = *alpha * y[jy];
          temp2 = *alpha * x[jx];
          ix = kx;
          iy = ky;
          i__2 = j;
          for (i = 1; i <= i__2; ++i)
          {
            a[i + j * a_dim1] = a[i + j * a_dim1] + x[ix] * temp1
                                + y[iy] * temp2;
            ix += *incx;
            iy += *incy;
            /* L30: */
          }
        }
        jx += *incx;
        jy += *incy;
        /* L40: */
      }
    }
  }
  else
  {

    /*        Form  A  when A is stored in the lower triangle. */

    if (*incx == 1 && *incy == 1)
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        if (x[j] != 0. || y[j] != 0.)
        {
          temp1 = *alpha * y[j];
          temp2 = *alpha * x[j];
          i__2 = *n;
          for (i = j; i <= i__2; ++i)
          {
            a[i + j * a_dim1] = a[i + j * a_dim1] + x[i] * temp1
                                + y[i] * temp2;
            /* L50: */
          }
        }
        /* L60: */
      }
    }
    else
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        if (x[jx] != 0. || y[jy] != 0.)
        {
          temp1 = *alpha * y[jy];
          temp2 = *alpha * x[jx];
          ix = jx;
          iy = jy;
          i__2 = *n;
          for (i = j; i <= i__2; ++i)
          {
            a[i + j * a_dim1] = a[i + j * a_dim1] + x[ix] * temp1
                                + y[iy] * temp2;
            ix += *incx;
            iy += *incy;
            /* L70: */
          }
        }
        jx += *incx;
        jy += *incy;
        /* L80: */
      }
    }
  }

  return 0;

  /*     End of DSYR2 . */

} /* dsyr2_ */

/* Subroutine */ int dsyr2k_(uplo, trans, n, k, alpha, a, lda, b, ldb, beta,
                             c, ldc, uplo_len, trans_len)
char *uplo, *trans;
integer *n, *k;
doublereal *alpha, *a;
integer *lda;
doublereal *b;
integer *ldb;
doublereal *beta, *c;
integer *ldc;
ftnlen uplo_len;
ftnlen trans_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
          i__3;

  /* Local variables */
  static integer info;
  static doublereal temp1, temp2;
  static integer i, j, l;
  extern logical lsame_();
  static integer nrowa;
  static logical upper;
  extern /* Subroutine */ int xerbla_();

  /*     .. Scalar Arguments .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */
  /*           lower triangular part of the array C must contain the lower
  */
  /*           triangular part  of the  symmetric matrix  and the strictly
  */
  /*           upper triangular part of C is not referenced.  On exit, the
  */
  /*           lower triangular part of the array  C is overwritten by the
  */
  /*           lower triangular part of the updated matrix. */

  /*  LDC    - INTEGER. */
  /*           On entry, LDC specifies the first dimension of C as declared
  */
  /*           in  the  calling  (sub)  program.   LDC  must  be  at  least
  */
  /*           max( 1, n ). */
  /*           Unchanged on exit. */


  /*  Level 3 Blas routine. */


  /*  -- Written on 8-February-1989. */
  /*     Jack Dongarra, Argonne National Laboratory. */
  /*     Iain Duff, AERE Harwell. */
  /*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
  /*     Sven Hammarling, Numerical Algorithms Group Ltd. */


  /*     .. External Functions .. */
  /*     .. External Subroutines .. */
  /*     .. Intrinsic Functions .. */
  /*     .. Local Scalars .. */
  /*     .. Parameters .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  c_dim1 = *ldc;
  c_offset = c_dim1 + 1;
  c -= c_offset;
  b_dim1 = *ldb;
  b_offset = b_dim1 + 1;
  b -= b_offset;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  if (lsame_(trans, "N", 1L, 1L))
  {
    nrowa = *n;
  }
  else
  {
    nrowa = *k;
  }
  upper = lsame_(uplo, "U", 1L, 1L);

  info = 0;
  if (! upper && ! lsame_(uplo, "L", 1L, 1L))
  {
    info = 1;
  }
  else if (! lsame_(trans, "N", 1L, 1L) && ! lsame_(trans, "T", 1L, 1L) &&
           ! lsame_(trans, "C", 1L, 1L))
  {
    info = 2;
  }
  else if (*n < 0)
  {
    info = 3;
  }
  else if (*k < 0)
  {
    info = 4;
  }
  else if (*lda < max(1, nrowa))
  {
    info = 7;
  }
  else if (*ldb < max(1, nrowa))
  {
    info = 9;
  }
  else if (*ldc < max(1, *n))
  {
    info = 12;
  }
  if (info != 0)
  {
    xerbla_("DSYR2K", &info, 6L);
    return 0;
  }

  /*     Quick return if possible. */

  if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.)
  {
    return 0;
  }

  /*     And when  alpha.eq.zero. */

  if (*alpha == 0.)
  {
    if (upper)
    {
      if (*beta == 0.)
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = j;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = 0.;
            /* L10: */
          }
          /* L20: */
        }
      }
      else
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = j;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = *beta * c[i + j * c_dim1];
            /* L30: */
          }
          /* L40: */
        }
      }
    }
    else
    {
      if (*beta == 0.)
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *n;
          for (i = j; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = 0.;
            /* L50: */
          }
          /* L60: */
        }
      }
      else
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *n;
          for (i = j; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = *beta * c[i + j * c_dim1];
            /* L70: */
          }
          /* L80: */
        }
      }
    }
    return 0;
  }

  /*     Start the operations. */

  if (lsame_(trans, "N", 1L, 1L))
  {

    /*        Form  C := alpha*A*B' + alpha*B*A' + C. */

    if (upper)
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        if (*beta == 0.)
        {
          i__2 = j;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = 0.;
            /* L90: */
          }
        }
        else if (*beta != 1.)
        {
          i__2 = j;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = *beta * c[i + j * c_dim1];
            /* L100: */
          }
        }
        i__2 = *k;
        for (l = 1; l <= i__2; ++l)
        {
          if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.)
          {
            temp1 = *alpha * b[j + l * b_dim1];
            temp2 = *alpha * a[j + l * a_dim1];
            i__3 = j;
            for (i = 1; i <= i__3; ++i)
            {
              c[i + j * c_dim1] = c[i + j * c_dim1] + a[i + l *
                                  a_dim1] * temp1 + b[i + l * b_dim1] *
                                  temp2;
              /* L110: */
            }
          }
          /* L120: */
        }
        /* L130: */
      }
    }
    else
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        if (*beta == 0.)
        {
          i__2 = *n;
          for (i = j; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = 0.;
            /* L140: */
          }
        }
        else if (*beta != 1.)
        {
          i__2 = *n;
          for (i = j; i <= i__2; ++i)
          {
            c[i + j * c_dim1] = *beta * c[i + j * c_dim1];
            /* L150: */
          }
        }
        i__2 = *k;
        for (l = 1; l <= i__2; ++l)
        {
          if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.)
          {
            temp1 = *alpha * b[j + l * b_dim1];
            temp2 = *alpha * a[j + l * a_dim1];
            i__3 = *n;
            for (i = j; i <= i__3; ++i)
            {
              c[i + j * c_dim1] = c[i + j * c_dim1] + a[i + l *
                                  a_dim1] * temp1 + b[i + l * b_dim1] *
                                  temp2;
              /* L160: */
            }
          }
          /* L170: */
        }
        /* L180: */
      }
    }
  }
  else
  {

    /*        Form  C := alpha*A'*B + alpha*B'*A + C. */

    if (upper)
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        i__2 = j;
        for (i = 1; i <= i__2; ++i)
        {
          temp1 = 0.;
          temp2 = 0.;
          i__3 = *k;
          for (l = 1; l <= i__3; ++l)
          {
            temp1 += a[l + i * a_dim1] * b[l + j * b_dim1];
            temp2 += b[l + i * b_dim1] * a[l + j * a_dim1];
            /* L190: */
          }
          if (*beta == 0.)
          {
            c[i + j * c_dim1] = *alpha * temp1 + *alpha * temp2;
          }
          else
          {
            c[i + j * c_dim1] = *beta * c[i + j * c_dim1] + *
                                alpha * temp1 + *alpha * temp2;
          }
          /* L200: */
        }
        /* L210: */
      }
    }
    else
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        i__2 = *n;
        for (i = j; i <= i__2; ++i)
        {
          temp1 = 0.;
          temp2 = 0.;
          i__3 = *k;
          for (l = 1; l <= i__3; ++l)
          {
            temp1 += a[l + i * a_dim1] * b[l + j * b_dim1];
            temp2 += b[l + i * b_dim1] * a[l + j * a_dim1];
            /* L220: */
          }
          if (*beta == 0.)
          {
            c[i + j * c_dim1] = *alpha * temp1 + *alpha * temp2;
          }
          else
          {
            c[i + j * c_dim1] = *beta * c[i + j * c_dim1] + *
                                alpha * temp1 + *alpha * temp2;
          }
          /* L230: */
        }
        /* L240: */
      }
    }
  }

  return 0;

  /*     End of DSYR2K. */

} /* dsyr2k_ */

/* Subroutine */ int dtrmm_(side, uplo, transa, diag, m, n, alpha, a, lda, b,
                            ldb, side_len, uplo_len, transa_len, diag_len)
char *side, *uplo, *transa, *diag;
integer *m, *n;
doublereal *alpha, *a;
integer *lda;
doublereal *b;
integer *ldb;
ftnlen side_len;
ftnlen uplo_len;
ftnlen transa_len;
ftnlen diag_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

  /* Local variables */
  static integer info;
  static doublereal temp;
  static integer i, j, k;
  static logical lside;
  extern logical lsame_();
  static integer nrowa;
  static logical upper;
  extern /* Subroutine */ int xerbla_();
  static logical nounit;

  /*     .. Scalar Arguments .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DTRMM  performs one of the matrix-matrix operations */

  /*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ), */

  /*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
  */
  /*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
  */

  /*     op( A ) = A   or   op( A ) = A'. */

  /*  Parameters */
  /*  ========== */

  /*  SIDE   - CHARACTER*1. */
  /*           On entry,  SIDE specifies whether  op( A ) multiplies B from
  */
  /*           the left or right as follows: */

  /*              SIDE = 'L' or 'l'   B := alpha*op( A )*B. */

  /*              SIDE = 'R' or 'r'   B := alpha*B*op( A ). */

  /*           Unchanged on exit. */

  /*  UPLO   - CHARACTER*1. */
  /*           On entry, UPLO specifies whether the matrix A is an upper or
  */
  /*           lower triangular matrix as follows: */

  /*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

  /*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

  /*           Unchanged on exit. */

  /*  TRANSA - CHARACTER*1. */
  /*           On entry, TRANSA specifies the form of op( A ) to be used in
  */
  /*           the matrix multiplication as follows: */

  /*              TRANSA = 'N' or 'n'   op( A ) = A. */

  /*              TRANSA = 'T' or 't'   op( A ) = A'. */

  /*              TRANSA = 'C' or 'c'   op( A ) = A'. */

  /*           Unchanged on exit. */

  /*  DIAG   - CHARACTER*1. */
  /*           On entry, DIAG specifies whether or not A is unit triangular
  */
  /*           as follows: */

  /*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

  /*              DIAG = 'N' or 'n'   A is not assumed to be unit */
  /*                                  triangular. */

  /*           Unchanged on exit. */

  /*  M      - INTEGER. */
  /*           On entry, M specifies the number of rows of B. M must be at
  */
  /*           least zero. */
  /*           Unchanged on exit. */

  /*  N      - INTEGER. */
  /*           On entry, N specifies the number of columns of B.  N must be
  */
  /*           at least zero. */
  /*           Unchanged on exit. */

  /*  ALPHA  - DOUBLE PRECISION. */
  /*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
  */
  /*           zero then  A is not referenced and  B need not be set before
  */
  /*           entry. */
  /*           Unchanged on exit. */

  /*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
  */
  /*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
  */
  /*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
  */
  /*           upper triangular part of the array  A must contain the upper
  */
  /*           triangular matrix  and the strictly lower triangular part of
  */
  /*           A is not referenced. */
  /*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
  */
  /*           lower triangular part of the array  A must contain the lower
  */
  /*           triangular matrix  and the strictly upper triangular part of
  */
  /*           A is not referenced. */
  /*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
  */
  /*           A  are not referenced either,  but are assumed to be  unity.
  */
  /*           Unchanged on exit. */

  /*  LDA    - INTEGER. */
  /*           On entry, LDA specifies the first dimension of A as declared
  */
  /*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
  */
  /*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
  */
  /*           then LDA must be at least max( 1, n ). */
  /*           Unchanged on exit. */

  /*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ). */
  /*           Before entry,  the leading  m by n part of the array  B must
  */
  /*           contain the matrix  B,  and  on exit  is overwritten  by the
  */
  /*           transformed matrix. */

  /*  LDB    - INTEGER. */
  /*           On entry, LDB specifies the first dimension of B as declared
  */
  /*           in  the  calling  (sub)  program.   LDB  must  be  at  least
  */
  /*           max( 1, m ). */
  /*           Unchanged on exit. */


  /*  Level 3 Blas routine. */

  /*  -- Written on 8-February-1989. */
  /*     Jack Dongarra, Argonne National Laboratory. */
  /*     Iain Duff, AERE Harwell. */
  /*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
  /*     Sven Hammarling, Numerical Algorithms Group Ltd. */


  /*     .. External Functions .. */
  /*     .. External Subroutines .. */
  /*     .. Intrinsic Functions .. */
  /*     .. Local Scalars .. */
  /*     .. Parameters .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  b_dim1 = *ldb;
  b_offset = b_dim1 + 1;
  b -= b_offset;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  lside = lsame_(side, "L", 1L, 1L);
  if (lside)
  {
    nrowa = *m;
  }
  else
  {
    nrowa = *n;
  }
  nounit = lsame_(diag, "N", 1L, 1L);
  upper = lsame_(uplo, "U", 1L, 1L);

  info = 0;
  if (! lside && ! lsame_(side, "R", 1L, 1L))
  {
    info = 1;
  }
  else if (! upper && ! lsame_(uplo, "L", 1L, 1L))
  {
    info = 2;
  }
  else if (! lsame_(transa, "N", 1L, 1L) && ! lsame_(transa, "T", 1L, 1L)
           && ! lsame_(transa, "C", 1L, 1L))
  {
    info = 3;
  }
  else if (! lsame_(diag, "U", 1L, 1L) && ! lsame_(diag, "N", 1L, 1L))
  {
    info = 4;
  }
  else if (*m < 0)
  {
    info = 5;
  }
  else if (*n < 0)
  {
    info = 6;
  }
  else if (*lda < max(1, nrowa))
  {
    info = 9;
  }
  else if (*ldb < max(1, *m))
  {
    info = 11;
  }
  if (info != 0)
  {
    xerbla_("DTRMM ", &info, 6L);
    return 0;
  }

  /*     Quick return if possible. */

  if (*n == 0)
  {
    return 0;
  }

  /*     And when  alpha.eq.zero. */

  if (*alpha == 0.)
  {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j)
    {
      i__2 = *m;
      for (i = 1; i <= i__2; ++i)
      {
        b[i + j * b_dim1] = 0.;
        /* L10: */
      }
      /* L20: */
    }
    return 0;
  }

  /*     Start the operations. */

  if (lside)
  {
    if (lsame_(transa, "N", 1L, 1L))
    {

      /*           Form  B := alpha*A*B. */

      if (upper)
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *m;
          for (k = 1; k <= i__2; ++k)
          {
            if (b[k + j * b_dim1] != 0.)
            {
              temp = *alpha * b[k + j * b_dim1];
              i__3 = k - 1;
              for (i = 1; i <= i__3; ++i)
              {
                b[i + j * b_dim1] += temp * a[i + k * a_dim1];
                /* L30: */
              }
              if (nounit)
              {
                temp *= a[k + k * a_dim1];
              }
              b[k + j * b_dim1] = temp;
            }
            /* L40: */
          }
          /* L50: */
        }
      }
      else
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          for (k = *m; k >= 1; --k)
          {
            if (b[k + j * b_dim1] != 0.)
            {
              temp = *alpha * b[k + j * b_dim1];
              b[k + j * b_dim1] = temp;
              if (nounit)
              {
                b[k + j * b_dim1] *= a[k + k * a_dim1];
              }
              i__2 = *m;
              for (i = k + 1; i <= i__2; ++i)
              {
                b[i + j * b_dim1] += temp * a[i + k * a_dim1];
                /* L60: */
              }
            }
            /* L70: */
          }
          /* L80: */
        }
      }
    }
    else
    {

      /*           Form  B := alpha*B*A'. */

      if (upper)
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          for (i = *m; i >= 1; --i)
          {
            temp = b[i + j * b_dim1];
            if (nounit)
            {
              temp *= a[i + i * a_dim1];
            }
            i__2 = i - 1;
            for (k = 1; k <= i__2; ++k)
            {
              temp += a[k + i * a_dim1] * b[k + j * b_dim1];
              /* L90: */
            }
            b[i + j * b_dim1] = *alpha * temp;
            /* L100: */
          }
          /* L110: */
        }
      }
      else
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            temp = b[i + j * b_dim1];
            if (nounit)
            {
              temp *= a[i + i * a_dim1];
            }
            i__3 = *m;
            for (k = i + 1; k <= i__3; ++k)
            {
              temp += a[k + i * a_dim1] * b[k + j * b_dim1];
              /* L120: */
            }
            b[i + j * b_dim1] = *alpha * temp;
            /* L130: */
          }
          /* L140: */
        }
      }
    }
  }
  else
  {
    if (lsame_(transa, "N", 1L, 1L))
    {

      /*           Form  B := alpha*B*A. */

      if (upper)
      {
        for (j = *n; j >= 1; --j)
        {
          temp = *alpha;
          if (nounit)
          {
            temp *= a[j + j * a_dim1];
          }
          i__1 = *m;
          for (i = 1; i <= i__1; ++i)
          {
            b[i + j * b_dim1] = temp * b[i + j * b_dim1];
            /* L150: */
          }
          i__1 = j - 1;
          for (k = 1; k <= i__1; ++k)
          {
            if (a[k + j * a_dim1] != 0.)
            {
              temp = *alpha * a[k + j * a_dim1];
              i__2 = *m;
              for (i = 1; i <= i__2; ++i)
              {
                b[i + j * b_dim1] += temp * b[i + k * b_dim1];
                /* L160: */
              }
            }
            /* L170: */
          }
          /* L180: */
        }
      }
      else
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          temp = *alpha;
          if (nounit)
          {
            temp *= a[j + j * a_dim1];
          }
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            b[i + j * b_dim1] = temp * b[i + j * b_dim1];
            /* L190: */
          }
          i__2 = *n;
          for (k = j + 1; k <= i__2; ++k)
          {
            if (a[k + j * a_dim1] != 0.)
            {
              temp = *alpha * a[k + j * a_dim1];
              i__3 = *m;
              for (i = 1; i <= i__3; ++i)
              {
                b[i + j * b_dim1] += temp * b[i + k * b_dim1];
                /* L200: */
              }
            }
            /* L210: */
          }
          /* L220: */
        }
      }
    }
    else
    {

      /*           Form  B := alpha*B*A'. */

      if (upper)
      {
        i__1 = *n;
        for (k = 1; k <= i__1; ++k)
        {
          i__2 = k - 1;
          for (j = 1; j <= i__2; ++j)
          {
            if (a[j + k * a_dim1] != 0.)
            {
              temp = *alpha * a[j + k * a_dim1];
              i__3 = *m;
              for (i = 1; i <= i__3; ++i)
              {
                b[i + j * b_dim1] += temp * b[i + k * b_dim1];
                /* L230: */
              }
            }
            /* L240: */
          }
          temp = *alpha;
          if (nounit)
          {
            temp *= a[k + k * a_dim1];
          }
          if (temp != 1.)
          {
            i__2 = *m;
            for (i = 1; i <= i__2; ++i)
            {
              b[i + k * b_dim1] = temp * b[i + k * b_dim1];
              /* L250: */
            }
          }
          /* L260: */
        }
      }
      else
      {
        for (k = *n; k >= 1; --k)
        {
          i__1 = *n;
          for (j = k + 1; j <= i__1; ++j)
          {
            if (a[j + k * a_dim1] != 0.)
            {
              temp = *alpha * a[j + k * a_dim1];
              i__2 = *m;
              for (i = 1; i <= i__2; ++i)
              {
                b[i + j * b_dim1] += temp * b[i + k * b_dim1];
                /* L270: */
              }
            }
            /* L280: */
          }
          temp = *alpha;
          if (nounit)
          {
            temp *= a[k + k * a_dim1];
          }
          if (temp != 1.)
          {
            i__1 = *m;
            for (i = 1; i <= i__1; ++i)
            {
              b[i + k * b_dim1] = temp * b[i + k * b_dim1];
              /* L290: */
            }
          }
          /* L300: */
        }
      }
    }
  }

  return 0;

  /*     End of DTRMM . */

} /* dtrmm_ */

/* Subroutine */ int dtrmv_(uplo, trans, diag, n, a, lda, x, incx, uplo_len,
                            trans_len, diag_len)
char *uplo, *trans, *diag;
integer *n;
doublereal *a;
integer *lda;
doublereal *x;
integer *incx;
ftnlen uplo_len;
ftnlen trans_len;
ftnlen diag_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;

  /* Local variables */
  static integer info;
  static doublereal temp;
  static integer i, j;
  extern logical lsame_();
  static integer ix, jx, kx;
  extern /* Subroutine */ int xerbla_();
  static logical nounit;

  /*     .. Scalar Arguments .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DTRMV  performs one of the matrix-vector operations */

  /*     x := A*x,   or   x := A'*x, */

  /*  where x is an n element vector and  A is an n by n unit, or non-unit,
  */
  /*  upper or lower triangular matrix. */

  /*  Parameters */
  /*  ========== */

  /*  UPLO   - CHARACTER*1. */
  /*           On entry, UPLO specifies whether the matrix is an upper or */
  /*           lower triangular matrix as follows: */

  /*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

  /*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

  /*           Unchanged on exit. */

  /*  TRANS  - CHARACTER*1. */
  /*           On entry, TRANS specifies the operation to be performed as */
  /*           follows: */

  /*              TRANS = 'N' or 'n'   x := A*x. */

  /*              TRANS = 'T' or 't'   x := A'*x. */

  /*              TRANS = 'C' or 'c'   x := A'*x. */

  /*           Unchanged on exit. */

  /*  DIAG   - CHARACTER*1. */
  /*           On entry, DIAG specifies whether or not A is unit */
  /*           triangular as follows: */

  /*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

  /*              DIAG = 'N' or 'n'   A is not assumed to be unit */
  /*                                  triangular. */

  /*           Unchanged on exit. */

  /*  N      - INTEGER. */
  /*           On entry, N specifies the order of the matrix A. */
  /*           N must be at least zero. */
  /*           Unchanged on exit. */

  /*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
  /*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
  /*           upper triangular part of the array A must contain the upper
  */
  /*           triangular matrix and the strictly lower triangular part of
  */
  /*           A is not referenced. */
  /*           Before entry with UPLO = 'L' or 'l', the leading n by n */
  /*           lower triangular part of the array A must contain the lower
  */
  /*           triangular matrix and the strictly upper triangular part of
  */
  /*           A is not referenced. */
  /*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
  */
  /*           A are not referenced either, but are assumed to be unity. */
  /*           Unchanged on exit. */

  /*  LDA    - INTEGER. */
  /*           On entry, LDA specifies the first dimension of A as declared
  */
  /*           in the calling (sub) program. LDA must be at least */
  /*           max( 1, n ). */
  /*           Unchanged on exit. */

  /*  X      - DOUBLE PRECISION array of dimension at least */
  /*           ( 1 + ( n - 1 )*abs( INCX ) ). */
  /*           Before entry, the incremented array X must contain the n */
  /*           element vector x. On exit, X is overwritten with the */
  /*           tranformed vector x. */

  /*  INCX   - INTEGER. */
  /*           On entry, INCX specifies the increment for the elements of */
  /*           X. INCX must not be zero. */
  /*           Unchanged on exit. */


  /*  Level 2 Blas routine. */

  /*  -- Written on 22-October-1986. */
  /*     Jack Dongarra, Argonne National Lab. */
  /*     Jeremy Du Croz, Nag Central Office. */
  /*     Sven Hammarling, Nag Central Office. */
  /*     Richard Hanson, Sandia National Labs. */


  /*     .. Parameters .. */
  /*     .. Local Scalars .. */
  /*     .. External Functions .. */
  /*     .. External Subroutines .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  --x;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  info = 0;
  if (! lsame_(uplo, "U", 1L, 1L) && ! lsame_(uplo, "L", 1L, 1L))
  {
    info = 1;
  }
  else if (! lsame_(trans, "N", 1L, 1L) && ! lsame_(trans, "T", 1L, 1L) &&
           ! lsame_(trans, "C", 1L, 1L))
  {
    info = 2;
  }
  else if (! lsame_(diag, "U", 1L, 1L) && ! lsame_(diag, "N", 1L, 1L))
  {
    info = 3;
  }
  else if (*n < 0)
  {
    info = 4;
  }
  else if (*lda < max(1, *n))
  {
    info = 6;
  }
  else if (*incx == 0)
  {
    info = 8;
  }
  if (info != 0)
  {
    xerbla_("DTRMV ", &info, 6L);
    return 0;
  }

  /*     Quick return if possible. */

  if (*n == 0)
  {
    return 0;
  }

  nounit = lsame_(diag, "N", 1L, 1L);

  /*     Set up the start point in X if the increment is not unity. This */
  /*     will be  ( N - 1 )*INCX  too small for descending loops. */

  if (*incx <= 0)
  {
    kx = 1 - (*n - 1) * *incx;
  }
  else if (*incx != 1)
  {
    kx = 1;
  }

  /*     Start the operations. In this version the elements of A are */
  /*     accessed sequentially with one pass through A. */

  if (lsame_(trans, "N", 1L, 1L))
  {

    /*        Form  x := A*x. */

    if (lsame_(uplo, "U", 1L, 1L))
    {
      if (*incx == 1)
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          if (x[j] != 0.)
          {
            temp = x[j];
            i__2 = j - 1;
            for (i = 1; i <= i__2; ++i)
            {
              x[i] += temp * a[i + j * a_dim1];
              /* L10: */
            }
            if (nounit)
            {
              x[j] *= a[j + j * a_dim1];
            }
          }
          /* L20: */
        }
      }
      else
      {
        jx = kx;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          if (x[jx] != 0.)
          {
            temp = x[jx];
            ix = kx;
            i__2 = j - 1;
            for (i = 1; i <= i__2; ++i)
            {
              x[ix] += temp * a[i + j * a_dim1];
              ix += *incx;
              /* L30: */
            }
            if (nounit)
            {
              x[jx] *= a[j + j * a_dim1];
            }
          }
          jx += *incx;
          /* L40: */
        }
      }
    }
    else
    {
      if (*incx == 1)
      {
        for (j = *n; j >= 1; --j)
        {
          if (x[j] != 0.)
          {
            temp = x[j];
            i__1 = j + 1;
            for (i = *n; i >= i__1; --i)
            {
              x[i] += temp * a[i + j * a_dim1];
              /* L50: */
            }
            if (nounit)
            {
              x[j] *= a[j + j * a_dim1];
            }
          }
          /* L60: */
        }
      }
      else
      {
        kx += (*n - 1) * *incx;
        jx = kx;
        for (j = *n; j >= 1; --j)
        {
          if (x[jx] != 0.)
          {
            temp = x[jx];
            ix = kx;
            i__1 = j + 1;
            for (i = *n; i >= i__1; --i)
            {
              x[ix] += temp * a[i + j * a_dim1];
              ix -= *incx;
              /* L70: */
            }
            if (nounit)
            {
              x[jx] *= a[j + j * a_dim1];
            }
          }
          jx -= *incx;
          /* L80: */
        }
      }
    }
  }
  else
  {

    /*        Form  x := A'*x. */

    if (lsame_(uplo, "U", 1L, 1L))
    {
      if (*incx == 1)
      {
        for (j = *n; j >= 1; --j)
        {
          temp = x[j];
          if (nounit)
          {
            temp *= a[j + j * a_dim1];
          }
          for (i = j - 1; i >= 1; --i)
          {
            temp += a[i + j * a_dim1] * x[i];
            /* L90: */
          }
          x[j] = temp;
          /* L100: */
        }
      }
      else
      {
        jx = kx + (*n - 1) * *incx;
        for (j = *n; j >= 1; --j)
        {
          temp = x[jx];
          ix = jx;
          if (nounit)
          {
            temp *= a[j + j * a_dim1];
          }
          for (i = j - 1; i >= 1; --i)
          {
            ix -= *incx;
            temp += a[i + j * a_dim1] * x[ix];
            /* L110: */
          }
          x[jx] = temp;
          jx -= *incx;
          /* L120: */
        }
      }
    }
    else
    {
      if (*incx == 1)
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          temp = x[j];
          if (nounit)
          {
            temp *= a[j + j * a_dim1];
          }
          i__2 = *n;
          for (i = j + 1; i <= i__2; ++i)
          {
            temp += a[i + j * a_dim1] * x[i];
            /* L130: */
          }
          x[j] = temp;
          /* L140: */
        }
      }
      else
      {
        jx = kx;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          temp = x[jx];
          ix = jx;
          if (nounit)
          {
            temp *= a[j + j * a_dim1];
          }
          i__2 = *n;
          for (i = j + 1; i <= i__2; ++i)
          {
            ix += *incx;
            temp += a[i + j * a_dim1] * x[ix];
            /* L150: */
          }
          x[jx] = temp;
          jx += *incx;
          /* L160: */
        }
      }
    }
  }

  return 0;

  /*     End of DTRMV . */

} /* dtrmv_ */

/* Subroutine */ int dtrsv_(uplo, trans, diag, n, a, lda, x, incx, uplo_len,
                            trans_len, diag_len)
char *uplo, *trans, *diag;
integer *n;
doublereal *a;
integer *lda;
doublereal *x;
integer *incx;
ftnlen uplo_len;
ftnlen trans_len;
ftnlen diag_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;

  /* Local variables */
  static integer info;
  static doublereal temp;
  static integer i, j;
  extern logical lsame_();
  static integer ix, jx, kx;
  extern /* Subroutine */ int xerbla_();
  static logical nounit;

  /*     .. Scalar Arguments .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DTRSV  solves one of the systems of equations */

  /*     A*x = b,   or   A'*x = b, */

  /*  where b and x are n element vectors and A is an n by n unit, or */
  /*  non-unit, upper or lower triangular matrix. */

  /*  No test for singularity or near-singularity is included in this */
  /*  routine. Such tests must be performed before calling this routine. */

  /*  Parameters */
  /*  ========== */

  /*  UPLO   - CHARACTER*1. */
  /*           On entry, UPLO specifies whether the matrix is an upper or */
  /*           lower triangular matrix as follows: */

  /*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

  /*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

  /*           Unchanged on exit. */

  /*  TRANS  - CHARACTER*1. */
  /*           On entry, TRANS specifies the equations to be solved as */
  /*           follows: */

  /*              TRANS = 'N' or 'n'   A*x = b. */

  /*              TRANS = 'T' or 't'   A'*x = b. */

  /*              TRANS = 'C' or 'c'   A'*x = b. */

  /*           Unchanged on exit. */

  /*  DIAG   - CHARACTER*1. */
  /*           On entry, DIAG specifies whether or not A is unit */
  /*           triangular as follows: */

  /*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

  /*              DIAG = 'N' or 'n'   A is not assumed to be unit */
  /*                                  triangular. */

  /*           Unchanged on exit. */

  /*  N      - INTEGER. */
  /*           On entry, N specifies the order of the matrix A. */
  /*           N must be at least zero. */
  /*           Unchanged on exit. */

  /*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
  /*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
  /*           upper triangular part of the array A must contain the upper
  */
  /*           triangular matrix and the strictly lower triangular part of
  */
  /*           A is not referenced. */
  /*           Before entry with UPLO = 'L' or 'l', the leading n by n */
  /*           lower triangular part of the array A must contain the lower
  */
  /*           triangular matrix and the strictly upper triangular part of
  */
  /*           A is not referenced. */
  /*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
  */
  /*           A are not referenced either, but are assumed to be unity. */
  /*           Unchanged on exit. */

  /*  LDA    - INTEGER. */
  /*           On entry, LDA specifies the first dimension of A as declared
  */
  /*           in the calling (sub) program. LDA must be at least */
  /*           max( 1, n ). */
  /*           Unchanged on exit. */

  /*  X      - DOUBLE PRECISION array of dimension at least */
  /*           ( 1 + ( n - 1 )*abs( INCX ) ). */
  /*           Before entry, the incremented array X must contain the n */
  /*           element right-hand side vector b. On exit, X is overwritten
  */
  /*           with the solution vector x. */

  /*  INCX   - INTEGER. */
  /*           On entry, INCX specifies the increment for the elements of */
  /*           X. INCX must not be zero. */
  /*           Unchanged on exit. */


  /*  Level 2 Blas routine. */

  /*  -- Written on 22-October-1986. */
  /*     Jack Dongarra, Argonne National Lab. */
  /*     Jeremy Du Croz, Nag Central Office. */
  /*     Sven Hammarling, Nag Central Office. */
  /*     Richard Hanson, Sandia National Labs. */


  /*     .. Parameters .. */
  /*     .. Local Scalars .. */
  /*     .. External Functions .. */
  /*     .. External Subroutines .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  --x;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  info = 0;
  if (! lsame_(uplo, "U", 1L, 1L) && ! lsame_(uplo, "L", 1L, 1L))
  {
    info = 1;
  }
  else if (! lsame_(trans, "N", 1L, 1L) && ! lsame_(trans, "T", 1L, 1L) &&
           ! lsame_(trans, "C", 1L, 1L))
  {
    info = 2;
  }
  else if (! lsame_(diag, "U", 1L, 1L) && ! lsame_(diag, "N", 1L, 1L))
  {
    info = 3;
  }
  else if (*n < 0)
  {
    info = 4;
  }
  else if (*lda < max(1, *n))
  {
    info = 6;
  }
  else if (*incx == 0)
  {
    info = 8;
  }
  if (info != 0)
  {
    xerbla_("DTRSV ", &info, 6L);
    return 0;
  }

  /*     Quick return if possible. */

  if (*n == 0)
  {
    return 0;
  }

  nounit = lsame_(diag, "N", 1L, 1L);

  /*     Set up the start point in X if the increment is not unity. This */
  /*     will be  ( N - 1 )*INCX  too small for descending loops. */

  if (*incx <= 0)
  {
    kx = 1 - (*n - 1) * *incx;
  }
  else if (*incx != 1)
  {
    kx = 1;
  }

  /*     Start the operations. In this version the elements of A are */
  /*     accessed sequentially with one pass through A. */

  if (lsame_(trans, "N", 1L, 1L))
  {

    /*        Form  x := inv( A )*x. */

    if (lsame_(uplo, "U", 1L, 1L))
    {
      if (*incx == 1)
      {
        for (j = *n; j >= 1; --j)
        {
          if (x[j] != 0.)
          {
            if (nounit)
            {
              x[j] /= a[j + j * a_dim1];
            }
            temp = x[j];
            for (i = j - 1; i >= 1; --i)
            {
              x[i] -= temp * a[i + j * a_dim1];
              /* L10: */
            }
          }
          /* L20: */
        }
      }
      else
      {
        jx = kx + (*n - 1) * *incx;
        for (j = *n; j >= 1; --j)
        {
          if (x[jx] != 0.)
          {
            if (nounit)
            {
              x[jx] /= a[j + j * a_dim1];
            }
            temp = x[jx];
            ix = jx;
            for (i = j - 1; i >= 1; --i)
            {
              ix -= *incx;
              x[ix] -= temp * a[i + j * a_dim1];
              /* L30: */
            }
          }
          jx -= *incx;
          /* L40: */
        }
      }
    }
    else
    {
      if (*incx == 1)
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          if (x[j] != 0.)
          {
            if (nounit)
            {
              x[j] /= a[j + j * a_dim1];
            }
            temp = x[j];
            i__2 = *n;
            for (i = j + 1; i <= i__2; ++i)
            {
              x[i] -= temp * a[i + j * a_dim1];
              /* L50: */
            }
          }
          /* L60: */
        }
      }
      else
      {
        jx = kx;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          if (x[jx] != 0.)
          {
            if (nounit)
            {
              x[jx] /= a[j + j * a_dim1];
            }
            temp = x[jx];
            ix = jx;
            i__2 = *n;
            for (i = j + 1; i <= i__2; ++i)
            {
              ix += *incx;
              x[ix] -= temp * a[i + j * a_dim1];
              /* L70: */
            }
          }
          jx += *incx;
          /* L80: */
        }
      }
    }
  }
  else
  {

    /*        Form  x := inv( A' )*x. */

    if (lsame_(uplo, "U", 1L, 1L))
    {
      if (*incx == 1)
      {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          temp = x[j];
          i__2 = j - 1;
          for (i = 1; i <= i__2; ++i)
          {
            temp -= a[i + j * a_dim1] * x[i];
            /* L90: */
          }
          if (nounit)
          {
            temp /= a[j + j * a_dim1];
          }
          x[j] = temp;
          /* L100: */
        }
      }
      else
      {
        jx = kx;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
          temp = x[jx];
          ix = kx;
          i__2 = j - 1;
          for (i = 1; i <= i__2; ++i)
          {
            temp -= a[i + j * a_dim1] * x[ix];
            ix += *incx;
            /* L110: */
          }
          if (nounit)
          {
            temp /= a[j + j * a_dim1];
          }
          x[jx] = temp;
          jx += *incx;
          /* L120: */
        }
      }
    }
    else
    {
      if (*incx == 1)
      {
        for (j = *n; j >= 1; --j)
        {
          temp = x[j];
          i__1 = j + 1;
          for (i = *n; i >= i__1; --i)
          {
            temp -= a[i + j * a_dim1] * x[i];
            /* L130: */
          }
          if (nounit)
          {
            temp /= a[j + j * a_dim1];
          }
          x[j] = temp;
          /* L140: */
        }
      }
      else
      {
        kx += (*n - 1) * *incx;
        jx = kx;
        for (j = *n; j >= 1; --j)
        {
          temp = x[jx];
          ix = kx;
          i__1 = j + 1;
          for (i = *n; i >= i__1; --i)
          {
            temp -= a[i + j * a_dim1] * x[ix];
            ix -= *incx;
            /* L150: */
          }
          if (nounit)
          {
            temp /= a[j + j * a_dim1];
          }
          x[jx] = temp;
          jx -= *incx;
          /* L160: */
        }
      }
    }
  }

  return 0;

  /*     End of DTRSV . */

} /* dtrsv_ */

integer idamax_(n, dx, incx)
integer *n;
doublereal *dx;
integer *incx;
{
  /* System generated locals */
  integer ret_val, i__1;
  doublereal d__1;

  /* Local variables */
  static doublereal dmax_;
  static integer i, ix;


  /*     finds the index of element having max. absolute value. */
  /*     jack dongarra, linpack, 3/11/78. */
  /*     modified 3/93 to return if incx .le. 0. */


  /* Parameter adjustments */
  --dx;

  /* Function Body */
  ret_val = 0;
  if (*n < 1 || *incx <= 0)
  {
    return ret_val;
  }
  ret_val = 1;
  if (*n == 1)
  {
    return ret_val;
  }
  if (*incx == 1)
  {
    goto L20;
  }

  /*        code for increment not equal to 1 */

  ix = 1;
  dmax_ = abs(dx[1]);
  ix += *incx;
  i__1 = *n;
  for (i = 2; i <= i__1; ++i)
  {
    if ((d__1 = dx[ix], abs(d__1)) <= dmax_)
    {
      goto L5;
    }
    ret_val = i;
    dmax_ = (d__1 = dx[ix], abs(d__1));
L5:
    ix += *incx;
    /* L10: */
  }
  return ret_val;

  /*        code for increment equal to 1 */

L20:
  dmax_ = abs(dx[1]);
  i__1 = *n;
  for (i = 2; i <= i__1; ++i)
  {
    if ((d__1 = dx[i], abs(d__1)) <= dmax_)
    {
      goto L30;
    }
    ret_val = i;
    dmax_ = (d__1 = dx[i], abs(d__1));
L30:
    ;
  }
  return ret_val;
} /* idamax_ */

