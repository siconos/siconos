#include "SiconosConfig.h"

#include "SiconosLapack.h"


/*  -- LAPACK routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRTRS solves a triangular system of the form */

/*     A * X = B  or  A**T * X = B, */

/*  where A is a triangular matrix of order N, and B is an N-by-NRHS */
/*  matrix.  A check is made to verify that A is nonsingular. */

/*  Arguments */
/*  ========= */

/*  UPLO    (input) CHARACTER*1 */
/*          = 'U':  A is upper triangular; */
/*          = 'L':  A is lower triangular. */

/*  TRANS   (input) CHARACTER*1 */
/*          Specifies the form of the system of equations: */
/*          = 'N':  A * X = B  (No transpose) */
/*          = 'T':  A**T * X = B  (Transpose) */
/*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) */

/*  DIAG    (input) CHARACTER*1 */
/*          = 'N':  A is non-unit triangular; */
/*          = 'U':  A is unit triangular. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  NRHS    (input) INTEGER */
/*          The number of right hand sides, i.e., the number of columns */
/*          of the matrix B.  NRHS >= 0. */

/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*          The triangular matrix A.  If UPLO = 'U', the leading N-by-N */
/*          upper triangular part of the array A contains the upper */
/*          triangular matrix, and the strictly lower triangular part of */
/*          A is not referenced.  If UPLO = 'L', the leading N-by-N lower */
/*          triangular part of the array A contains the lower triangular */
/*          matrix, and the strictly upper triangular part of A is not */
/*          referenced.  If DIAG = 'U', the diagonal elements of A are */
/*          also not referenced and are assumed to be 1. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS) */
/*          On entry, the right hand side matrix B. */
/*          On exit, if INFO = 0, the solution matrix X. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B.  LDB >= max(1,N). */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0: if INFO = -i, the i-th argument had an illegal value */
/*          > 0: if INFO = i, the i-th diagonal element of A is zero, */
/*               indicating that the matrix is singular and the solutions */
/*               X have not been computed. */

/*  ===================================================================== */

#if defined(HAS_ATLAS_LAPACK)

int clapack_dtrtrs(const enum ATLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum ATLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const enum CBLAS_DIAG Diag, int n, int nrhs, double *a, int lda, double *b, int ldb)
{
  /* Table of constant values */
  static double c_b12 = 1.;

  /* System generated locals */
  int a_dim1, a_offset = 0, /* b_dim1,*/ b_offset = 0, i__1, info;
  /* Local variables */
  int nounit;


  a_dim1 = lda;

  /* ?
     a_dim1 = lda;
     a_offset = 1 + a_dim1 * 1;
     a -= a_offset;
     b_dim1 = ldb;
     b_offset = 1 + b_dim1 * 1;
     b -= b_offset;
  */

  /* Function Body */
  info = 0;
  nounit = (Diag == CblasNonUnit);
  if (!(Uplo == CblasUpper) && !(Uplo == CblasLower))
  {
    info = -1;
  }
  else if (!(Trans == CblasNoTrans) && !(Trans == CblasTrans) && !(Trans == CblasConjTrans))
  {
    info = -2;
  }
  else if (! nounit && !(Uplo == CblasUpper))
  {
    info = -3;
  }
  else if (n < 0)
  {
    info = -4;
  }
  else if (nrhs < 0)
  {
    info = -5;
  }
  else if (lda < max(1, n))
  {
    info = -7;
  }
  else if (ldb < max(1, n))
  {
    info = -9;
  }
  if (info != 0)
  {
    i__1 = -(info);
    return info;
  }

  /*     Quick return if possible */

  if (n == 0)
  {
    return 0;
  }

  /*     Check for singularity. */

  if (nounit)
  {
    i__1 = n;
    for (info = 0; info < i__1; ++info)
    {
      if (a[info + info * a_dim1] == 0.)
      {
        return ++info;
      }
      /* L10: */
    }
  }
  info = 0;

  /*     Solve A * x = b  or  A' * x = b. */

  cblas_dtrsm(Order, Side, Uplo, Trans, Diag, n, nrhs, c_b12, &a[a_offset], lda, &b[b_offset], ldb);

  return info;
}
#else

#endif
/* clapack_dtrtrs */
