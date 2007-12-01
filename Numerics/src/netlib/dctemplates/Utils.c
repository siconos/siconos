/* Utils.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct
{
  doublereal a[40000], m[200];
} system_;

#define system_1 system_

struct
{
  integer n, lda;
} matdim_;

#define matdim_1 matdim_

struct
{
  char curpform[5];
} forms_;

#define forms_1 forms_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;


/*  This file contains routines used by Jacobi, SOR, and Chebyshev: */

/*  Jacobi/SOR: */

/*     MATSPLIT  calls specific matrix splitting routine */
/*     JACSPLIT */
/*     SORSPLIT */
/*     BACKSOLVE */

/*  Chebyshev: */

/*     GETEIG    computes eigenvalue of iteration matrix. */

/*     =========================================================== */
/* Subroutine */
int matsplit_(omega, b, work, ldw, method, flag_, method_len,
              flag_len)
doublereal *omega, *b, *work;
integer *ldw;
char *method, *flag_;
ftnlen method_len;
ftnlen flag_len;
{
  /* System generated locals */
  integer work_dim1, work_offset;

  /* Builtin functions */
  integer s_wsle(), do_lio(), e_wsle();
  /* Subroutine */
  int s_stop();

  /* Local variables */
  extern /* Subroutine */ int jacsplit_(), sorsplit_();
  extern logical lsame_();

  /* Fortran I/O blocks */
  static cilist io___1 = { 0, 6, 0, 0, 0 };



  /*     .. */
  /*     .. Parameters .. */

  /*     MAXDIM2 = MAXDIM*MAXDIM. */


  /*     .. Common Blocks .. */


  /* Parameter adjustments */
  work_dim1 = *ldw;
  work_offset = work_dim1 + 1;
  work -= work_offset;
  --b;

  /* Function Body */
  if (lsame_(method, "JACOBI", 1L, 6L))
  {
    jacsplit_(&matdim_1.n, system_1.a, &matdim_1.lda, &work[work_offset],
              ldw, flag_, 1L);
  }
  else if (lsame_(method, "SOR", 1L, 3L))
  {
    sorsplit_(omega, &matdim_1.n, system_1.a, &matdim_1.lda, &b[1], &work[
                work_offset], ldw, flag_, 1L);
  }
  else
  {
    s_wsle(&io___1);
    do_lio(&c__9, &c__1, "ERROR: UNKNOW METHOD. QUITTING...", 33L);
    e_wsle();
    s_stop("", 0L);
  }

  return 0;

} /* matsplit_ */


/*     =========================================================== */
/* Subroutine */ int jacsplit_(n, a, lda, work, ldw, flag_, flag_len)
integer *n;
doublereal *a;
integer *lda;
doublereal *work;
integer *ldw;
char *flag_;
ftnlen flag_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2;

  /* Builtin functions */
  integer s_wsle(), do_lio(), e_wsle();
  /* Subroutine */
  int s_stop();

  /* Local variables */
  static integer i, j;
  extern logical lsame_();

  /* Fortran I/O blocks */
  static cilist io___4 = { 0, 6, 0, 0, 0 };




  /* Parameter adjustments */
  work_dim1 = *ldw;
  work_offset = work_dim1 + 1;
  work -= work_offset;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  if (lsame_(flag_, "SPLIT", 1L, 5L))
  {
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
      work[i + work_dim1] = 1. / a[i + i * a_dim1];
      a[i + i * a_dim1] = 0.;
      i__2 = *n;
      for (j = 1; j <= i__2; ++j)
      {
        a[i + j * a_dim1] = -a[i + j * a_dim1];
        /* L10: */
      }
      /* L20: */
    }
  }
  else if (lsame_(flag_, "RECONSTRUCT", 1L, 11L))
  {
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
      i__2 = *n;
      for (j = 1; j <= i__2; ++j)
      {
        a[i + j * a_dim1] = -a[i + j * a_dim1];
        /* L30: */
      }
      a[i + i * a_dim1] = 1. / work[i + work_dim1];
      /* L40: */
    }
  }
  else
  {
    s_wsle(&io___4);
    do_lio(&c__9, &c__1, "UNKNOWN SPLITTING OPTION. QUITTING...", 37L);
    e_wsle();
    s_stop("", 0L);
  }

  return 0;

} /* jacsplit_ */


/*     =========================================================== */
/* Subroutine */ int sorsplit_(omega, n, a, lda, b, work, ldw, flag_,
                               flag_len)
doublereal *omega;
integer *n;
doublereal *a;
integer *lda;
doublereal *b, *work;
integer *ldw;
char *flag_;
ftnlen flag_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2;

  /* Builtin functions */
  integer s_wsle(), do_lio(), e_wsle();
  /* Subroutine */
  int s_stop();

  /* Local variables */
  static integer i, j;
  extern logical lsame_();

  /* Fortran I/O blocks */
  static cilist io___7 = { 0, 6, 0, 0, 0 };




  /* Parameter adjustments */
  work_dim1 = *ldw;
  work_offset = work_dim1 + 1;
  work -= work_offset;
  --b;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  if (lsame_(flag_, "SPLIT", 3L, 5L))
  {

    /*        Set M. */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
      work[i + i * work_dim1] = a[i + i * a_dim1];
      i__2 = i - 1;
      for (j = 1; j <= i__2; ++j)
      {
        work[i + j * work_dim1] = *omega * a[i + j * a_dim1];
        /* L10: */
      }
      /* L20: */
    }

    /*        Set NN and B. */

    /*        Temporarily store the matrix A in order to reconstruct */
    /*        the original matrix. Because the lower triangular portion */
    /*        of A must be zeroed, this is the easiest way to deal with it
    . */
    /*        This causes the requirement that WORK be N x (2N+3). */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
      i__2 = *n;
      for (j = 1; j <= i__2; ++j)
      {
        work[i + (j + *n + 3) * work_dim1] = a[i + j * a_dim1];
        /* L30: */
      }
      /* L40: */
    }

    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
      b[i] = *omega * b[i];
      a[i + i * a_dim1] = (1. - *omega) * a[i + i * a_dim1];
      i__2 = *n;
      for (j = i + 1; j <= i__2; ++j)
      {
        a[i + j * a_dim1] = -(*omega) * a[i + j * a_dim1];
        /* L50: */
      }
      /* L60: */
    }

    i__1 = *n;
    for (i = 2; i <= i__1; ++i)
    {
      i__2 = i - 1;
      for (j = 1; j <= i__2; ++j)
      {
        a[i + j * a_dim1] = 0.;
        /* L70: */
      }
      /* L80: */
    }

  }
  else if (lsame_(flag_, "RECONSTRUCT", 3L, 11L))
  {
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
      b[i] /= *omega;
      i__2 = *n;
      for (j = 1; j <= i__2; ++j)
      {
        a[i + j * a_dim1] = work[i + (j + *n + 3) * work_dim1];
        /* L90: */
      }
      /* L100: */
    }

  }
  else
  {
    s_wsle(&io___7);
    do_lio(&c__9, &c__1, "UNKNOWN SPLITTING OPTION. QUITTING...", 37L);
    e_wsle();
    s_stop("", 0L);
  }
  return 0;

} /* sorsplit_ */


/*     ======================================================== */
/* Subroutine */ int backsolve_(n, a, lda, x)
integer *n;
doublereal *a;
integer *lda;
doublereal *x;
{
  /* System generated locals */
  integer a_dim1, a_offset;

  /* Local variables */
  extern /* Subroutine */ int dtrsv_();


  /*     .. Argument Declarations .. */

  /*     Mask to BLAS routine. X overwritten with inv(A)*X */

  /* Parameter adjustments */
  --x;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  dtrsv_("LOWER", "NOTRANS", "NONUNIT", n, &a[a_offset], lda, &x[1], &c__1,
         5L, 7L, 7L);

  return 0;

} /* backsolve_ */


/*     =========================================================== */
/* Subroutine */ int geteig_(work, ldw, eigmax, eigmin)
doublereal *work;
integer *ldw;
doublereal *eigmax, *eigmin;
{
  /* System generated locals */
  integer work_dim1, work_offset, i__1;

  /* Builtin functions */
  integer s_wsle(), e_wsle(), do_lio();

  /* Local variables */
  static integer info, i;
  extern /* Subroutine */ int dsyev_();
  extern doublereal dlamch_();
  extern logical lsamen_();
  extern /* Subroutine */ int matcopy_();
  extern doublereal matnorm_();

  /* Fortran I/O blocks */
  static cilist io___10 = { 0, 6, 0, 0, 0 };
  static cilist io___11 = { 0, 6, 0, 0, 0 };
  static cilist io___12 = { 0, 6, 0, 0, 0 };
  static cilist io___13 = { 0, 6, 0, 0, 0 };
  static cilist io___14 = { 0, 6, 0, 0, 0 };
  static cilist io___15 = { 0, 6, 0, 0, 0 };



  /*     .. Argument Declarations .. */

  /*     This routine using an LAPACK routine for computing all the */
  /*     eigenvalues of the matrix A. This is for testing purposes only, */
  /*     as this is more expensive than a direct dense solver for the */
  /*     linear system. This is for testing purposes only. */
  /*     .. */
  /*     .. Parameters .. */

  /*     MAXDIM2 = MAXDIM*MAXDIM. */


  /*     .. Common Blocks .. */

  /*     .. */
  /*     .. Local Scalars .. */

  /*     .. Executable Statements .. */

  /*     As the matrix A is overwritten in the following routine, we */
  /*     copy it to temporary workspace. */

  /* Parameter adjustments */
  work_dim1 = *ldw;
  work_offset = work_dim1 + 1;
  work -= work_offset;

  /* Function Body */
  matcopy_(&matdim_1.n, system_1.a, &matdim_1.lda, &work[work_offset], ldw);
  if (lsamen_(&c__3, forms_1.curpform, "JACBI", 5L, 5L))
  {
    i__1 = matdim_1.n;
    for (i = 1; i <= i__1; ++i)
    {
      work[i + i * work_dim1] /= system_1.m[i - 1];
      /* L30: */
    }
  }

  /*     Call LAPACK eigenvalue routine. */

  i__1 = matdim_1.n * 3 - 1;
  dsyev_("NO_VEC", "UPPER", &matdim_1.n, &work[work_offset], ldw, &work[(
           matdim_1.n + 1) * work_dim1 + 1], &work[(matdim_1.n + 2) *
               work_dim1 + 1], &i__1, &info, 6L, 5L);

  if (info != 0)
  {
    s_wsle(&io___10);
    e_wsle();
    s_wsle(&io___11);
    do_lio(&c__9, &c__1, "CHEBYSHEV WARNING: DSYEV COULD NOT COMPUTE ALL\
 EIGENVALUES.", 59L);
    e_wsle();
    s_wsle(&io___12);
    do_lio(&c__9, &c__1, "SETTING EIGMIN/MAX TO DEFAULT VALUES EPS AND |\
A|", 48L);
    e_wsle();
    s_wsle(&io___13);
    e_wsle();
    *eigmin = dlamch_("EPS", 3L);
    *eigmax = matnorm_(&matdim_1.n, system_1.a, &matdim_1.lda);
    return 0;
  }
  else
  {
    *eigmin = work[(matdim_1.n + 1) * work_dim1 + 1];
    *eigmax = work[matdim_1.n + (matdim_1.n + 1) * work_dim1];
  }

  /*     Eigenvalues should be positive. */

  if (*eigmin < 0.)
  {
    s_wsle(&io___14);
    do_lio(&c__9, &c__1, "CHEBYSHEV WARNING: COMPUTED MIN EIGENVALUE <= \
0: SET TO EPSILON", 63L);
    e_wsle();
    *eigmin = dlamch_("EPS", 3L);
  }
  if (*eigmax < *eigmin)
  {
    s_wsle(&io___15);
    do_lio(&c__9, &c__1, "CHEBYSHEV WARNING: MAX EIGENVALUE < MIN: SET T\
O |A|", 51L);
    e_wsle();
    *eigmax = matnorm_(&matdim_1.n, system_1.a, &matdim_1.lda);
  }

  return 0;

} /* geteig_ */


/*     ================================================================ */
/* Subroutine */ int matcopy_(n, a, lda, b, ldb)
integer *n;
doublereal *a;
integer *lda;
doublereal *b;
integer *ldb;
{
  /* System generated locals */
  integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

  /* Local variables */
  static integer i, j;




  /* Parameter adjustments */
  b_dim1 = *ldb;
  b_offset = b_dim1 + 1;
  b -= b_offset;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  i__1 = *n;
  for (j = 1; j <= i__1; ++j)
  {
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
      b[i + j * b_dim1] = a[i + j * a_dim1];
      /* L10: */
    }
    /* L20: */
  }

  return 0;

} /* matcopy_ */

