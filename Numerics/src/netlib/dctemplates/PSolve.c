/* PSolve.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct
{
  char curpform[5];
} forms_;

#define forms_1 forms_

struct
{
  integer n, lda;
} matdim_;

#define matdim_1 matdim_

struct
{
  doublereal a[40000], m[200];
} system_;

#define system_1 system_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;


/*  This file contains the preconditioner solve routines: */

/*  PSOLVE and PSOLVETRANS call the appropriate solver: */

/*     PSOLVENONE and PSOLVENONETRANS for using no preconditioning. */

/*     PSOLVEJAC and PSOLVEJACTRANS for Jacobi preconditioning. */

/*     Also included are the solvers for QMR which require left and right */
/*     preconditioning: PSOLVEQ and PSOLVETRANSQ */

/*     ================================================================ */
int psolve_(x, b)

doublereal *x, *b;
{
  /* Builtin functions */
  integer s_wsle(), do_lio(), e_wsle();
  /* Subroutine */
  int s_stop();

  /* Local variables */
  extern logical lsame_();
  extern /* Subroutine */ int psolvejac_(), psolvenone_();

  /* Fortran I/O blocks */
  static cilist io___1 = { 0, 6, 0, 0, 0 };



  /*     .. Array Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Common Blocks .. */

  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --b;
  --x;

  /* Function Body */
  if (lsame_(forms_1.curpform, "IDENT", 5L, 5L))
  {
    psolvenone_(&x[1], &b[1]);
  }
  else if (lsame_(forms_1.curpform, "JACBI", 5L, 5L))
  {
    psolvejac_(&x[1], &b[1]);
  }
  else
  {
    s_wsle(&io___1);
    do_lio(&c__9, &c__1, "IN PSOLVE: UNKNOWN PRECONDITIONER", 33L);
    do_lio(&c__9, &c__1, forms_1.curpform, 5L);
    do_lio(&c__9, &c__1, " QUITTING", 9L);
    e_wsle();
    s_stop("", 0L);
  }

  return 0;

} /* psolve_ */

/*     ================================================================ */

int psolvetrans_(x, b)
doublereal *x, *b;
{
  /* Builtin functions */
  integer s_wsle(), do_lio(), e_wsle();
  /* Subroutine */
  int s_stop();

  /* Local variables */
  extern /* Subroutine */ int psolvenonetrans_();
  extern logical lsame_();
  extern /* Subroutine */ int psolvejactrans_();

  /* Fortran I/O blocks */
  static cilist io___2 = { 0, 6, 0, 0, 0 };

  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --b;
  --x;

  /* Function Body */
  if (lsame_(forms_1.curpform, "IDENT", 5L, 5L))
  {
    psolvenonetrans_(&x[1], &b[1]);
  }
  else if (lsame_(forms_1.curpform, "JACBI", 5L, 5L))
  {
    psolvejactrans_(&x[1], &b[1]);
  }
  else
  {
    s_wsle(&io___2);
    do_lio(&c__9, &c__1, "IN PSOLVE: UNKNOWN PRECONDITIONER", 33L);
    do_lio(&c__9, &c__1, forms_1.curpform, 5L);
    do_lio(&c__9, &c__1, " QUITTING", 9L);
    e_wsle();
    s_stop("", 0L);
  }

  return 0;

} /* psolvetrans_ */

/*     ================================================================ */

int psolvenone_(x, b)
doublereal *x, *b;
{
  extern  int dcopy_();

  /*  Purpose */
  /*  ======= */

  /*  This PSOLVE is for the unpreconditioned version, i.e. just does */
  /*  a vector copy ( B to X ) then returns. */

  /*  Arguments */
  /*  ========= */

  /*  B       (input) DOUBLE PRECISION array, dimension N. */
  /*          On entry, right hand side vector B. */
  /*          Unchanged on exit. */

  /*  X       (output) DOUBLE PRECISION array, dimension N. */
  /*          Set to solution on output. */

  /*  BLAS:  DCOPY */
  /*  ============================================================ */

  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --b;
  --x;

  /* Function Body */
  dcopy_(&matdim_1.n, &b[1], &c__1, &x[1], &c__1);

  return 0;

  /*     End of PSolveNone */

}


/*     ===================================================== */
int psolvenonetrans_(x, b)
doublereal *x, *b;
{
  extern /* Subroutine */ int dcopy_();

  /*  Purpose */
  /*  ======= */

  /*  This PSOLVE is for the unpreconditioned version, i.e. just does */
  /*  a vector copy ( B to X ) then returns. */

  /*  Arguments */
  /*  ========= */

  /*  B       (input) DOUBLE PRECISION array, dimension N. */
  /*          On entry, right hand side vector B. */
  /*          Unchanged on exit. */

  /*  X       (output) DOUBLE PRECISION array, dimension N. */
  /*          Set to solution on output. */

  /*  BLAS:  DCOPY */
  /*  ============================================================ */

  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --b;
  --x;

  /* Executable Statements */
  dcopy_(&matdim_1.n, &b[1], &c__1, &x[1], &c__1);

  return 0;

  /*     End of PSolve */

}


/*     =========================================================== */
int psolvejac_(x, b)
doublereal *x, *b;
{
  /* System generated locals */
  integer i__1;

  /* Local variables */
  static integer i;

  /*  Purpose */
  /*  ======= */

  /*  PSOLVE solves the linear system Mx = b where matrix M has */
  /*  is diagonal. */

  /*  Arguments */
  /*  ========= */

  /*  B       (input) DOUBLE PRECISION array, dimension N. */
  /*          On entry, right hand side vector B. */
  /*          Unchanged on exit. */

  /*  X       (output) DOUBLE PRECISION array, dimension N. */
  /*          Set to solution on output. */
  /*  ============================================================ */

  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --b;
  --x;

  /* Function Body */
  i__1 = matdim_1.n;
  for (i = 1; i <= i__1; ++i)
  {
    x[i] = b[i] / system_1.m[i - 1];
    /* L10: */
  }

  return 0;

  /*     End of PSolveJac */

}


/*     ========================================================= */
int psolvejactrans_(x, b)

doublereal *x, *b;
{
  extern /* Subroutine */ int psolvejac_();

  /*  Purpose */
  /*  ======= */

  /*  PSOLVETRANS solves the linear system Mx = b where matrix M has */
  /*  is diagonal. Since this is the same as the non-transpose version, */
  /*  this routine is actual just a mask to PSOLVEJAC. */

  /*  Arguments */
  /*  ========= */

  /*  B       (input) DOUBLE PRECISION array, dimension N. */
  /*          On entry, right hand side vector B. */
  /*          Unchanged on exit. */

  /*  X       (output) DOUBLE PRECISION array, dimension N. */
  /*          Set to solution on output. */
  /*  ============================================================ */

  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --b;
  --x;

  /* Function Body */
  psolvejac_(&x[1], &b[1]);

  return 0;

  /*     End of PSolveJacTrans */

}

/*     ================================================================ */
/*     Following are the solvers for QMR, allowing left and right */
/*     preconditioning. */
/*     ================================================================ */
int psolveq_(x, b, which, which_len)
doublereal *x, *b;
char *which;
ftnlen which_len;
{
  /* Builtin functions */
  integer s_wsle(), do_lio(), e_wsle();
  /* Subroutine */
  int s_stop();

  /* Local variables */
  extern logical lsame_();
  extern /* Subroutine */ int psolvejac_(), psolvenone_();

  /* Fortran I/O blocks */
  static cilist io___4 = { 0, 6, 0, 0, 0 };

  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --b;
  --x;

  /* Function Body */
  if (lsame_(forms_1.curpform, "IDENT", 5L, 5L))
  {
    psolvenone_(&x[1], &b[1]);
  }
  else if (lsame_(forms_1.curpform, "JACBI", 5L, 5L))
  {
    if (lsame_(which, "LEFT", 4L, 4L))
    {
      psolvejac_(&x[1], &b[1]);
    }
    else
    {
      psolvenone_(&x[1], &b[1]);
    }
  }
  else
  {
    s_wsle(&io___4);
    do_lio(&c__9, &c__1, "IN PSOLVEQ: UNKNOWN PRECONDITIONER", 34L);
    do_lio(&c__9, &c__1, forms_1.curpform, 5L);
    do_lio(&c__9, &c__1, " QUITTING", 9L);
    e_wsle();
    s_stop("", 0L);
  }

  return 0;

}

/*     ================================================================ */

int psolvetransq_(x, b, which, which_len)
doublereal *x, *b;
char *which;
ftnlen which_len;
{
  /* Builtin functions */
  integer s_wsle(), do_lio(), e_wsle();
  /* Subroutine */
  int s_stop();

  /* Local variables */
  extern logical lsame_();
  extern /* Subroutine */ int psolvejac_(), psolvenone_();

  /* Fortran I/O blocks */
  static cilist io___5 = { 0, 6, 0, 0, 0 };

  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --b;
  --x;

  /* Executable Statements */
  if (lsame_(forms_1.curpform, "IDENT", 5L, 5L))
  {
    psolvenone_(&x[1], &b[1]);
  }
  else if (lsame_(forms_1.curpform, "JACBI", 5L, 5L))
  {
    if (lsame_(which, "LEFT", 4L, 4L))
    {
      psolvejac_(&x[1], &b[1]);
    }
    else
    {
      psolvenone_(&x[1], &b[1]);
    }
  }
  else
  {
    s_wsle(&io___5);
    do_lio(&c__9, &c__1, "IN PSOLVEQ: UNKNOWN PRECONDITIONER", 34L);
    do_lio(&c__9, &c__1, forms_1.curpform, 5L);
    do_lio(&c__9, &c__1, " QUITTING", 9L);
    e_wsle();
    s_stop("", 0L);
  }

  return 0;

}
