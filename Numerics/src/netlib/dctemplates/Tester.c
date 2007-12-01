/* Tester.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct
{
  integer n, lda;
} matdim_;

#define matdim_1 matdim_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__6 = 6;
static integer c_b23 = 100000;
static integer c__9 = 9;
static integer c__3 = 3;
static doublereal c_b125 = -1.;
static doublereal c_b127 = 1.;
static integer c__4 = 4;
static integer c__2 = 2;


/* Main program */
MAIN__()
{
  /* Initialized data */

  static char snames[6 * 9 + 1] = "CG    CHEBY SOR   BICG  CGS   BICGS GMRES Q\
MR   JACOB ";

  /* Format strings */
  static char fmt_998[] = "(a6,l2)";
  static char fmt_999[] = "(\002 SUBPROGRAM NAME \002,a6,\002 NOT RECOGNIZ\
ED\002,/\002 ******* T\002,\002ESTS ABANDONED *******\002)";
  static char fmt_900[] = "(\002   SYMMETRIC POSITIVE DEFINITE ROUTINES PA\
SSED. (\002,i3,\002 TESTS )\002)";
  static char fmt_901[] = "(\002   NONSYMMETRIC ROUTINES PASSED. (\002,i3\
,\002 TESTS )\002)";
  static char fmt_910[] = "(\002 PASSED FOR SYMMETRIC POSITIVE DEFINITE MA\
TRICES (\002,i3,\002 TESTS )\002)";
  static char fmt_911[] = "(\002 SYMMETRIC POSITIVE DEFINITE MATRICES: \
(\002,i3,\002 TESTS )\002)";
  static char fmt_990[] = "(\002 THERE ARE\002,i3,\002 CRITICAL ERRORS A\
ND \002,i3,\002 SUSPICIOUS RESULTS.\002)";
  static char fmt_991[] = "(\002 SEE FILE test.results FOR ADDITIONAL INFO\
RMATION\002)";
  static char fmt_920[] = "(\002 PASSED FOR NONSYMMETRIC MATRICES (\002,\
i3,\002 TESTS )\002)";
  static char fmt_921[] = "(\002 NONSYMMETRIC MATRICES: (\002,i3,\002 TEST\
S )\002)";

  /* System generated locals */
  integer i__1;
  olist o__1;
  cllist cl__1;

  /* Builtin functions */
  /* Subroutine */
  int s_copy();
  integer f_open(), s_rsle(), do_lio(), e_rsle(), s_rsfe(), do_fio(),
          e_rsfe(), s_wsfe(), e_wsfe();
  /* Subroutine */
  int s_stop();
  integer f_clos(), s_wsle(), e_wsle();

  /* Local variables */
  static integer suspnsy;
  extern /* Subroutine */ int psolvetrans_();
  static doublereal work[100000];
  static integer spdtests;
  static doublereal b[200];
  static integer i, nsytests;
  static doublereal x[1800] /* was [200][9] */;
  extern /* Subroutine */ int psolvetransq_();
  static char pform[5 * 2];
  extern /* Subroutine */ int backsolve_();
  static doublereal scaledtol;
  static logical ltest[9];
  static doublereal x0[200];
  extern /* Subroutine */ int header_();
  extern logical lsamen_();
  extern /* Subroutine */ int matvec_();
  static char snamet[6];
  extern /* Subroutine */ int footer_();
  static logical spdres;
  extern /* Subroutine */ int psolve_();
  static logical ltestt, nsyres;
  static integer ldw, ldx;
  static doublereal tol;
  extern /* Subroutine */ int dspdchk_(), dnsychk_();
  static integer critspd;
  extern /* Subroutine */ int psolveq_();
  static integer critnsy, suspspd;
  extern /* Subroutine */ int matvectrans_();

  /* Fortran I/O blocks */
  static cilist io___7 = { 0, 9, 0, 0, 0 };
  static cilist io___9 = { 0, 9, 0, 0, 0 };
  static cilist io___13 = { 0, 9, 0, fmt_998, 0 };
  static cilist io___16 = { 0, 6, 0, fmt_999, 0 };
  static cilist io___24 = { 0, 9, 1, fmt_998, 0 };
  static cilist io___25 = { 0, 6, 0, fmt_999, 0 };
  static cilist io___29 = { 0, 6, 0, 0, 0 };
  static cilist io___30 = { 0, 6, 0, 0, 0 };
  static cilist io___31 = { 0, 6, 0, 0, 0 };
  static cilist io___32 = { 0, 6, 0, fmt_900, 0 };
  static cilist io___33 = { 0, 6, 0, fmt_901, 0 };
  static cilist io___34 = { 0, 6, 0, 0, 0 };
  static cilist io___35 = { 0, 6, 0, 0, 0 };
  static cilist io___36 = { 0, 6, 0, fmt_910, 0 };
  static cilist io___37 = { 0, 6, 0, 0, 0 };
  static cilist io___38 = { 0, 6, 0, 0, 0 };
  static cilist io___39 = { 0, 6, 0, fmt_911, 0 };
  static cilist io___40 = { 0, 6, 0, fmt_990, 0 };
  static cilist io___41 = { 0, 6, 0, fmt_991, 0 };
  static cilist io___42 = { 0, 6, 0, 0, 0 };
  static cilist io___43 = { 0, 6, 0, 0, 0 };
  static cilist io___44 = { 0, 6, 0, 0, 0 };
  static cilist io___45 = { 0, 6, 0, fmt_920, 0 };
  static cilist io___46 = { 0, 6, 0, 0, 0 };
  static cilist io___47 = { 0, 6, 0, 0, 0 };
  static cilist io___48 = { 0, 6, 0, 0, 0 };
  static cilist io___49 = { 0, 6, 0, fmt_921, 0 };
  static cilist io___50 = { 0, 6, 0, fmt_990, 0 };
  static cilist io___51 = { 0, 6, 0, fmt_991, 0 };
  static cilist io___52 = { 0, 6, 0, 0, 0 };



  /*  Test program for the DOUBLE PRECISION iterative templates. */

  /*  The program must be driven by a short data file. An annotated example
  */
  /*  of a data file as follows (save to file test.data, delete the first */
  /*  3 columns): */

  /*  1.0D-15                        CONVERGENCE TOLERANCE */
  /*  10                             SCALED RESIDUAL TOLERANCE */
  /*  CG     T PUT F FOR NO TEST.    ALGORITHMS TO BE TESTED */
  /*  CHEBY  T PUT F FOR NO TEST. */
  /*  SOR    T PUT F FOR NO TEST. */
  /*  BICG   T PUT F FOR NO TEST. */
  /*  CGS    T PUT F FOR NO TEST. */
  /*  BICGS  T PUT F FOR NO TEST. */
  /*  GMRES  T PUT F FOR NO TEST. */
  /*  QMR    T PUT F FOR NO TEST. */
  /*  JACOB  T PUT F FOR NO TEST. */
  /*  3                              NUMBER OF SPD MATRICES TO BE GENERATED
  */
  /*  WATH  2, 2, 1, ONES, ZERO      MATRIX, NX, NY, NZ, RHS, INITIAL GUESS
  */
  /*  F2SH  6, 6, 1, SUMR, ZERO */
  /*  F3SH  3, 3, 3, ONES, ZERO */
  /*  BICG   T  PUT F FOR NO TEST.   ALGORITHMS TO BE TESTED */
  /*  CGS    T  PUT F FOR NO TEST. */
  /*  BICGS  T  PUT F FOR NO TEST. */
  /*  GMRES  T  PUT F FOR NO TEST. */
  /*  QMR    T  PUT F FOR NO TEST. */
  /*  4                              NUMBER OF MATRICES TO BE GENERATED */
  /*  PDE1, 5, 5, 5, SUMR , ZERO     MATRIX, NX, NY, NZ, RHS, INITIAL GUESS
  */
  /*  PDE2, 5, 5, 5, SUMR , ZERO */
  /*  PDE3, 5, 5, 5, ONES , ZERO */
  /*  PDE4, 6, 6, 1, ONES , ZERO */

  /*  See: */

  /*     Barrett, Berry, Chan, Demmel, Donato, Dongarra, */
  /*     Eijkhout, Pozo, Romine, and van der Vorst. */
  /*     Templates for the Solution of Linear Systems: Building Blocks */
  /*     for Iterative Methods, SIAM Publications, 1993. */
  /*     (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps). */

  /*  -- Written on 1-November-1993. */
  /*     Richard Barrett, University of Tennessee */
  /*     Jack Dongarra, Univ. of Tennessee and Oak Ridge National */
  /*     Laboratory */

  /*     .. Parameters .. */
  /*     MAXLEN must be greater than or equal to (2N**2)+3N, i.e. WORK */
  /*     must have dimension N x (2N+3). This is for SOR (see StatUtils */
  /*     for details). Chebyshev requires N*2. For workspace requirements */
  /*     of the algorithms, see the individial template. */

  /*     .. */
  /*     .. Scalar Declarations .. */
  /*     .. */
  /*     .. Array Declarations .. */
  /*     .. */
  /*     .. Common Blocks .. */
  /*     .. */
  /*     .. External Routines .. */

  /*     .. */
  /*     .. Executable Statements .. */

  /*     Initializations. */

  matdim_1.lda = 200;
  ldx = 200;
  ldw = 200;

  spdres = TRUE_;
  nsyres = TRUE_;

  s_copy(pform, "IDENT", 5L, 5L);
  s_copy(pform + 5, "JACBI", 5L, 5L);

  o__1.oerr = 0;
  o__1.ounit = 9;
  o__1.ofnmlen = 9;
  o__1.ofnm = "test.data";
  o__1.orl = 0;
  o__1.osta = 0;
  o__1.oacc = 0;
  o__1.ofm = 0;
  o__1.oblnk = 0;
  f_open(&o__1);
  o__1.oerr = 0;
  o__1.ounit = 10;
  o__1.ofnmlen = 12;
  o__1.ofnm = "test.results";
  o__1.orl = 0;
  o__1.osta = 0;
  o__1.oacc = 0;
  o__1.ofm = 0;
  o__1.oblnk = 0;
  f_open(&o__1);

  /*     Get the convergence tolerance, the tolerance for the normalized */
  /*     scaled residual, and the number of systems to be generated. */
  /*     and the algorithms to be tested. */

  s_rsle(&io___7);
  do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
  e_rsle();
  s_rsle(&io___9);
  do_lio(&c__5, &c__1, (char *)&scaledtol, (ftnlen)sizeof(doublereal));
  e_rsle();

  /*     Get input data for SPD testing: */
  /*     Read names of subroutines and flags which indicate whether */
  /*     they are to be tested. */

  for (i = 1; i <= 9; ++i)
  {
    ltest[i - 1] = FALSE_;
    /* L10: */
  }
L20:
  s_rsfe(&io___13);
  do_fio(&c__1, snamet, 6L);
  do_fio(&c__1, (char *)&ltestt, (ftnlen)sizeof(logical));
  e_rsfe();
  for (i = 1; i <= 9; ++i)
  {
    if (lsamen_(&c__6, snamet, snames + (i - 1) * 6, 6L, 6L))
    {
      goto L40;
    }
    /* L30: */
  }
  s_wsfe(&io___16);
  do_fio(&c__1, snamet, 6L);
  e_wsfe();
  s_stop("", 0L);
L40:
  ltest[i - 1] = ltestt;
  if (i < 9)
  {
    goto L20;
  }

  /* L50: */

  /*     Begin testing. */

  header_(&tol);

  /*     Symmetric Positive Definite Routine Tester. */

  dspdchk_(x, &ldx, b, x0, work, &ldw, pform, matvec_, matvectrans_,
           psolve_, psolvetrans_, psolveq_, psolvetransq_, backsolve_, &tol,
           &scaledtol, ltest, &c_b23, &spdres, &spdtests, &suspspd, &critspd,
           5L);

  /*     Get input data for Nonsymmetric testing: */
  /*     Read names of subroutines and flags which indicate whether */
  /*     they are to be tested. */

  for (i = 1; i <= 9; ++i)
  {
    ltest[i - 1] = FALSE_;
    /* L60: */
  }
L70:
  i__1 = s_rsfe(&io___24);
  if (i__1 != 0)
  {
    goto L100;
  }
  i__1 = do_fio(&c__1, snamet, 6L);
  if (i__1 != 0)
  {
    goto L100;
  }
  i__1 = do_fio(&c__1, (char *)&ltestt, (ftnlen)sizeof(logical));
  if (i__1 != 0)
  {
    goto L100;
  }
  i__1 = e_rsfe();
  if (i__1 != 0)
  {
    goto L100;
  }
  for (i = 4; i <= 8; ++i)
  {
    if (lsamen_(&c__6, snamet, snames + (i - 1) * 6, 6L, 6L))
    {
      goto L90;
    }
    /* L80: */
  }
  s_wsfe(&io___25);
  do_fio(&c__1, snamet, 6L);
  e_wsfe();
  s_stop("", 0L);
L90:
  ltest[i - 1] = ltestt;
  if (i < 8)
  {
    goto L70;
  }

L100:

  /*     Nonsymmetric Routine Tester. */

  dnsychk_(x, &ldx, b, x0, work, &ldw, pform, matvec_, matvectrans_,
           psolve_, psolvetrans_, psolveq_, psolvetransq_, backsolve_, &tol,
           &scaledtol, ltest, &c_b23, &nsyres, &nsytests, &suspnsy, &critnsy,
           5L);

  /*     End of testing. */

  footer_();

  cl__1.cerr = 0;
  cl__1.cunit = 9;
  cl__1.csta = 0;
  f_clos(&cl__1);
  cl__1.cerr = 0;
  cl__1.cunit = 10;
  cl__1.csta = 0;
  f_clos(&cl__1);

  /*     Print overall results to screen. */

  s_wsle(&io___29);
  e_wsle();
  if (spdres && nsyres)
  {

    /*        All tests passed. */

    s_wsle(&io___30);
    do_lio(&c__9, &c__1, "TESTS COMPLETE:", 15L);
    e_wsle();
    s_wsle(&io___31);
    e_wsle();
    if (spdtests > 0)
    {
      s_wsfe(&io___32);
      do_fio(&c__1, (char *)&spdtests, (ftnlen)sizeof(integer));
      e_wsfe();
    }
    if (nsytests > 0)
    {
      s_wsfe(&io___33);
      do_fio(&c__1, (char *)&nsytests, (ftnlen)sizeof(integer));
      e_wsfe();
    }
  }
  else
  {
    if (spdres)
    {

      if (spdtests > 0)
      {

        /*              SPD tests passed. */

        s_wsle(&io___34);
        do_lio(&c__9, &c__1, "TESTS COMPLETE:", 15L);
        e_wsle();
        s_wsle(&io___35);
        e_wsle();
        s_wsfe(&io___36);
        do_fio(&c__1, (char *)&spdtests, (ftnlen)sizeof(integer));
        e_wsfe();
      }
      else
      {
        s_wsle(&io___37);
        e_wsle();
        s_wsle(&io___38);
        do_lio(&c__9, &c__1, "SPD TESTING NOT PERFORMED.", 26L);
        e_wsle();
      }
    }
    else
    {

      /*           SPD testing failed. */

      s_wsfe(&io___39);
      do_fio(&c__1, (char *)&spdtests, (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___40);
      do_fio(&c__1, (char *)&critspd, (ftnlen)sizeof(integer));
      do_fio(&c__1, (char *)&suspspd, (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___41);
      e_wsfe();
    }
    s_wsle(&io___42);
    e_wsle();
    if (nsyres)
    {

      if (spdtests > 0)
      {

        /*              Nonsymmetric tests passed. */

        s_wsle(&io___43);
        do_lio(&c__9, &c__1, "TESTS COMPLETE:", 15L);
        e_wsle();
        s_wsle(&io___44);
        e_wsle();
        s_wsfe(&io___45);
        do_fio(&c__1, (char *)&nsytests, (ftnlen)sizeof(integer));
        e_wsfe();
      }
      else
      {
        s_wsle(&io___46);
        e_wsle();
        s_wsle(&io___47);
        do_lio(&c__9, &c__1, "NONSYMMETRIC TESTING NOT PERFORMED.",
               35L);
        e_wsle();
      }
    }
    else
    {

      /*           Nonsymmetric testing failed. */

      s_wsle(&io___48);
      e_wsle();
      s_wsfe(&io___49);
      do_fio(&c__1, (char *)&nsytests, (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___50);
      do_fio(&c__1, (char *)&critnsy, (ftnlen)sizeof(integer));
      do_fio(&c__1, (char *)&suspnsy, (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___51);
      e_wsfe();
    }
  }
  s_wsle(&io___52);
  e_wsle();

  /*     Format statements for screen output of general test results. */






  s_stop("", 0L);

  /*     End of Driver for Testing the Iterative Templates */

} /* MAIN__ */


/*     =============================================================== */
/* Subroutine */ int vecgen_(form, n, a, lda, b, info, form_len)
char *form;
integer *n;
doublereal *a;
integer *lda;
doublereal *b;
integer *info;
ftnlen form_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;

  /* Local variables */
  static integer i, j;
  extern logical lsamen_();
  static doublereal tmp;




  /* Parameter adjustments */
  --b;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  *info = 0;

  if (lsamen_(&c__3, form, "ONES", 4L, 4L))
  {
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
      b[i] = 1.;
      /* L10: */
    }
  }
  else if (lsamen_(&c__3, form, "ZEROS", 4L, 5L))
  {
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
      b[i] = 0.;
      /* L20: */
    }
  }
  else if (lsamen_(&c__3, form, "SUMROW", 4L, 6L))
  {
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
      tmp = 0.;
      i__2 = *n;
      for (j = 1; j <= i__2; ++j)
      {
        tmp += a[i + j * a_dim1];
        /* L30: */
      }
      b[i] = tmp;
      /* L40: */
    }
  }
  else
  {
    *info = -1;
  }

  return 0;

} /* vecgen_ */


/*     =============================================================== */
/* Subroutine */ int precon_(n, a, lda, pform, m, info, pform_len)
integer *n;
doublereal *a;
integer *lda;
char *pform;
doublereal *m;
integer *info;
ftnlen pform_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1;

  /* Builtin functions */
  integer s_wsle(), e_wsle(), do_lio();

  /* Local variables */
  static integer i;
  extern logical lsamen_();

  /* Fortran I/O blocks */
  static cilist io___57 = { 0, 6, 0, 0, 0 };
  static cilist io___58 = { 0, 6, 0, 0, 0 };
  static cilist io___59 = { 0, 6, 0, 0, 0 };
  static cilist io___60 = { 0, 6, 0, 0, 0 };
  static cilist io___61 = { 0, 6, 0, 0, 0 };
  static cilist io___62 = { 0, 6, 0, 0, 0 };



  /*     .. Scalar and Array Declarations .. */


  /*  Purpose: */
  /*  ======= */

  /*  PRECON forms a preconditioner matrix of type PROFRM for */
  /*  iterative solvers of the linear system Ax = b. */

  /*  PFORM: */

  /*      IDENT        identity matrix (for testing) */

  /*      JACBI        diagonal scaling */

  /*     ============================================== */

  /*     .. Local Scalars .. */


  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --m;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  if (lsamen_(&c__5, pform, "IDENT", 5L, 5L))
  {

    /*         Identity matrix need not be formed, since the solve involvi
    ng */
    /*         the preconditioner (PSolve) merely copies the right hand si
    de */
    /*         to the solution vector. */

    return 0;

  }
  else if (lsamen_(&c__5, pform, "JACBI", 5L, 5L))
  {

    /*        Diagonal Scaling: diag(A). Note that we actually form inv(M)
     so that*/
    /*         solver can use multiplication. */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
      m[i] = a[i + i * a_dim1];
      /* L10: */
    }

  }
  else
  {

    /*        Selected preconditioner not implemented */

    s_wsle(&io___57);
    e_wsle();
    s_wsle(&io___58);
    e_wsle();
    s_wsle(&io___59);
    do_lio(&c__9, &c__1, "PRECONDITIONER ", 15L);
    do_lio(&c__9, &c__1, pform, 5L);
    do_lio(&c__9, &c__1, " NOT YET IMPLEMENTED", 20L);
    e_wsle();
    s_wsle(&io___60);
    e_wsle();
    s_wsle(&io___61);
    e_wsle();
    s_wsle(&io___62);
    e_wsle();
    *info = -1;

  }

  return 0;
} /* precon_ */


/*     ================================================================ */
doublereal getbreak_()
{
  /* System generated locals */
  doublereal ret_val, d__1;

  /* Local variables */
  extern doublereal dlamch_();
  static doublereal eps;


  /*     Get breakdown parameter tolerance; for the test routine, */
  /*     set to machine precision. */


  eps = dlamch_("EPS", 3L);
  /* Computing 2nd power */
  d__1 = eps;
  ret_val = d__1 * d__1;

  return ret_val;

} /* getbreak_ */

/*     =============================================================== */
doublereal scaledresid_(anorm, n, x, rk, tol)
doublereal *anorm;
integer *n;
doublereal *x, *rk, *tol;
{
  /* System generated locals */
  doublereal ret_val, d__1;

  /* Local variables */
  static doublereal xnorm;
  extern integer idamax_();
  static doublereal resnorm;


  /*     Returns |B-A*X| / ( |A||X|*N*TOL ), using the infinity norm. */


  /* Parameter adjustments */
  --rk;
  --x;

  /* Function Body */
  xnorm = (d__1 = x[idamax_(n, &x[1], &c__1)], abs(d__1));
  resnorm = (d__1 = rk[idamax_(n, &rk[1], &c__1)], abs(d__1));

  ret_val = resnorm / (*tol * *n * *anorm * xnorm);

  return ret_val;

} /* scaledresid_ */

/*     =========================================================== */
doublereal matnorm_(n, a, lda)
integer *n;
doublereal *a;
integer *lda;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;
  doublereal ret_val, d__1;

  /* Local variables */
  static doublereal temp;
  static integer i, j;
  static doublereal rowsum;


  /*     Compute infinity norm of matrix A. */


  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  temp = 0.;
  i__1 = *n;
  for (i = 1; i <= i__1; ++i)
  {
    rowsum = 0.;
    i__2 = *n;
    for (j = 1; j <= i__2; ++j)
    {
      rowsum += (d__1 = a[i + j * a_dim1], abs(d__1));
      /* L10: */
    }
    temp = max(rowsum, temp);
    /* L20: */
  }

  ret_val = temp;

  return ret_val;

} /* matnorm_ */


/*     =============================================================== */
/* Subroutine */ int result_(n, a, lda, x, ldx, b, rk, mattype, pform, iter,
                             resid, tol, info, aform, anorm, ltest, scaledtol, testpassed, criterr,
                             mattype_len, pform_len, aform_len)
integer *n;
doublereal *a;
integer *lda;
doublereal *x;
integer *ldx;
doublereal *b, *rk;
char *mattype, *pform;
integer *iter;
doublereal *resid, *tol;
integer *info;
char *aform;
doublereal *anorm;
logical *ltest;
doublereal *scaledtol;
logical *testpassed;
integer *criterr;
ftnlen mattype_len;
ftnlen pform_len;
ftnlen aform_len;
{
  /* Format strings */
  static char fmt_900[] = "(\002Order\002,i4,\002 SPD 2-d Poisson matrix (\
no preconditioning)\002)";
  static char fmt_901[] = "(\002Order\002,i4,\002 SPD 2-d Poisson matrix (\
Jacobi preconditioning)\002)";
  static char fmt_902[] = "(\002Order\002,i4,\002 SPD 3-d Poisson matrix (\
no preconditioning)\002)";
  static char fmt_903[] = "(\002Order\002,i4,\002 SPD 3-d Poisson matrix (\
Jacobi preconditioning)\002)";
  static char fmt_910[] = "(\002Order \002,i4,\002 SPD Wathen matrix (no p\
reconditioning)\002)";
  static char fmt_911[] = "(\002Order \002,i4,\002 SPD Wathen matrix (Jaco\
bi preconditioning)\002)";
  static char fmt_920[] = "(\002Order \002,i4,\002 \002,a4,\002 nonsymmetr\
ic matrix (no preconditioning)\002)";
  static char fmt_921[] = "(\002Order \002,i4,\002 \002,a4,\002 nonsymmetr\
ic matrix (Jacobi preconditioning)\002)";
  static char fmt_991[] = "(\002  \002,a9,\002 \002,1pe8.2,\002    \002,1p\
e8.2,\002   \002,i5)";
  static char fmt_992[] = "(\002  \002,a9,\002 \002,1pe8.2,\002    \002,1p\
e8.2,\002   \002,i5,\002            X\002)";
  static char fmt_993[] = "(\002  \002,a9,\002 \002,1pe8.2,\002    \002,1p\
e8.2,\002   \002,i5,\002    \002,i3)";

  /* System generated locals */
  integer a_dim1, a_offset, x_dim1, x_offset, i__1;

  /* Builtin functions */
  /* Subroutine */
  int s_copy();
  integer s_wsfe(), do_fio(), e_wsfe(), s_wsle(), e_wsle(), do_lio();

  /* Local variables */
  static integer firstalg, i;
  extern logical lsame_();
  extern /* Subroutine */ int dgemv_(), dcopy_();
  extern logical lsamen_();
  static char method[9 * 9];
  static integer numalg;
  static doublereal sresid[9];
  extern doublereal scaledresid_();

  /* Fortran I/O blocks */
  static cilist io___75 = { 0, 10, 0, fmt_900, 0 };
  static cilist io___76 = { 0, 6, 0, fmt_900, 0 };
  static cilist io___77 = { 0, 10, 0, fmt_901, 0 };
  static cilist io___78 = { 0, 6, 0, fmt_901, 0 };
  static cilist io___79 = { 0, 10, 0, fmt_902, 0 };
  static cilist io___80 = { 0, 6, 0, fmt_902, 0 };
  static cilist io___81 = { 0, 10, 0, fmt_903, 0 };
  static cilist io___82 = { 0, 6, 0, fmt_903, 0 };
  static cilist io___83 = { 0, 10, 0, fmt_910, 0 };
  static cilist io___84 = { 0, 6, 0, fmt_910, 0 };
  static cilist io___85 = { 0, 10, 0, fmt_911, 0 };
  static cilist io___86 = { 0, 6, 0, fmt_911, 0 };
  static cilist io___87 = { 0, 10, 0, fmt_920, 0 };
  static cilist io___88 = { 0, 6, 0, fmt_920, 0 };
  static cilist io___89 = { 0, 10, 0, fmt_921, 0 };
  static cilist io___90 = { 0, 6, 0, fmt_921, 0 };
  static cilist io___91 = { 0, 10, 0, 0, 0 };
  static cilist io___92 = { 0, 10, 0, fmt_991, 0 };
  static cilist io___93 = { 0, 6, 0, fmt_991, 0 };
  static cilist io___94 = { 0, 10, 0, fmt_992, 0 };
  static cilist io___95 = { 0, 6, 0, fmt_992, 0 };
  static cilist io___96 = { 0, 10, 0, fmt_993, 0 };
  static cilist io___97 = { 0, 6, 0, fmt_993, 0 };
  static cilist io___98 = { 0, 10, 0, 0, 0 };



  /*     .. Argument Declaractions .. */


  /*  Purpose */
  /*  ======= */

  /*  Report results of METHOD on matrix type MATTYPE. If the residual */
  /*  is not directly computed by the algorithm, then the residual RESID */
  /*  as returned by the algorithm is compared with the residual as */
  /*  computed using the solution returned, i.e. || B-AX ||. */
  /*  ======================================================= */

  /*     .. Local Declarations .. */




  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --ltest;
  --info;
  --resid;
  --iter;
  --rk;
  --b;
  x_dim1 = *ldx;
  x_offset = x_dim1 + 1;
  x -= x_offset;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  s_copy(method, "CG       ", 9L, 9L);
  s_copy(method + 9, "Chebyshev", 9L, 9L);
  s_copy(method + 18, "SOR      ", 9L, 9L);
  s_copy(method + 27, "BiCG     ", 9L, 9L);
  s_copy(method + 36, "CGS      ", 9L, 9L);
  s_copy(method + 45, "BiCGSTAB ", 9L, 9L);
  s_copy(method + 54, "GMRESm   ", 9L, 9L);
  s_copy(method + 63, "QMR      ", 9L, 9L);
  s_copy(method + 72, "Jacobi   ", 9L, 9L);

  /*     Compare algorithm reported residual with |b-AX|/(|A||x|n*TOL) */

  if (lsame_(mattype, "SPD", 3L, 3L))
  {
    firstalg = 1;
    numalg = 9;
  }
  else
  {
    firstalg = 4;
    numalg = 8;
  }
  i__1 = numalg;
  for (i = firstalg; i <= i__1; ++i)
  {
    if (resid[i] != 0.)
    {
      dcopy_(n, &b[1], &c__1, &rk[1], &c__1);
      dgemv_("N", n, n, &c_b125, &a[a_offset], lda, &x[i * x_dim1 + 1],
             &c__1, &c_b127, &rk[1], &c__1, 1L);
      sresid[i - 1] = scaledresid_(anorm, n, &x[i * x_dim1 + 1], &rk[1],
                                   tol);
    }
    /* L10: */
  }

  if (lsamen_(&c__4, aform, "F2SH", 4L, 4L))
  {
    if (lsamen_(&c__2, pform, "IDENT", 5L, 5L))
    {
      s_wsfe(&io___75);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___76);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
    }
    else if (lsamen_(&c__2, pform, "JACBI", 5L, 5L))
    {
      s_wsfe(&io___77);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___78);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
    }
  }
  else if (lsamen_(&c__4, aform, "F3SH", 4L, 4L))
  {
    if (lsamen_(&c__2, pform, "IDENT", 5L, 5L))
    {
      s_wsfe(&io___79);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___80);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
    }
    else if (lsamen_(&c__2, pform, "JACBI", 5L, 5L))
    {
      s_wsfe(&io___81);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___82);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
    }
  }
  else if (lsamen_(&c__4, aform, "WATH", 4L, 4L))
  {
    if (lsamen_(&c__2, pform, "IDENT", 5L, 5L))
    {
      s_wsfe(&io___83);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___84);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
    }
    else if (lsamen_(&c__2, pform, "JACBI", 5L, 5L))
    {
      s_wsfe(&io___85);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___86);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      e_wsfe();
    }
  }
  else if (lsamen_(&c__3, aform, "PDE", 4L, 3L))
  {
    if (lsamen_(&c__2, pform, "IDENT", 5L, 5L))
    {
      s_wsfe(&io___87);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      do_fio(&c__1, aform, 4L);
      e_wsfe();
      s_wsfe(&io___88);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      do_fio(&c__1, aform, 4L);
      e_wsfe();
    }
    else if (lsamen_(&c__2, pform, "JACBI", 5L, 5L))
    {
      s_wsfe(&io___89);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      do_fio(&c__1, aform, 4L);
      e_wsfe();
      s_wsfe(&io___90);
      do_fio(&c__1, (char *) & (*n), (ftnlen)sizeof(integer));
      do_fio(&c__1, aform, 4L);
      e_wsfe();
    }
  }
  s_wsle(&io___91);
  e_wsle();

  /*     Loop over the algorithms, with a final error check. */

  i__1 = numalg;
  for (i = firstalg; i <= i__1; ++i)
  {

    /*        Check updated residual vs. scaled residual. */

    if (ltest[i])
    {
      if (info[i] == 0)
      {

        /*              Method claims to have found solution. */
        /*              Check scaled residual. */

        if (sresid[i - 1] <= *scaledtol)
        {

          /*                 Scaled residual check passed. */

          s_wsfe(&io___92);
          do_fio(&c__1, method + (i - 1) * 9, 9L);
          do_fio(&c__1, (char *)&resid[i], (ftnlen)sizeof(
                   doublereal));
          do_fio(&c__1, (char *)&sresid[i - 1], (ftnlen)sizeof(
                   doublereal));
          do_fio(&c__1, (char *)&iter[i], (ftnlen)sizeof(integer));
          e_wsfe();
          s_wsfe(&io___93);
          do_fio(&c__1, method + (i - 1) * 9, 9L);
          do_fio(&c__1, (char *)&resid[i], (ftnlen)sizeof(
                   doublereal));
          do_fio(&c__1, (char *)&sresid[i - 1], (ftnlen)sizeof(
                   doublereal));
          do_fio(&c__1, (char *)&iter[i], (ftnlen)sizeof(integer));
          e_wsfe();
        }
        else
        {
          ++(*criterr);
          *testpassed = FALSE_;
          s_wsfe(&io___94);
          do_fio(&c__1, method + (i - 1) * 9, 9L);
          do_fio(&c__1, (char *)&resid[i], (ftnlen)sizeof(
                   doublereal));
          do_fio(&c__1, (char *)&sresid[i - 1], (ftnlen)sizeof(
                   doublereal));
          do_fio(&c__1, (char *)&iter[i], (ftnlen)sizeof(integer));
          e_wsfe();
          s_wsfe(&io___95);
          do_fio(&c__1, method + (i - 1) * 9, 9L);
          do_fio(&c__1, (char *)&resid[i], (ftnlen)sizeof(
                   doublereal));
          do_fio(&c__1, (char *)&sresid[i - 1], (ftnlen)sizeof(
                   doublereal));
          do_fio(&c__1, (char *)&iter[i], (ftnlen)sizeof(integer));
          e_wsfe();
        }
      }
      else if (info[i] == 100)
      {
        goto L30;
      }
      else
      {
        *testpassed = FALSE_;

        /*              Method claims to have not found solution to to
        lerance, */
        /*              either because the maximum number of iteration
        s were */
        /*              performed, or breakdown occured. */

        s_wsfe(&io___96);
        do_fio(&c__1, method + (i - 1) * 9, 9L);
        do_fio(&c__1, (char *)&resid[i], (ftnlen)sizeof(doublereal));
        do_fio(&c__1, (char *)&sresid[i - 1], (ftnlen)sizeof(
                 doublereal));
        do_fio(&c__1, (char *)&iter[i], (ftnlen)sizeof(integer));
        do_fio(&c__1, (char *)&info[i], (ftnlen)sizeof(integer));
        e_wsfe();
        s_wsfe(&io___97);
        do_fio(&c__1, method + (i - 1) * 9, 9L);
        do_fio(&c__1, (char *)&resid[i], (ftnlen)sizeof(doublereal));
        do_fio(&c__1, (char *)&sresid[i - 1], (ftnlen)sizeof(
                 doublereal));
        do_fio(&c__1, (char *)&iter[i], (ftnlen)sizeof(integer));
        do_fio(&c__1, (char *)&info[i], (ftnlen)sizeof(integer));
        e_wsfe();
      }

    }
    else
    {

      /*           Method was not involved in test */

      goto L30;

    }

L30:
    ;
  }

  s_wsle(&io___98);
  do_lio(&c__9, &c__1, "--------------------------------------------------\
-----", 55L);
  e_wsle();

  /*     Header for each system. */


  /*     Reporting of results. */


  return 0;

  /*     End of Result.f */

} /* result_ */


/*     ================================================================== */
/* Subroutine */ int header_(tol)
doublereal *tol;
{
  /* Format strings */
  static char fmt_21[] = "(\002MACHINE PRECISION = \002,1pe8.2)";
  static char fmt_22[] = "(\002CONVERGENCE TEST TOLERANCE = \002,1pe8.2)";

  /* Builtin functions */
  integer s_wsle(), e_wsle(), do_lio(), s_wsfe(), do_fio(), e_wsfe();

  /* Local variables */
  extern doublereal dlamch_();
  static doublereal eps;

  /* Fortran I/O blocks */
  static cilist io___100 = { 0, 10, 0, 0, 0 };
  static cilist io___101 = { 0, 10, 0, 0, 0 };
  static cilist io___102 = { 0, 10, 0, 0, 0 };
  static cilist io___103 = { 0, 10, 0, 0, 0 };
  static cilist io___104 = { 0, 10, 0, 0, 0 };
  static cilist io___105 = { 0, 10, 0, 0, 0 };
  static cilist io___106 = { 0, 10, 0, 0, 0 };
  static cilist io___107 = { 0, 10, 0, 0, 0 };
  static cilist io___108 = { 0, 10, 0, 0, 0 };
  static cilist io___109 = { 0, 10, 0, 0, 0 };
  static cilist io___110 = { 0, 10, 0, 0, 0 };
  static cilist io___111 = { 0, 10, 0, 0, 0 };
  static cilist io___112 = { 0, 10, 0, 0, 0 };
  static cilist io___113 = { 0, 10, 0, fmt_21, 0 };
  static cilist io___114 = { 0, 10, 0, fmt_22, 0 };
  static cilist io___115 = { 0, 10, 0, 0, 0 };
  static cilist io___116 = { 0, 10, 0, 0, 0 };
  static cilist io___117 = { 0, 10, 0, 0, 0 };
  static cilist io___118 = { 0, 10, 0, 0, 0 };
  static cilist io___119 = { 0, 10, 0, 0, 0 };
  static cilist io___120 = { 0, 10, 0, 0, 0 };
  static cilist io___121 = { 0, 10, 0, 0, 0 };
  static cilist io___122 = { 0, 10, 0, 0, 0 };
  static cilist io___123 = { 0, 10, 0, 0, 0 };
  static cilist io___124 = { 0, 10, 0, 0, 0 };
  static cilist io___125 = { 0, 6, 0, 0, 0 };
  static cilist io___126 = { 0, 6, 0, 0, 0 };
  static cilist io___127 = { 0, 6, 0, 0, 0 };
  static cilist io___128 = { 0, 6, 0, 0, 0 };
  static cilist io___129 = { 0, 6, 0, 0, 0 };
  static cilist io___130 = { 0, 6, 0, 0, 0 };
  static cilist io___131 = { 0, 6, 0, 0, 0 };
  static cilist io___132 = { 0, 6, 0, 0, 0 };
  static cilist io___133 = { 0, 6, 0, 0, 0 };
  static cilist io___134 = { 0, 6, 0, 0, 0 };
  static cilist io___135 = { 0, 6, 0, 0, 0 };
  static cilist io___136 = { 0, 6, 0, 0, 0 };
  static cilist io___137 = { 0, 6, 0, 0, 0 };
  static cilist io___138 = { 0, 6, 0, fmt_21, 0 };
  static cilist io___139 = { 0, 6, 0, fmt_22, 0 };
  static cilist io___140 = { 0, 6, 0, 0, 0 };
  static cilist io___141 = { 0, 6, 0, 0, 0 };
  static cilist io___142 = { 0, 6, 0, 0, 0 };
  static cilist io___143 = { 0, 6, 0, 0, 0 };
  static cilist io___144 = { 0, 6, 0, 0, 0 };
  static cilist io___145 = { 0, 6, 0, 0, 0 };
  static cilist io___146 = { 0, 6, 0, 0, 0 };
  static cilist io___147 = { 0, 6, 0, 0, 0 };
  static cilist io___148 = { 0, 6, 0, 0, 0 };
  static cilist io___149 = { 0, 6, 0, 0, 0 };




  eps = dlamch_("E", 1L);

  /*     Print header to file. */

  s_wsle(&io___100);
  e_wsle();
  s_wsle(&io___101);
  do_lio(&c__9, &c__1, "DETAILS OF ITERATIVE TEMPLATES TEST:", 36L);
  e_wsle();
  s_wsle(&io___102);
  e_wsle();
  s_wsle(&io___103);
  do_lio(&c__9, &c__1, "   Univ. of Tennessee and Oak Ridge National Labor\
atory", 55L);
  e_wsle();
  s_wsle(&io___104);
  do_lio(&c__9, &c__1, "   October 1, 1993", 18L);
  e_wsle();
  s_wsle(&io___105);
  do_lio(&c__9, &c__1, "   Details of these algorithms are described in \"\
Templates", 58L);
  e_wsle();
  s_wsle(&io___106);
  do_lio(&c__9, &c__1, "   for the Solution of Linear Systems: Building Bl\
ocks for", 58L);
  e_wsle();
  s_wsle(&io___107);
  do_lio(&c__9, &c__1, "   Iterative Methods\", Barrett, Berry, Chan, Demm\
el, Donato,", 60L);
  e_wsle();
  s_wsle(&io___108);
  do_lio(&c__9, &c__1, "   Dongarra, Eijkhout, Pozo, Romine, and van der V\
orst,", 55L);
  e_wsle();
  s_wsle(&io___109);
  do_lio(&c__9, &c__1, "   SIAM Publications, 1993.", 27L);
  e_wsle();
  s_wsle(&io___110);
  do_lio(&c__9, &c__1, "   (ftp netlib2.cs.utk.edu; cd linalg; get templat\
es.ps).", 57L);
  e_wsle();
  s_wsle(&io___111);
  e_wsle();
  s_wsle(&io___112);
  e_wsle();
  s_wsfe(&io___113);
  do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
  e_wsfe();
  s_wsfe(&io___114);
  do_fio(&c__1, (char *) & (*tol), (ftnlen)sizeof(doublereal));
  e_wsfe();
  s_wsle(&io___115);
  e_wsle();
  s_wsle(&io___116);
  e_wsle();
  s_wsle(&io___117);
  do_lio(&c__9, &c__1, " For a detailed description of the following infor\
mation,", 57L);
  e_wsle();
  s_wsle(&io___118);
  do_lio(&c__9, &c__1, " see the end of this file.", 26L);
  e_wsle();
  s_wsle(&io___119);
  e_wsle();
  s_wsle(&io___120);
  do_lio(&c__9, &c__1, "==================================================\
====", 54L);
  e_wsle();
  s_wsle(&io___121);
  do_lio(&c__9, &c__1, "           CONVERGENCE  NORMALIZED  NUM", 39L);
  e_wsle();
  s_wsle(&io___122);
  do_lio(&c__9, &c__1, "  METHOD    CRITERION    RESIDUAL   ITER  INFO  FL\
AG", 52L);
  e_wsle();
  s_wsle(&io___123);
  do_lio(&c__9, &c__1, "==================================================\
====", 54L);
  e_wsle();
  s_wsle(&io___124);
  e_wsle();

  /*     Print header to screen. */

  s_wsle(&io___125);
  e_wsle();
  s_wsle(&io___126);
  do_lio(&c__9, &c__1, "DETAILS OF ITERATIVE TEMPLATES TEST:", 36L);
  e_wsle();
  s_wsle(&io___127);
  e_wsle();
  s_wsle(&io___128);
  do_lio(&c__9, &c__1, "   Univ. of Tennessee and Oak Ridge National Labor\
atory", 55L);
  e_wsle();
  s_wsle(&io___129);
  do_lio(&c__9, &c__1, "   October 1, 1993", 18L);
  e_wsle();
  s_wsle(&io___130);
  do_lio(&c__9, &c__1, "   Details of these algorithms are described in \"\
Templates", 58L);
  e_wsle();
  s_wsle(&io___131);
  do_lio(&c__9, &c__1, "   for the Solution of Linear Systems: Building Bl\
ocks for", 58L);
  e_wsle();
  s_wsle(&io___132);
  do_lio(&c__9, &c__1, "   Iterative Methods\", Barrett, Berry, Chan, Demm\
el, Donato,", 60L);
  e_wsle();
  s_wsle(&io___133);
  do_lio(&c__9, &c__1, "   Dongarra, Eijkhout, Pozo, Romine, and van der V\
orst,", 55L);
  e_wsle();
  s_wsle(&io___134);
  do_lio(&c__9, &c__1, "   SIAM Publications, 1993.", 27L);
  e_wsle();
  s_wsle(&io___135);
  do_lio(&c__9, &c__1, "   (ftp netlib2.cs.utk.edu; cd linalg; get templat\
es.ps).", 57L);
  e_wsle();
  s_wsle(&io___136);
  e_wsle();
  s_wsle(&io___137);
  e_wsle();
  s_wsfe(&io___138);
  do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
  e_wsfe();
  s_wsfe(&io___139);
  do_fio(&c__1, (char *) & (*tol), (ftnlen)sizeof(doublereal));
  e_wsfe();
  s_wsle(&io___140);
  e_wsle();
  s_wsle(&io___141);
  e_wsle();
  s_wsle(&io___142);
  do_lio(&c__9, &c__1, " For a detailed description of the following infor\
mation,", 57L);
  e_wsle();
  s_wsle(&io___143);
  do_lio(&c__9, &c__1, " see the end of this file.", 26L);
  e_wsle();
  s_wsle(&io___144);
  e_wsle();
  s_wsle(&io___145);
  do_lio(&c__9, &c__1, "==================================================\
===", 53L);
  e_wsle();
  s_wsle(&io___146);
  do_lio(&c__9, &c__1, "           CONVERGENCE  NORMALIZED  NUM", 39L);
  e_wsle();
  s_wsle(&io___147);
  do_lio(&c__9, &c__1, "  METHOD    CRITERION    RESIDUAL   ITER  INFO FLAG"
         , 51L);
  e_wsle();
  s_wsle(&io___148);
  do_lio(&c__9, &c__1, "==================================================\
===", 53L);
  e_wsle();
  s_wsle(&io___149);
  e_wsle();


  return 0;

} /* header_ */


/*     ================================================================== */
/* Subroutine */ int footer_()
{
  /* Builtin functions */
  integer s_wsle(), e_wsle(), do_lio();

  /* Fortran I/O blocks */
  static cilist io___150 = { 0, 10, 0, 0, 0 };
  static cilist io___151 = { 0, 10, 0, 0, 0 };
  static cilist io___152 = { 0, 10, 0, 0, 0 };
  static cilist io___153 = { 0, 10, 0, 0, 0 };
  static cilist io___154 = { 0, 10, 0, 0, 0 };
  static cilist io___155 = { 0, 10, 0, 0, 0 };
  static cilist io___156 = { 0, 10, 0, 0, 0 };
  static cilist io___157 = { 0, 10, 0, 0, 0 };
  static cilist io___158 = { 0, 10, 0, 0, 0 };
  static cilist io___159 = { 0, 10, 0, 0, 0 };
  static cilist io___160 = { 0, 10, 0, 0, 0 };
  static cilist io___161 = { 0, 10, 0, 0, 0 };
  static cilist io___162 = { 0, 10, 0, 0, 0 };
  static cilist io___163 = { 0, 10, 0, 0, 0 };
  static cilist io___164 = { 0, 10, 0, 0, 0 };
  static cilist io___165 = { 0, 10, 0, 0, 0 };
  static cilist io___166 = { 0, 10, 0, 0, 0 };
  static cilist io___167 = { 0, 10, 0, 0, 0 };
  static cilist io___168 = { 0, 10, 0, 0, 0 };
  static cilist io___169 = { 0, 10, 0, 0, 0 };
  static cilist io___170 = { 0, 10, 0, 0, 0 };
  static cilist io___171 = { 0, 10, 0, 0, 0 };
  static cilist io___172 = { 0, 10, 0, 0, 0 };
  static cilist io___173 = { 0, 10, 0, 0, 0 };
  static cilist io___174 = { 0, 10, 0, 0, 0 };
  static cilist io___175 = { 0, 10, 0, 0, 0 };
  static cilist io___176 = { 0, 10, 0, 0, 0 };
  static cilist io___177 = { 0, 10, 0, 0, 0 };
  static cilist io___178 = { 0, 10, 0, 0, 0 };
  static cilist io___179 = { 0, 10, 0, 0, 0 };
  static cilist io___180 = { 0, 10, 0, 0, 0 };
  static cilist io___181 = { 0, 10, 0, 0, 0 };
  static cilist io___182 = { 0, 10, 0, 0, 0 };
  static cilist io___183 = { 0, 10, 0, 0, 0 };
  static cilist io___184 = { 0, 10, 0, 0, 0 };
  static cilist io___185 = { 0, 10, 0, 0, 0 };
  static cilist io___186 = { 0, 10, 0, 0, 0 };
  static cilist io___187 = { 0, 10, 0, 0, 0 };
  static cilist io___188 = { 0, 10, 0, 0, 0 };
  static cilist io___189 = { 0, 10, 0, 0, 0 };
  static cilist io___190 = { 0, 10, 0, 0, 0 };
  static cilist io___191 = { 0, 10, 0, 0, 0 };
  static cilist io___192 = { 0, 10, 0, 0, 0 };
  static cilist io___193 = { 0, 10, 0, 0, 0 };
  static cilist io___194 = { 0, 10, 0, 0, 0 };
  static cilist io___195 = { 0, 10, 0, 0, 0 };
  static cilist io___196 = { 0, 10, 0, 0, 0 };
  static cilist io___197 = { 0, 10, 0, 0, 0 };
  static cilist io___198 = { 0, 10, 0, 0, 0 };
  static cilist io___199 = { 0, 10, 0, 0, 0 };
  static cilist io___200 = { 0, 10, 0, 0, 0 };
  static cilist io___201 = { 0, 10, 0, 0, 0 };
  static cilist io___202 = { 0, 10, 0, 0, 0 };
  static cilist io___203 = { 0, 10, 0, 0, 0 };
  static cilist io___204 = { 0, 10, 0, 0, 0 };
  static cilist io___205 = { 0, 10, 0, 0, 0 };
  static cilist io___206 = { 0, 10, 0, 0, 0 };
  static cilist io___207 = { 0, 10, 0, 0, 0 };
  static cilist io___208 = { 0, 10, 0, 0, 0 };
  static cilist io___209 = { 0, 10, 0, 0, 0 };
  static cilist io___210 = { 0, 10, 0, 0, 0 };
  static cilist io___211 = { 0, 10, 0, 0, 0 };
  static cilist io___212 = { 0, 10, 0, 0, 0 };
  static cilist io___213 = { 0, 10, 0, 0, 0 };
  static cilist io___214 = { 0, 10, 0, 0, 0 };
  static cilist io___215 = { 0, 10, 0, 0, 0 };
  static cilist io___216 = { 0, 10, 0, 0, 0 };
  static cilist io___217 = { 0, 10, 0, 0, 0 };
  static cilist io___218 = { 0, 10, 0, 0, 0 };
  static cilist io___219 = { 0, 10, 0, 0, 0 };
  static cilist io___220 = { 0, 10, 0, 0, 0 };
  static cilist io___221 = { 0, 10, 0, 0, 0 };
  static cilist io___222 = { 0, 10, 0, 0, 0 };
  static cilist io___223 = { 0, 10, 0, 0, 0 };
  static cilist io___224 = { 0, 10, 0, 0, 0 };
  static cilist io___225 = { 0, 10, 0, 0, 0 };
  static cilist io___226 = { 0, 10, 0, 0, 0 };
  static cilist io___227 = { 0, 10, 0, 0, 0 };
  static cilist io___228 = { 0, 10, 0, 0, 0 };
  static cilist io___229 = { 0, 10, 0, 0, 0 };
  static cilist io___230 = { 0, 10, 0, 0, 0 };
  static cilist io___231 = { 0, 10, 0, 0, 0 };
  static cilist io___232 = { 0, 10, 0, 0, 0 };
  static cilist io___233 = { 0, 10, 0, 0, 0 };
  static cilist io___234 = { 0, 10, 0, 0, 0 };
  static cilist io___235 = { 0, 10, 0, 0, 0 };
  static cilist io___236 = { 0, 10, 0, 0, 0 };
  static cilist io___237 = { 0, 10, 0, 0, 0 };
  static cilist io___238 = { 0, 10, 0, 0, 0 };
  static cilist io___239 = { 0, 10, 0, 0, 0 };
  static cilist io___240 = { 0, 10, 0, 0, 0 };
  static cilist io___241 = { 0, 10, 0, 0, 0 };
  static cilist io___242 = { 0, 10, 0, 0, 0 };
  static cilist io___243 = { 0, 10, 0, 0, 0 };
  static cilist io___244 = { 0, 10, 0, 0, 0 };
  static cilist io___245 = { 0, 10, 0, 0, 0 };



  /*     Puts descriptive information at bottom of results file */

  s_wsle(&io___150);
  e_wsle();
  s_wsle(&io___151);
  do_lio(&c__9, &c__1, "======", 6L);
  e_wsle();
  s_wsle(&io___152);
  do_lio(&c__9, &c__1, "LEGEND:", 7L);
  e_wsle();
  s_wsle(&io___153);
  do_lio(&c__9, &c__1, "======", 6L);
  e_wsle();
  s_wsle(&io___154);
  e_wsle();
  s_wsle(&io___155);
  do_lio(&c__9, &c__1, "   ==================", 21L);
  e_wsle();
  s_wsle(&io___156);
  do_lio(&c__9, &c__1, "   SYSTEM DESCRIPTION", 21L);
  e_wsle();
  s_wsle(&io___157);
  do_lio(&c__9, &c__1, "   ==================", 21L);
  e_wsle();
  s_wsle(&io___158);
  e_wsle();
  s_wsle(&io___159);
  do_lio(&c__9, &c__1, "   SPD matrices:", 16L);
  e_wsle();
  s_wsle(&io___160);
  e_wsle();
  s_wsle(&io___161);
  do_lio(&c__9, &c__1, "      WATH: \"Wathen Matrix\": consistent mass mat\
rix", 51L);
  e_wsle();
  s_wsle(&io___162);
  do_lio(&c__9, &c__1, "      F2SH: 2-d Poisson problem", 31L);
  e_wsle();
  s_wsle(&io___163);
  do_lio(&c__9, &c__1, "      F3SH: 3-d Poisson problem", 31L);
  e_wsle();
  s_wsle(&io___164);
  e_wsle();
  s_wsle(&io___165);
  do_lio(&c__9, &c__1, "   Nonsymmetric matrices:", 25L);
  e_wsle();
  s_wsle(&io___166);
  e_wsle();
  s_wsle(&io___167);
  do_lio(&c__9, &c__1, "      PDE1: u_xx+u_yy+au_x+(a_x/2)u", 35L);
  e_wsle();
  s_wsle(&io___168);
  do_lio(&c__9, &c__1, "            for a = 20exp[3.5(x**2+y**2 )]", 42L);
  e_wsle();
  s_wsle(&io___169);
  do_lio(&c__9, &c__1, "      PDE2: u_xx+u_yy+u_zz+1000u_x", 34L);
  e_wsle();
  s_wsle(&io___170);
  do_lio(&c__9, &c__1, "      PDE3  u_xx+u_yy+u_zz-10**5x**2(u_x+u_y+u_z )",
         50L);
  e_wsle();
  s_wsle(&io___171);
  do_lio(&c__9, &c__1, "      PDE4: u_xx+u_yy+u_zz+1000exp(xyz)(u_x+u_y-u_\
z)", 52L);
  e_wsle();
  s_wsle(&io___172);
  e_wsle();
  s_wsle(&io___173);
  do_lio(&c__9, &c__1, "   =====================", 24L);
  e_wsle();
  s_wsle(&io___174);
  do_lio(&c__9, &c__1, "   CONVERGENCE CRITERION", 24L);
  e_wsle();
  s_wsle(&io___175);
  do_lio(&c__9, &c__1, "   =====================", 24L);
  e_wsle();
  s_wsle(&io___176);
  e_wsle();
  s_wsle(&io___177);
  do_lio(&c__9, &c__1, "   Convergence criteria: residual as reported by t\
he", 52L);
  e_wsle();
  s_wsle(&io___178);
  do_lio(&c__9, &c__1, "   algorithm: ||AX - B|| / ||B||. Note that NaN ma\
y signify", 59L);
  e_wsle();
  s_wsle(&io___179);
  do_lio(&c__9, &c__1, "   divergence of the residual to the point of nume\
rical overflow.", 65L);
  e_wsle();
  s_wsle(&io___180);
  e_wsle();
  s_wsle(&io___181);
  do_lio(&c__9, &c__1, "   ===================", 22L);
  e_wsle();
  s_wsle(&io___182);
  do_lio(&c__9, &c__1, "   NORMALIZED RESIDUAL", 22L);
  e_wsle();
  s_wsle(&io___183);
  do_lio(&c__9, &c__1, "   ===================", 22L);
  e_wsle();
  s_wsle(&io___184);
  e_wsle();
  s_wsle(&io___185);
  do_lio(&c__9, &c__1, "   Normalized Residual: ||AX - B|| / (||A||||X||*N\
*TOL).", 56L);
  e_wsle();
  s_wsle(&io___186);
  do_lio(&c__9, &c__1, "   This is an apostiori check of the iterated solu\
tion. Note that if a", 70L);
  e_wsle();
  s_wsle(&io___187);
  do_lio(&c__9, &c__1, "   NaN occurs in the convergence criterion (discus\
sed above,", 60L);
  e_wsle();
  s_wsle(&io___188);
  do_lio(&c__9, &c__1, "   a NaN will be seen here also.", 32L);
  e_wsle();
  s_wsle(&io___189);
  e_wsle();
  s_wsle(&io___190);
  do_lio(&c__9, &c__1, "   ====", 7L);
  e_wsle();
  s_wsle(&io___191);
  do_lio(&c__9, &c__1, "   INFO", 7L);
  e_wsle();
  s_wsle(&io___192);
  do_lio(&c__9, &c__1, "   ====", 7L);
  e_wsle();
  s_wsle(&io___193);
  e_wsle();
  s_wsle(&io___194);
  do_lio(&c__9, &c__1, "   If this column is blank, then the algorithm cla\
ims to have", 61L);
  e_wsle();
  s_wsle(&io___195);
  do_lio(&c__9, &c__1, "   found the solution to tolerance (i.e. INFO = 0)."
         , 51L);
  e_wsle();
  s_wsle(&io___196);
  do_lio(&c__9, &c__1, "   This should be verified by checking the normali\
zed residual.", 63L);
  e_wsle();
  s_wsle(&io___197);
  e_wsle();
  s_wsle(&io___198);
  do_lio(&c__9, &c__1, "   Otherwise:", 13L);
  e_wsle();
  s_wsle(&io___199);
  e_wsle();
  s_wsle(&io___200);
  do_lio(&c__9, &c__1, "      = 1: Convergence not achieved given the maxi\
mum number of iterations.", 75L);
  e_wsle();
  s_wsle(&io___201);
  e_wsle();
  s_wsle(&io___202);
  do_lio(&c__9, &c__1, "      Input parameter errors:", 29L);
  e_wsle();
  s_wsle(&io___203);
  e_wsle();
  s_wsle(&io___204);
  do_lio(&c__9, &c__1, "      = -1: matrix dimension N < 0", 34L);
  e_wsle();
  s_wsle(&io___205);
  do_lio(&c__9, &c__1, "      = -2: LDW < N", 19L);
  e_wsle();
  s_wsle(&io___206);
  do_lio(&c__9, &c__1, "      = -3: Maximum number of iterations <= 0.",
         46L);
  e_wsle();
  s_wsle(&io___207);
  do_lio(&c__9, &c__1, "      = -4: For SOR: OMEGA not in interval (0,2)",
         48L);
  e_wsle();
  s_wsle(&io___208);
  do_lio(&c__9, &c__1, "            For GMRES: LDW2 < 2*RESTRT", 38L);
  e_wsle();
  s_wsle(&io___209);
  e_wsle();
  s_wsle(&io___210);
  do_lio(&c__9, &c__1, "      <= -10: Algorithm was terminated due to brea\
kdown.", 56L);
  e_wsle();
  s_wsle(&io___211);
  do_lio(&c__9, &c__1, "              See algorithm documentation for deta\
ils.", 54L);
  e_wsle();
  s_wsle(&io___212);
  e_wsle();
  s_wsle(&io___213);
  do_lio(&c__9, &c__1, "   ====", 7L);
  e_wsle();
  s_wsle(&io___214);
  do_lio(&c__9, &c__1, "   FLAG", 7L);
  e_wsle();
  s_wsle(&io___215);
  do_lio(&c__9, &c__1, "   ====", 7L);
  e_wsle();
  s_wsle(&io___216);
  e_wsle();
  s_wsle(&io___217);
  do_lio(&c__9, &c__1, "      X: Algorithm has reported convergence, but",
         48L);
  e_wsle();
  s_wsle(&io___218);
  do_lio(&c__9, &c__1, "         approximate solution fails scaled", 42L);
  e_wsle();
  s_wsle(&io___219);
  do_lio(&c__9, &c__1, "         residual check.", 24L);
  e_wsle();
  s_wsle(&io___220);
  e_wsle();
  s_wsle(&io___221);
  do_lio(&c__9, &c__1, "   =====", 8L);
  e_wsle();
  s_wsle(&io___222);
  do_lio(&c__9, &c__1, "   NOTES", 8L);
  e_wsle();
  s_wsle(&io___223);
  do_lio(&c__9, &c__1, "   =====", 8L);
  e_wsle();
  s_wsle(&io___224);
  e_wsle();
  s_wsle(&io___225);
  do_lio(&c__9, &c__1, "   GMRES: For the symmetric test matrices, the res\
tart parameter is", 67L);
  e_wsle();
  s_wsle(&io___226);
  do_lio(&c__9, &c__1, "   set to N. This should, theoretically, result in\
 no restarting. For", 69L);
  e_wsle();
  s_wsle(&io___227);
  do_lio(&c__9, &c__1, "   nonsymmetric testing the restart parameter is s\
et to N / 2.", 62L);
  e_wsle();
  s_wsle(&io___228);
  e_wsle();
  s_wsle(&io___229);
  do_lio(&c__9, &c__1, "   Stationary methods:", 22L);
  e_wsle();
  s_wsle(&io___230);
  e_wsle();
  s_wsle(&io___231);
  do_lio(&c__9, &c__1, "   - Since the residual norm ||b-Ax|| is not avail\
able as part of", 65L);
  e_wsle();
  s_wsle(&io___232);
  do_lio(&c__9, &c__1, "     the algorithm, the convergence criteria is di\
fferent from the", 66L);
  e_wsle();
  s_wsle(&io___233);
  do_lio(&c__9, &c__1, "     nonstationary methods. Here we use", 39L);
  e_wsle();
  s_wsle(&io___234);
  e_wsle();
  s_wsle(&io___235);
  do_lio(&c__9, &c__1, "        || X - X1 || / || X ||.", 31L);
  e_wsle();
  s_wsle(&io___236);
  e_wsle();
  s_wsle(&io___237);
  do_lio(&c__9, &c__1, "     That is, we compare the current approximated \
solution with the", 67L);
  e_wsle();
  s_wsle(&io___238);
  do_lio(&c__9, &c__1, "     approximation from the previous step.", 42L);
  e_wsle();
  s_wsle(&io___239);
  e_wsle();
  s_wsle(&io___240);
  do_lio(&c__9, &c__1, "   - Since Jacobi and SOR do not use preconditioni\
ng,", 53L);
  e_wsle();
  s_wsle(&io___241);
  do_lio(&c__9, &c__1, "     Jacobi is only iterated once per system, and \
SOR loops over", 64L);
  e_wsle();
  s_wsle(&io___242);
  do_lio(&c__9, &c__1, "     different values for OMEGA (the first time th\
rough OMEGA = 1,", 66L);
  e_wsle();
  s_wsle(&io___243);
  do_lio(&c__9, &c__1, "     i.e. the algorithm defaults to Gauss-Siedel).\
 This explains the ", 69L);
  e_wsle();
  s_wsle(&io___244);
  do_lio(&c__9, &c__1, "     different residual norms for SOR with the sam\
e matrix.", 59L);
  e_wsle();
  s_wsle(&io___245);
  e_wsle();

  return 0;

} /* footer_ */

/* Main program alias */ int templatestester_()
{
  MAIN__();
}
