/* dspdchk.f -- translated by f2c (version of 20 August 1993  13:15:44).
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

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__9 = 9;
static integer c__4 = 4;


/* Subroutine */
int dspdchk_(x, ldx, b, x0, work, ldw, pform, matvec,
             matvectrans, psolve, psolvetrans, psolveq, psolvetransq, backsolve,
             tol, scaledtol, ltest, maxws, spdres, numtests, numsusp, criterr,
             pform_len)
doublereal *x;
integer *ldx;
doublereal *b, *x0, *work;
integer *ldw;
char *pform;
/* Subroutine */
int (*matvec)(), (*matvectrans)(), (*psolve)(), (*
    psolvetrans)(), (*psolveq)(), (*psolvetransq)(), (*backsolve)();
doublereal *tol, *scaledtol;
logical *ltest;
integer *maxws;
logical *spdres;
integer *numtests, *numsusp, *criterr;
ftnlen pform_len;
{
  /* Format strings */
  static char fmt_81[] = "(\002WARNING: COULD NOT FORM 2-D POISSON (N=\002\
,i4,\002),INFO=\002,i2)";
  static char fmt_82[] = "(\002WARNING: COULD NOT FORM 3-D POISSON (N=\002\
,i4,\002),INFO=\002,i2)";
  static char fmt_83[] = "(\002WARNING: COULD NOT FORM ORDER\002,i4,\002 W\
ATHEN MATRIX\002)";
  static char fmt_92[] = "(\002ERROR: RHS FOR TEST MATRIX\002,a4,\002 NOT \
FORMED\002)";
  static char fmt_93[] = "(\002ERROR: INITIAL GUESS FOR TEST MATRIX\002,\
a4,\002 NOT FORMED\002)";
  static char fmt_94[] = "(\002ERROR: PRECONDITIONER\002,a5,\002 NOT FOR\
MED\002)";

  /* System generated locals */
  integer x_dim1, x_offset, i__1, i__2;

  /* Builtin functions */
  integer s_rsle(), do_lio(), e_rsle(), s_wsle(), e_wsle(), s_wsfe(),
          do_fio(), e_wsfe();
  /* Subroutine */
  int s_copy();

  /* Local variables */
  extern /* Subroutine */ int bicg_();
  static logical flag_;
  static integer ival, info[9];
  extern /* Subroutine */ int bicgstab_();
  static integer iter[9], nsys;
  static logical tstbicgs, tstcheby, tstgmres;
  static char initialguess[4];
  static integer i, j, k;
  extern /* Subroutine */ int cheby_();
  static char aform[4];
  static doublereal resid[9];
  static integer iindx;
  static doublereal anorm;
  extern /* Subroutine */ int gmres_(), dcopy_();
  static integer maxit;
  static logical tstcg, tstjacobi;
  extern doublereal negonefun_();
  extern /* Subroutine */ int cg_(), comp2dense_(), jacobi_();
  static integer nx, ny, nz;
  extern /* Subroutine */ int vecgen_();
  extern logical lsamen_();
  static integer needws;
  extern /* Subroutine */ int wathen_(), precon_(), gen57pt_();
  static logical tstcgs;
  extern /* Subroutine */ int result_();
  static integer restrt;
  static logical tstqmr, tstsor;
  extern /* Subroutine */ int cgs_();
  static char rhs[4];
  extern /* Subroutine */ int qmr_(), sor_();
  static logical tstbicg;
  static integer matform;
  extern doublereal matnorm_();
  static integer npforms, ipointr;
  extern doublereal zerofun_();

  /* Fortran I/O blocks */
  static cilist io___2 = { 0, 9, 0, 0, 0 };
  static cilist io___4 = { 0, 6, 0, 0, 0 };
  static cilist io___17 = { 0, 9, 0, 0, 0 };
  static cilist io___28 = { 0, 6, 0, fmt_81, 0 };
  static cilist io___29 = { 0, 10, 0, fmt_81, 0 };
  static cilist io___30 = { 0, 6, 0, fmt_82, 0 };
  static cilist io___31 = { 0, 10, 0, fmt_82, 0 };
  static cilist io___33 = { 0, 6, 0, fmt_83, 0 };
  static cilist io___34 = { 0, 10, 0, fmt_83, 0 };
  static cilist io___35 = { 0, 6, 0, 0, 0 };
  static cilist io___38 = { 0, 6, 0, fmt_92, 0 };
  static cilist io___39 = { 0, 10, 0, fmt_92, 0 };
  static cilist io___40 = { 0, 6, 0, fmt_93, 0 };
  static cilist io___41 = { 0, 10, 0, fmt_93, 0 };
  static cilist io___45 = { 0, 6, 0, fmt_94, 0 };
  static cilist io___46 = { 0, 10, 0, fmt_94, 0 };
  static cilist io___49 = { 0, 6, 0, 0, 0 };
  static cilist io___50 = { 0, 10, 0, 0, 0 };



  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Common Blocks .. */
  /*     MAXDIM2 = MAXDIM*MAXDIM. */




  /*  Purpose */
  /*  ======= */

  /*  Subroutine to test the performance of the template kernels */
  /*  on symmetric positivie definite matrices. */

  /*  Generates, solves, and checks accuracy of linear systems. */

  /*  Algorithms tested: */

  /*     1. CG */
  /*     2. Cheby */
  /*     3. SOR */
  /*     4. BiCG */
  /*     5. CGS */
  /*     6. BiCGSTAB */
  /*     7. GMRES */
  /*     8. QMR */
  /*     9. Jacobi ( for diagonally dominant Poisson matrices only ) */

  /*  Various systems are generated. Each method attempts to solve the */
  /*  system to the input TOL in MAXIT iterations. Each method iterates */
  /*  using various preconditioners. */

  /*  The result is compared with the normalized scaled residual */

  /*     || b - A*x || / ( ||A||||x||*N*TOL ). */

  /*  In order to do this, the solution vectors are stored in  matrix */
  /*  X( LDX,* ). Column j contains the solution vector for algorithm j, */
  /*  j as defined above. */

  /*  ================================================================= */

  /*     .. Parameters .. */

  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Local Arrays .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     PDE Coefficient functions. */
  /*     .. */
  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --ltest;
  pform -= 5;
  --work;
  --x0;
  --b;
  x_dim1 = *ldx;
  x_offset = x_dim1 + 1;
  x -= x_offset;

  /* Function Body */
  npforms = 2;
  *numtests = 0;
  *numsusp = 0;
  *criterr = 0;

  s_rsle(&io___2);
  do_lio(&c__3, &c__1, (char *)&nsys, (ftnlen)sizeof(integer));
  e_rsle();

  /*     Check for quick return. */

  if (nsys < 0)
  {
    s_wsle(&io___4);
    do_lio(&c__9, &c__1, "ERROR IN NONSYMMETRIC TESTER: NUMBER OF SYSTEM\
S TO BE GENERATED IS LESS THAN 0", 78L);
    e_wsle();
    return 0;
  }

  flag_ = FALSE_;
  for (i = 1; i <= 9; ++i)
  {
    if (ltest[i])
    {
      flag_ = TRUE_;
    }
    /* L5: */
  }
  if (! flag_)
  {
    return 0;
  }

  tstcg = ltest[1];
  tstcheby = ltest[2];
  tstsor = ltest[3];
  tstbicg = ltest[4];
  tstcgs = ltest[5];
  tstbicgs = ltest[6];
  tstgmres = ltest[7];
  tstqmr = ltest[8];
  tstjacobi = ltest[9];

  /* L10: */

  i__1 = nsys;
  for (matform = 1; matform <= i__1; ++matform)
  {

    s_rsle(&io___17);
    do_lio(&c__9, &c__1, aform, 4L);
    do_lio(&c__3, &c__1, (char *)&nx, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&ny, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nz, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, rhs, 4L);
    do_lio(&c__9, &c__1, initialguess, 4L);
    e_rsle();

    /*        The following two matrices are generated using a 5- or 7-poi
    nt */
    /*        stencil using centered differences on a 1d, 2d, or 3d grid,
    */
    /*        with Dirichlet boundary conditions. */

    /*        The last 7 arguments to this routine are the coefficient */
    /*        functions for the PDE: */

    /*           delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
     */
    /*           delx ( d u ) + dely (e u) + delz( f u ) + g u */

    if (lsamen_(&c__4, aform, "F2SH", 4L, 4L))
    {

      /*           -u_xx - u_yy */

      matdim_1.n = nx * ny * nz;
      ipointr = 1;
      iindx = ipointr + matdim_1.n + 1;
      ival = iindx + matdim_1.n * 5;
      gen57pt_(&nx, &ny, &nz, &work[ival], &work[iindx], &work[ipointr],
               negonefun_, negonefun_, zerofun_, zerofun_, zerofun_,
               zerofun_, zerofun_);
      comp2dense_(&work[ival], &work[ipointr], &work[iindx], &
                  matdim_1.n, system_1.a, &matdim_1.lda, "ROW", info, 3L);
      if (info[0] != 0)
      {
        s_wsfe(&io___28);
        do_fio(&c__1, (char *)&matdim_1.n, (ftnlen)sizeof(integer));
        do_fio(&c__1, (char *)&info[0], (ftnlen)sizeof(integer));
        e_wsfe();
        s_wsfe(&io___29);
        do_fio(&c__1, (char *)&matdim_1.n, (ftnlen)sizeof(integer));
        do_fio(&c__1, (char *)&info[0], (ftnlen)sizeof(integer));
        e_wsfe();
        goto L60;
      }

    }
    else if (lsamen_(&c__4, aform, "F3SH", 4L, 4L))
    {

      /*           -u_xx - u_yy - u_zz */

      matdim_1.n = nx * ny * nz;
      ipointr = 1;
      iindx = ipointr + matdim_1.n + 1;
      ival = iindx + matdim_1.n * 10;
      gen57pt_(&nx, &ny, &nz, &work[ival], &work[iindx], &work[ipointr],
               negonefun_, negonefun_, negonefun_, zerofun_, zerofun_,
               zerofun_, zerofun_);
      comp2dense_(&work[ival], &work[ipointr], &work[iindx], &
                  matdim_1.n, system_1.a, &matdim_1.lda, "ROW", info, 3L);
      if (info[0] != 0)
      {
        s_wsfe(&io___30);
        do_fio(&c__1, (char *)&matdim_1.n, (ftnlen)sizeof(integer));
        do_fio(&c__1, (char *)&info[0], (ftnlen)sizeof(integer));
        e_wsfe();
        s_wsfe(&io___31);
        do_fio(&c__1, (char *)&matdim_1.n, (ftnlen)sizeof(integer));
        do_fio(&c__1, (char *)&info[0], (ftnlen)sizeof(integer));
        e_wsfe();
        goto L60;
      }
    }
    else if (lsamen_(&c__4, aform, "WATHEN", 4L, 6L))
    {

      /*           Form Wathen matrix. */
      /*           ( Matrix order will be 3*NX*NY + 2*NX + 2*NY + 1 ) */

      k = 0;
      wathen_(&nx, &ny, &k, &matdim_1.n, system_1.a, &matdim_1.lda, &
              work[1], ldw, info);
      if (info[0] != 0)
      {
        s_wsfe(&io___33);
        i__2 = nx * 3 * ny + (nx << 1) + (ny << 1) + 1;
        do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
        do_fio(&c__1, (char *)&info[0], (ftnlen)sizeof(integer));
        e_wsfe();
        s_wsfe(&io___34);
        i__2 = nx * 3 * ny + (nx << 1) + (ny << 1) + 1;
        do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
        do_fio(&c__1, (char *)&info[0], (ftnlen)sizeof(integer));
        e_wsfe();
        goto L60;
      }
    }
    else
    {
      s_wsle(&io___35);
      do_lio(&c__9, &c__1, aform, 4L);
      do_lio(&c__9, &c__1, "IS AN UNKNOWM MATRIX TYPE", 25L);
      e_wsle();
      goto L60;
    }

    maxit = matdim_1.n << 2;
    anorm = matnorm_(&matdim_1.n, system_1.a, &matdim_1.lda);

    /*        Form RHS and initial guess vectors. */

    vecgen_(rhs, &matdim_1.n, system_1.a, &matdim_1.lda, &b[1], &info[1],
            4L);
    vecgen_(initialguess, &matdim_1.n, system_1.a, &matdim_1.lda, &x0[1],
            &info[2], 4L);

    if (info[1] != 0)
    {
      s_wsfe(&io___38);
      do_fio(&c__1, rhs, 4L);
      e_wsfe();
      s_wsfe(&io___39);
      do_fio(&c__1, rhs, 4L);
      e_wsfe();
      goto L60;
    }
    else if (info[2] != 0)
    {
      s_wsfe(&io___40);
      do_fio(&c__1, initialguess, 4L);
      e_wsfe();
      s_wsfe(&io___41);
      do_fio(&c__1, initialguess, 4L);
      e_wsfe();
      goto L60;
    }

    /* L20: */

    /*        Solve system using the various algorithms, using no */
    /*        preconditioning, then diagonal preconditioning. */

    i__2 = npforms;
    for (i = 1; i <= i__2; ++i)
    {

      for (j = 1; j <= 9; ++j)
      {
        info[j - 1] = 0;
        iter[j - 1] = maxit;
        resid[j - 1] = *tol;
        dcopy_(&matdim_1.n, &x0[1], &c__1, &x[j * x_dim1 + 1], &c__1);
        /* L30: */
      }

      precon_(&matdim_1.n, system_1.a, &matdim_1.lda, pform + i * 5,
              system_1.m, info, 5L);
      if (info[0] != 0)
      {
        s_wsfe(&io___45);
        do_fio(&c__1, pform + i * 5, 5L);
        e_wsfe();
        s_wsfe(&io___46);
        do_fio(&c__1, pform + i * 5, 5L);
        e_wsfe();
        goto L50;
      }
      s_copy(forms_1.curpform, pform + i * 5, 5L, 5L);

      if (tstcg)
      {
        cg_(&matdim_1.n, &b[1], &x[x_dim1 + 1], &work[1], ldw, iter,
            resid, matvec, psolve, info);
      }

      if (tstcheby)
      {
        cheby_(&matdim_1.n, &b[1], &x[(x_dim1 << 1) + 1], &work[1],
               ldw, &iter[1], &resid[1], matvec, psolve, &info[1]);
      }

      if (tstsor)
      {

        /*              Set OMEGA */

        if (matform == 1)
        {

          /*                 Guass-Seidel */

          work[1] = 1.;
        }
        else
        {
          work[1] = 1.2;
        }
        sor_(&matdim_1.n, &b[1], &x[x_dim1 * 3 + 1], &work[1], ldw, &
             iter[2], &resid[2], matvec, backsolve, &info[2]);
      }

      if (tstbicg)
      {
        bicg_(&matdim_1.n, &b[1], &x[(x_dim1 << 2) + 1], &work[1],
              ldw, &iter[3], &resid[3], matvec, matvectrans, psolve,
              psolvetrans, &info[3]);
      }

      if (tstcgs)
      {
        cgs_(&matdim_1.n, &b[1], &x[x_dim1 * 5 + 1], &work[1], ldw, &
             iter[4], &resid[4], matvec, psolve, &info[4]);
      }

      if (tstbicgs)
      {
        bicgstab_(&matdim_1.n, &b[1], &x[x_dim1 * 6 + 1], &work[1],
                  ldw, &iter[5], &resid[5], matvec, psolve, &info[5]);
      }

      if (tstgmres)
      {

        /*              For the symmetric case, restarts = N. */

        restrt = matdim_1.n;
        needws = matdim_1.n * (restrt + 4) + (restrt + 1) * (restrt +
                 2);
        if (needws > *maxws)
        {
          info[6] = 100;
          s_wsle(&io___49);
          do_lio(&c__9, &c__1, "WARNING: NOT ENOUGH WORKSPACE FOR \
GMRES, (TEST MATRIX", 53L);
          do_lio(&c__3, &c__1, (char *)&matform, (ftnlen)sizeof(
                   integer));
          do_lio(&c__9, &c__1, ")", 1L);
          e_wsle();
          s_wsle(&io___50);
          do_lio(&c__9, &c__1, "NOT ENOUGH WORKSPACE FOR GMRES, (T\
EST MATRIX", 44L);
          do_lio(&c__3, &c__1, (char *)&matform, (ftnlen)sizeof(
                   integer));
          do_lio(&c__9, &c__1, " REQUIRES MAXLEN=", 17L);
          do_lio(&c__3, &c__1, (char *)&needws, (ftnlen)sizeof(
                   integer));
          do_lio(&c__9, &c__1, ")", 1L);
          e_wsle();
        }
        else
        {
          gmres_(&matdim_1.n, &b[1], &x[x_dim1 * 7 + 1], &restrt, &
                 work[1], ldw, &work[matdim_1.n * (restrt + 5) + 1]
                 , ldw, &iter[6], &resid[6], matvec, psolve, &info[
                   6]);
        }
      }

      if (tstqmr)
      {
        qmr_(&matdim_1.n, &b[1], &x[(x_dim1 << 3) + 1], &work[1], ldw,
             &iter[7], &resid[7], matvec, matvectrans, psolveq,
             psolvetransq, &info[7]);
      }

      if (tstjacobi)
      {
        if (! lsamen_(&c__3, aform, "WATH", 4L, 4L) && i == 1)
        {

          /*                 Since preconditioning does not apply to
           Jacobi, it is */
          /*                 only called once per test matrix. (Wath
          en matrix is */
          /*                 not diagonally dominant.) */

          jacobi_(&matdim_1.n, &b[1], &x[x_dim1 * 9 + 1], &work[1],
                  ldw, &iter[8], &resid[8], matvec, &info[8]);
          if (info[8] != 0)
          {
            ++(*numsusp);
          }
        }
        else
        {

          /*                 Flag not to check accuracy. */

          info[8] = 100;
        }
      }

      /*           Check for convergence. */

      for (j = 1; j <= 8; ++j)
      {
        if (info[j - 1] != 0)
        {
          ++(*numsusp);
        }
        /* L40: */
      }

      /*           Write results to file. */

      result_(&matdim_1.n, system_1.a, &matdim_1.lda, &x[x_offset], ldx,
              &b[1], &work[1], "SPD", pform + i * 5, iter, resid, tol,
              info, aform, &anorm, &ltest[1], scaledtol, spdres,
              criterr, 3L, 5L, 4L);
      ++(*numtests);

L50:
      ;
    }

L60:
    ;
  }



  /* L200: */

  return 0;

  /*     -- End of DSPDCHK */

} /* dspdchk_ */

