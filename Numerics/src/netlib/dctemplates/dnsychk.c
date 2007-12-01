/* dnsychk.f -- translated by f2c (version of 20 August 1993  13:15:44).
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
int dnsychk_(x, ldx, b, x0, work, ldw, pform, matvec,
             matvectrans, psolve, psolvetrans, psolveq, psolvetransq, backsolve,
             tol, scaledtol, ltest, maxws, nsyres, numtests, numsusp, criterr,
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
logical *nsyres;
integer *numtests, *numsusp, *criterr;
ftnlen pform_len;
{
  /* Format strings */
  static char fmt_81[] = "(\002WARNING: COULD NOT FORM ORDER \002,i4,\002\
 \002,a4,\002MATRIX; INFO=\002,i2)";
  static char fmt_92[] = "(\002ERROR: RHS\002,a4,\002 NOT FORMED\002)";
  static char fmt_93[] = "(\002ERROR: INITIAL GUESS\002,a4,\002 NOT FORME\
D\002)";
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
  static integer iter[9];
  extern doublereal henkdfun_();
  static integer nsys;
  static logical tstbicgs;
  extern doublereal ten5x2fun_();
  static logical tstgmres;
  extern doublereal thousfun_();
  static char initialguess[4];
  static integer i, j;
  extern doublereal negthousxfun_();
  static char aform[4];
  static doublereal resid[9];
  static integer iindx;
  static doublereal anorm;
  extern /* Subroutine */ int gmres_(), dcopy_();
  static integer maxit;
  extern /* Subroutine */ int comp2dense_();
  extern doublereal thousxfun_();
  static integer nx, ny, nz;
  extern /* Subroutine */ int vecgen_();
  extern logical lsamen_();
  static integer needws;
  extern /* Subroutine */ int precon_(), gen57pt_();
  extern doublereal onefun_();
  static logical tstcgs;
  extern /* Subroutine */ int result_();
  static integer restrt;
  static logical tstqmr;
  extern /* Subroutine */ int cgs_();
  static char rhs[4];
  extern /* Subroutine */ int qmr_();
  extern doublereal henkfun_();
  static logical tstbicg;
  static integer matform;
  extern doublereal matnorm_();
  static integer npforms, ipointr;
  extern doublereal zerofun_();

  /* Fortran I/O blocks */
  static cilist io___2 = { 0, 9, 0, 0, 0 };
  static cilist io___4 = { 0, 6, 0, 0, 0 };
  static cilist io___13 = { 0, 9, 0, 0, 0 };
  static cilist io___23 = { 0, 6, 0, 0, 0 };
  static cilist io___25 = { 0, 6, 0, fmt_81, 0 };
  static cilist io___26 = { 0, 10, 0, fmt_81, 0 };
  static cilist io___27 = { 0, 6, 0, fmt_92, 0 };
  static cilist io___28 = { 0, 10, 0, fmt_92, 0 };
  static cilist io___29 = { 0, 6, 0, fmt_93, 0 };
  static cilist io___30 = { 0, 10, 0, fmt_93, 0 };
  static cilist io___36 = { 0, 6, 0, fmt_94, 0 };
  static cilist io___37 = { 0, 10, 0, fmt_94, 0 };
  static cilist io___40 = { 0, 6, 0, 0, 0 };
  static cilist io___41 = { 0, 10, 0, 0, 0 };



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

  /*  Subroutine to test the performance of the nonstationary template */
  /*  kernels on nonsymmetric matrices. */

  /*  Generates, solves, and check accuracy of linear systems. */

  /*  Algorithms tested: */

  /*     4. BiCG */
  /*     5. CGS */
  /*     6. BiCGSTAB */
  /*     7. GMRESm */
  /*     8. QMR */

  /*  Various systems are generated. Each method attempts to solve the */
  /*  system to the input TOL in MAXIT iterations. Each method iterates */
  /*  using various preconditioners. */

  /*  The result is compared with the normalized scaled residual */

  /*     || b - A*x || / ( ||A||||x||*N*TOL ). */

  /*  In order to do this, the solution vectors are stored in  matrix */
  /*  X( LDX,* ). Column j contains the solution vector for algorithm j, */
  /*  j as defined above. */

  /*  ================================================================= */

  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Local Arrays .. */

  /*     PDE Coefficient functions. */


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
  for (i = 4; i <= 8; ++i)
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

  tstbicg = ltest[4];
  tstcgs = ltest[5];
  tstbicgs = ltest[6];
  tstgmres = ltest[7];
  tstqmr = ltest[8];

L10:

  i__1 = nsys;
  for (matform = 1; matform <= i__1; ++matform)
  {

    s_rsle(&io___13);
    do_lio(&c__9, &c__1, aform, 4L);
    do_lio(&c__3, &c__1, (char *)&nx, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&ny, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nz, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, rhs, 4L);
    do_lio(&c__9, &c__1, initialguess, 4L);
    e_rsle();

    matdim_1.n = nx * ny * nz;
    ipointr = 1;
    iindx = ipointr + matdim_1.n + 1;
    ival = iindx + matdim_1.n * 5;

    /*        The following matrices are generated using a 5- or 7-point
    */
    /*        stencil using centered differences on a 1d, 2d, or 3d grid,
    */
    /*        with Dirichlet boundary conditions. */

    /*        The last 7 arguments to this routine are the coefficient */
    /*        functions for the PDE: */

    /*           delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
     */
    /*           delx ( d u ) + dely (e u) + delz( f u ) + g u */

    if (lsamen_(&c__4, aform, "PDE1", 4L, 4L))
    {

      /*           u_xx + u_yy + au_x + (a_x/2)u for a = 20exp[3.5(x^2 +
       y^2 )] */

      gen57pt_(&nx, &ny, &nz, &work[ival], &work[iindx], &work[ipointr],
               onefun_, onefun_, zerofun_, henkfun_, zerofun_, zerofun_,
               henkdfun_);

      /*        The following three PDE are from Yang, "PCG-like methods
       for */
      /*        nonsymmetric linear systems". */

    }
    else if (lsamen_(&c__4, aform, "PDE2", 4L, 4L))
    {

      /*           u_xx + u_yy + u_zz + 1000u_x */

      gen57pt_(&nx, &ny, &nz, &work[ival], &work[iindx], &work[ipointr],
               onefun_, onefun_, onefun_, thousfun_, zerofun_, zerofun_,
               zerofun_);
    }
    else if (lsamen_(&c__4, aform, "PDE3", 4L, 4L))
    {

      /*           u_xx + u_yy + u_zz - 10^5x^2(u_x + u_y + u_z ) */

      gen57pt_(&nx, &ny, &nz, &work[ival], &work[iindx], &work[ipointr],
               onefun_, onefun_, onefun_, ten5x2fun_, ten5x2fun_,
               ten5x2fun_, zerofun_);
    }
    else if (lsamen_(&c__4, aform, "PDE4", 4L, 4L))
    {

      /*           u_xx + u_yy + u_zz + 1000exp(xyz)( u_x + u_y - u_z )
      */

      gen57pt_(&nx, &ny, &nz, &work[ival], &work[iindx], &work[ipointr],
               onefun_, onefun_, onefun_, thousxfun_, thousxfun_,
               negthousxfun_, zerofun_);
    }
    else
    {
      s_wsle(&io___23);
      do_lio(&c__9, &c__1, aform, 4L);
      do_lio(&c__9, &c__1, "IS AN UNKNOWM MATRIX TYPE", 25L);
      e_wsle();
      goto L60;
    }

    /*        Convert to dense form. */

    comp2dense_(&work[ival], &work[ipointr], &work[iindx], &matdim_1.n,
                system_1.a, &matdim_1.lda, "ROW", info, 3L);
    if (info[0] != 0)
    {
      s_wsfe(&io___25);
      do_fio(&c__1, (char *)&matdim_1.n, (ftnlen)sizeof(integer));
      do_fio(&c__1, aform, 4L);
      do_fio(&c__1, (char *)&info[0], (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___26);
      do_fio(&c__1, (char *)&matdim_1.n, (ftnlen)sizeof(integer));
      do_fio(&c__1, aform, 4L);
      do_fio(&c__1, (char *)&info[0], (ftnlen)sizeof(integer));
      e_wsfe();
      goto L60;
    }

    /* L15: */

    vecgen_(rhs, &matdim_1.n, system_1.a, &matdim_1.lda, &b[1], &info[1],
            4L);
    vecgen_(initialguess, &matdim_1.n, system_1.a, &matdim_1.lda, &x0[1],
            &info[2], 4L);

    if (info[1] != 0)
    {
      s_wsfe(&io___27);
      do_fio(&c__1, (char *)&matform, (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___28);
      do_fio(&c__1, (char *)&matform, (ftnlen)sizeof(integer));
      e_wsfe();
      goto L10;
    }
    else if (info[2] != 0)
    {
      s_wsfe(&io___29);
      do_fio(&c__1, (char *)&matform, (ftnlen)sizeof(integer));
      e_wsfe();
      s_wsfe(&io___30);
      do_fio(&c__1, (char *)&matform, (ftnlen)sizeof(integer));
      e_wsfe();
      goto L10;
    }

    /* L20: */

    maxit = matdim_1.n << 2;
    anorm = matnorm_(&matdim_1.n, system_1.a, &matdim_1.lda);

    /*        Solve system using the various algorithms, using no */
    /*        preconditioning, then diagonal preconditioning. */

    i__2 = npforms;
    for (i = 1; i <= i__2; ++i)
    {

      /*           Initializations. */

      for (j = 3; j <= 9; ++j)
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
        s_wsfe(&io___36);
        do_fio(&c__1, pform + i * 5, 5L);
        e_wsfe();
        s_wsfe(&io___37);
        do_fio(&c__1, pform + i * 5, 5L);
        e_wsfe();
        goto L50;
      }
      s_copy(forms_1.curpform, pform + i * 5, 5L, 5L);

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
        restrt = matdim_1.n / 2;
        if (restrt == 0)
        {
          restrt = matdim_1.n;
        }
        needws = matdim_1.n * (restrt + 4) + (restrt + 1) * (restrt +
                 2);
        if (needws > *maxws)
        {
          info[6] = 100;
          s_wsle(&io___40);
          do_lio(&c__9, &c__1, "WARNING: NOT ENOUGH WORKSPACE FOR \
GMRES, (TEST MATRIX", 53L);
          do_lio(&c__3, &c__1, (char *)&matform, (ftnlen)sizeof(
                   integer));
          do_lio(&c__9, &c__1, ")", 1L);
          e_wsle();
          s_wsle(&io___41);
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

      /*           Check for convergence. */

      for (j = 4; j <= 8; ++j)
      {
        if (info[j - 1] != 0)
        {
          ++(*numsusp);
        }
        /* L40: */
      }

      /*           Write results to file. */

      result_(&matdim_1.n, system_1.a, &matdim_1.lda, &x[x_offset], ldx,
              &b[1], &work[1], "NSY", pform + i * 5, iter, resid, tol,
              info, aform, &anorm, &ltest[1], scaledtol, nsyres,
              criterr, 3L, 5L, 4L);
      ++(*numtests);

L50:
      ;
    }

L60:
    ;
  }


  return 0;

  /*     -- End of DNSYCHK */

} /* dnsychk_ */

