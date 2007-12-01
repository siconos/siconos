/* MatGen.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b14 = 2.7182817459106;
static integer c__3 = 3;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__8 = 8;


/*     This file contains routines for forming matrices that result from */
/*     a 5- or 7-point discretization of elliptic PDEs with Dirichlet */
/*     boundary conditions, and a consistent mass matrix "Wathen". */

/*        GEN57PT and GETSTEN are from SPARSEKIT. They actually form the */
/*        row compressed matrix. */

/*        COEFF provides the functions for computing the coefficients */
/*        of the PDE. */

/*       Finally, for testing the iterative templates, COMP2DENSE converts*/
/*        the row compressed matrix to dense form. */

/*     ================================================================= */
/* Subroutine */
int gen57pt_(nx, ny, nz, a, indx, pointr, afun, bfun, cfun,
             dfun, efun, ffun, gfun)
integer *nx, *ny, *nz;
doublereal *a;
integer *indx, *pointr;
/* Subroutine */
int (*afun)(), (*bfun)(), (*cfun)(), (*dfun)(), (*efun)(
), (*ffun)(), (*gfun)();
{
  /* System generated locals */
  integer i__1, i__2, i__3;

  /* Local variables */
  static integer node;
  static doublereal h;
  static integer iedge, ix, iy, kx, ky, kz, iz;
  static doublereal stencil[7];
  extern /* Subroutine */ int getsten_();


  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */
  /*     .. External Functions .. */

  /*  Purpose */
  /*  ======= */

  /*  Adapted/altered from SPARSEKIT */

  /*  This subroutine computes the sparse matrix in row compressed */
  /*  format for the elliptic operator */

  /*  L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) + */
  /*        delx ( d u ) + dely (e u) + delz( f u ) + g u */

  /*  with Dirichlet Boundary conditions, on a rectangular 1-D, */
  /*  2-D or 3-D grid using centered difference schemes. */

  /*  The functions a, b, ..., g are known through the */
  /*  subroutines  afun, bfun, ..., gfun. */
  /*  Note that to obtain the correct matrix, any function that is not */
  /*  needed should be set to zero. For example for two-dimensional */
  /*  problems, nz should be set to 1 and the functions cfun and ffun */
  /*  should be zero functions. */

  /*  Uses natural ordering, first x direction, then y, then z */
  /*  mesh size h is uniform and determined by grid points */
  /*  in the x-direction. */

  /*  Arguments */
  /*  ========= */

  /*  NX     (input) INTEGER */
  /*         Number of points in X direction. */

  /*  NY     (input) INTEGER */
  /*         Number of points in Y direction. */

  /*  NZ     (input) INTEGER */
  /*         Number of points in Z direction. */

  /*  A,     (output) DOUBLE PRECISION array. */
  /*         Nonzero elements of the matrix. Stored in row compressed form.
  */

  /*  INDX   (output) INTEGER array. */
  /*         Column index of matrix element. */

  /*  POINTR (output) INTEGER array. */
  /*         Each element = P+1, where P is the number of nonzero elements
  */
  /*         in the preceding rows of the matrix. */

  /*  AFUN, */
  /*  BFUN, */
  /*  CFUN, */
  /*  DFUN, */
  /*  EFUN, */
  /*  FFUN, */
  /*  GFUN   (external subroutine) */
  /*        The user must supply the functions for computing the coefficient
  s*/
  /*         of the PDE. */

  /*  Description of the STENCIL: */

  /*     stencil [1:7] has the following meaning: */

  /*        center point = stencil(1) */
  /*        west point   = stencil(2) */
  /*        east point   = stencil(3) */
  /*        south point  = stencil(4) */
  /*        north point  = stencil(5) */
  /*        front point  = stencil(6) */
  /*        back point   = stencil(7) */

  /*                           st(5) */
  /*                            | */
  /*                            | */
  /*                            | */
  /*                            |          .st(7) */
  /*                            |     . */
  /*                            | . */
  /*         st(2) ----------- st(1) ---------- st(3) */
  /*                       .    | */
  /*                   .        | */
  /*               .            | */
  /*            st(6)           | */
  /*                            | */
  /*                            | */
  /*                           st(4) */

  /*     =============================================================== */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */

  /*     .. Executable Statements .. */

  /*     Initializations */

  /* Parameter adjustments */
  --pointr;
  --indx;
  --a;

  /* Function Body */
  h = 1. / (*nx + 1);
  kx = 1;
  ky = *nx;
  kz = *nx * *ny;
  iedge = 1;
  node = 1;

  i__1 = *nz;
  for (iz = 1; iz <= i__1; ++iz)
  {
    i__2 = *ny;
    for (iy = 1; iy <= i__2; ++iy)
    {
      i__3 = *nx;
      for (ix = 1; ix <= i__3; ++ix)
      {
        pointr[node] = iedge;

        /*              Get stencil. */

        getsten_(ny, nz, &ix, &iy, &iz, stencil, &h, afun, bfun, cfun,
                 dfun, efun, ffun, gfun);

        /*              West. */

        if (ix > 1)
        {
          indx[iedge] = node - kx;
          a[iedge] = stencil[1];
          ++iedge;
        }

        /*              South. */

        if (iy > 1)
        {
          indx[iedge] = node - ky;
          a[iedge] = stencil[3];
          ++iedge;
        }

        /*              Front Plane. */

        if (iz > 1)
        {
          indx[iedge] = node - kz;
          a[iedge] = stencil[5];
          ++iedge;
        }

        /*              Center node. */

        indx[iedge] = node;
        a[iedge] = stencil[0];
        ++iedge;

        /*              Upper part. */

        /*              East. */

        if (ix < *nx)
        {
          indx[iedge] = node + kx;
          a[iedge] = stencil[2];
          ++iedge;
        }
        if (iy < *ny)
        {
          indx[iedge] = node + ky;
          a[iedge] = stencil[4];
          ++iedge;
        }

        /*              Back plane. */

        if (iz < *nz)
        {
          indx[iedge] = node + kz;
          a[iedge] = stencil[6];
          ++iedge;
        }

        /*              Next node. */

        ++node;

        /* L10: */
      }
      /* L20: */
    }
    /* L30: */
  }

  pointr[node] = iedge;

  return 0;

  /*     -- End of GEN57PT */

} /* gen57pt_ */

/*     =============================================================== */
/* Subroutine */ int getsten_(ny, nz, kx, ky, kz, stencil, h, afun, bfun,
                              cfun, dfun, efun, ffun, gfun)
integer *ny, *nz, *kx, *ky, *kz;
doublereal *stencil, *h;
doublereal(*afun)(), (*bfun)(), (*cfun)(), (*dfun)(), (*efun)(), (*ffun)
(), (*gfun)();
{
  /* System generated locals */
  doublereal d__1;

  /* Local variables */
  static doublereal cntr;
  static integer k;
  static doublereal hhalf, coeff, x, y, z;


  /*     .. Argument Declarations .. */
  /*     .. */
  /*     .. External Functions .. */

  /*  Purpose */
  /*  ======= */

  /*  This subroutine calcultes the correct stencil values for */
  /*  elliptic operator */

  /*     L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) + */
  /*           delx ( d u ) + dely (e u) + delz( f u ) + g u. */

  /*  For 2-D problems the discretization formula that is used is: */

  /*  h**2 * Lu == a(i+1/2,j)*{u(i+1,j) - u(i,j)} + */
  /*               a(i-1/2,j)*{u(i-1,j) - u(i,j)} + */
  /*               b(i,j+1/2)*{u(i,j+1) - u(i,j)} + */
  /*               b(i,j-1/2)*{u(i,j-1) - u(i,j)} + */
  /*              (h/2)*d(i,j)*{u(i+1,j) - u(i-1,j)} + */
  /*              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} + */
  /*              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} + */
  /*              (h**2)*g(i,j)*u(i,j) */

  /*  =================================================================== */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */

  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --stencil;

  /* Function Body */
  for (k = 1; k <= 7; ++k)
  {
    stencil[k] = 0.;
    /* L10: */
  }

  hhalf = *h * .5;
  x = *h * *kx;
  y = *h * *ky;
  z = *h * *kz;
  cntr = 0.;

  /*     Differentiation w.r.t. X. */

  d__1 = x + hhalf;
  coeff = (*afun)(&d__1, &y, &z);
  stencil[3] += coeff;
  cntr += coeff;

  d__1 = x - hhalf;
  coeff = (*afun)(&d__1, &y, &z);
  stencil[2] += coeff;
  cntr += coeff;

  coeff = (*dfun)(&x, &y, &z) * hhalf;
  stencil[3] += coeff;
  stencil[2] -= coeff;
  if (*ny <= 1)
  {
    goto L99;
  }

  /*     Differentiation w.r.t. Y. */

  d__1 = y + hhalf;
  coeff = (*bfun)(&x, &d__1, &z);
  stencil[5] += coeff;
  cntr += coeff;

  d__1 = y - hhalf;
  coeff = (*bfun)(&x, &d__1, &z);
  stencil[4] += coeff;
  cntr += coeff;

  coeff = (*efun)(&x, &y, &z) * hhalf;
  stencil[5] += coeff;
  stencil[4] -= coeff;
  if (*nz <= 1)
  {
    goto L99;
  }

  /*     Differentiation w.r.t. Z. */

  d__1 = z + hhalf;
  coeff = (*cfun)(&x, &y, &d__1);
  stencil[7] += coeff;
  cntr += coeff;

  d__1 = z - hhalf;
  coeff = (*cfun)(&x, &y, &d__1);
  stencil[6] += coeff;
  cntr += coeff;

  coeff = (*ffun)(&x, &y, &z) * hhalf;
  stencil[7] += coeff;
  stencil[6] -= coeff;

  /*     Discretization of product by G. */

L99:
  coeff = (*gfun)(&x, &y, &z);
  stencil[1] = *h * *h * coeff - cntr;

  return 0;

} /* getsten_ */

/*     ============================================================= */
/*     Below are some functions for computing the value of the */
/*     coefficients. */
/*     ============================================================= */
doublereal zerofun_(x, y, z)
doublereal *x, *y, *z;
{
  /* System generated locals */
  doublereal ret_val;


  /*     .. Argument Declarations .. */

  /*     Purpose: Function to return ZERO. */
  /*     ======= */

  /*     .. Parameters .. */

  /*     .. Executable Statements .. */

  ret_val = 0.;

  /*     RETURN */

  return ret_val;
} /* zerofun_ */


/*     ============================================================= */
doublereal onefun_(x, y, z)
doublereal *x, *y, *z;
{
  /* System generated locals */
  doublereal ret_val;


  /*     .. Argument Declarations .. */

  /*     Purpose: Function to return ONE. */
  /*     ======= */

  /*     .. Parameters .. */

  /*     .. Executable Statements .. */

  ret_val = 1.;

  /*     RETURN */

  return ret_val;
} /* onefun_ */


/*     ============================================================= */
doublereal negonefun_(x, y, z)
doublereal *x, *y, *z;
{
  /* System generated locals */
  doublereal ret_val;


  /*     .. Argument Declarations .. */

  /*     Purpose: Function to return -ONE. */
  /*     ======= */

  /*     .. Parameters .. */

  /*     .. Executable Statements .. */

  ret_val = -1.;

  /*     RETURN */

  return ret_val;
} /* negonefun_ */


/*     ============================================================= */
doublereal thousfun_(x, y, z)
doublereal *x, *y, *z;
{
  /* System generated locals */
  doublereal ret_val;


  /*     .. Argument Declarations .. */

  /*     Purpose: Function to return 1000. */
  /*     ======= */

  /*     .. Executable Statements .. */

  ret_val = 1e3;

  /*     RETURN */

  return ret_val;
} /* thousfun_ */


/*     ============================================================= */
doublereal ten5x2fun_(x, y, z)
doublereal *x, *y, *z;
{
  /* System generated locals */
  doublereal ret_val;


  /*     .. Argument Declarations .. */

  /*     Purpose: Function to return 10 * X^2. */
  /*     ======= */

  /*     .. Executable Statements .. */

  ret_val = *x * 1e6 * *x;

  /*     RETURN */

  return ret_val;
} /* ten5x2fun_ */


/*     ============================================================= */
doublereal thousxfun_(x, y, z)
doublereal *x, *y, *z;
{
  /* System generated locals */
  doublereal ret_val, d__1;

  /* Builtin functions */
  double pow_dd();


  /*     .. Argument Declarations .. */

  /*     Purpose: Evaluates the coefficient function */
  /*     ======= */

  /*     .. Parameter .. */

  /*     .. Executable Statements .. */

  d__1 = *x * *y * *z;
  ret_val = pow_dd(&c_b14, &d__1) * 1e3;

  /*     RETURN */

  return ret_val;
} /* thousxfun_ */


/*     ============================================================= */
doublereal negthousxfun_(x, y, z)
doublereal *x, *y, *z;
{
  /* System generated locals */
  doublereal ret_val, d__1;

  /* Builtin functions */
  double pow_dd();


  /*     .. Argument Declarations .. */

  /*     Purpose: Evaluates the coefficient function */
  /*     ======= */

  /*     .. Parameter .. */

  /*     .. Executable Statements .. */

  d__1 = *x * *y * *z;
  ret_val = pow_dd(&c_b14, &d__1) * -1e3;

  /*     RETURN */

  return ret_val;
} /* negthousxfun_ */


/*     ============================================================= */
doublereal henkfun_(x, y, z)
doublereal *x, *y, *z;
{
  /* System generated locals */
  doublereal ret_val, d__1, d__2, d__3;

  /* Builtin functions */
  double pow_dd();


  /*     .. Argument Declarations .. */

  /*     Purpose: Evaluates the derivative of the above coefficient function
   */
  /*     ======= */

  /*     .. Parameter .. */

  /*     .. Executable Statements .. */

  /* Computing 2nd power */
  d__2 = *x;
  /* Computing 2nd power */
  d__3 = *y;
  d__1 = (d__2 * d__2 + d__3 * d__3) * 3.5;
  ret_val = pow_dd(&c_b14, &d__1) * 20;

  /*     RETURN */

  return ret_val;
} /* henkfun_ */


/*     ============================================================= */
doublereal henkdfun_(x, y, z)
doublereal *x, *y, *z;
{
  /* System generated locals */
  doublereal ret_val, d__1, d__2, d__3;

  /* Builtin functions */
  double pow_dd();


  /*     .. Argument Declarations .. */

  /*    Purpose: Evaluates the derivative of the above coefficient function.
  */
  /*     ======= */

  /*     .. Parameter .. */

  /*     .. Executable Statements .. */

  /* Computing 2nd power */
  d__2 = *x;
  /* Computing 2nd power */
  d__3 = *y;
  d__1 = (d__2 * d__2 + d__3 * d__3) * 3.5;
  ret_val = *x * 70 * pow_dd(&c_b14, &d__1);

  /*     RETURN */

  return ret_val;
} /* henkdfun_ */


/*     =============================================================== */
/* Subroutine */ int comp2dense_(asparse, pointr, indx, n, adense, lda, flag_,
                                 info, flag_len)
doublereal *asparse;
integer *pointr, *indx, *n;
doublereal *adense;
integer *lda;
char *flag_;
integer *info;
ftnlen flag_len;
{
  /* System generated locals */
  integer adense_dim1, adense_offset, i__1, i__2;

  /* Local variables */
  static integer i, j;
  extern logical lsame_(), lsamen_();



  /*     Convert sparse matrix storage to dense. */


  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  adense_dim1 = *lda;
  adense_offset = adense_dim1 + 1;
  adense -= adense_offset;
  --indx;
  --pointr;
  --asparse;

  /* Function Body */
  *info = 0;
  if (*n <= 0)
  {
    *info = -1;
  }
  else if (*n > *lda)
  {
    *info = -2;
  }
  else if (! lsamen_(&c__3, flag_, "ROW", 3L, 3L) || ! lsamen_(&c__3,
           flag_, "ROW", 3L, 3L))
  {
    *info = -3;
  }
  if (*info != 0)
  {
    return 0;
  }

  i__1 = *n;
  for (j = 1; j <= i__1; ++j)
  {
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
      adense[i + j * adense_dim1] = 0.;
      /* L10: */
    }
    /* L20: */
  }

  if (lsame_(flag_, "ROW", 3L, 3L))
  {
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
      i__2 = pointr[i + 1] - 1;
      for (j = pointr[i]; j <= i__2; ++j)
      {
        adense[i + indx[j] * adense_dim1] = asparse[j];
        /* L30: */
      }
      /* L40: */
    }
  }
  else if (lsame_(flag_, "COL", 3L, 3L))
  {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j)
    {
      i__2 = pointr[j + 1] - 1;
      for (i = pointr[j]; i <= i__2; ++i)
      {
        adense[indx[i] + j * adense_dim1] = asparse[i];
        /* L50: */
      }
      /* L60: */
    }
  }

  return 0;

} /* comp2dense_ */


/*     ================================================================ */
/* Subroutine */ int wathen_(nx, ny, kk, n, a, lda, work, ldw, info)
integer *nx, *ny, *kk, *n;
doublereal *a;
integer *lda;
doublereal *work;
integer *ldw, *info;
{
  /* System generated locals */
  integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2;

  /* Builtin functions */
  integer s_wsle(), do_lio(), e_wsle();

  /* Local variables */
  static integer kcol, krow, e, i, j, iseed[4];
  extern /* Subroutine */ int set_e__();
  static doublereal rhoit;
  static integer em, nn[8];
  extern doublereal dlaran_();
  static integer rho;

  /* Fortran I/O blocks */
  static cilist io___20 = { 0, 6, 0, 0, 0 };




  /*     Translated from the matlab version found on netlib. */

  /*     A is a random N-by-N finite element matrix where */
  /*     N = 3*NX*NY + 2*NX + 2*NY + 1. A is precisely the "consistent */
  /*     mass matrix" for a regular NX-by-NY grid of 8-node (serendipity) */
  /*     elements in 2 space dimensions. A is symmetric positive definite */
  /*     for any (positive) values of the "density", RHO(NX,NY), which is */
  /*     chosen randomly in this routine. In particular, if D=DIAG(DIAG(A)),
   */
  /*    then 0.25 <= EIG(INV(D)*A) <= 4.5 for any positive integers NX and N
  Y*/
  /*     and any densities RHO(NX,NY). This diagonally scaled matrix is */
  /*     returned by WATHEN(NX,NY,1). */

  /*     Reference: A.J.Wathen, DOUBLE PRECISIONistic eigenvalue bounds for
  */
  /*    the Galerkin mass matrix, IMA J. Numer. Anal., 7 (1987), pp. 449-457
  .*/

  /*     BEWARE - this is a sparse matrix, stored in -dense- form, and */
  /*              it quickly gets large! */

  /*     .. Local Scalars .. */


  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  work_dim1 = *ldw;
  work_offset = work_dim1 + 1;
  work -= work_offset;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  *info = 0;
  *n = *nx * 3 * *ny + (*nx << 1) + (*ny << 1) + 1;
  if (*n > *lda)
  {
    s_wsle(&io___20);
    do_lio(&c__9, &c__1, "NOT ENOUGH ROOM ALLOCATED FOR WATHEN MATRIX",
           43L);
    e_wsle();
    *info = -1;
    return 0;
  }
  else if (*nx < 1)
  {
    *info = -2;
  }
  else if (*ny < 1)
  {
    *info = -3;
  }
  else if (max(*nx, *ny) > *ldw)
  {
    *info = -4;
  }
  if (*info != 0)
  {
    return 0;
  }

  /*     Alias workspace columns. */

  e = 1;
  em = e + 8;
  rho = em + 8;

  set_e__(&work[e * work_dim1 + 1], ldw);

  i__1 = *n;
  for (j = 1; j <= i__1; ++j)
  {
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
      a[i + j * a_dim1] = 0.;
      /* L10: */
    }
    /* L20: */
  }

  iseed[0] = 304;
  iseed[1] = 152;
  iseed[2] = 2042;
  iseed[3] = 77;
  i__1 = *ny;
  for (j = 1; j <= i__1; ++j)
  {
    i__2 = *nx;
    for (i = 1; i <= i__2; ++i)
    {
      work[i + (rho + j - 1) * work_dim1] = dlaran_(iseed) * 100;
      /* L30: */
    }
    /* L40: */
  }

  i__1 = *ny;
  for (j = 1; j <= i__1; ++j)
  {
    i__2 = *nx;
    for (i = 1; i <= i__2; ++i)
    {

      nn[0] = j * 3 * *nx + (i << 1) + (j << 1) + 1;
      nn[1] = nn[0] - 1;
      nn[2] = nn[1] - 1;
      nn[3] = (j * 3 - 1) * *nx + (j << 1) + i - 1;
      nn[4] = (j - 1) * 3 * *nx + (i << 1) + (j << 1) - 3;
      nn[5] = nn[4] + 1;
      nn[6] = nn[5] + 1;
      nn[7] = nn[3] + 1;

      rhoit = work[i + (rho + j - 1) * work_dim1];
      for (krow = 1; krow <= 8; ++krow)
      {
        for (kcol = 1; kcol <= 8; ++kcol)
        {
          work[krow + (em + kcol - 1) * work_dim1] = rhoit * work[
                krow + (e + kcol - 1) * work_dim1];
          /* L50: */
        }
        /* L60: */
      }

      for (krow = 1; krow <= 8; ++krow)
      {
        for (kcol = 1; kcol <= 8; ++kcol)
        {
          a[nn[krow - 1] + nn[kcol - 1] * a_dim1] += work[krow + (
                em + kcol - 1) * work_dim1];
          /* L70: */
        }
        /* L80: */
      }

      /* L90: */
    }
    /* L100: */
  }

  if (*kk == 1)
  {

    /*        A = diag(diag(A)) \ A (the result being unit diagonal); */

    i__1 = j;
    for (i = 1; i <= i__1; ++i)
    {
      a[i + i * a_dim1] = 1.;
      /* L110: */
    }
  }

  return 0;

} /* wathen_ */

/*     ======================================================== */
/* Subroutine */ int set_e__(e, lde)
doublereal *e;
integer *lde;
{
  /* System generated locals */
  integer e_dim1, e_offset;

  /* Local variables */
  static integer i, j;
  extern /* Subroutine */ int dscal_();
  static doublereal scale;



  /* Parameter adjustments */
  e_dim1 = *lde;
  e_offset = e_dim1 + 1;
  e -= e_offset;

  /* Function Body */
  e[e_dim1 + 1] = 6.;
  e[e_dim1 + 2] = -6.;
  e[e_dim1 + 3] = 2.;
  e[e_dim1 + 4] = -8.;
  e[(e_dim1 << 1) + 1] = -6.;
  e[(e_dim1 << 1) + 2] = 32.;
  e[(e_dim1 << 1) + 3] = -6.;
  e[(e_dim1 << 1) + 4] = 20.;
  e[e_dim1 * 3 + 1] = 2.;
  e[e_dim1 * 3 + 2] = -6.;
  e[e_dim1 * 3 + 3] = 6.;
  e[e_dim1 * 3 + 4] = -6.;
  e[(e_dim1 << 2) + 1] = -8.;
  e[(e_dim1 << 2) + 2] = 20.;
  e[(e_dim1 << 2) + 3] = -6.;
  e[(e_dim1 << 2) + 4] = 32.;

  e[e_dim1 * 5 + 1] = 3.;
  e[e_dim1 * 5 + 2] = -8.;
  e[e_dim1 * 5 + 3] = 2.;
  e[e_dim1 * 5 + 4] = -6.;
  e[e_dim1 * 6 + 1] = -8.;
  e[e_dim1 * 6 + 2] = 16.;
  e[e_dim1 * 6 + 3] = -8.;
  e[e_dim1 * 6 + 4] = 20.;
  e[e_dim1 * 7 + 1] = 2.;
  e[e_dim1 * 7 + 2] = -8.;
  e[e_dim1 * 7 + 3] = 3.;
  e[e_dim1 * 7 + 4] = -8.;
  e[(e_dim1 << 3) + 1] = -6.;
  e[(e_dim1 << 3) + 2] = 20.;
  e[(e_dim1 << 3) + 3] = -8.;
  e[(e_dim1 << 3) + 4] = 16.;

  for (j = 1; j <= 4; ++j)
  {
    for (i = 5; i <= 8; ++i)
    {
      e[i + j * e_dim1] = e[j + i * e_dim1];
      /* L10: */
    }
    /* L20: */
  }

  for (j = 5; j <= 8; ++j)
  {
    for (i = 5; i <= 8; ++i)
    {
      e[i + j * e_dim1] = e[i - 4 + (j - 4) * e_dim1];
      /* L30: */
    }
    /* L40: */
  }

  scale = .022222222222222223;
  for (i = 1; i <= 8; ++i)
  {
    dscal_(&c__8, &scale, &e[i * e_dim1 + 1], &c__1);
    /* L50: */
  }

  return 0;

} /* set_e__ */

