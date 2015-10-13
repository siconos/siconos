/* dlapack.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
  -lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b21 = 0.;
static doublereal c_b22 = 1.;
static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__3 = 3;
static doublereal c_b211 = -1.;

/* Subroutine */
int dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info,
           jobz_len, uplo_len)
char *jobz, *uplo;
integer *n;
doublereal *a;
integer *lda;
doublereal *w, *work;
integer *lwork, *info;
ftnlen jobz_len;
ftnlen uplo_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;
  doublereal d__1;

  /* Builtin functions */
  double sqrt();

  /* Local variables */
  static integer inde;
  static doublereal anrm;
  static integer imax;
  static doublereal rmin, rmax;
  static integer lopt, j;
  extern /* Subroutine */ int dscal_();
  static doublereal sigma;
  extern logical lsame_();
  static integer iinfo;
  static logical lower, wantz;
  extern doublereal dlamch_();
  static integer iscale;
  static doublereal safmin;
  extern /* Subroutine */ int xerbla_();
  static doublereal bignum;
  static integer indtau;
  extern /* Subroutine */ int dsterf_();
  extern doublereal dlansy_();
  static integer indwrk;
  extern /* Subroutine */ int dorgtr_(), dsteqr_(), dsytrd_();
  static integer llwork;
  static doublereal smlnum, eps;


  /*  -- LAPACK driver routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     March 31, 1993 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DSYEV computes all eigenvalues and, optionally, eigenvectors of a */
  /*  real symmetric matrix A. */

  /*  Arguments */
  /*  ========= */

  /*  JOBZ    (input) CHARACTER*1 */
  /*          = 'N':  Compute eigenvalues only; */
  /*          = 'V':  Compute eigenvalues and eigenvectors. */

  /*  UPLO    (input) CHARACTER*1 */
  /*          = 'U':  Upper triangle of A is stored; */
  /*          = 'L':  Lower triangle of A is stored. */

  /*  N       (input) INTEGER */
  /*          The order of the matrix A.  N >= 0. */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N) */
  /*          On entry, the symmetric matrix A.  If UPLO = 'U', the */
  /*          leading N-by-N upper triangular part of A contains the */
  /*          upper triangular part of the matrix A.  If UPLO = 'L', */
  /*          the leading N-by-N lower triangular part of A contains */
  /*          the lower triangular part of the matrix A. */
  /*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the */
  /*          orthonormal eigenvectors of the matrix A. */
  /*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
  */
  /*          or the upper triangle (if UPLO='U') of A, including the */
  /*          diagonal, is destroyed. */

  /*  LDA     (input) INTEGER */
  /*          The leading dimension of the array A.  LDA >= max(1,N). */

  /*  W       (output) DOUBLE PRECISION array, dimension (N) */
  /*          If INFO = 0, the eigenvalues in ascending order. */

  /*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK) */
  /*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

  /*  LWORK   (input) INTEGER */
  /*          The length of the array WORK.  LWORK >= max(1,3*N-1). */
  /*          For optimal efficiency, LWORK >= (NB+2)*N, */
  /*          where NB is the blocksize for DSYTRD returned by ILAENV. */

  /*  INFO    (output) INTEGER */
  /*          = 0:  successful exit */
  /*          < 0:  if INFO = -i, the i-th argument had an illegal value */
  /*          > 0:  if INFO = i, the algorithm failed to converge; i */
  /*                off-diagonal elements of an intermediate tridiagonal */
  /*                form did not converge to zero. */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  --work;
  --w;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  wantz = lsame_(jobz, "V", 1L, 1L);
  lower = lsame_(uplo, "L", 1L, 1L);

  *info = 0;
  if (!(wantz || lsame_(jobz, "N", 1L, 1L)))
  {
    *info = -1;
  }
  else if (!(lower || lsame_(uplo, "U", 1L, 1L)))
  {
    *info = -2;
  }
  else if (*n < 0)
  {
    *info = -3;
  }
  else if (*lda < max(1, *n))
  {
    *info = -5;
  }
  else /* if(complicated condition) */
  {
    /* Computing MAX */
    i__1 = 1, i__2 = *n * 3 - 1;
    if (*lwork < max(i__1, i__2))
    {
      *info = -8;
    }
  }

  if (*info != 0)
  {
    i__1 = -(*info);
    xerbla_("DSYEV ", &i__1, 6L);
    return 0;
  }

  /*     Quick return if possible */

  if (*n == 0)
  {
    work[1] = 1.;
    return 0;
  }

  if (*n == 1)
  {
    w[1] = a[a_dim1 + 1];
    work[1] = 3.;
    if (wantz)
    {
      a[a_dim1 + 1] = 1.;
    }
    return 0;
  }

  /*     Get machine constants. */

  safmin = dlamch_("Safe minimum", 12L);
  eps = dlamch_("Precision", 9L);
  smlnum = safmin / eps;
  bignum = 1. / smlnum;
  rmin = sqrt(smlnum);
  rmax = sqrt(bignum);

  /*     Scale matrix to allowable range, if necessary. */

  anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1], 1L, 1L);
  iscale = 0;
  if (anrm > 0. && anrm < rmin)
  {
    iscale = 1;
    sigma = rmin / anrm;
  }
  else if (anrm > rmax)
  {
    iscale = 1;
    sigma = rmax / anrm;
  }
  if (iscale == 1)
  {
    if (lower)
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        i__2 = *n - j + 1;
        dscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
        /* L10: */
      }
    }
    else
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        dscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
        /* L20: */
      }
    }
  }

  /*     Call DSYTRD to reduce symmetric matrix to tridiagonal form. */

  inde = 1;
  indtau = inde + *n;
  indwrk = indtau + *n;
  llwork = *lwork - indwrk + 1;
  dsytrd_(uplo, n, &a[a_offset], lda, &w[1], &work[inde], &work[indtau], &
          work[indwrk], &llwork, &iinfo, 1L);
  lopt = (integer)((*n << 1) + work[indwrk]);

  /*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
  /*     DORGTR to generate the orthogonal matrix, then call DSTEQR. */

  if (! wantz)
  {
    dsterf_(n, &w[1], &work[inde], info);
  }
  else
  {
    dorgtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
            llwork, &iinfo, 1L);
    dsteqr_(jobz, n, &w[1], &work[inde], &a[a_offset], lda, &work[indtau],
            info, 1L);
  }

  /*     If matrix was scaled, then rescale eigenvalues appropriately. */

  if (iscale == 1)
  {
    if (*info == 0)
    {
      imax = *n;
    }
    else
    {
      imax = *info - 1;
    }
    d__1 = 1. / sigma;
    dscal_(&imax, &d__1, &w[1], &c__1);
  }

  /*     Set WORK(1) to optimal workspace size. */

  /* Computing MAX */
  i__1 = *n * 3 - 1;
  work[1] = (doublereal) max(i__1, lopt);

  return 0;

  /*     End of DSYEV */

} /* dsyev_ */

/* Subroutine */ int dsteqr_(compz, n, d, e, z, ldz, work, info, compz_len)
char *compz;
integer *n;
doublereal *d, *e, *z;
integer *ldz;
doublereal *work;
integer *info;
ftnlen compz_len;
{
  /* System generated locals */
  integer z_dim1, z_offset, i__1, i__2;
  doublereal d__1, d__2;

  /* Builtin functions */
  double d_sign();

  /* Local variables */
  static integer lend, jtot;
  extern /* Subroutine */ int dlae2_();
  static doublereal b, c, f, g;
  static integer i, j, k, l, m;
  static doublereal p, r, s;
  extern logical lsame_();
  extern /* Subroutine */ int dlasr_(), dswap_();
  static integer l1;
  extern /* Subroutine */ int dlaev2_();
  static integer lendm1, lendp1;
  extern doublereal dlapy2_();
  static integer ii;
  extern doublereal dlamch_();
  static integer mm;
  extern /* Subroutine */ int dlartg_(), xerbla_(), dlazro_();
  static integer nmaxit, icompz, lm1, mm1, nm1;
  static doublereal rt1, rt2, eps, tst;


  /*  -- LAPACK routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     March 31, 1993 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a */
  /*  symmetric tridiagonal matrix using the implicit QL or QR method. */
  /*  The eigenvectors of a full or band symmetric matrix can also be found
  */
  /*  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
  */
  /*  tridiagonal form. */

  /*  Arguments */
  /*  ========= */

  /*  COMPZ   (input) CHARACTER*1 */
  /*          = 'N':  Compute eigenvalues only. */
  /*          = 'V':  Compute eigenvalues and eigenvectors of the original
  */
  /*                  symmetric matrix.  On entry, Z must contain the */
  /*                  orthogonal matrix used to reduce the original matrix
  */
  /*                  to tridiagonal form. */
  /*          = 'I':  Compute eigenvalues and eigenvectors of the */
  /*                  tridiagonal matrix.  Z is initialized to the identity
  */
  /*                  matrix. */

  /*  N       (input) INTEGER */
  /*          The order of the matrix.  N >= 0. */

  /*  D       (input/output) DOUBLE PRECISION array, dimension (N) */
  /*          On entry, the diagonal elements of the tridiagonal matrix. */
  /*          On exit, if INFO = 0, the eigenvalues in ascending order. */

  /*  E       (input/output) DOUBLE PRECISION array, dimension (N-1) */
  /*          On entry, the (n-1) subdiagonal elements of the tridiagonal */
  /*          matrix. */
  /*          On exit, E has been destroyed. */

  /*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N) */
  /*          On entry, if  COMPZ = 'V', then Z contains the orthogonal */
  /*          matrix used in the reduction to tridiagonal form. */
  /*          On exit, if  COMPZ = 'V', Z contains the orthonormal */
  /*          eigenvectors of the original symmetric matrix, and if */
  /*          COMPZ = 'I', Z contains the orthonormal eigenvectors of */
  /*          the symmetric tridiagonal matrix.  If an error exit is */
  /*          made, Z contains the eigenvectors associated with the */
  /*          stored eigenvalues. */
  /*          If COMPZ = 'N', then Z is not referenced. */

  /*  LDZ     (input) INTEGER */
  /*          The leading dimension of the array Z.  LDZ >= 1, and if */
  /*          eigenvectors are desired, then  LDZ >= max(1,N). */

  /*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
  */
  /*          If COMPZ = 'N', then WORK is not referenced. */

  /*  INFO    (output) INTEGER */
  /*          = 0:  successful exit */
  /*          < 0:  if INFO = -i, the i-th argument had an illegal value */
  /*          > 0:  the algorithm has failed to find all the eigenvalues in
  */
  /*                a total of 30*N iterations; if INFO = i, then i */
  /*                elements of E have not converged to zero; on exit, D */
  /*                and E contain the elements of a symmetric tridiagonal */
  /*                matrix which is orthogonally similar to the original */
  /*                matrix. */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  --work;
  z_dim1 = *ldz;
  z_offset = z_dim1 + 1;
  z -= z_offset;
  --e;
  --d;

  /* Function Body */
  *info = 0;

  if (lsame_(compz, "N", 1L, 1L))
  {
    icompz = 0;
  }
  else if (lsame_(compz, "V", 1L, 1L))
  {
    icompz = 1;
  }
  else if (lsame_(compz, "I", 1L, 1L))
  {
    icompz = 2;
  }
  else
  {
    icompz = -1;
  }
  if (icompz < 0)
  {
    *info = -1;
  }
  else if (*n < 0)
  {
    *info = -2;
  }
  else if (*ldz < 1 || icompz > 0 && *ldz < max(1, *n))
  {
    *info = -6;
  }
  if (*info != 0)
  {
    i__1 = -(*info);
    xerbla_("DSTEQR", &i__1, 6L);
    return 0;
  }

  /*     Quick return if possible */

  if (*n == 0)
  {
    return 0;
  }

  if (*n == 1)
  {
    if (icompz > 0)
    {
      z[z_dim1 + 1] = 1.;
    }
    return 0;
  }

  /*     Determine the unit roundoff for this environment. */

  eps = dlamch_("E", 1L);

  /*     Compute the eigenvalues and eigenvectors of the tridiagonal */
  /*     matrix. */

  if (icompz == 2)
  {
    dlazro_(n, n, &c_b21, &c_b22, &z[z_offset], ldz);
  }

  nmaxit = *n * 30;
  jtot = 0;

  /*     Determine where the matrix splits and choose QL or QR iteration */
  /*     for each block, according to whether top or bottom diagonal */
  /*     element is smaller. */

  l1 = 1;
  nm1 = *n - 1;

L10:
  if (l1 > *n)
  {
    goto L160;
  }
  if (l1 > 1)
  {
    e[l1 - 1] = 0.;
  }
  if (l1 <= nm1)
  {
    i__1 = nm1;
    for (m = l1; m <= i__1; ++m)
    {
      tst = (d__1 = e[m], abs(d__1));
      if (tst <= eps * ((d__1 = d[m], abs(d__1)) + (d__2 = d[m + 1],
                        abs(d__2))))
      {
        goto L30;
      }
      /* L20: */
    }
  }
  m = *n;

L30:
  l = l1;
  lend = m;
  if ((d__1 = d[lend], abs(d__1)) < (d__2 = d[l], abs(d__2)))
  {
    l = lend;
    lend = l1;
  }
  l1 = m + 1;

  if (lend >= l)
  {

    /*        QL Iteration */

    /*        Look for small subdiagonal element. */

L40:
    if (l != lend)
    {
      lendm1 = lend - 1;
      i__1 = lendm1;
      for (m = l; m <= i__1; ++m)
      {
        tst = (d__1 = e[m], abs(d__1));
        if (tst <= eps * ((d__1 = d[m], abs(d__1)) + (d__2 = d[m + 1],
                          abs(d__2))))
        {
          goto L60;
        }
        /* L50: */
      }
    }

    m = lend;

L60:
    if (m < lend)
    {
      e[m] = 0.;
    }
    p = d[l];
    if (m == l)
    {
      goto L80;
    }

    /*        If remaining matrix is 2-by-2, use DLAE2 or DLAEV2 */
    /*        to compute its eigensystem. */

    if (m == l + 1)
    {
      if (icompz > 0)
      {
        dlaev2_(&d[l], &e[l], &d[l + 1], &rt1, &rt2, &c, &s);
        work[l] = c;
        work[*n - 1 + l] = s;
        dlasr_("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
               z[l * z_dim1 + 1], ldz, 1L, 1L, 1L);
      }
      else
      {
        dlae2_(&d[l], &e[l], &d[l + 1], &rt1, &rt2);
      }
      d[l] = rt1;
      d[l + 1] = rt2;
      e[l] = 0.;
      l += 2;
      if (l <= lend)
      {
        goto L40;
      }
      goto L10;
    }

    if (jtot == nmaxit)
    {
      goto L140;
    }
    ++jtot;

    /*        Form shift. */

    g = (d[l + 1] - p) / (e[l] * 2.);
    r = dlapy2_(&g, &c_b22);
    g = d[m] - p + e[l] / (g + d_sign(&r, &g));

    s = 1.;
    c = 1.;
    p = 0.;

    /*        Inner loop */

    mm1 = m - 1;
    i__1 = l;
    for (i = mm1; i >= i__1; --i)
    {
      f = s * e[i];
      b = c * e[i];
      dlartg_(&g, &f, &c, &s, &r);
      if (i != m - 1)
      {
        e[i + 1] = r;
      }
      g = d[i + 1] - p;
      r = (d[i] - g) * s + c * 2. * b;
      p = s * r;
      d[i + 1] = g + p;
      g = c * r - b;

      /*           If eigenvectors are desired, then save rotations. */

      if (icompz > 0)
      {
        work[i] = c;
        work[*n - 1 + i] = -s;
      }

      /* L70: */
    }

    /*        If eigenvectors are desired, then apply saved rotations. */

    if (icompz > 0)
    {
      mm = m - l + 1;
      dlasr_("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &z[l *
             z_dim1 + 1], ldz, 1L, 1L, 1L);
    }

    d[l] -= p;
    e[l] = g;
    goto L40;

    /*        Eigenvalue found. */

L80:
    d[l] = p;

    ++l;
    if (l <= lend)
    {
      goto L40;
    }
    goto L10;

  }
  else
  {

    /*        QR Iteration */

    /*        Look for small superdiagonal element. */

L90:
    if (l != lend)
    {
      lendp1 = lend + 1;
      i__1 = lendp1;
      for (m = l; m >= i__1; --m)
      {
        tst = (d__1 = e[m - 1], abs(d__1));
        if (tst <= eps * ((d__1 = d[m], abs(d__1)) + (d__2 = d[m - 1],
                          abs(d__2))))
        {
          goto L110;
        }
        /* L100: */
      }
    }

    m = lend;

L110:
    if (m > lend)
    {
      e[m - 1] = 0.;
    }
    p = d[l];
    if (m == l)
    {
      goto L130;
    }

    /*        If remaining matrix is 2-by-2, use DLAE2 or DLAEV2 */
    /*        to compute its eigensystem. */

    if (m == l - 1)
    {
      if (icompz > 0)
      {
        dlaev2_(&d[l - 1], &e[l - 1], &d[l], &rt1, &rt2, &c, &s);
        work[m] = c;
        work[*n - 1 + m] = s;
        dlasr_("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
               z[(l - 1) * z_dim1 + 1], ldz, 1L, 1L, 1L);
      }
      else
      {
        dlae2_(&d[l - 1], &e[l - 1], &d[l], &rt1, &rt2);
      }
      d[l - 1] = rt1;
      d[l] = rt2;
      e[l - 1] = 0.;
      l += -2;
      if (l >= lend)
      {
        goto L90;
      }
      goto L10;
    }

    if (jtot == nmaxit)
    {
      goto L140;
    }
    ++jtot;

    /*        Form shift. */

    g = (d[l - 1] - p) / (e[l - 1] * 2.);
    r = dlapy2_(&g, &c_b22);
    g = d[m] - p + e[l - 1] / (g + d_sign(&r, &g));

    s = 1.;
    c = 1.;
    p = 0.;

    /*        Inner loop */

    lm1 = l - 1;
    i__1 = lm1;
    for (i = m; i <= i__1; ++i)
    {
      f = s * e[i];
      b = c * e[i];
      dlartg_(&g, &f, &c, &s, &r);
      if (i != m)
      {
        e[i - 1] = r;
      }
      g = d[i] - p;
      r = (d[i + 1] - g) * s + c * 2. * b;
      p = s * r;
      d[i] = g + p;
      g = c * r - b;

      /*           If eigenvectors are desired, then save rotations. */

      if (icompz > 0)
      {
        work[i] = c;
        work[*n - 1 + i] = s;
      }

      /* L120: */
    }

    /*        If eigenvectors are desired, then apply saved rotations. */

    if (icompz > 0)
    {
      mm = l - m + 1;
      dlasr_("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &z[m *
             z_dim1 + 1], ldz, 1L, 1L, 1L);
    }

    d[l] -= p;
    e[lm1] = g;
    goto L90;

    /*        Eigenvalue found. */

L130:
    d[l] = p;

    --l;
    if (l >= lend)
    {
      goto L90;
    }
    goto L10;

  }

  /*     Set error -- no convergence to an eigenvalue after a total */
  /*     of N*MAXIT iterations. */

L140:
  i__1 = *n - 1;
  for (i = 1; i <= i__1; ++i)
  {
    if (e[i] != 0.)
    {
      ++(*info);
    }
    /* L150: */
  }
  return 0;

  /*     Order eigenvalues and eigenvectors. */

L160:
  i__1 = *n;
  for (ii = 2; ii <= i__1; ++ii)
  {
    i = ii - 1;
    k = i;
    p = d[i];
    i__2 = *n;
    for (j = ii; j <= i__2; ++j)
    {
      if (d[j] < p)
      {
        k = j;
        p = d[j];
      }
      /* L170: */
    }
    if (k != i)
    {
      d[k] = d[i];
      d[i] = p;
      if (icompz > 0)
      {
        dswap_(n, &z[i * z_dim1 + 1], &c__1, &z[k * z_dim1 + 1], &
               c__1);
      }
    }
    /* L180: */
  }

  return 0;

  /*     End of DSTEQR */

} /* dsteqr_ */

/* Subroutine */ int dlartg_(f, g, cs, sn, r)
doublereal *f, *g, *cs, *sn, *r;
{
  /* Builtin functions */
  double sqrt();

  /* Local variables */
  static doublereal t, tt;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLARTG generate a plane rotation so that */

  /*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1. */
  /*     [ -SN  CS  ]     [ G ]     [ 0 ] */

  /*  This is a faster version of the BLAS1 routine DROTG, except for */
  /*  the following differences: */
  /*     F and G are unchanged on return. */
  /*     If G=0, then CS=1 and SN=0. */
  /*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any */
  /*        floating point operations (saves work in DBDSQR when */
  /*        there are zeros on the diagonal). */

  /*  Arguments */
  /*  ========= */

  /*  F       (input) DOUBLE PRECISION */
  /*          The first component of vector to be rotated. */

  /*  G       (input) DOUBLE PRECISION */
  /*          The second component of vector to be rotated. */

  /*  CS      (output) DOUBLE PRECISION */
  /*          The cosine of the rotation. */

  /*  SN      (output) DOUBLE PRECISION */
  /*          The sine of the rotation. */

  /*  R       (output) DOUBLE PRECISION */
  /*          The nonzero component of the rotated vector. */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  if (*g == 0.)
  {
    *cs = 1.;
    *sn = 0.;
    *r = *f;
  }
  else if (*f == 0.)
  {
    *cs = 0.;
    *sn = 1.;
    *r = *g;
  }
  else
  {
    if (abs(*f) > abs(*g))
    {
      t = *g / *f;
      tt = sqrt(t * t + 1.);
      *cs = 1. / tt;
      *sn = t * *cs;
      *r = *f * tt;
    }
    else
    {
      t = *f / *g;
      tt = sqrt(t * t + 1.);
      *sn = 1. / tt;
      *cs = t * *sn;
      *r = *g * tt;
    }
  }
  return 0;

  /*     End of DLARTG */

} /* dlartg_ */

/* Subroutine */ int dlasr_(side, pivot, direct, m, n, c, s, a, lda, side_len,
                            pivot_len, direct_len)
char *side, *pivot, *direct;
integer *m, *n;
doublereal *c, *s, *a;
integer *lda;
ftnlen side_len;
ftnlen pivot_len;
ftnlen direct_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;

  /* Local variables */
  static integer info;
  static doublereal temp;
  static integer i, j;
  extern logical lsame_();
  static doublereal ctemp, stemp;
  extern /* Subroutine */ int xerbla_();


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLASR   performs the transformation */

  /*     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side ) */

  /*     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side ) */

  /*  where A is an m by n real matrix and P is an orthogonal matrix, */
  /*  consisting of a sequence of plane rotations determined by the */
  /*  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
  */
  /*  and z = n when SIDE = 'R' or 'r' ): */

  /*  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then */

  /*     P = P( z - 1 )*...*P( 2 )*P( 1 ), */

  /*  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then */

  /*     P = P( 1 )*P( 2 )*...*P( z - 1 ), */

  /*  where  P( k ) is a plane rotation matrix for the following planes: */

  /*     when  PIVOT = 'V' or 'v'  ( Variable pivot ), */
  /*        the plane ( k, k + 1 ) */

  /*     when  PIVOT = 'T' or 't'  ( Top pivot ), */
  /*        the plane ( 1, k + 1 ) */

  /*     when  PIVOT = 'B' or 'b'  ( Bottom pivot ), */
  /*        the plane ( k, z ) */

  /*  c( k ) and s( k )  must contain the  cosine and sine that define the
  */
  /*  matrix  P( k ).  The two by two plane rotation part of the matrix */
  /*  P( k ), R( k ), is assumed to be of the form */

  /*     R( k ) = (  c( k )  s( k ) ). */
  /*              ( -s( k )  c( k ) ) */

  /*  This version vectorises across rows of the array A when SIDE = 'L'. */

  /*  Arguments */
  /*  ========= */

  /*  SIDE    (input) CHARACTER*1 */
  /*          Specifies whether the plane rotation matrix P is applied to */
  /*          A on the left or the right. */
  /*          = 'L':  Left, compute A := P*A */
  /*          = 'R':  Right, compute A:= A*P' */

  /*  DIRECT  (input) CHARACTER*1 */
  /*          Specifies whether P is a forward or backward sequence of */
  /*          plane rotations. */
  /*          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 ) */
  /*          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 ) */

  /*  PIVOT   (input) CHARACTER*1 */
  /*          Specifies the plane for which P(k) is a plane rotation */
  /*          matrix. */
  /*          = 'V':  Variable pivot, the plane (k,k+1) */
  /*          = 'T':  Top pivot, the plane (1,k+1) */
  /*          = 'B':  Bottom pivot, the plane (k,z) */

  /*  M       (input) INTEGER */
  /*          The number of rows of the matrix A.  If m <= 1, an immediate
  */
  /*          return is effected. */

  /*  N       (input) INTEGER */
  /*          The number of columns of the matrix A.  If n <= 1, an */
  /*          immediate return is effected. */

  /*  C, S    (input) DOUBLE PRECISION arrays, dimension */
  /*                  (M-1) if SIDE = 'L' */
  /*                  (N-1) if SIDE = 'R' */
  /*          c(k) and s(k) contain the cosine and sine that define the */
  /*          matrix P(k).  The two by two plane rotation part of the */
  /*          matrix P(k), R(k), is assumed to be of the form */
  /*          R( k ) = (  c( k )  s( k ) ). */
  /*                   ( -s( k )  c( k ) ) */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          The m by n matrix A.  On exit, A is overwritten by P*A if */
  /*          SIDE = 'R' or by A*P' if SIDE = 'L'. */

  /*  LDA     (input) INTEGER */
  /*          The leading dimension of the array A.  LDA >= max(1,M). */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters */

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;
  --s;
  --c;

  /* Function Body */
  info = 0;
  if (!(lsame_(side, "L", 1L, 1L) || lsame_(side, "R", 1L, 1L)))
  {
    info = 1;
  }
  else if (!(lsame_(pivot, "V", 1L, 1L) || lsame_(pivot, "T", 1L, 1L) ||
             lsame_(pivot, "B", 1L, 1L)))
  {
    info = 2;
  }
  else if (!(lsame_(direct, "F", 1L, 1L) || lsame_(direct, "B", 1L, 1L)))
  {
    info = 3;
  }
  else if (*m < 0)
  {
    info = 4;
  }
  else if (*n < 0)
  {
    info = 5;
  }
  else if (*lda < max(1, *m))
  {
    info = 9;
  }
  if (info != 0)
  {
    xerbla_("DLASR ", &info, 6L);
    return 0;
  }

  /*     Quick return if possible */

  if (*m == 0 || *n == 0)
  {
    return 0;
  }
  if (lsame_(side, "L", 1L, 1L))
  {

    /*        Form  P * A */

    if (lsame_(pivot, "V", 1L, 1L))
    {
      if (lsame_(direct, "F", 1L, 1L))
      {
        i__1 = *m - 1;
        for (j = 1; j <= i__1; ++j)
        {
          ctemp = c[j];
          stemp = s[j];
          if (ctemp != 1. || stemp != 0.)
          {
            i__2 = *n;
            for (i = 1; i <= i__2; ++i)
            {
              temp = a[j + 1 + i * a_dim1];
              a[j + 1 + i * a_dim1] = ctemp * temp - stemp * a[
                                        j + i * a_dim1];
              a[j + i * a_dim1] = stemp * temp + ctemp * a[j +
                                  i * a_dim1];
              /* L10: */
            }
          }
          /* L20: */
        }
      }
      else if (lsame_(direct, "B", 1L, 1L))
      {
        for (j = *m - 1; j >= 1; --j)
        {
          ctemp = c[j];
          stemp = s[j];
          if (ctemp != 1. || stemp != 0.)
          {
            i__1 = *n;
            for (i = 1; i <= i__1; ++i)
            {
              temp = a[j + 1 + i * a_dim1];
              a[j + 1 + i * a_dim1] = ctemp * temp - stemp * a[
                                        j + i * a_dim1];
              a[j + i * a_dim1] = stemp * temp + ctemp * a[j +
                                  i * a_dim1];
              /* L30: */
            }
          }
          /* L40: */
        }
      }
    }
    else if (lsame_(pivot, "T", 1L, 1L))
    {
      if (lsame_(direct, "F", 1L, 1L))
      {
        i__1 = *m;
        for (j = 2; j <= i__1; ++j)
        {
          ctemp = c[j - 1];
          stemp = s[j - 1];
          if (ctemp != 1. || stemp != 0.)
          {
            i__2 = *n;
            for (i = 1; i <= i__2; ++i)
            {
              temp = a[j + i * a_dim1];
              a[j + i * a_dim1] = ctemp * temp - stemp * a[i *
                                  a_dim1 + 1];
              a[i * a_dim1 + 1] = stemp * temp + ctemp * a[i *
                                  a_dim1 + 1];
              /* L50: */
            }
          }
          /* L60: */
        }
      }
      else if (lsame_(direct, "B", 1L, 1L))
      {
        for (j = *m; j >= 2; --j)
        {
          ctemp = c[j - 1];
          stemp = s[j - 1];
          if (ctemp != 1. || stemp != 0.)
          {
            i__1 = *n;
            for (i = 1; i <= i__1; ++i)
            {
              temp = a[j + i * a_dim1];
              a[j + i * a_dim1] = ctemp * temp - stemp * a[i *
                                  a_dim1 + 1];
              a[i * a_dim1 + 1] = stemp * temp + ctemp * a[i *
                                  a_dim1 + 1];
              /* L70: */
            }
          }
          /* L80: */
        }
      }
    }
    else if (lsame_(pivot, "B", 1L, 1L))
    {
      if (lsame_(direct, "F", 1L, 1L))
      {
        i__1 = *m - 1;
        for (j = 1; j <= i__1; ++j)
        {
          ctemp = c[j];
          stemp = s[j];
          if (ctemp != 1. || stemp != 0.)
          {
            i__2 = *n;
            for (i = 1; i <= i__2; ++i)
            {
              temp = a[j + i * a_dim1];
              a[j + i * a_dim1] = stemp * a[*m + i * a_dim1] +
                                  ctemp * temp;
              a[*m + i * a_dim1] = ctemp * a[*m + i * a_dim1] -
                                   stemp * temp;
              /* L90: */
            }
          }
          /* L100: */
        }
      }
      else if (lsame_(direct, "B", 1L, 1L))
      {
        for (j = *m - 1; j >= 1; --j)
        {
          ctemp = c[j];
          stemp = s[j];
          if (ctemp != 1. || stemp != 0.)
          {
            i__1 = *n;
            for (i = 1; i <= i__1; ++i)
            {
              temp = a[j + i * a_dim1];
              a[j + i * a_dim1] = stemp * a[*m + i * a_dim1] +
                                  ctemp * temp;
              a[*m + i * a_dim1] = ctemp * a[*m + i * a_dim1] -
                                   stemp * temp;
              /* L110: */
            }
          }
          /* L120: */
        }
      }
    }
  }
  else if (lsame_(side, "R", 1L, 1L))
  {

    /*        Form A * P' */

    if (lsame_(pivot, "V", 1L, 1L))
    {
      if (lsame_(direct, "F", 1L, 1L))
      {
        i__1 = *n - 1;
        for (j = 1; j <= i__1; ++j)
        {
          ctemp = c[j];
          stemp = s[j];
          if (ctemp != 1. || stemp != 0.)
          {
            i__2 = *m;
            for (i = 1; i <= i__2; ++i)
            {
              temp = a[i + (j + 1) * a_dim1];
              a[i + (j + 1) * a_dim1] = ctemp * temp - stemp *
                                        a[i + j * a_dim1];
              a[i + j * a_dim1] = stemp * temp + ctemp * a[i +
                                  j * a_dim1];
              /* L130: */
            }
          }
          /* L140: */
        }
      }
      else if (lsame_(direct, "B", 1L, 1L))
      {
        for (j = *n - 1; j >= 1; --j)
        {
          ctemp = c[j];
          stemp = s[j];
          if (ctemp != 1. || stemp != 0.)
          {
            i__1 = *m;
            for (i = 1; i <= i__1; ++i)
            {
              temp = a[i + (j + 1) * a_dim1];
              a[i + (j + 1) * a_dim1] = ctemp * temp - stemp *
                                        a[i + j * a_dim1];
              a[i + j * a_dim1] = stemp * temp + ctemp * a[i +
                                  j * a_dim1];
              /* L150: */
            }
          }
          /* L160: */
        }
      }
    }
    else if (lsame_(pivot, "T", 1L, 1L))
    {
      if (lsame_(direct, "F", 1L, 1L))
      {
        i__1 = *n;
        for (j = 2; j <= i__1; ++j)
        {
          ctemp = c[j - 1];
          stemp = s[j - 1];
          if (ctemp != 1. || stemp != 0.)
          {
            i__2 = *m;
            for (i = 1; i <= i__2; ++i)
            {
              temp = a[i + j * a_dim1];
              a[i + j * a_dim1] = ctemp * temp - stemp * a[i +
                                  a_dim1];
              a[i + a_dim1] = stemp * temp + ctemp * a[i +
                              a_dim1];
              /* L170: */
            }
          }
          /* L180: */
        }
      }
      else if (lsame_(direct, "B", 1L, 1L))
      {
        for (j = *n; j >= 2; --j)
        {
          ctemp = c[j - 1];
          stemp = s[j - 1];
          if (ctemp != 1. || stemp != 0.)
          {
            i__1 = *m;
            for (i = 1; i <= i__1; ++i)
            {
              temp = a[i + j * a_dim1];
              a[i + j * a_dim1] = ctemp * temp - stemp * a[i +
                                  a_dim1];
              a[i + a_dim1] = stemp * temp + ctemp * a[i +
                              a_dim1];
              /* L190: */
            }
          }
          /* L200: */
        }
      }
    }
    else if (lsame_(pivot, "B", 1L, 1L))
    {
      if (lsame_(direct, "F", 1L, 1L))
      {
        i__1 = *n - 1;
        for (j = 1; j <= i__1; ++j)
        {
          ctemp = c[j];
          stemp = s[j];
          if (ctemp != 1. || stemp != 0.)
          {
            i__2 = *m;
            for (i = 1; i <= i__2; ++i)
            {
              temp = a[i + j * a_dim1];
              a[i + j * a_dim1] = stemp * a[i + *n * a_dim1] +
                                  ctemp * temp;
              a[i + *n * a_dim1] = ctemp * a[i + *n * a_dim1] -
                                   stemp * temp;
              /* L210: */
            }
          }
          /* L220: */
        }
      }
      else if (lsame_(direct, "B", 1L, 1L))
      {
        for (j = *n - 1; j >= 1; --j)
        {
          ctemp = c[j];
          stemp = s[j];
          if (ctemp != 1. || stemp != 0.)
          {
            i__1 = *m;
            for (i = 1; i <= i__1; ++i)
            {
              temp = a[i + j * a_dim1];
              a[i + j * a_dim1] = stemp * a[i + *n * a_dim1] +
                                  ctemp * temp;
              a[i + *n * a_dim1] = ctemp * a[i + *n * a_dim1] -
                                   stemp * temp;
              /* L230: */
            }
          }
          /* L240: */
        }
      }
    }
  }

  return 0;

  /*     End of DLASR */

} /* dlasr_ */

/* Subroutine */ int dlaev2_(a, b, c, rt1, rt2, cs1, sn1)
doublereal *a, *b, *c, *rt1, *rt2, *cs1, *sn1;
{
  /* System generated locals */
  doublereal d__1;

  /* Builtin functions */
  double sqrt();

  /* Local variables */
  static doublereal acmn, acmx, ab, df, cs, ct, tb, sm, tn, rt, adf, acs;
  static integer sgn1, sgn2;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix */
  /*     [  A   B  ] */
  /*     [  B   C  ]. */
  /*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
  */
  /*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
  */
  /*  eigenvector for RT1, giving the decomposition */

  /*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ] */
  /*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ]. */

  /*  Arguments */
  /*  ========= */

  /*  A       (input) DOUBLE PRECISION */
  /*          The (1,1) entry of the 2-by-2 matrix. */

  /*  B       (input) DOUBLE PRECISION */
  /*          The (1,2) entry and the conjugate of the (2,1) entry of the */
  /*          2-by-2 matrix. */

  /*  C       (input) DOUBLE PRECISION */
  /*          The (2,2) entry of the 2-by-2 matrix. */

  /*  RT1     (output) DOUBLE PRECISION */
  /*          The eigenvalue of larger absolute value. */

  /*  RT2     (output) DOUBLE PRECISION */
  /*          The eigenvalue of smaller absolute value. */

  /*  CS1     (output) DOUBLE PRECISION */
  /*  SN1     (output) DOUBLE PRECISION */
  /*          The vector (CS1, SN1) is a unit right eigenvector for RT1. */

  /*  Further Details */
  /*  =============== */

  /*  RT1 is accurate to a few ulps barring over/underflow. */

  /*  RT2 may be inaccurate if there is massive cancellation in the */
  /*  determinant A*C-B*B; higher precision or correctly rounded or */
  /*  correctly truncated arithmetic would be needed to compute RT2 */
  /*  accurately in all cases. */

  /*  CS1 and SN1 are accurate to a few ulps barring over/underflow. */

  /*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
  */
  /*  Underflow is harmless if the input data is 0 or exceeds */
  /*     underflow_threshold / macheps. */

  /* =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Compute the eigenvalues */

  sm = *a + *c;
  df = *a - *c;
  adf = abs(df);
  tb = *b + *b;
  ab = abs(tb);
  if (abs(*a) > abs(*c))
  {
    acmx = *a;
    acmn = *c;
  }
  else
  {
    acmx = *c;
    acmn = *a;
  }
  if (adf > ab)
  {
    /* Computing 2nd power */
    d__1 = ab / adf;
    rt = adf * sqrt(d__1 * d__1 + 1.);
  }
  else if (adf < ab)
  {
    /* Computing 2nd power */
    d__1 = adf / ab;
    rt = ab * sqrt(d__1 * d__1 + 1.);
  }
  else
  {

    /*        Includes case AB=ADF=0 */

    rt = ab * sqrt(2.);
  }
  if (sm < 0.)
  {
    *rt1 = (sm - rt) * .5;
    sgn1 = -1;

    /*        Order of execution important. */
    /*        To get fully accurate smaller eigenvalue, */
    /*        next line needs to be executed in higher precision. */

    *rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
  }
  else if (sm > 0.)
  {
    *rt1 = (sm + rt) * .5;
    sgn1 = 1;

    /*        Order of execution important. */
    /*        To get fully accurate smaller eigenvalue, */
    /*        next line needs to be executed in higher precision. */

    *rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
  }
  else
  {

    /*        Includes case RT1 = RT2 = 0 */

    *rt1 = rt * .5;
    *rt2 = rt * -.5;
    sgn1 = 1;
  }

  /*     Compute the eigenvector */

  if (df >= 0.)
  {
    cs = df + rt;
    sgn2 = 1;
  }
  else
  {
    cs = df - rt;
    sgn2 = -1;
  }
  acs = abs(cs);
  if (acs > ab)
  {
    ct = -tb / cs;
    *sn1 = 1. / sqrt(ct * ct + 1.);
    *cs1 = ct * *sn1;
  }
  else
  {
    if (ab == 0.)
    {
      *cs1 = 1.;
      *sn1 = 0.;
    }
    else
    {
      tn = -cs / tb;
      *cs1 = 1. / sqrt(tn * tn + 1.);
      *sn1 = tn * *cs1;
    }
  }
  if (sgn1 == sgn2)
  {
    tn = *cs1;
    *cs1 = -(*sn1);
    *sn1 = tn;
  }
  return 0;

  /*     End of DLAEV2 */

} /* dlaev2_ */

/* Subroutine */ int dlazro_(m, n, alpha, beta, a, lda)
integer *m, *n;
doublereal *alpha, *beta, *a;
integer *lda;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;

  /* Local variables */
  static integer i, j;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAZRO initializes a 2-D array A to BETA on the diagonal and */
  /*  ALPHA on the offdiagonals. */

  /*  Arguments */
  /*  ========= */

  /*  M       (input) INTEGER */
  /*          The number of rows of the matrix A.  M >= 0. */

  /*  N       (input) INTEGER */
  /*          The number of columns of the matrix A.  N >= 0. */

  /*  ALPHA   (input) DOUBLE PRECISION */
  /*          The constant to which the offdiagonal elements are to be set.
  */

  /*  BETA    (input) DOUBLE PRECISION */
  /*          The constant to which the diagonal elements are to be set. */

  /*  A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On exit, the leading m by n submatrix of A is set such that */
  /*             A(i,j) = ALPHA,  1 <= i <= m, 1 <= j <= n, i <> j */
  /*             A(i,i) = BETA,   1 <= i <= min(m,n). */

  /*  LDA     (input) INTEGER */
  /*          The leading dimension of the array A.  LDA >= max(1,M). */

  /*  =====================================================================
  */

  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  i__1 = *n;
  for (j = 1; j <= i__1; ++j)
  {
    i__2 = *m;
    for (i = 1; i <= i__2; ++i)
    {
      a[i + j * a_dim1] = *alpha;
      /* L10: */
    }
    /* L20: */
  }

  i__1 = min(*m, *n);
  for (i = 1; i <= i__1; ++i)
  {
    a[i + i * a_dim1] = *beta;
    /* L30: */
  }

  return 0;

  /*     End of DLAZRO */

} /* dlazro_ */

/* Subroutine */ int dorgtr_(uplo, n, a, lda, tau, work, lwork, info,
                             uplo_len)
char *uplo;
integer *n;
doublereal *a;
integer *lda;
doublereal *tau, *work;
integer *lwork, *info;
ftnlen uplo_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2, i__3;

  /* Local variables */
  static integer i, j;
  extern logical lsame_();
  static integer iinfo;
  static logical upper;
  extern /* Subroutine */ int xerbla_(), dorgql_(), dorgqr_();


  /*  -- LAPACK routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     March 31, 1993 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DORGTR generates a real orthogonal matrix Q which is defined as the */
  /*  product of n-1 elementary reflectors of order N, as returned by */
  /*  DSYTRD: */

  /*  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1), */

  /*  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1). */

  /*  Arguments */
  /*  ========= */

  /*  UPLO    (input) CHARACTER*1 */
  /*          = 'U': Upper triangle of A contains elementary reflectors */
  /*                 from DSYTRD; */
  /*          = 'L': Lower triangle of A contains elementary reflectors */
  /*                 from DSYTRD. */

  /*  N       (input) INTEGER */
  /*          The order of the matrix Q. N >= 0. */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On entry, the vectors which define the elementary reflectors,
  */
  /*          as returned by DSYTRD. */
  /*          On exit, the N-by-N orthogonal matrix Q. */

  /*  LDA     (input) INTEGER */
  /*          The leading dimension of the array A. LDA >= max(1,N). */

  /*  TAU     (input) DOUBLE PRECISION array, dimension (N-1) */
  /*          TAU(i) must contain the scalar factor of the elementary */
  /*          reflector H(i), as returned by DSYTRD. */

  /*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK) */
  /*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

  /*  LWORK   (input) INTEGER */
  /*          The dimension of the array WORK. LWORK >= max(1,N-1). */
  /*          For optimum performance LWORK >= (N-1)*NB, where NB is */
  /*          the optimal blocksize. */

  /*  INFO    (output) INTEGER */
  /*          = 0:  successful exit */
  /*          < 0:  if INFO = -i, the i-th argument had an illegal value */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input arguments */

  /* Parameter adjustments */
  --work;
  --tau;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  *info = 0;
  upper = lsame_(uplo, "U", 1L, 1L);
  if (! upper && ! lsame_(uplo, "L", 1L, 1L))
  {
    *info = -1;
  }
  else if (*n < 0)
  {
    *info = -2;
  }
  else if (*lda < max(1, *n))
  {
    *info = -4;
  }
  else /* if(complicated condition) */
  {
    /* Computing MAX */
    i__1 = 1, i__2 = *n - 1;
    if (*lwork < max(i__1, i__2))
    {
      *info = -7;
    }
  }
  if (*info != 0)
  {
    i__1 = -(*info);
    xerbla_("DORGTR", &i__1, 6L);
    return 0;
  }

  /*     Quick return if possible */

  if (*n == 0)
  {
    work[1] = 1.;
    return 0;
  }

  if (upper)
  {

    /*        Q was determined by a call to DSYTRD with UPLO = 'U' */

    /*        Shift the vectors which define the elementary reflectors one
     */
    /*        column to the left, and set the last row and column of Q to
    */
    /*        those of the unit matrix */

    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j)
    {
      i__2 = j - 1;
      for (i = 1; i <= i__2; ++i)
      {
        a[i + j * a_dim1] = a[i + (j + 1) * a_dim1];
        /* L10: */
      }
      a[*n + j * a_dim1] = 0.;
      /* L20: */
    }
    i__1 = *n - 1;
    for (i = 1; i <= i__1; ++i)
    {
      a[i + *n * a_dim1] = 0.;
      /* L30: */
    }
    a[*n + *n * a_dim1] = 1.;

    /*        Generate Q(1:n-1,1:n-1) */

    i__1 = *n - 1;
    i__2 = *n - 1;
    i__3 = *n - 1;
    dorgql_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1],
            lwork, &iinfo);

  }
  else
  {

    /*        Q was determined by a call to DSYTRD with UPLO = 'L'. */

    /*        Shift the vectors which define the elementary reflectors one
     */
    /*        column to the right, and set the first row and column of Q t
    o */
    /*        those of the unit matrix */

    for (j = *n; j >= 2; --j)
    {
      a[j * a_dim1 + 1] = 0.;
      i__1 = *n;
      for (i = j + 1; i <= i__1; ++i)
      {
        a[i + j * a_dim1] = a[i + (j - 1) * a_dim1];
        /* L40: */
      }
      /* L50: */
    }
    a[a_dim1 + 1] = 1.;
    i__1 = *n;
    for (i = 2; i <= i__1; ++i)
    {
      a[i + a_dim1] = 0.;
      /* L60: */
    }
    if (*n > 1)
    {

      /*           Generate Q(2:n,2:n) */

      i__1 = *n - 1;
      i__2 = *n - 1;
      i__3 = *n - 1;
      dorgqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[1],
              &work[1], lwork, &iinfo);
    }
  }
  return 0;

  /*     End of DORGTR */

} /* dorgtr_ */

/* Subroutine */ int dorgqr_(m, n, k, a, lda, tau, work, lwork, info)
integer *m, *n, *k;
doublereal *a;
integer *lda;
doublereal *tau, *work;
integer *lwork, *info;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2, i__3;

  /* Local variables */
  static integer i, j, l, nbmin, iinfo;
  extern /* Subroutine */ int dorg2r_();
  static integer ib, nb, ki, kk;
  extern /* Subroutine */ int dlarfb_();
  static integer nx;
  extern /* Subroutine */ int dlarft_(), xerbla_();
  extern integer ilaenv_();
  static integer ldwork, iws;


  /*  -- LAPACK routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     March 31, 1993 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DORGQR generates an M-by-N real matrix Q with orthonormal columns, */
  /*  which is defined as the first N columns of a product of K elementary
  */
  /*  reflectors of order M */

  /*        Q  =  H(1) H(2) . . . H(k) */

  /*  as returned by DGEQRF. */

  /*  Arguments */
  /*  ========= */

  /*  M       (input) INTEGER */
  /*          The number of rows of the matrix Q. M >= 0. */

  /*  N       (input) INTEGER */
  /*          The number of columns of the matrix Q. M >= N >= 0. */

  /*  K       (input) INTEGER */
  /*          The number of elementary reflectors whose product defines the
  */
  /*          matrix Q. N >= K >= 0. */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On entry, the i-th column must contain the vector which */
  /*          defines the elementary reflector H(i), for i = 1,2,...,k, as
  */
  /*          returned by DGEQRF in the first k columns of its array */
  /*          argument A. */
  /*          On exit, the M-by-N matrix Q. */

  /*  LDA     (input) INTEGER */
  /*          The first dimension of the array A. LDA >= max(1,M). */

  /*  TAU     (input) DOUBLE PRECISION array, dimension (K) */
  /*          TAU(i) must contain the scalar factor of the elementary */
  /*          reflector H(i), as returned by DGEQRF. */

  /*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK) */
  /*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

  /*  LWORK   (input) INTEGER */
  /*          The dimension of the array WORK. LWORK >= max(1,N). */
  /*          For optimum performance LWORK >= N*NB, where NB is the */
  /*          optimal blocksize. */

  /*  INFO    (output) INTEGER */
  /*          = 0:  successful exit */
  /*          < 0:  if INFO = -i, the i-th argument has an illegal value */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input arguments */

  /* Parameter adjustments */
  --work;
  --tau;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  *info = 0;
  if (*m < 0)
  {
    *info = -1;
  }
  else if (*n < 0 || *n > *m)
  {
    *info = -2;
  }
  else if (*k < 0 || *k > *n)
  {
    *info = -3;
  }
  else if (*lda < max(1, *m))
  {
    *info = -5;
  }
  else if (*lwork < max(1, *n))
  {
    *info = -8;
  }
  if (*info != 0)
  {
    i__1 = -(*info);
    xerbla_("DORGQR", &i__1, 6L);
    return 0;
  }

  /*     Quick return if possible */

  if (*n <= 0)
  {
    work[1] = 1.;
    return 0;
  }

  /*     Determine the block size. */

  nb = ilaenv_(&c__1, "DORGQR", " ", m, n, k, &c_n1, 6L, 1L);
  nbmin = 2;
  nx = 0;
  iws = *n;
  if (nb > 1 && nb < *k)
  {

    /*        Determine when to cross over from blocked to unblocked code.
     */

    /* Computing MAX */
    i__1 = 0, i__2 = ilaenv_(&c__3, "DORGQR", " ", m, n, k, &c_n1, 6L, 1L)
                     ;
    nx = max(i__1, i__2);
    if (nx < *k)
    {

      /*           Determine if workspace is large enough for blocked co
      de. */

      ldwork = *n;
      iws = ldwork * nb;
      if (*lwork < iws)
      {

        /*              Not enough workspace to use optimal NB:  reduc
        e NB and */
        /*              determine the minimum value of NB. */

        nb = *lwork / ldwork;
        /* Computing MAX */
        i__1 = 2, i__2 = ilaenv_(&c__2, "DORGQR", " ", m, n, k, &c_n1,
                                 6L, 1L);
        nbmin = max(i__1, i__2);
      }
    }
  }

  if (nb >= nbmin && nb < *k && nx < *k)
  {

    /*        Use blocked code after the last block. */
    /*        The first kk columns are handled by the block method. */

    ki = (*k - nx - 1) / nb * nb;
    /* Computing MIN */
    i__1 = *k, i__2 = ki + nb;
    kk = min(i__1, i__2);

    /*        Set A(1:kk,kk+1:n) to zero. */

    i__1 = *n;
    for (j = kk + 1; j <= i__1; ++j)
    {
      i__2 = kk;
      for (i = 1; i <= i__2; ++i)
      {
        a[i + j * a_dim1] = 0.;
        /* L10: */
      }
      /* L20: */
    }
  }
  else
  {
    kk = 0;
  }

  /*     Use unblocked code for the last or only block. */

  if (kk < *n)
  {
    i__1 = *m - kk;
    i__2 = *n - kk;
    i__3 = *k - kk;
    dorg2r_(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
            tau[kk + 1], &work[1], &iinfo);
  }

  if (kk > 0)
  {

    /*        Use blocked code */

    i__1 = -nb;
    for (i = ki + 1; i__1 < 0 ? i >= 1 : i <= 1; i += i__1)
    {
      /* Computing MIN */
      i__2 = nb, i__3 = *k - i + 1;
      ib = min(i__2, i__3);
      if (i + ib <= *n)
      {

        /*              Form the triangular factor of the block reflec
        tor */
        /*              H = H(i) H(i+1) . . . H(i+ib-1) */

        i__2 = *m - i + 1;
        dlarft_("Forward", "Columnwise", &i__2, &ib, &a[i + i *
                a_dim1], lda, &tau[i], &work[1], &ldwork, 7L, 10L);

        /*              Apply H to A(i:m,i+ib:n) from the left */

        i__2 = *m - i + 1;
        i__3 = *n - i - ib + 1;
        dlarfb_("Left", "No transpose", "Forward", "Columnwise", &
                i__2, &i__3, &ib, &a[i + i * a_dim1], lda, &work[1], &
                ldwork, &a[i + (i + ib) * a_dim1], lda, &work[ib + 1],
                &ldwork, 4L, 12L, 7L, 10L);
      }

      /*           Apply H to rows i:m of current block */

      i__2 = *m - i + 1;
      dorg2r_(&i__2, &ib, &ib, &a[i + i * a_dim1], lda, &tau[i], &work[
                1], &iinfo);

      /*           Set rows 1:i-1 of current block to zero */

      i__2 = i + ib - 1;
      for (j = i; j <= i__2; ++j)
      {
        i__3 = i - 1;
        for (l = 1; l <= i__3; ++l)
        {
          a[l + j * a_dim1] = 0.;
          /* L30: */
        }
        /* L40: */
      }
      /* L50: */
    }
  }

  work[1] = (doublereal) iws;
  return 0;

  /*     End of DORGQR */

} /* dorgqr_ */

/* Subroutine */ int dorg2r_(m, n, k, a, lda, tau, work, info)
integer *m, *n, *k;
doublereal *a;
integer *lda;
doublereal *tau, *work;
integer *info;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;
  doublereal d__1;

  /* Local variables */
  static integer i, j, l;
  extern /* Subroutine */ int dscal_(), dlarf_(), xerbla_();


  /*  -- LAPACK routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     February 29, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DORG2R generates an m by n real matrix Q with orthonormal columns, */
  /*  which is defined as the first n columns of a product of k elementary
  */
  /*  reflectors of order m */

  /*        Q  =  H(1) H(2) . . . H(k) */

  /*  as returned by DGEQRF. */

  /*  Arguments */
  /*  ========= */

  /*  M       (input) INTEGER */
  /*          The number of rows of the matrix Q. M >= 0. */

  /*  N       (input) INTEGER */
  /*          The number of columns of the matrix Q. M >= N >= 0. */

  /*  K       (input) INTEGER */
  /*          The number of elementary reflectors whose product defines the
  */
  /*          matrix Q. N >= K >= 0. */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On entry, the i-th column must contain the vector which */
  /*          defines the elementary reflector H(i), for i = 1,2,...,k, as
  */
  /*          returned by DGEQRF in the first k columns of its array */
  /*          argument A. */
  /*          On exit, the m-by-n matrix Q. */

  /*  LDA     (input) INTEGER */
  /*          The first dimension of the array A. LDA >= max(1,M). */

  /*  TAU     (input) DOUBLE PRECISION array, dimension (K) */
  /*          TAU(i) must contain the scalar factor of the elementary */
  /*          reflector H(i), as returned by DGEQRF. */

  /*  WORK    (workspace) DOUBLE PRECISION array, dimension (N) */

  /*  INFO    (output) INTEGER */
  /*          = 0: successful exit */
  /*          < 0: if INFO = -i, the i-th argument has an illegal value */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input arguments */

  /* Parameter adjustments */
  --work;
  --tau;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  *info = 0;
  if (*m < 0)
  {
    *info = -1;
  }
  else if (*n < 0 || *n > *m)
  {
    *info = -2;
  }
  else if (*k < 0 || *k > *n)
  {
    *info = -3;
  }
  else if (*lda < max(1, *m))
  {
    *info = -5;
  }
  if (*info != 0)
  {
    i__1 = -(*info);
    xerbla_("DORG2R", &i__1, 6L);
    return 0;
  }

  /*     Quick return if possible */

  if (*n <= 0)
  {
    return 0;
  }

  /*     Initialise columns k+1:n to columns of the unit matrix */

  i__1 = *n;
  for (j = *k + 1; j <= i__1; ++j)
  {
    i__2 = *m;
    for (l = 1; l <= i__2; ++l)
    {
      a[l + j * a_dim1] = 0.;
      /* L10: */
    }
    a[j + j * a_dim1] = 1.;
    /* L20: */
  }

  for (i = *k; i >= 1; --i)
  {

    /*        Apply H(i) to A(i:m,i:n) from the left */

    if (i < *n)
    {
      a[i + i * a_dim1] = 1.;
      i__1 = *m - i + 1;
      i__2 = *n - i;
      dlarf_("Left", &i__1, &i__2, &a[i + i * a_dim1], &c__1, &tau[i], &
             a[i + (i + 1) * a_dim1], lda, &work[1], 4L);
    }
    if (i < *m)
    {
      i__1 = *m - i;
      d__1 = -tau[i];
      dscal_(&i__1, &d__1, &a[i + 1 + i * a_dim1], &c__1);
    }
    a[i + i * a_dim1] = 1. - tau[i];

    /*        Set A(1:i-1,i) to zero */

    i__1 = i - 1;
    for (l = 1; l <= i__1; ++l)
    {
      a[l + i * a_dim1] = 0.;
      /* L30: */
    }
    /* L40: */
  }
  return 0;

  /*     End of DORG2R */

} /* dorg2r_ */

/* Subroutine */ int dorgql_(m, n, k, a, lda, tau, work, lwork, info)
integer *m, *n, *k;
doublereal *a;
integer *lda;
doublereal *tau, *work;
integer *lwork, *info;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

  /* Local variables */
  static integer i, j, l, nbmin, iinfo;
  extern /* Subroutine */ int dorg2l_();
  static integer ib, nb, kk;
  extern /* Subroutine */ int dlarfb_();
  static integer nx;
  extern /* Subroutine */ int dlarft_(), xerbla_();
  extern integer ilaenv_();
  static integer ldwork, iws;


  /*  -- LAPACK routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     March 31, 1993 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DORGQL generates an M-by-N real matrix Q with orthonormal columns, */
  /*  which is defined as the last N columns of a product of K elementary */
  /*  reflectors of order M */

  /*        Q  =  H(k) . . . H(2) H(1) */

  /*  as returned by DGEQLF. */

  /*  Arguments */
  /*  ========= */

  /*  M       (input) INTEGER */
  /*          The number of rows of the matrix Q. M >= 0. */

  /*  N       (input) INTEGER */
  /*          The number of columns of the matrix Q. M >= N >= 0. */

  /*  K       (input) INTEGER */
  /*          The number of elementary reflectors whose product defines the
  */
  /*          matrix Q. N >= K >= 0. */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On entry, the (n-k+i)-th column must contain the vector which
  */
  /*          defines the elementary reflector H(i), for i = 1,2,...,k, as
  */
  /*          returned by DGEQLF in the last k columns of its array */
  /*          argument A. */
  /*          On exit, the M-by-N matrix Q. */

  /*  LDA     (input) INTEGER */
  /*          The first dimension of the array A. LDA >= max(1,M). */

  /*  TAU     (input) DOUBLE PRECISION array, dimension (K) */
  /*          TAU(i) must contain the scalar factor of the elementary */
  /*          reflector H(i), as returned by DGEQLF. */

  /*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK) */
  /*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

  /*  LWORK   (input) INTEGER */
  /*          The dimension of the array WORK. LWORK >= max(1,N). */
  /*          For optimum performance LWORK >= N*NB, where NB is the */
  /*          optimal blocksize. */

  /*  INFO    (output) INTEGER */
  /*          = 0:  successful exit */
  /*          < 0:  if INFO = -i, the i-th argument has an illegal value */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input arguments */

  /* Parameter adjustments */
  --work;
  --tau;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  *info = 0;
  if (*m < 0)
  {
    *info = -1;
  }
  else if (*n < 0 || *n > *m)
  {
    *info = -2;
  }
  else if (*k < 0 || *k > *n)
  {
    *info = -3;
  }
  else if (*lda < max(1, *m))
  {
    *info = -5;
  }
  else if (*lwork < max(1, *n))
  {
    *info = -8;
  }
  if (*info != 0)
  {
    i__1 = -(*info);
    xerbla_("DORGQL", &i__1, 6L);
    return 0;
  }

  /*     Quick return if possible */

  if (*n <= 0)
  {
    work[1] = 1.;
    return 0;
  }

  /*     Determine the block size. */

  nb = ilaenv_(&c__1, "DORGQL", " ", m, n, k, &c_n1, 6L, 1L);
  nbmin = 2;
  nx = 0;
  iws = *n;
  if (nb > 1 && nb < *k)
  {

    /*        Determine when to cross over from blocked to unblocked code.
     */

    /* Computing MAX */
    i__1 = 0, i__2 = ilaenv_(&c__3, "DORGQL", " ", m, n, k, &c_n1, 6L, 1L)
                     ;
    nx = max(i__1, i__2);
    if (nx < *k)
    {

      /*           Determine if workspace is large enough for blocked co
      de. */

      ldwork = *n;
      iws = ldwork * nb;
      if (*lwork < iws)
      {

        /*              Not enough workspace to use optimal NB:  reduc
        e NB and */
        /*              determine the minimum value of NB. */

        nb = *lwork / ldwork;
        /* Computing MAX */
        i__1 = 2, i__2 = ilaenv_(&c__2, "DORGQL", " ", m, n, k, &c_n1,
                                 6L, 1L);
        nbmin = max(i__1, i__2);
      }
    }
  }

  if (nb >= nbmin && nb < *k && nx < *k)
  {

    /*        Use blocked code after the first block. */
    /*        The last kk columns are handled by the block method. */

    /* Computing MIN */
    i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
    kk = min(i__1, i__2);

    /*        Set A(m-kk+1:m,1:n-kk) to zero. */

    i__1 = *n - kk;
    for (j = 1; j <= i__1; ++j)
    {
      i__2 = *m;
      for (i = *m - kk + 1; i <= i__2; ++i)
      {
        a[i + j * a_dim1] = 0.;
        /* L10: */
      }
      /* L20: */
    }
  }
  else
  {
    kk = 0;
  }

  /*     Use unblocked code for the first or only block. */

  i__1 = *m - kk;
  i__2 = *n - kk;
  i__3 = *k - kk;
  dorg2l_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &iinfo)
  ;

  if (kk > 0)
  {

    /*        Use blocked code */

    i__1 = *k;
    i__2 = nb;
    for (i = *k - kk + 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2)
    {
      /* Computing MIN */
      i__3 = nb, i__4 = *k - i + 1;
      ib = min(i__3, i__4);
      if (*n - *k + i > 1)
      {

        /*              Form the triangular factor of the block reflec
        tor */
        /*              H = H(i+ib-1) . . . H(i+1) H(i) */

        i__3 = *m - *k + i + ib - 1;
        dlarft_("Backward", "Columnwise", &i__3, &ib, &a[(*n - *k + i)
                * a_dim1 + 1], lda, &tau[i], &work[1], &ldwork, 8L,
                10L);

        /*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the
        left */

        i__3 = *m - *k + i + ib - 1;
        i__4 = *n - *k + i - 1;
        dlarfb_("Left", "No transpose", "Backward", "Columnwise", &
                i__3, &i__4, &ib, &a[(*n - *k + i) * a_dim1 + 1], lda,
                &work[1], &ldwork, &a[a_offset], lda, &work[ib + 1],
                &ldwork, 4L, 12L, 8L, 10L);
      }

      /*           Apply H to rows 1:m-k+i+ib-1 of current block */

      i__3 = *m - *k + i + ib - 1;
      dorg2l_(&i__3, &ib, &ib, &a[(*n - *k + i) * a_dim1 + 1], lda, &
              tau[i], &work[1], &iinfo);

      /*           Set rows m-k+i+ib:m of current block to zero */

      i__3 = *n - *k + i + ib - 1;
      for (j = *n - *k + i; j <= i__3; ++j)
      {
        i__4 = *m;
        for (l = *m - *k + i + ib; l <= i__4; ++l)
        {
          a[l + j * a_dim1] = 0.;
          /* L30: */
        }
        /* L40: */
      }
      /* L50: */
    }
  }

  work[1] = (doublereal) iws;
  return 0;

  /*     End of DORGQL */

} /* dorgql_ */

/* Subroutine */ int dlarfb_(side, trans, direct, storev, m, n, k, v, ldv, t,
                             ldt, c, ldc, work, ldwork, side_len, trans_len, direct_len,
                             storev_len)
char *side, *trans, *direct, *storev;
integer *m, *n, *k;
doublereal *v;
integer *ldv;
doublereal *t;
integer *ldt;
doublereal *c;
integer *ldc;
doublereal *work;
integer *ldwork;
ftnlen side_len;
ftnlen trans_len;
ftnlen direct_len;
ftnlen storev_len;
{
  /* System generated locals */
  integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1,
          work_offset, i__1, i__2;

  /* Local variables */
  static integer i, j;
  extern /* Subroutine */ int dgemm_();
  extern logical lsame_();
  extern /* Subroutine */ int dcopy_(), dtrmm_();
  static char transt[1];


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     February 29, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLARFB applies a real block reflector H or its transpose H' to a */
  /*  real m by n matrix C, from either the left or the right. */

  /*  Arguments */
  /*  ========= */

  /*  SIDE    (input) CHARACTER*1 */
  /*          = 'L': apply H or H' from the Left */
  /*          = 'R': apply H or H' from the Right */

  /*  TRANS   (input) CHARACTER*1 */
  /*          = 'N': apply H (No transpose) */
  /*          = 'T': apply H' (Transpose) */

  /*  DIRECT  (input) CHARACTER*1 */
  /*          Indicates how H is formed from a product of elementary */
  /*          reflectors */
  /*          = 'F': H = H(1) H(2) . . . H(k) (Forward) */
  /*          = 'B': H = H(k) . . . H(2) H(1) (Backward) */

  /*  STOREV  (input) CHARACTER*1 */
  /*          Indicates how the vectors which define the elementary */
  /*          reflectors are stored: */
  /*          = 'C': Columnwise */
  /*          = 'R': Rowwise */

  /*  M       (input) INTEGER */
  /*          The number of rows of the matrix C. */

  /*  N       (input) INTEGER */
  /*          The number of columns of the matrix C. */

  /*  K       (input) INTEGER */
  /*          The order of the matrix T (= the number of elementary */
  /*          reflectors whose product defines the block reflector). */

  /*  V       (input) DOUBLE PRECISION array, dimension */
  /*                                (LDV,K) if STOREV = 'C' */
  /*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
  */
  /*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
  */
  /*          The matrix V. See further details. */

  /*  LDV     (input) INTEGER */
  /*          The leading dimension of the array V. */
  /*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M); */
  /*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N); */
  /*          if STOREV = 'R', LDV >= K. */

  /*  T       (input) DOUBLE PRECISION array, dimension (LDT,K) */
  /*          The triangular k by k matrix T in the representation of the */
  /*          block reflector. */

  /*  LDT     (input) INTEGER */
  /*          The leading dimension of the array T. LDT >= K. */

  /*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
  /*          On entry, the m by n matrix C. */
  /*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'. */

  /*  LDC     (input) INTEGER */
  /*          The leading dimension of the array C. LDA >= max(1,M). */

  /*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K) */

  /*  LDWORK  (input) INTEGER */
  /*          The leading dimension of the array WORK. */
  /*          If SIDE = 'L', LDWORK >= max(1,N); */
  /*          if SIDE = 'R', LDWORK >= max(1,M). */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Quick return if possible */

  /* Parameter adjustments */
  work_dim1 = *ldwork;
  work_offset = work_dim1 + 1;
  work -= work_offset;
  c_dim1 = *ldc;
  c_offset = c_dim1 + 1;
  c -= c_offset;
  t_dim1 = *ldt;
  t_offset = t_dim1 + 1;
  t -= t_offset;
  v_dim1 = *ldv;
  v_offset = v_dim1 + 1;
  v -= v_offset;

  /* Function Body */
  if (*m <= 0 || *n <= 0)
  {
    return 0;
  }

  if (lsame_(trans, "N", 1L, 1L))
  {
    *transt = 'T';
  }
  else
  {
    *transt = 'N';
  }

  if (lsame_(storev, "C", 1L, 1L))
  {

    if (lsame_(direct, "F", 1L, 1L))
    {

      /*           Let  V =  ( V1 )    (first K rows) */
      /*                     ( V2 ) */
      /*           where  V1  is unit lower triangular. */

      if (lsame_(side, "L", 1L, 1L))
      {

        /*              Form  H * C  or  H' * C  where  C = ( C1 ) */
        /*                                                  ( C2 ) */

        /*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in
        WORK) */

        /*              W := C1' */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          dcopy_(n, &c[j + c_dim1], ldc, &work[j * work_dim1 + 1], &
                 c__1);
          /* L10: */
        }

        /*              W := W * V1 */

        dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b22,
               &v[v_offset], ldv, &work[work_offset], ldwork, 5L,
               5L, 12L, 4L);
        if (*m > *k)
        {

          /*                 W := W + C2'*V2 */

          i__1 = *m - *k;
          dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b22, &
                 c[*k + 1 + c_dim1], ldc, &v[*k + 1 + v_dim1], ldv,
                 &c_b22, &work[work_offset], ldwork, 9L, 12L);
        }

        /*              W := W * T'  or  W * T */

        dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b22, &t[
                 t_offset], ldt, &work[work_offset], ldwork, 5L, 5L,
               1L, 8L);

        /*              C := C - V * W' */

        if (*m > *k)
        {

          /*                 C2 := C2 - V2 * W' */

          i__1 = *m - *k;
          dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b211,
                 &v[*k + 1 + v_dim1], ldv, &work[work_offset],
                 ldwork, &c_b22, &c[*k + 1 + c_dim1], ldc, 12L, 9L)
          ;
        }

        /*              W := W * V1' */

        dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b22, &
               v[v_offset], ldv, &work[work_offset], ldwork, 5L, 5L,
               9L, 4L);

        /*              C1 := C1 - W' */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *n;
          for (i = 1; i <= i__2; ++i)
          {
            c[j + i * c_dim1] -= work[i + j * work_dim1];
            /* L20: */
          }
          /* L30: */
        }

      }
      else if (lsame_(side, "R", 1L, 1L))
      {

        /*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
         */

        /*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WOR
        K) */

        /*              W := C1 */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          dcopy_(m, &c[j * c_dim1 + 1], &c__1, &work[j * work_dim1
                 + 1], &c__1);
          /* L40: */
        }

        /*              W := W * V1 */

        dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b22,
               &v[v_offset], ldv, &work[work_offset], ldwork, 5L,
               5L, 12L, 4L);
        if (*n > *k)
        {

          /*                 W := W + C2 * V2 */

          i__1 = *n - *k;
          dgemm_("No transpose", "No transpose", m, k, &i__1, &
                 c_b22, &c[(*k + 1) * c_dim1 + 1], ldc, &v[*k + 1
                     + v_dim1], ldv, &c_b22, &work[work_offset],
                 ldwork, 12L, 12L);
        }

        /*              W := W * T  or  W * T' */

        dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b22, &t[
                 t_offset], ldt, &work[work_offset], ldwork, 5L, 5L,
               1L, 8L);

        /*              C := C - W * V' */

        if (*n > *k)
        {

          /*                 C2 := C2 - W * V2' */

          i__1 = *n - *k;
          dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b211,
                 &work[work_offset], ldwork, &v[*k + 1 + v_dim1],
                 ldv, &c_b22, &c[(*k + 1) * c_dim1 + 1], ldc, 12L,
                 9L);
        }

        /*              W := W * V1' */

        dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b22, &
               v[v_offset], ldv, &work[work_offset], ldwork, 5L, 5L,
               9L, 4L);

        /*              C1 := C1 - W */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + j * c_dim1] -= work[i + j * work_dim1];
            /* L50: */
          }
          /* L60: */
        }
      }

    }
    else
    {

      /*           Let  V =  ( V1 ) */
      /*                     ( V2 )    (last K rows) */
      /*           where  V2  is unit upper triangular. */

      if (lsame_(side, "L", 1L, 1L))
      {

        /*              Form  H * C  or  H' * C  where  C = ( C1 ) */
        /*                                                  ( C2 ) */

        /*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in
        WORK) */

        /*              W := C2' */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          dcopy_(n, &c[*m - *k + j + c_dim1], ldc, &work[j *
                 work_dim1 + 1], &c__1);
          /* L70: */
        }

        /*              W := W * V2 */

        dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b22,
               &v[*m - *k + 1 + v_dim1], ldv, &work[work_offset],
               ldwork, 5L, 5L, 12L, 4L);
        if (*m > *k)
        {

          /*                 W := W + C1'*V1 */

          i__1 = *m - *k;
          dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b22, &
                 c[c_offset], ldc, &v[v_offset], ldv, &c_b22, &
                 work[work_offset], ldwork, 9L, 12L);
        }

        /*              W := W * T'  or  W * T */

        dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b22, &t[
                 t_offset], ldt, &work[work_offset], ldwork, 5L, 5L,
               1L, 8L);

        /*              C := C - V * W' */

        if (*m > *k)
        {

          /*                 C1 := C1 - V1 * W' */

          i__1 = *m - *k;
          dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b211,
                 &v[v_offset], ldv, &work[work_offset], ldwork, &
                 c_b22, &c[c_offset], ldc, 12L, 9L);
        }

        /*              W := W * V2' */

        dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &
               v[*m - *k + 1 + v_dim1], ldv, &work[work_offset],
               ldwork, 5L, 5L, 9L, 4L);

        /*              C2 := C2 - W' */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *n;
          for (i = 1; i <= i__2; ++i)
          {
            c[*m - *k + j + i * c_dim1] -= work[i + j * work_dim1]
                                           ;
            /* L80: */
          }
          /* L90: */
        }

      }
      else if (lsame_(side, "R", 1L, 1L))
      {

        /*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
         */

        /*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WOR
        K) */

        /*              W := C2 */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          dcopy_(m, &c[(*n - *k + j) * c_dim1 + 1], &c__1, &work[j *
                 work_dim1 + 1], &c__1);
          /* L100: */
        }

        /*              W := W * V2 */

        dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b22,
               &v[*n - *k + 1 + v_dim1], ldv, &work[work_offset],
               ldwork, 5L, 5L, 12L, 4L);
        if (*n > *k)
        {

          /*                 W := W + C1 * V1 */

          i__1 = *n - *k;
          dgemm_("No transpose", "No transpose", m, k, &i__1, &
                 c_b22, &c[c_offset], ldc, &v[v_offset], ldv, &
                 c_b22, &work[work_offset], ldwork, 12L, 12L);
        }

        /*              W := W * T  or  W * T' */

        dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b22, &t[
                 t_offset], ldt, &work[work_offset], ldwork, 5L, 5L,
               1L, 8L);

        /*              C := C - W * V' */

        if (*n > *k)
        {

          /*                 C1 := C1 - W * V1' */

          i__1 = *n - *k;
          dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b211,
                 &work[work_offset], ldwork, &v[v_offset], ldv, &
                 c_b22, &c[c_offset], ldc, 12L, 9L);
        }

        /*              W := W * V2' */

        dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b22, &
               v[*n - *k + 1 + v_dim1], ldv, &work[work_offset],
               ldwork, 5L, 5L, 9L, 4L);

        /*              C2 := C2 - W */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + (*n - *k + j) * c_dim1] -= work[i + j *
                                                  work_dim1];
            /* L110: */
          }
          /* L120: */
        }
      }
    }

  }
  else if (lsame_(storev, "R", 1L, 1L))
  {

    if (lsame_(direct, "F", 1L, 1L))
    {

      /*           Let  V =  ( V1  V2 )    (V1: first K columns) */
      /*           where  V1  is unit upper triangular. */

      if (lsame_(side, "L", 1L, 1L))
      {

        /*              Form  H * C  or  H' * C  where  C = ( C1 ) */
        /*                                                  ( C2 ) */

        /*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored i
        n WORK) */

        /*              W := C1' */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          dcopy_(n, &c[j + c_dim1], ldc, &work[j * work_dim1 + 1], &
                 c__1);
          /* L130: */
        }

        /*              W := W * V1' */

        dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &
               v[v_offset], ldv, &work[work_offset], ldwork, 5L, 5L,
               9L, 4L);
        if (*m > *k)
        {

          /*                 W := W + C2'*V2' */

          i__1 = *m - *k;
          dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b22, &c[*
                 k + 1 + c_dim1], ldc, &v[(*k + 1) * v_dim1 + 1],
                 ldv, &c_b22, &work[work_offset], ldwork, 9L, 9L);
        }

        /*              W := W * T'  or  W * T */

        dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b22, &t[
                 t_offset], ldt, &work[work_offset], ldwork, 5L, 5L,
               1L, 8L);

        /*              C := C - V' * W' */

        if (*m > *k)
        {

          /*                 C2 := C2 - V2' * W' */

          i__1 = *m - *k;
          dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b211, &v[
                   (*k + 1) * v_dim1 + 1], ldv, &work[work_offset],
                 ldwork, &c_b22, &c[*k + 1 + c_dim1], ldc, 9L, 9L);
        }

        /*              W := W * V1 */

        dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b22,
               &v[v_offset], ldv, &work[work_offset], ldwork, 5L,
               5L, 12L, 4L);

        /*              C1 := C1 - W' */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *n;
          for (i = 1; i <= i__2; ++i)
          {
            c[j + i * c_dim1] -= work[i + j * work_dim1];
            /* L140: */
          }
          /* L150: */
        }

      }
      else if (lsame_(side, "R", 1L, 1L))
      {

        /*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
         */

        /*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in
        WORK) */

        /*              W := C1 */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          dcopy_(m, &c[j * c_dim1 + 1], &c__1, &work[j * work_dim1
                 + 1], &c__1);
          /* L160: */
        }

        /*              W := W * V1' */

        dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b22, &
               v[v_offset], ldv, &work[work_offset], ldwork, 5L, 5L,
               9L, 4L);
        if (*n > *k)
        {

          /*                 W := W + C2 * V2' */

          i__1 = *n - *k;
          dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b22, &
                 c[(*k + 1) * c_dim1 + 1], ldc, &v[(*k + 1) *
                     v_dim1 + 1], ldv, &c_b22, &work[work_offset],
                 ldwork, 12L, 9L);
        }

        /*              W := W * T  or  W * T' */

        dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b22, &t[
                 t_offset], ldt, &work[work_offset], ldwork, 5L, 5L,
               1L, 8L);

        /*              C := C - W * V */

        if (*n > *k)
        {

          /*                 C2 := C2 - W * V2 */

          i__1 = *n - *k;
          dgemm_("No transpose", "No transpose", m, &i__1, k, &
                 c_b211, &work[work_offset], ldwork, &v[(*k + 1) *
                     v_dim1 + 1], ldv, &c_b22, &c[(*k + 1) * c_dim1 +
                                                  1], ldc, 12L, 12L);
        }

        /*              W := W * V1 */

        dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b22,
               &v[v_offset], ldv, &work[work_offset], ldwork, 5L,
               5L, 12L, 4L);

        /*              C1 := C1 - W */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + j * c_dim1] -= work[i + j * work_dim1];
            /* L170: */
          }
          /* L180: */
        }

      }

    }
    else
    {

      /*           Let  V =  ( V1  V2 )    (V2: last K columns) */
      /*           where  V2  is unit lower triangular. */

      if (lsame_(side, "L", 1L, 1L))
      {

        /*              Form  H * C  or  H' * C  where  C = ( C1 ) */
        /*                                                  ( C2 ) */

        /*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored i
        n WORK) */

        /*              W := C2' */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          dcopy_(n, &c[*m - *k + j + c_dim1], ldc, &work[j *
                 work_dim1 + 1], &c__1);
          /* L190: */
        }

        /*              W := W * V2' */

        dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b22, &
               v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
               , ldwork, 5L, 5L, 9L, 4L);
        if (*m > *k)
        {

          /*                 W := W + C1'*V1' */

          i__1 = *m - *k;
          dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b22, &c[
                   c_offset], ldc, &v[v_offset], ldv, &c_b22, &work[
                   work_offset], ldwork, 9L, 9L);
        }

        /*              W := W * T'  or  W * T */

        dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b22, &t[
                 t_offset], ldt, &work[work_offset], ldwork, 5L, 5L,
               1L, 8L);

        /*              C := C - V' * W' */

        if (*m > *k)
        {

          /*                 C1 := C1 - V1' * W' */

          i__1 = *m - *k;
          dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b211, &v[
                   v_offset], ldv, &work[work_offset], ldwork, &
                 c_b22, &c[c_offset], ldc, 9L, 9L);
        }

        /*              W := W * V2 */

        dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b22,
               &v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
                 work_offset], ldwork, 5L, 5L, 12L, 4L);

        /*              C2 := C2 - W' */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *n;
          for (i = 1; i <= i__2; ++i)
          {
            c[*m - *k + j + i * c_dim1] -= work[i + j * work_dim1]
                                           ;
            /* L200: */
          }
          /* L210: */
        }

      }
      else if (lsame_(side, "R", 1L, 1L))
      {

        /*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
         */

        /*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in
        WORK) */

        /*              W := C2 */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          dcopy_(m, &c[(*n - *k + j) * c_dim1 + 1], &c__1, &work[j *
                 work_dim1 + 1], &c__1);
          /* L220: */
        }

        /*              W := W * V2' */

        dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b22, &
               v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
               , ldwork, 5L, 5L, 9L, 4L);
        if (*n > *k)
        {

          /*                 W := W + C1 * V1' */

          i__1 = *n - *k;
          dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b22, &
                 c[c_offset], ldc, &v[v_offset], ldv, &c_b22, &
                 work[work_offset], ldwork, 12L, 9L);
        }

        /*              W := W * T  or  W * T' */

        dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b22, &t[
                 t_offset], ldt, &work[work_offset], ldwork, 5L, 5L,
               1L, 8L);

        /*              C := C - W * V */

        if (*n > *k)
        {

          /*                 C1 := C1 - W * V1 */

          i__1 = *n - *k;
          dgemm_("No transpose", "No transpose", m, &i__1, k, &
                 c_b211, &work[work_offset], ldwork, &v[v_offset],
                 ldv, &c_b22, &c[c_offset], ldc, 12L, 12L);
        }

        /*              W := W * V2 */

        dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b22,
               &v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
                 work_offset], ldwork, 5L, 5L, 12L, 4L);

        /*              C1 := C1 - W */

        i__1 = *k;
        for (j = 1; j <= i__1; ++j)
        {
          i__2 = *m;
          for (i = 1; i <= i__2; ++i)
          {
            c[i + (*n - *k + j) * c_dim1] -= work[i + j *
                                                  work_dim1];
            /* L230: */
          }
          /* L240: */
        }

      }

    }
  }

  return 0;

  /*     End of DLARFB */

} /* dlarfb_ */

/* Subroutine */ int dlarft_(direct, storev, n, k, v, ldv, tau, t, ldt,
                             direct_len, storev_len)
char *direct, *storev;
integer *n, *k;
doublereal *v;
integer *ldv;
doublereal *tau, *t;
integer *ldt;
ftnlen direct_len;
ftnlen storev_len;
{
  /* System generated locals */
  integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
  doublereal d__1;

  /* Local variables */
  static integer i, j;
  extern logical lsame_();
  extern /* Subroutine */ int dgemv_(), dtrmv_();
  static doublereal vii;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     February 29, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLARFT forms the triangular factor T of a real block reflector H */
  /*  of order n, which is defined as a product of k elementary reflectors.
  */

  /*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
  */

  /*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
  */

  /*  If STOREV = 'C', the vector which defines the elementary reflector */
  /*  H(i) is stored in the i-th column of the array V, and */

  /*     H  =  I - V * T * V' */

  /*  If STOREV = 'R', the vector which defines the elementary reflector */
  /*  H(i) is stored in the i-th row of the array V, and */

  /*     H  =  I - V' * T * V */

  /*  Arguments */
  /*  ========= */

  /*  DIRECT  (input) CHARACTER*1 */
  /*          Specifies the order in which the elementary reflectors are */
  /*          multiplied to form the block reflector: */
  /*          = 'F': H = H(1) H(2) . . . H(k) (Forward) */
  /*          = 'B': H = H(k) . . . H(2) H(1) (Backward) */

  /*  STOREV  (input) CHARACTER*1 */
  /*          Specifies how the vectors which define the elementary */
  /*          reflectors are stored (see also Further Details): */
  /*          = 'C': columnwise */
  /*          = 'R': rowwise */

  /*  N       (input) INTEGER */
  /*          The order of the block reflector H. N >= 0. */

  /*  K       (input) INTEGER */
  /*          The order of the triangular factor T (= the number of */
  /*          elementary reflectors). K >= 1. */

  /*  V       (input/output) DOUBLE PRECISION array, dimension */
  /*                               (LDV,K) if STOREV = 'C' */
  /*                               (LDV,N) if STOREV = 'R' */
  /*          The matrix V. See further details. */

  /*  LDV     (input) INTEGER */
  /*          The leading dimension of the array V. */
  /*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
  */

  /*  TAU     (input) DOUBLE PRECISION array, dimension (K) */
  /*          TAU(i) must contain the scalar factor of the elementary */
  /*          reflector H(i). */

  /*  T       (output) DOUBLE PRECISION array, dimension (LDT,K) */
  /*          The k by k triangular factor T of the block reflector. */
  /*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
  */
  /*          lower triangular. The rest of the array is not used. */

  /*  LDT     (input) INTEGER */
  /*          The leading dimension of the array T. LDT >= K. */

  /*  Further Details */
  /*  =============== */

  /*  The shape of the matrix V and the storage of the vectors which define
  */
  /*  the H(i) is best illustrated by the following example with n = 5 and
  */
  /*  k = 3. The elements equal to 1 are not stored; the corresponding */
  /*  array elements are modified but restored on exit. The rest of the */
  /*  array is not used. */

  /*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
  */

  /*               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
  */
  /*                   ( v1  1    )                     (     1 v2 v2 v2 )
  */
  /*                   ( v1 v2  1 )                     (        1 v3 v3 )
  */
  /*                   ( v1 v2 v3 ) */
  /*                   ( v1 v2 v3 ) */

  /*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
  */

  /*               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
  */
  /*                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
  */
  /*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
  */
  /*                   (     1 v3 ) */
  /*                   (        1 ) */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Quick return if possible */

  /* Parameter adjustments */
  t_dim1 = *ldt;
  t_offset = t_dim1 + 1;
  t -= t_offset;
  --tau;
  v_dim1 = *ldv;
  v_offset = v_dim1 + 1;
  v -= v_offset;

  /* Function Body */
  if (*n == 0)
  {
    return 0;
  }

  if (lsame_(direct, "F", 1L, 1L))
  {
    i__1 = *k;
    for (i = 1; i <= i__1; ++i)
    {
      if (tau[i] == 0.)
      {

        /*              H(i)  =  I */

        i__2 = i;
        for (j = 1; j <= i__2; ++j)
        {
          t[j + i * t_dim1] = 0.;
          /* L10: */
        }
      }
      else
      {

        /*              general case */

        vii = v[i + i * v_dim1];
        v[i + i * v_dim1] = 1.;
        if (lsame_(storev, "C", 1L, 1L))
        {

          /*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)'
          * V(i:n,i) */

          i__2 = *n - i + 1;
          i__3 = i - 1;
          d__1 = -tau[i];
          dgemv_("Transpose", &i__2, &i__3, &d__1, &v[i + v_dim1],
                 ldv, &v[i + i * v_dim1], &c__1, &c_b21, &t[i *
                     t_dim1 + 1], &c__1, 9L);
        }
        else
        {

          /*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) *
           V(i,i:n)' */

          i__2 = i - 1;
          i__3 = *n - i + 1;
          d__1 = -tau[i];
          dgemv_("No transpose", &i__2, &i__3, &d__1, &v[i * v_dim1
                 + 1], ldv, &v[i + i * v_dim1], ldv, &c_b21, &t[i *
                     t_dim1 + 1], &c__1, 12L);
        }
        v[i + i * v_dim1] = vii;

        /*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i) */

        i__2 = i - 1;
        dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[
                 t_offset], ldt, &t[i * t_dim1 + 1], &c__1, 5L, 12L,
               8L);
        t[i + i * t_dim1] = tau[i];
      }
      /* L20: */
    }
  }
  else
  {
    for (i = *k; i >= 1; --i)
    {
      if (tau[i] == 0.)
      {

        /*              H(i)  =  I */

        i__1 = *k;
        for (j = i; j <= i__1; ++j)
        {
          t[j + i * t_dim1] = 0.;
          /* L30: */
        }
      }
      else
      {

        /*              general case */

        if (i < *k)
        {
          if (lsame_(storev, "C", 1L, 1L))
          {
            vii = v[*n - *k + i + i * v_dim1];
            v[*n - *k + i + i * v_dim1] = 1.;

            /*                    T(i+1:k,i) := */
            /*                            - tau(i) * V(1:n-k+i,i+1
            :k)' * V(1:n-k+i,i) */

            i__1 = *n - *k + i;
            i__2 = *k - i;
            d__1 = -tau[i];
            dgemv_("Transpose", &i__1, &i__2, &d__1, &v[(i + 1) *
                   v_dim1 + 1], ldv, &v[i * v_dim1 + 1], &c__1, &
                   c_b21, &t[i + 1 + i * t_dim1], &c__1, 9L);
            v[*n - *k + i + i * v_dim1] = vii;
          }
          else
          {
            vii = v[i + (*n - *k + i) * v_dim1];
            v[i + (*n - *k + i) * v_dim1] = 1.;

            /*                    T(i+1:k,i) := */
            /*                            - tau(i) * V(i+1:k,1:n-k
            +i) * V(i,1:n-k+i)' */

            i__1 = *k - i;
            i__2 = *n - *k + i;
            d__1 = -tau[i];
            dgemv_("No transpose", &i__1, &i__2, &d__1, &v[i + 1
                   + v_dim1], ldv, &v[i + v_dim1], ldv, &c_b21, &
                   t[i + 1 + i * t_dim1], &c__1, 12L);
            v[i + (*n - *k + i) * v_dim1] = vii;
          }

          /*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,
          i) */

          i__1 = *k - i;
          dtrmv_("Lower", "No transpose", "Non-unit", &i__1, &t[i +
                 1 + (i + 1) * t_dim1], ldt, &t[i + 1 + i * t_dim1]
                 , &c__1, 5L, 12L, 8L);
        }
        t[i + i * t_dim1] = tau[i];
      }
      /* L40: */
    }
  }
  return 0;

  /*     End of DLARFT */

} /* dlarft_ */

/* Subroutine */ int dorg2l_(m, n, k, a, lda, tau, work, info)
integer *m, *n, *k;
doublereal *a;
integer *lda;
doublereal *tau, *work;
integer *info;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2, i__3;
  doublereal d__1;

  /* Local variables */
  static integer i, j, l;
  extern /* Subroutine */ int dscal_(), dlarf_();
  static integer ii;
  extern /* Subroutine */ int xerbla_();


  /*  -- LAPACK routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     February 29, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DORG2L generates an m by n real matrix Q with orthonormal columns, */
  /*  which is defined as the last n columns of a product of k elementary */
  /*  reflectors of order m */

  /*        Q  =  H(k) . . . H(2) H(1) */

  /*  as returned by DGEQLF. */

  /*  Arguments */
  /*  ========= */

  /*  M       (input) INTEGER */
  /*          The number of rows of the matrix Q. M >= 0. */

  /*  N       (input) INTEGER */
  /*          The number of columns of the matrix Q. M >= N >= 0. */

  /*  K       (input) INTEGER */
  /*          The number of elementary reflectors whose product defines the
  */
  /*          matrix Q. N >= K >= 0. */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On entry, the (n-k+i)-th column must contain the vector which
  */
  /*          defines the elementary reflector H(i), for i = 1,2,...,k, as
  */
  /*          returned by DGEQLF in the last k columns of its array */
  /*          argument A. */
  /*          On exit, the m by n matrix Q. */

  /*  LDA     (input) INTEGER */
  /*          The first dimension of the array A. LDA >= max(1,M). */

  /*  TAU     (input) DOUBLE PRECISION array, dimension (K) */
  /*          TAU(i) must contain the scalar factor of the elementary */
  /*          reflector H(i), as returned by DGEQLF. */

  /*  WORK    (workspace) DOUBLE PRECISION array, dimension (N) */

  /*  INFO    (output) INTEGER */
  /*          = 0: successful exit */
  /*          < 0: if INFO = -i, the i-th argument has an illegal value */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input arguments */

  /* Parameter adjustments */
  --work;
  --tau;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  *info = 0;
  if (*m < 0)
  {
    *info = -1;
  }
  else if (*n < 0 || *n > *m)
  {
    *info = -2;
  }
  else if (*k < 0 || *k > *n)
  {
    *info = -3;
  }
  else if (*lda < max(1, *m))
  {
    *info = -5;
  }
  if (*info != 0)
  {
    i__1 = -(*info);
    xerbla_("DORG2L", &i__1, 6L);
    return 0;
  }

  /*     Quick return if possible */

  if (*n <= 0)
  {
    return 0;
  }

  /*     Initialise columns 1:n-k to columns of the unit matrix */

  i__1 = *n - *k;
  for (j = 1; j <= i__1; ++j)
  {
    i__2 = *m;
    for (l = 1; l <= i__2; ++l)
    {
      a[l + j * a_dim1] = 0.;
      /* L10: */
    }
    a[*m - *n + j + j * a_dim1] = 1.;
    /* L20: */
  }

  i__1 = *k;
  for (i = 1; i <= i__1; ++i)
  {
    ii = *n - *k + i;

    /*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left */

    a[*m - *n + ii + ii * a_dim1] = 1.;
    i__2 = *m - *n + ii;
    i__3 = ii - 1;
    dlarf_("Left", &i__2, &i__3, &a[ii * a_dim1 + 1], &c__1, &tau[i], &a[
             a_offset], lda, &work[1], 4L);
    i__2 = *m - *n + ii - 1;
    d__1 = -tau[i];
    dscal_(&i__2, &d__1, &a[ii * a_dim1 + 1], &c__1);
    a[*m - *n + ii + ii * a_dim1] = 1. - tau[i];

    /*        Set A(m-k+i+1:m,n-k+i) to zero */

    i__2 = *m;
    for (l = *m - *n + ii + 1; l <= i__2; ++l)
    {
      a[l + ii * a_dim1] = 0.;
      /* L30: */
    }
    /* L40: */
  }
  return 0;

  /*     End of DORG2L */

} /* dorg2l_ */

/* Subroutine */ int dlarf_(side, m, n, v, incv, tau, c, ldc, work, side_len)
char *side;
integer *m, *n;
doublereal *v;
integer *incv;
doublereal *tau, *c;
integer *ldc;
doublereal *work;
ftnlen side_len;
{
  /* System generated locals */
  integer c_dim1, c_offset;
  doublereal d__1;

  /* Local variables */
  extern /* Subroutine */ int dger_();
  extern logical lsame_();
  extern /* Subroutine */ int dgemv_();


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     February 29, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLARF applies a real elementary reflector H to a real m by n matrix */
  /*  C, from either the left or the right. H is represented in the form */

  /*        H = I - tau * v * v' */

  /*  where tau is a real scalar and v is a real vector. */

  /*  If tau = 0, then H is taken to be the unit matrix. */

  /*  Arguments */
  /*  ========= */

  /*  SIDE    (input) CHARACTER*1 */
  /*          = 'L': form  H * C */
  /*          = 'R': form  C * H */

  /*  M       (input) INTEGER */
  /*          The number of rows of the matrix C. */

  /*  N       (input) INTEGER */
  /*          The number of columns of the matrix C. */

  /*  V       (input) DOUBLE PRECISION array, dimension */
  /*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L' */
  /*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R' */
  /*          The vector v in the representation of H. V is not used if */
  /*          TAU = 0. */

  /*  INCV    (input) INTEGER */
  /*          The increment between elements of v. INCV <> 0. */

  /*  TAU     (input) DOUBLE PRECISION */
  /*          The value tau in the representation of H. */

  /*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
  /*          On entry, the m by n matrix C. */
  /*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
  */
  /*          or C * H if SIDE = 'R'. */

  /*  LDC     (input) INTEGER */
  /*          The leading dimension of the array C. LDC >= max(1,M). */

  /*  WORK    (workspace) DOUBLE PRECISION array, dimension */
  /*                         (N) if SIDE = 'L' */
  /*                      or (M) if SIDE = 'R' */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --work;
  c_dim1 = *ldc;
  c_offset = c_dim1 + 1;
  c -= c_offset;
  --v;

  /* Function Body */
  if (lsame_(side, "L", 1L, 1L))
  {

    /*        Form  H * C */

    if (*tau != 0.)
    {

      /*           w := C' * v */

      dgemv_("Transpose", m, n, &c_b22, &c[c_offset], ldc, &v[1], incv,
             &c_b21, &work[1], &c__1, 9L);

      /*           C := C - v * w' */

      d__1 = -(*tau);
      dger_(m, n, &d__1, &v[1], incv, &work[1], &c__1, &c[c_offset],
            ldc);
    }
  }
  else
  {

    /*        Form  C * H */

    if (*tau != 0.)
    {

      /*           w := C * v */

      dgemv_("No transpose", m, n, &c_b22, &c[c_offset], ldc, &v[1],
             incv, &c_b21, &work[1], &c__1, 12L);

      /*           C := C - w * v' */

      d__1 = -(*tau);
      dger_(m, n, &d__1, &work[1], &c__1, &v[1], incv, &c[c_offset],
            ldc);
    }
  }
  return 0;

  /*     End of DLARF */

} /* dlarf_ */

/* Subroutine */ int dsterf_(n, d, e, info)
integer *n;
doublereal *d, *e;
integer *info;
{
  /* System generated locals */
  integer i__1, i__2;
  doublereal d__1, d__2;

  /* Builtin functions */
  double sqrt(), d_sign();

  /* Local variables */
  static doublereal oldc;
  static integer lend, jtot;
  extern /* Subroutine */ int dlae2_();
  static doublereal c;
  static integer i, j, k, l, m;
  static doublereal p, gamma, r, s, alpha, sigma;
  static integer l1, lendm1, lendp1;
  extern doublereal dlapy2_();
  static doublereal bb;
  static integer ii;
  extern doublereal dlamch_();
  static doublereal oldgam;
  extern /* Subroutine */ int xerbla_();
  static integer nmaxit, lm1, mm1, nm1;
  static doublereal rt1, rt2, eps, rte, tst;


  /*  -- LAPACK routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     March 31, 1993 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix */
  /*  using the Pal-Walker-Kahan variant of the QL or QR algorithm. */

  /*  Arguments */
  /*  ========= */

  /*  N       (input) INTEGER */
  /*          The order of the matrix.  N >= 0. */

  /*  D       (input/output) DOUBLE PRECISION array, dimension (N) */
  /*          On entry, the n diagonal elements of the tridiagonal matrix.
  */
  /*          On exit, if INFO = 0, the eigenvalues in ascending order. */

  /*  E       (input/output) DOUBLE PRECISION array, dimension (N-1) */
  /*          On entry, the (n-1) subdiagonal elements of the tridiagonal */
  /*          matrix. */
  /*          On exit, E has been destroyed. */

  /*  INFO    (output) INTEGER */
  /*          = 0:  successful exit */
  /*          < 0:  if INFO = -i, the i-th argument had an illegal value */
  /*          > 0:  the algorithm failed to find all of the eigenvalues in
  */
  /*                a total of 30*N iterations; if INFO = i, then i */
  /*                elements of E have not converged to zero. */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters. */

  /* Parameter adjustments */
  --e;
  --d;

  /* Function Body */
  *info = 0;

  /*     Quick return if possible */

  if (*n < 0)
  {
    *info = -1;
    i__1 = -(*info);
    xerbla_("DSTERF", &i__1, 6L);
    return 0;
  }
  if (*n <= 1)
  {
    return 0;
  }

  /*     Determine the unit roundoff for this environment. */

  eps = dlamch_("E", 1L);

  /*     Compute the eigenvalues of the tridiagonal matrix. */

  i__1 = *n - 1;
  for (i = 1; i <= i__1; ++i)
  {
    /* Computing 2nd power */
    d__1 = e[i];
    e[i] = d__1 * d__1;
    /* L10: */
  }

  nmaxit = *n * 30;
  sigma = 0.;
  jtot = 0;

  /*     Determine where the matrix splits and choose QL or QR iteration */
  /*     for each block, according to whether top or bottom diagonal */
  /*     element is smaller. */

  l1 = 1;
  nm1 = *n - 1;

L20:
  if (l1 > *n)
  {
    goto L170;
  }
  if (l1 > 1)
  {
    e[l1 - 1] = 0.;
  }
  if (l1 <= nm1)
  {
    i__1 = nm1;
    for (m = l1; m <= i__1; ++m)
    {
      tst = sqrt((d__1 = e[m], abs(d__1)));
      if (tst <= eps * ((d__1 = d[m], abs(d__1)) + (d__2 = d[m + 1],
                        abs(d__2))))
      {
        goto L40;
      }
      /* L30: */
    }
  }
  m = *n;

L40:
  l = l1;
  lend = m;
  if ((d__1 = d[lend], abs(d__1)) < (d__2 = d[l], abs(d__2)))
  {
    l = lend;
    lend = l1;
  }
  l1 = m + 1;

  if (lend >= l)
  {

    /*        QL Iteration */

    /*        Look for small subdiagonal element. */

L50:
    if (l != lend)
    {
      lendm1 = lend - 1;
      i__1 = lendm1;
      for (m = l; m <= i__1; ++m)
      {
        tst = sqrt((d__1 = e[m], abs(d__1)));
        if (tst <= eps * ((d__1 = d[m], abs(d__1)) + (d__2 = d[m + 1],
                          abs(d__2))))
        {
          goto L70;
        }
        /* L60: */
      }
    }

    m = lend;

L70:
    if (m < lend)
    {
      e[m] = 0.;
    }
    p = d[l];
    if (m == l)
    {
      goto L90;
    }

    /*        If remaining matrix is 2 by 2, use DLAE2 to compute its */
    /*        eigenvalues. */

    if (m == l + 1)
    {
      rte = sqrt(e[l]);
      dlae2_(&d[l], &rte, &d[l + 1], &rt1, &rt2);
      d[l] = rt1;
      d[l + 1] = rt2;
      e[l] = 0.;
      l += 2;
      if (l <= lend)
      {
        goto L50;
      }
      goto L20;
    }

    if (jtot == nmaxit)
    {
      goto L150;
    }
    ++jtot;

    /*        Form shift. */

    rte = sqrt(e[l]);
    sigma = (d[l + 1] - p) / (rte * 2.);
    r = dlapy2_(&sigma, &c_b22);
    sigma = p - rte / (sigma + d_sign(&r, &sigma));

    c = 1.;
    s = 0.;
    gamma = d[m] - sigma;
    p = gamma * gamma;

    /*        Inner loop */

    mm1 = m - 1;
    i__1 = l;
    for (i = mm1; i >= i__1; --i)
    {
      bb = e[i];
      r = p + bb;
      if (i != m - 1)
      {
        e[i + 1] = s * r;
      }
      oldc = c;
      c = p / r;
      s = bb / r;
      oldgam = gamma;
      alpha = d[i];
      gamma = c * (alpha - sigma) - s * oldgam;
      d[i + 1] = oldgam + (alpha - gamma);
      if (c != 0.)
      {
        p = gamma * gamma / c;
      }
      else
      {
        p = oldc * bb;
      }
      /* L80: */
    }

    e[l] = s * p;
    d[l] = sigma + gamma;
    goto L50;

    /*        Eigenvalue found. */

L90:
    d[l] = p;

    ++l;
    if (l <= lend)
    {
      goto L50;
    }
    goto L20;

  }
  else
  {

    /*        QR Iteration */

    /*        Look for small superdiagonal element. */

L100:
    if (l != lend)
    {
      lendp1 = lend + 1;
      i__1 = lendp1;
      for (m = l; m >= i__1; --m)
      {
        tst = sqrt((d__1 = e[m - 1], abs(d__1)));
        if (tst <= eps * ((d__1 = d[m], abs(d__1)) + (d__2 = d[m - 1],
                          abs(d__2))))
        {
          goto L120;
        }
        /* L110: */
      }
    }

    m = lend;

L120:
    if (m > lend)
    {
      e[m - 1] = 0.;
    }
    p = d[l];
    if (m == l)
    {
      goto L140;
    }

    /*        If remaining matrix is 2 by 2, use DLAE2 to compute its */
    /*        eigenvalues. */

    if (m == l - 1)
    {
      rte = sqrt(e[l - 1]);
      dlae2_(&d[l], &rte, &d[l - 1], &rt1, &rt2);
      d[l] = rt1;
      d[l - 1] = rt2;
      e[l - 1] = 0.;
      l += -2;
      if (l >= lend)
      {
        goto L100;
      }
      goto L20;
    }

    if (jtot == nmaxit)
    {
      goto L150;
    }
    ++jtot;

    /*        Form shift. */

    rte = sqrt(e[l - 1]);
    sigma = (d[l - 1] - p) / (rte * 2.);
    r = dlapy2_(&sigma, &c_b22);
    sigma = p - rte / (sigma + d_sign(&r, &sigma));

    c = 1.;
    s = 0.;
    gamma = d[m] - sigma;
    p = gamma * gamma;

    /*        Inner loop */

    lm1 = l - 1;
    i__1 = lm1;
    for (i = m; i <= i__1; ++i)
    {
      bb = e[i];
      r = p + bb;
      if (i != m)
      {
        e[i - 1] = s * r;
      }
      oldc = c;
      c = p / r;
      s = bb / r;
      oldgam = gamma;
      alpha = d[i + 1];
      gamma = c * (alpha - sigma) - s * oldgam;
      d[i] = oldgam + (alpha - gamma);
      if (c != 0.)
      {
        p = gamma * gamma / c;
      }
      else
      {
        p = oldc * bb;
      }
      /* L130: */
    }

    e[lm1] = s * p;
    d[l] = sigma + gamma;
    goto L100;

    /*        Eigenvalue found. */

L140:
    d[l] = p;

    --l;
    if (l >= lend)
    {
      goto L100;
    }
    goto L20;

  }

  /*     Set error -- no convergence to an eigenvalue after a total */
  /*     of N*MAXIT iterations. */

L150:
  i__1 = *n - 1;
  for (i = 1; i <= i__1; ++i)
  {
    if (e[i] != 0.)
    {
      ++(*info);
    }
    /* L160: */
  }
  return 0;

  /*     Sort eigenvalues in increasing order. */

L170:
  i__1 = *n;
  for (ii = 2; ii <= i__1; ++ii)
  {
    i = ii - 1;
    k = i;
    p = d[i];
    i__2 = *n;
    for (j = ii; j <= i__2; ++j)
    {
      if (d[j] < p)
      {
        k = j;
        p = d[j];
      }
      /* L180: */
    }
    if (k != i)
    {
      d[k] = d[i];
      d[i] = p;
    }
    /* L190: */
  }

  return 0;

  /*     End of DSTERF */

} /* dsterf_ */

/* Subroutine */ int dlae2_(a, b, c, rt1, rt2)
doublereal *a, *b, *c, *rt1, *rt2;
{
  /* System generated locals */
  doublereal d__1;

  /* Builtin functions */
  double sqrt();

  /* Local variables */
  static doublereal acmn, acmx, ab, df, tb, sm, rt, adf;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix */
  /*     [  A   B  ] */
  /*     [  B   C  ]. */
  /*  On return, RT1 is the eigenvalue of larger absolute value, and RT2 */
  /*  is the eigenvalue of smaller absolute value. */

  /*  Arguments */
  /*  ========= */

  /*  A       (input) DOUBLE PRECISION */
  /*          The (1,1) entry of the 2-by-2 matrix. */

  /*  B       (input) DOUBLE PRECISION */
  /*          The (1,2) and (2,1) entries of the 2-by-2 matrix. */

  /*  C       (input) DOUBLE PRECISION */
  /*          The (2,2) entry of the 2-by-2 matrix. */

  /*  RT1     (output) DOUBLE PRECISION */
  /*          The eigenvalue of larger absolute value. */

  /*  RT2     (output) DOUBLE PRECISION */
  /*          The eigenvalue of smaller absolute value. */

  /*  Further Details */
  /*  =============== */

  /*  RT1 is accurate to a few ulps barring over/underflow. */

  /*  RT2 may be inaccurate if there is massive cancellation in the */
  /*  determinant A*C-B*B; higher precision or correctly rounded or */
  /*  correctly truncated arithmetic would be needed to compute RT2 */
  /*  accurately in all cases. */

  /*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
  */
  /*  Underflow is harmless if the input data is 0 or exceeds */
  /*     underflow_threshold / macheps. */

  /* =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Compute the eigenvalues */

  sm = *a + *c;
  df = *a - *c;
  adf = abs(df);
  tb = *b + *b;
  ab = abs(tb);
  if (abs(*a) > abs(*c))
  {
    acmx = *a;
    acmn = *c;
  }
  else
  {
    acmx = *c;
    acmn = *a;
  }
  if (adf > ab)
  {
    /* Computing 2nd power */
    d__1 = ab / adf;
    rt = adf * sqrt(d__1 * d__1 + 1.);
  }
  else if (adf < ab)
  {
    /* Computing 2nd power */
    d__1 = adf / ab;
    rt = ab * sqrt(d__1 * d__1 + 1.);
  }
  else
  {

    /*        Includes case AB=ADF=0 */

    rt = ab * sqrt(2.);
  }
  if (sm < 0.)
  {
    *rt1 = (sm - rt) * .5;

    /*        Order of execution important. */
    /*        To get fully accurate smaller eigenvalue, */
    /*        next line needs to be executed in higher precision. */

    *rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
  }
  else if (sm > 0.)
  {
    *rt1 = (sm + rt) * .5;

    /*        Order of execution important. */
    /*        To get fully accurate smaller eigenvalue, */
    /*        next line needs to be executed in higher precision. */

    *rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
  }
  else
  {

    /*        Includes case RT1 = RT2 = 0 */

    *rt1 = rt * .5;
    *rt2 = rt * -.5;
  }
  return 0;

  /*     End of DLAE2 */

} /* dlae2_ */

/* Subroutine */ int dsytrd_(uplo, n, a, lda, d, e, tau, work, lwork, info,
                             uplo_len)
char *uplo;
integer *n;
doublereal *a;
integer *lda;
doublereal *d, *e, *tau, *work;
integer *lwork, *info;
ftnlen uplo_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2, i__3;

  /* Local variables */
  static integer i, j;
  extern logical lsame_();
  static integer nbmin, iinfo;
  static logical upper;
  extern /* Subroutine */ int dsytd2_(), dsyr2k_();
  static integer nb, kk, nx;
  extern /* Subroutine */ int dlatrd_(), xerbla_();
  extern integer ilaenv_();
  static integer ldwork, iws;


  /*  -- LAPACK routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     March 31, 1993 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DSYTRD reduces a real symmetric matrix A to real symmetric */
  /*  tridiagonal form T by an orthogonal similarity transformation: */
  /*  Q**T * A * Q = T. */

  /*  Arguments */
  /*  ========= */

  /*  UPLO    (input) CHARACTER*1 */
  /*          = 'U':  Upper triangle of A is stored; */
  /*          = 'L':  Lower triangle of A is stored. */

  /*  N       (input) INTEGER */
  /*          The order of the matrix A.  N >= 0. */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  */
  /*          N-by-N upper triangular part of A contains the upper */
  /*          triangular part of the matrix A, and the strictly lower */
  /*          triangular part of A is not referenced.  If UPLO = 'L', the */
  /*          leading N-by-N lower triangular part of A contains the lower
  */
  /*          triangular part of the matrix A, and the strictly upper */
  /*          triangular part of A is not referenced. */
  /*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
  */
  /*          of A are overwritten by the corresponding elements of the */
  /*          tridiagonal matrix T, and the elements above the first */
  /*          superdiagonal, with the array TAU, represent the orthogonal */
  /*          matrix Q as a product of elementary reflectors; if UPLO */
  /*          = 'L', the diagonal and first subdiagonal of A are over- */
  /*          written by the corresponding elements of the tridiagonal */
  /*          matrix T, and the elements below the first subdiagonal, with
  */
  /*          the array TAU, represent the orthogonal matrix Q as a product
  */
  /*          of elementary reflectors. See Further Details. */

  /*  LDA     (input) INTEGER */
  /*          The leading dimension of the array A.  LDA >= max(1,N). */

  /*  D       (output) DOUBLE PRECISION array, dimension (N) */
  /*          The diagonal elements of the tridiagonal matrix T: */
  /*          D(i) = A(i,i). */

  /*  E       (output) DOUBLE PRECISION array, dimension (N-1) */
  /*          The off-diagonal elements of the tridiagonal matrix T: */
  /*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
  */

  /*  TAU     (output) DOUBLE PRECISION array, dimension (N-1) */
  /*          The scalar factors of the elementary reflectors (see Further
  */
  /*          Details). */

  /*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK) */
  /*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

  /*  LWORK   (input) INTEGER */
  /*          The dimension of the array WORK.  LWORK >= 1. */
  /*          For optimum performance LWORK >= N*NB, where NB is the */
  /*          optimal blocksize. */

  /*  INFO    (output) INTEGER */
  /*          = 0:  successful exit */
  /*          < 0:  if INFO = -i, the i-th argument had an illegal value */

  /*  Further Details */
  /*  =============== */

  /*  If UPLO = 'U', the matrix Q is represented as a product of elementary
  */
  /*  reflectors */

  /*     Q = H(n-1) . . . H(2) H(1). */

  /*  Each H(i) has the form */

  /*     H(i) = I - tau * v * v' */

  /*  where tau is a real scalar, and v is a real vector with */
  /*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in */
  /*  A(1:i-1,i+1), and tau in TAU(i). */

  /*  If UPLO = 'L', the matrix Q is represented as a product of elementary
  */
  /*  reflectors */

  /*     Q = H(1) H(2) . . . H(n-1). */

  /*  Each H(i) has the form */

  /*     H(i) = I - tau * v * v' */

  /*  where tau is a real scalar, and v is a real vector with */
  /*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
  */
  /*  and tau in TAU(i). */

  /*  The contents of A on exit are illustrated by the following examples */
  /*  with n = 5: */

  /*  if UPLO = 'U':                       if UPLO = 'L': */

  /*    (  d   e   v2  v3  v4 )              (  d                  ) */
  /*    (      d   e   v3  v4 )              (  e   d              ) */
  /*    (          d   e   v4 )              (  v1  e   d          ) */
  /*    (              d   e  )              (  v1  v2  e   d      ) */
  /*    (                  d  )              (  v1  v2  v3  e   d  ) */

  /*  where d and e denote diagonal and off-diagonal elements of T, and vi
  */
  /*  denotes an element of the vector defining H(i). */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters */

  /* Parameter adjustments */
  --work;
  --tau;
  --e;
  --d;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  *info = 0;
  upper = lsame_(uplo, "U", 1L, 1L);
  if (! upper && ! lsame_(uplo, "L", 1L, 1L))
  {
    *info = -1;
  }
  else if (*n < 0)
  {
    *info = -2;
  }
  else if (*lda < max(1, *n))
  {
    *info = -4;
  }
  else if (*lwork < 1)
  {
    *info = -9;
  }
  if (*info != 0)
  {
    i__1 = -(*info);
    xerbla_("DSYTRD", &i__1, 6L);
    return 0;
  }

  /*     Quick return if possible */

  if (*n == 0)
  {
    work[1] = 1.;
    return 0;
  }

  /*     Determine the block size. */

  nb = ilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, 6L, 1L);
  nx = *n;
  iws = 1;
  if (nb > 1 && nb < *n)
  {

    /*        Determine when to cross over from blocked to unblocked code
    */
    /*        (last block is always handled by unblocked code). */

    /* Computing MAX */
    i__1 = nb, i__2 = ilaenv_(&c__3, "DSYTRD", uplo, n, &c_n1, &c_n1, &
                              c_n1, 6L, 1L);
    nx = max(i__1, i__2);
    if (nx < *n)
    {

      /*           Determine if workspace is large enough for blocked co
      de. */

      ldwork = *n;
      iws = ldwork * nb;
      if (*lwork < iws)
      {

        /*              Not enough workspace to use optimal NB:  deter
        mine the */
        /*              minimum value of NB, and reduce NB or force us
        e of */
        /*              unblocked code by setting NX = N. */

        nb = *lwork / ldwork;
        nbmin = ilaenv_(&c__2, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1,
                        6L, 1L);
        if (nb < nbmin)
        {
          nx = *n;
        }
      }
    }
    else
    {
      nx = *n;
    }
  }
  else
  {
    nb = 1;
  }

  if (upper)
  {

    /*        Reduce the upper triangle of A. */
    /*        Columns 1:kk are handled by the unblocked method. */

    kk = *n - (*n - nx + nb - 1) / nb * nb;
    i__1 = kk + 1;
    i__2 = -nb;
    for (i = *n - nb + 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2)
    {

      /*           Reduce columns i:i+nb-1 to tridiagonal form and form
      the */
      /*           matrix W which is needed to update the unreduced part
       of */
      /*           the matrix */

      i__3 = i + nb - 1;
      dlatrd_(uplo, &i__3, &nb, &a[a_offset], lda, &e[1], &tau[1], &
              work[1], &ldwork, 1L);

      /*           Update the unreduced submatrix A(1:i-1,1:i-1), using
      an */
      /*           update of the form:  A := A - V*W' - W*V' */

      i__3 = i - 1;
      dsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b211, &a[i * a_dim1
              + 1], lda, &work[1], &ldwork, &c_b22, &a[a_offset], lda,
              1L, 12L);

      /*           Copy superdiagonal elements back into A, and diagonal
       */
      /*           elements into D */

      i__3 = i + nb - 1;
      for (j = i; j <= i__3; ++j)
      {
        a[j - 1 + j * a_dim1] = e[j - 1];
        d[j] = a[j + j * a_dim1];
        /* L10: */
      }
      /* L20: */
    }

    /*        Use unblocked code to reduce the last or only block */

    dsytd2_(uplo, &kk, &a[a_offset], lda, &d[1], &e[1], &tau[1], &iinfo,
            1L);
  }
  else
  {

    /*        Reduce the lower triangle of A */

    i__2 = *n - nx;
    i__1 = nb;
    for (i = 1; i__1 < 0 ? i >= i__2 : i <= i__2; i += i__1)
    {

      /*           Reduce columns i:i+nb-1 to tridiagonal form and form
      the */
      /*           matrix W which is needed to update the unreduced part
       of */
      /*           the matrix */

      i__3 = *n - i + 1;
      dlatrd_(uplo, &i__3, &nb, &a[i + i * a_dim1], lda, &e[i], &tau[i],
              &work[1], &ldwork, 1L);

      /*           Update the unreduced submatrix A(i+ib:n,i+ib:n), usin
      g */
      /*           an update of the form:  A := A - V*W' - W*V' */

      i__3 = *n - i - nb + 1;
      dsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b211, &a[i + nb + i *
              a_dim1], lda, &work[nb + 1], &ldwork, &c_b22, &a[i + nb
                  + (i + nb) * a_dim1], lda, 1L, 12L);

      /*           Copy subdiagonal elements back into A, and diagonal
      */
      /*           elements into D */

      i__3 = i + nb - 1;
      for (j = i; j <= i__3; ++j)
      {
        a[j + 1 + j * a_dim1] = e[j];
        d[j] = a[j + j * a_dim1];
        /* L30: */
      }
      /* L40: */
    }

    /*        Use unblocked code to reduce the last or only block */

    i__1 = *n - i + 1;
    dsytd2_(uplo, &i__1, &a[i + i * a_dim1], lda, &d[i], &e[i], &tau[i], &
            iinfo, 1L);
  }

  work[1] = (doublereal) iws;
  return 0;

  /*     End of DSYTRD */

} /* dsytrd_ */

/* Subroutine */ int dsytd2_(uplo, n, a, lda, d, e, tau, info, uplo_len)
char *uplo;
integer *n;
doublereal *a;
integer *lda;
doublereal *d, *e, *tau;
integer *info;
ftnlen uplo_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2, i__3;

  /* Local variables */
  extern doublereal ddot_();
  static doublereal taui;
  extern /* Subroutine */ int dsyr2_();
  static integer i;
  static doublereal alpha;
  extern logical lsame_();
  extern /* Subroutine */ int daxpy_();
  static logical upper;
  extern /* Subroutine */ int dsymv_(), dlarfg_(), xerbla_();


  /*  -- LAPACK routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal */
  /*  form T by an orthogonal similarity transformation: Q' * A * Q = T. */

  /*  Arguments */
  /*  ========= */

  /*  UPLO    (input) CHARACTER*1 */
  /*          Specifies whether the upper or lower triangular part of the */
  /*          symmetric matrix A is stored: */
  /*          = 'U':  Upper triangular */
  /*          = 'L':  Lower triangular */

  /*  N       (input) INTEGER */
  /*          The order of the matrix A.  N >= 0. */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  */
  /*          n-by-n upper triangular part of A contains the upper */
  /*          triangular part of the matrix A, and the strictly lower */
  /*          triangular part of A is not referenced.  If UPLO = 'L', the */
  /*          leading n-by-n lower triangular part of A contains the lower
  */
  /*          triangular part of the matrix A, and the strictly upper */
  /*          triangular part of A is not referenced. */
  /*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
  */
  /*          of A are overwritten by the corresponding elements of the */
  /*          tridiagonal matrix T, and the elements above the first */
  /*          superdiagonal, with the array TAU, represent the orthogonal */
  /*          matrix Q as a product of elementary reflectors; if UPLO */
  /*          = 'L', the diagonal and first subdiagonal of A are over- */
  /*          written by the corresponding elements of the tridiagonal */
  /*          matrix T, and the elements below the first subdiagonal, with
  */
  /*          the array TAU, represent the orthogonal matrix Q as a product
  */
  /*          of elementary reflectors. See Further Details. */

  /*  LDA     (input) INTEGER */
  /*          The leading dimension of the array A.  LDA >= max(1,N). */

  /*  D       (output) DOUBLE PRECISION array, dimension (N) */
  /*          The diagonal elements of the tridiagonal matrix T: */
  /*          D(i) = A(i,i). */

  /*  E       (output) DOUBLE PRECISION array, dimension (N-1) */
  /*          The off-diagonal elements of the tridiagonal matrix T: */
  /*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
  */

  /*  TAU     (output) DOUBLE PRECISION array, dimension (N-1) */
  /*          The scalar factors of the elementary reflectors (see Further
  */
  /*          Details). */

  /*  INFO    (output) INTEGER */
  /*          = 0:  successful exit */
  /*          < 0:  if INFO = -i, the i-th argument had an illegal value. */

  /*  Further Details */
  /*  =============== */

  /*  If UPLO = 'U', the matrix Q is represented as a product of elementary
  */
  /*  reflectors */

  /*     Q = H(n-1) . . . H(2) H(1). */

  /*  Each H(i) has the form */

  /*     H(i) = I - tau * v * v' */

  /*  where tau is a real scalar, and v is a real vector with */
  /*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in */
  /*  A(1:i-1,i+1), and tau in TAU(i). */

  /*  If UPLO = 'L', the matrix Q is represented as a product of elementary
  */
  /*  reflectors */

  /*     Q = H(1) H(2) . . . H(n-1). */

  /*  Each H(i) has the form */

  /*     H(i) = I - tau * v * v' */

  /*  where tau is a real scalar, and v is a real vector with */
  /*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
  */
  /*  and tau in TAU(i). */

  /*  The contents of A on exit are illustrated by the following examples */
  /*  with n = 5: */

  /*  if UPLO = 'U':                       if UPLO = 'L': */

  /*    (  d   e   v2  v3  v4 )              (  d                  ) */
  /*    (      d   e   v3  v4 )              (  e   d              ) */
  /*    (          d   e   v4 )              (  v1  e   d          ) */
  /*    (              d   e  )              (  v1  v2  e   d      ) */
  /*    (                  d  )              (  v1  v2  v3  e   d  ) */

  /*  where d and e denote diagonal and off-diagonal elements of T, and vi
  */
  /*  denotes an element of the vector defining H(i). */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test the input parameters */

  /* Parameter adjustments */
  --tau;
  --e;
  --d;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  *info = 0;
  upper = lsame_(uplo, "U", 1L, 1L);
  if (! upper && ! lsame_(uplo, "L", 1L, 1L))
  {
    *info = -1;
  }
  else if (*n < 0)
  {
    *info = -2;
  }
  else if (*lda < max(1, *n))
  {
    *info = -4;
  }
  if (*info != 0)
  {
    i__1 = -(*info);
    xerbla_("DSYTD2", &i__1, 6L);
    return 0;
  }

  /*     Quick return if possible */

  if (*n <= 0)
  {
    return 0;
  }

  if (upper)
  {

    /*        Reduce the upper triangle of A */

    for (i = *n - 1; i >= 1; --i)
    {

      /*           Generate elementary reflector H(i) = I - tau * v * v'
       */
      /*           to annihilate A(1:i-1,i+1) */

      dlarfg_(&i, &a[i + (i + 1) * a_dim1], &a[(i + 1) * a_dim1 + 1], &
              c__1, &taui);
      e[i] = a[i + (i + 1) * a_dim1];

      if (taui != 0.)
      {

        /*              Apply H(i) from both sides to A(1:i,1:i) */

        a[i + (i + 1) * a_dim1] = 1.;

        /*              Compute  x := tau * A * v  storing x in TAU(1:
        i) */

        dsymv_(uplo, &i, &taui, &a[a_offset], lda, &a[(i + 1) *
               a_dim1 + 1], &c__1, &c_b21, &tau[1], &c__1, 1L);

        /*              Compute  w := x - 1/2 * tau * (x'*v) * v */

        alpha = taui * -.5 * ddot_(&i, &tau[1], &c__1, &a[(i + 1) *
                                   a_dim1 + 1], &c__1);
        daxpy_(&i, &alpha, &a[(i + 1) * a_dim1 + 1], &c__1, &tau[1], &
               c__1);

        /*              Apply the transformation as a rank-2 update:
        */
        /*                 A := A - v * w' - w * v' */

        dsyr2_(uplo, &i, &c_b211, &a[(i + 1) * a_dim1 + 1], &c__1, &
               tau[1], &c__1, &a[a_offset], lda, 1L);

        a[i + (i + 1) * a_dim1] = e[i];
      }
      d[i + 1] = a[i + 1 + (i + 1) * a_dim1];
      tau[i] = taui;
      /* L10: */
    }
    d[1] = a[a_dim1 + 1];
  }
  else
  {

    /*        Reduce the lower triangle of A */

    i__1 = *n - 1;
    for (i = 1; i <= i__1; ++i)
    {

      /*           Generate elementary reflector H(i) = I - tau * v * v'
       */
      /*           to annihilate A(i+2:n,i) */

      i__2 = *n - i;
      /* Computing MIN */
      i__3 = i + 2;
      dlarfg_(&i__2, &a[i + 1 + i * a_dim1], &a[min(i__3, *n) + i *
              a_dim1], &c__1, &taui);
      e[i] = a[i + 1 + i * a_dim1];

      if (taui != 0.)
      {

        /*              Apply H(i) from both sides to A(i+1:n,i+1:n)
        */

        a[i + 1 + i * a_dim1] = 1.;

        /*              Compute  x := tau * A * v  storing y in TAU(i:
        n-1) */

        i__2 = *n - i;
        dsymv_(uplo, &i__2, &taui, &a[i + 1 + (i + 1) * a_dim1], lda,
               &a[i + 1 + i * a_dim1], &c__1, &c_b21, &tau[i], &c__1,
               1L);

        /*              Compute  w := x - 1/2 * tau * (x'*v) * v */

        i__2 = *n - i;
        alpha = taui * -.5 * ddot_(&i__2, &tau[i], &c__1, &a[i + 1 +
                                   i * a_dim1], &c__1);
        i__2 = *n - i;
        daxpy_(&i__2, &alpha, &a[i + 1 + i * a_dim1], &c__1, &tau[i],
               &c__1);

        /*              Apply the transformation as a rank-2 update:
        */
        /*                 A := A - v * w' - w * v' */

        i__2 = *n - i;
        dsyr2_(uplo, &i__2, &c_b211, &a[i + 1 + i * a_dim1], &c__1, &
               tau[i], &c__1, &a[i + 1 + (i + 1) * a_dim1], lda, 1L);

        a[i + 1 + i * a_dim1] = e[i];
      }
      d[i] = a[i + i * a_dim1];
      tau[i] = taui;
      /* L20: */
    }
    d[*n] = a[*n + *n * a_dim1];
  }

  return 0;

  /*     End of DSYTD2 */

} /* dsytd2_ */

/* Subroutine */ int xerbla_(srname, info, srname_len)
char *srname;
integer *info;
ftnlen srname_len;
{
  /* Format strings */
  static char fmt_9999[] = "(\002 ** On entry to \002,a6,\002 parameter nu\
mber \002,i2,\002 had \002,\002an illegal value\002)";

  /* Builtin functions */
  integer s_wsfe(), do_fio(), e_wsfe();
  /* Subroutine */
  int s_stop();

  /* Fortran I/O blocks */
  static cilist io___164 = { 0, 6, 0, fmt_9999, 0 };



  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     February 29, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  XERBLA  is an error handler for the LAPACK routines. */
  /*  It is called by an LAPACK routine if an input parameter has an */
  /*  invalid value.  A message is printed and execution stops. */

  /*  Installers may consider modifying the STOP statement in order to */
  /*  call system-specific exception-handling facilities. */

  /*  Arguments */
  /*  ========= */

  /*  SRNAME  (input) CHARACTER*6 */
  /*          The name of the routine which called XERBLA. */

  /*  INFO    (input) INTEGER */
  /*          The position of the invalid parameter in the parameter list */
  /*          of the calling routine. */

  /*     .. Executable Statements .. */

  s_wsfe(&io___164);
  do_fio(&c__1, srname, 6L);
  do_fio(&c__1, (char *) & (*info), (ftnlen)sizeof(integer));
  e_wsfe();

  s_stop("", 0L);


  /*     End of XERBLA */

} /* xerbla_ */

/* Subroutine */ int dlatrd_(uplo, n, nb, a, lda, e, tau, w, ldw, uplo_len)
char *uplo;
integer *n, *nb;
doublereal *a;
integer *lda;
doublereal *e, *tau, *w;
integer *ldw;
ftnlen uplo_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3;

  /* Local variables */
  extern doublereal ddot_();
  static integer i;
  static doublereal alpha;
  extern /* Subroutine */ int dscal_();
  extern logical lsame_();
  extern /* Subroutine */ int dgemv_(), daxpy_(), dsymv_(), dlarfg_();
  static integer iw;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLATRD reduces NB rows and columns of a real symmetric matrix A to */
  /*  symmetric tridiagonal form by an orthogonal similarity */
  /*  transformation Q' * A * Q, and returns the matrices V and W which are
  */
  /*  needed to apply the transformation to the unreduced part of A. */

  /*  If UPLO = 'U', DLATRD reduces the last NB rows and columns of a */
  /*  matrix, of which the upper triangle is supplied; */
  /*  if UPLO = 'L', DLATRD reduces the first NB rows and columns of a */
  /*  matrix, of which the lower triangle is supplied. */

  /*  This is an auxiliary routine called by DSYTRD. */

  /*  Arguments */
  /*  ========= */

  /*  UPLO    (input) CHARACTER */
  /*          Specifies whether the upper or lower triangular part of the */
  /*          symmetric matrix A is stored: */
  /*          = 'U': Upper triangular */
  /*          = 'L': Lower triangular */

  /*  N       (input) INTEGER */
  /*          The order of the matrix A. */

  /*  NB      (input) INTEGER */
  /*          The number of rows and columns to be reduced. */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  */
  /*          n-by-n upper triangular part of A contains the upper */
  /*          triangular part of the matrix A, and the strictly lower */
  /*          triangular part of A is not referenced.  If UPLO = 'L', the */
  /*          leading n-by-n lower triangular part of A contains the lower
  */
  /*          triangular part of the matrix A, and the strictly upper */
  /*          triangular part of A is not referenced. */
  /*          On exit: */
  /*          if UPLO = 'U', the last NB columns have been reduced to */
  /*            tridiagonal form, with the diagonal elements overwriting */
  /*            the diagonal elements of A; the elements above the diagonal
  */
  /*            with the array TAU, represent the orthogonal matrix Q as a
  */
  /*            product of elementary reflectors; */
  /*          if UPLO = 'L', the first NB columns have been reduced to */
  /*            tridiagonal form, with the diagonal elements overwriting */
  /*            the diagonal elements of A; the elements below the diagonal
  */
  /*            with the array TAU, represent the  orthogonal matrix Q as a
  */
  /*            product of elementary reflectors. */
  /*          See Further Details. */

  /*  LDA     (input) INTEGER */
  /*          The leading dimension of the array A.  LDA >= (1,N). */

  /*  E       (output) DOUBLE PRECISION array, dimension (N-1) */
  /*          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal */
  /*          elements of the last NB columns of the reduced matrix; */
  /*          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of */
  /*          the first NB columns of the reduced matrix. */

  /*  TAU     (output) DOUBLE PRECISION array, dimension (N-1) */
  /*          The scalar factors of the elementary reflectors, stored in */
  /*          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
  */
  /*          See Further Details. */

  /*  W       (output) DOUBLE PRECISION array, dimension (LDW,NB) */
  /*          The n-by-nb matrix W required to update the unreduced part */
  /*          of A. */

  /*  LDW     (input) INTEGER */
  /*          The leading dimension of the array W. LDW >= max(1,N). */

  /*  Further Details */
  /*  =============== */

  /*  If UPLO = 'U', the matrix Q is represented as a product of elementary
  */
  /*  reflectors */

  /*     Q = H(n) H(n-1) . . . H(n-nb+1). */

  /*  Each H(i) has the form */

  /*     H(i) = I - tau * v * v' */

  /*  where tau is a real scalar, and v is a real vector with */
  /*  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
  */
  /*  and tau in TAU(i-1). */

  /*  If UPLO = 'L', the matrix Q is represented as a product of elementary
  */
  /*  reflectors */

  /*     Q = H(1) H(2) . . . H(nb). */

  /*  Each H(i) has the form */

  /*     H(i) = I - tau * v * v' */

  /*  where tau is a real scalar, and v is a real vector with */
  /*  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
  */
  /*  and tau in TAU(i). */

  /*  The elements of the vectors v together form the n-by-nb matrix V */
  /*  which is needed, with W, to apply the transformation to the unreduced
  */
  /*  part of the matrix, using a symmetric rank-2k update of the form: */
  /*  A := A - V*W' - W*V'. */

  /*  The contents of A on exit are illustrated by the following examples */
  /*  with n = 5 and nb = 2: */

  /*  if UPLO = 'U':                       if UPLO = 'L': */

  /*    (  a   a   a   v4  v5 )              (  d                  ) */
  /*    (      a   a   v4  v5 )              (  1   d              ) */
  /*    (          a   1   v5 )              (  v1  1   a          ) */
  /*    (              d   1  )              (  v1  v2  a   a      ) */
  /*    (                  d  )              (  v1  v2  a   a   a  ) */

  /*  where d denotes a diagonal element of the reduced matrix, a denotes */
  /*  an element of the original matrix that is unchanged, and vi denotes */
  /*  an element of the vector defining H(i). */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Quick return if possible */

  /* Parameter adjustments */
  w_dim1 = *ldw;
  w_offset = w_dim1 + 1;
  w -= w_offset;
  --tau;
  --e;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  if (*n <= 0)
  {
    return 0;
  }

  if (lsame_(uplo, "U", 1L, 1L))
  {

    /*        Reduce last NB columns of upper triangle */

    i__1 = *n - *nb + 1;
    for (i = *n; i >= i__1; --i)
    {
      iw = i - *n + *nb;
      if (i < *n)
      {

        /*              Update A(1:i,i) */

        i__2 = *n - i;
        dgemv_("No transpose", &i, &i__2, &c_b211, &a[(i + 1) *
               a_dim1 + 1], lda, &w[i + (iw + 1) * w_dim1], ldw, &
               c_b22, &a[i * a_dim1 + 1], &c__1, 12L);
        i__2 = *n - i;
        dgemv_("No transpose", &i, &i__2, &c_b211, &w[(iw + 1) *
               w_dim1 + 1], ldw, &a[i + (i + 1) * a_dim1], lda, &
               c_b22, &a[i * a_dim1 + 1], &c__1, 12L);
      }
      if (i > 1)
      {

        /*              Generate elementary reflector H(i) to annihila
        te */
        /*              A(1:i-2,i) */

        i__2 = i - 1;
        dlarfg_(&i__2, &a[i - 1 + i * a_dim1], &a[i * a_dim1 + 1], &
                c__1, &tau[i - 1]);
        e[i - 1] = a[i - 1 + i * a_dim1];
        a[i - 1 + i * a_dim1] = 1.;

        /*              Compute W(1:i-1,i) */

        i__2 = i - 1;
        dsymv_("Upper", &i__2, &c_b22, &a[a_offset], lda, &a[i *
               a_dim1 + 1], &c__1, &c_b21, &w[iw * w_dim1 + 1], &
               c__1, 5L);
        if (i < *n)
        {
          i__2 = i - 1;
          i__3 = *n - i;
          dgemv_("Transpose", &i__2, &i__3, &c_b22, &w[(iw + 1) *
                 w_dim1 + 1], ldw, &a[i * a_dim1 + 1], &c__1, &
                 c_b21, &w[i + 1 + iw * w_dim1], &c__1, 9L);
          i__2 = i - 1;
          i__3 = *n - i;
          dgemv_("No transpose", &i__2, &i__3, &c_b211, &a[(i + 1) *
                 a_dim1 + 1], lda, &w[i + 1 + iw * w_dim1], &c__1,
                 &c_b22, &w[iw * w_dim1 + 1], &c__1, 12L);
          i__2 = i - 1;
          i__3 = *n - i;
          dgemv_("Transpose", &i__2, &i__3, &c_b22, &a[(i + 1) *
                 a_dim1 + 1], lda, &a[i * a_dim1 + 1], &c__1, &
                 c_b21, &w[i + 1 + iw * w_dim1], &c__1, 9L);
          i__2 = i - 1;
          i__3 = *n - i;
          dgemv_("No transpose", &i__2, &i__3, &c_b211, &w[(iw + 1)
                 * w_dim1 + 1], ldw, &w[i + 1 + iw * w_dim1], &
                 c__1, &c_b22, &w[iw * w_dim1 + 1], &c__1, 12L);
        }
        i__2 = i - 1;
        dscal_(&i__2, &tau[i - 1], &w[iw * w_dim1 + 1], &c__1);
        i__2 = i - 1;
        alpha = tau[i - 1] * -.5 * ddot_(&i__2, &w[iw * w_dim1 + 1], &
                                         c__1, &a[i * a_dim1 + 1], &c__1);
        i__2 = i - 1;
        daxpy_(&i__2, &alpha, &a[i * a_dim1 + 1], &c__1, &w[iw *
               w_dim1 + 1], &c__1);
      }

      /* L10: */
    }
  }
  else
  {

    /*        Reduce first NB columns of lower triangle */

    i__1 = *nb;
    for (i = 1; i <= i__1; ++i)
    {

      /*           Update A(i:n,i) */

      i__2 = *n - i + 1;
      i__3 = i - 1;
      dgemv_("No transpose", &i__2, &i__3, &c_b211, &a[i + a_dim1], lda,
             &w[i + w_dim1], ldw, &c_b22, &a[i + i * a_dim1], &c__1,
             12L);
      i__2 = *n - i + 1;
      i__3 = i - 1;
      dgemv_("No transpose", &i__2, &i__3, &c_b211, &w[i + w_dim1], ldw,
             &a[i + a_dim1], lda, &c_b22, &a[i + i * a_dim1], &c__1,
             12L);
      if (i < *n)
      {

        /*              Generate elementary reflector H(i) to annihila
        te */
        /*              A(i+2:n,i) */

        i__2 = *n - i;
        /* Computing MIN */
        i__3 = i + 2;
        dlarfg_(&i__2, &a[i + 1 + i * a_dim1], &a[min(i__3, *n) + i *
                a_dim1], &c__1, &tau[i]);
        e[i] = a[i + 1 + i * a_dim1];
        a[i + 1 + i * a_dim1] = 1.;

        /*              Compute W(i+1:n,i) */

        i__2 = *n - i;
        dsymv_("Lower", &i__2, &c_b22, &a[i + 1 + (i + 1) * a_dim1],
               lda, &a[i + 1 + i * a_dim1], &c__1, &c_b21, &w[i + 1
                   + i * w_dim1], &c__1, 5L);
        i__2 = *n - i;
        i__3 = i - 1;
        dgemv_("Transpose", &i__2, &i__3, &c_b22, &w[i + 1 + w_dim1],
               ldw, &a[i + 1 + i * a_dim1], &c__1, &c_b21, &w[i *
                   w_dim1 + 1], &c__1, 9L);
        i__2 = *n - i;
        i__3 = i - 1;
        dgemv_("No transpose", &i__2, &i__3, &c_b211, &a[i + 1 +
               a_dim1], lda, &w[i * w_dim1 + 1], &c__1, &c_b22, &w[i
                   + 1 + i * w_dim1], &c__1, 12L);
        i__2 = *n - i;
        i__3 = i - 1;
        dgemv_("Transpose", &i__2, &i__3, &c_b22, &a[i + 1 + a_dim1],
               lda, &a[i + 1 + i * a_dim1], &c__1, &c_b21, &w[i *
                   w_dim1 + 1], &c__1, 9L);
        i__2 = *n - i;
        i__3 = i - 1;
        dgemv_("No transpose", &i__2, &i__3, &c_b211, &w[i + 1 +
               w_dim1], ldw, &w[i * w_dim1 + 1], &c__1, &c_b22, &w[i
                   + 1 + i * w_dim1], &c__1, 12L);
        i__2 = *n - i;
        dscal_(&i__2, &tau[i], &w[i + 1 + i * w_dim1], &c__1);
        i__2 = *n - i;
        alpha = tau[i] * -.5 * ddot_(&i__2, &w[i + 1 + i * w_dim1], &
                                     c__1, &a[i + 1 + i * a_dim1], &c__1);
        i__2 = *n - i;
        daxpy_(&i__2, &alpha, &a[i + 1 + i * a_dim1], &c__1, &w[i + 1
               + i * w_dim1], &c__1);
      }

      /* L20: */
    }
  }

  return 0;

  /*     End of DLATRD */

} /* dlatrd_ */

/* Subroutine */ int dlarfg_(n, alpha, x, incx, tau)
integer *n;
doublereal *alpha, *x;
integer *incx;
doublereal *tau;
{
  /* System generated locals */
  integer i__1;
  doublereal d__1;

  /* Builtin functions */
  double d_sign();

  /* Local variables */
  static doublereal beta;
  extern doublereal dnrm2_();
  static integer j;
  extern /* Subroutine */ int dscal_();
  static doublereal xnorm;
  extern doublereal dlapy2_(), dlamch_();
  static doublereal safmin, rsafmn;
  static integer knt;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     February 29, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLARFG generates a real elementary reflector H of order n, such */
  /*  that */

  /*        H * ( alpha ) = ( beta ),   H' * H = I. */
  /*            (   x   )   (   0  ) */

  /*  where alpha and beta are scalars, and x is an (n-1)-element real */
  /*  vector. H is represented in the form */

  /*        H = I - tau * ( 1 ) * ( 1 v' ) , */
  /*                      ( v ) */

  /*  where tau is a real scalar and v is a real (n-1)-element */
  /*  vector. */

  /*  If the elements of x are all zero, then tau = 0 and H is taken to be
  */
  /*  the unit matrix. */

  /*  Otherwise  1 <= tau <= 2. */

  /*  Arguments */
  /*  ========= */

  /*  N       (input) INTEGER */
  /*          The order of the elementary reflector. */

  /*  ALPHA   (input/output) DOUBLE PRECISION */
  /*          On entry, the value alpha. */
  /*          On exit, it is overwritten with the value beta. */

  /*  X       (input/output) DOUBLE PRECISION array, dimension */
  /*                         (1+(N-2)*abs(INCX)) */
  /*          On entry, the vector x. */
  /*          On exit, it is overwritten with the vector v. */

  /*  INCX    (input) INTEGER */
  /*          The increment between elements of X. INCX <> 0. */

  /*  TAU     (output) DOUBLE PRECISION */
  /*          The value tau. */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --x;

  /* Function Body */
  if (*n <= 1)
  {
    *tau = 0.;
    return 0;
  }

  i__1 = *n - 1;
  xnorm = dnrm2_(&i__1, &x[1], incx);

  if (xnorm == 0.)
  {

    /*        H  =  I */

    *tau = 0.;
  }
  else
  {

    /*        general case */

    d__1 = dlapy2_(alpha, &xnorm);
    beta = -d_sign(&d__1, alpha);
    safmin = dlamch_("S", 1L);
    if (abs(beta) < safmin)
    {

      /*           XNORM, BETA may be inaccurate; scale X and recompute
      them */

      rsafmn = 1. / safmin;
      knt = 0;
L10:
      ++knt;
      i__1 = *n - 1;
      dscal_(&i__1, &rsafmn, &x[1], incx);
      beta *= rsafmn;
      *alpha *= rsafmn;
      if (abs(beta) < safmin)
      {
        goto L10;
      }

      /*           New BETA is at most 1, at least SAFMIN */

      i__1 = *n - 1;
      xnorm = dnrm2_(&i__1, &x[1], incx);
      d__1 = dlapy2_(alpha, &xnorm);
      beta = -d_sign(&d__1, alpha);
      *tau = (beta - *alpha) / beta;
      i__1 = *n - 1;
      d__1 = 1. / (*alpha - beta);
      dscal_(&i__1, &d__1, &x[1], incx);

      /*           If ALPHA is subnormal, it may lose relative accuracy
      */

      *alpha = beta;
      i__1 = knt;
      for (j = 1; j <= i__1; ++j)
      {
        *alpha *= safmin;
        /* L20: */
      }
    }
    else
    {
      *tau = (beta - *alpha) / beta;
      i__1 = *n - 1;
      d__1 = 1. / (*alpha - beta);
      dscal_(&i__1, &d__1, &x[1], incx);
      *alpha = beta;
    }
  }

  return 0;

  /*     End of DLARFG */

} /* dlarfg_ */

doublereal dlamch_(cmach, cmach_len)
char *cmach;
ftnlen cmach_len;
{
  /* Initialized data */

  static logical first = TRUE_;

  /* System generated locals */
  integer i__1;
  doublereal ret_val;

  /* Builtin functions */
  double pow_di();

  /* Local variables */
  static doublereal base;
  static integer beta;
  static doublereal emin, prec, emax;
  static integer imin, imax;
  static logical lrnd;
  static doublereal rmin, rmax, t, rmach;
  extern logical lsame_();
  static doublereal small, sfmin;
  extern /* Subroutine */ int dlamc2_();
  static integer it;
  static doublereal rnd, eps;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMCH determines double precision machine parameters. */

  /*  Arguments */
  /*  ========= */

  /*  CMACH   (input) CHARACTER*1 */
  /*          Specifies the value to be returned by DLAMCH: */
  /*          = 'E' or 'e',   DLAMCH := eps */
  /*          = 'S' or 's ,   DLAMCH := sfmin */
  /*          = 'B' or 'b',   DLAMCH := base */
  /*          = 'P' or 'p',   DLAMCH := eps*base */
  /*          = 'N' or 'n',   DLAMCH := t */
  /*          = 'R' or 'r',   DLAMCH := rnd */
  /*          = 'M' or 'm',   DLAMCH := emin */
  /*          = 'U' or 'u',   DLAMCH := rmin */
  /*          = 'L' or 'l',   DLAMCH := emax */
  /*          = 'O' or 'o',   DLAMCH := rmax */

  /*          where */

  /*          eps   = relative machine precision */
  /*          sfmin = safe minimum, such that 1/sfmin does not overflow */
  /*          base  = base of the machine */
  /*          prec  = eps*base */
  /*          t     = number of (base) digits in the mantissa */
  /*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise */
  /*          emin  = minimum exponent before (gradual) underflow */
  /*          rmin  = underflow threshold - base**(emin-1) */
  /*          emax  = largest exponent before overflow */
  /*          rmax  = overflow threshold  - (base**emax)*(1-eps) */

  /* =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Save statement .. */
  /*     .. */
  /*     .. Data statements .. */
  /*     .. */
  /*     .. Executable Statements .. */

  if (first)
  {
    first = FALSE_;
    dlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
    base = (doublereal) beta;
    t = (doublereal) it;
    if (lrnd)
    {
      rnd = 1.;
      i__1 = 1 - it;
      eps = pow_di(&base, &i__1) / 2;
    }
    else
    {
      rnd = 0.;
      i__1 = 1 - it;
      eps = pow_di(&base, &i__1);
    }
    prec = eps * base;
    emin = (doublereal) imin;
    emax = (doublereal) imax;
    sfmin = rmin;
    small = 1. / rmax;
    if (small >= sfmin)
    {

      /*           Use SMALL plus a bit, to avoid the possibility of rou
      nding */
      /*           causing overflow when computing  1/sfmin. */

      sfmin = small * (eps + 1.);
    }
  }

  if (lsame_(cmach, "E", 1L, 1L))
  {
    rmach = eps;
  }
  else if (lsame_(cmach, "S", 1L, 1L))
  {
    rmach = sfmin;
  }
  else if (lsame_(cmach, "B", 1L, 1L))
  {
    rmach = base;
  }
  else if (lsame_(cmach, "P", 1L, 1L))
  {
    rmach = prec;
  }
  else if (lsame_(cmach, "N", 1L, 1L))
  {
    rmach = t;
  }
  else if (lsame_(cmach, "R", 1L, 1L))
  {
    rmach = rnd;
  }
  else if (lsame_(cmach, "M", 1L, 1L))
  {
    rmach = emin;
  }
  else if (lsame_(cmach, "U", 1L, 1L))
  {
    rmach = rmin;
  }
  else if (lsame_(cmach, "L", 1L, 1L))
  {
    rmach = emax;
  }
  else if (lsame_(cmach, "O", 1L, 1L))
  {
    rmach = rmax;
  }

  ret_val = rmach;
  return ret_val;

  /*     End of DLAMCH */

} /* dlamch_ */


/* *********************************************************************** */

/* Subroutine */ int dlamc1_(beta, t, rnd, ieee1)
integer *beta, *t;
logical *rnd, *ieee1;
{
  /* Initialized data */

  static logical first = TRUE_;

  /* System generated locals */
  doublereal d__1, d__2;

  /* Local variables */
  static logical lrnd;
  static doublereal a, b, c, f;
  static integer lbeta;
  static doublereal savec;
  extern doublereal dlamc3_();
  static logical lieee1;
  static doublereal t1, t2;
  static integer lt;
  static doublereal one, qtr;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMC1 determines the machine parameters given by BETA, T, RND, and */
  /*  IEEE1. */

  /*  Arguments */
  /*  ========= */

  /*  BETA    (output) INTEGER */
  /*          The base of the machine. */

  /*  T       (output) INTEGER */
  /*          The number of ( BETA ) digits in the mantissa. */

  /*  RND     (output) LOGICAL */
  /*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
  /*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
  */
  /*          be a reliable guide to the way in which the machine performs
  */
  /*          its arithmetic. */

  /*  IEEE1   (output) LOGICAL */
  /*          Specifies whether rounding appears to be done in the IEEE */
  /*          'round to nearest' style. */

  /*  Further Details */
  /*  =============== */

  /*  The routine is based on the routine  ENVRON  by Malcolm and */
  /*  incorporates suggestions by Gentleman and Marovich. See */

  /*     Malcolm M. A. (1972) Algorithms to reveal properties of */
  /*        floating-point arithmetic. Comms. of the ACM, 15, 949-951. */

  /*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms */
  /*        that reveal properties of floating point arithmetic units. */
  /*        Comms. of the ACM, 17, 276-277. */

  /* =====================================================================
  */

  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Save statement .. */
  /*     .. */
  /*     .. Data statements .. */
  /*     .. */
  /*     .. Executable Statements .. */

  if (first)
  {
    first = FALSE_;
    one = 1.;

    /*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BE
    TA, */
    /*        IEEE1, T and RND. */

    /*        Throughout this routine  we use the function  DLAMC3  to ens
    ure */
    /*        that relevant values are  stored and not held in registers,
     or */
    /*        are not affected by optimizers. */

    /*        Compute  a = 2.0**m  with the  smallest positive integer m s
    uch */
    /*        that */

    /*           fl( a + 1.0 ) = a. */

    a = 1.;
    c = 1.;

    /* +       WHILE( C.EQ.ONE )LOOP */
L10:
    if (c == one)
    {
      a *= 2;
      c = dlamc3_(&a, &one);
      d__1 = -a;
      c = dlamc3_(&c, &d__1);
      goto L10;
    }
    /* +       END WHILE */

    /*        Now compute  b = 2.0**m  with the smallest positive integer
    m */
    /*        such that */

    /*           fl( a + b ) .gt. a. */

    b = 1.;
    c = dlamc3_(&a, &b);

    /* +       WHILE( C.EQ.A )LOOP */
L20:
    if (c == a)
    {
      b *= 2;
      c = dlamc3_(&a, &b);
      goto L20;
    }
    /* +       END WHILE */

    /*        Now compute the base.  a and c  are neighbouring floating po
    int */
    /*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and
     so */
    /*        their difference is beta. Adding 0.25 to c is to ensure that
     it */
    /*        is truncated to beta and not ( beta - 1 ). */

    qtr = one / 4;
    savec = c;
    d__1 = -a;
    c = dlamc3_(&c, &d__1);
    lbeta = (integer)(c + qtr);

    /*        Now determine whether rounding or chopping occurs,  by addin
    g a */
    /*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to
     a. */

    b = (doublereal) lbeta;
    d__1 = b / 2;
    d__2 = -b / 100;
    f = dlamc3_(&d__1, &d__2);
    c = dlamc3_(&f, &a);
    if (c == a)
    {
      lrnd = TRUE_;
    }
    else
    {
      lrnd = FALSE_;
    }
    d__1 = b / 2;
    d__2 = b / 100;
    f = dlamc3_(&d__1, &d__2);
    c = dlamc3_(&f, &a);
    if (lrnd && c == a)
    {
      lrnd = FALSE_;
    }

    /*        Try and decide whether rounding is done in the  IEEE  'round
     to */
    /*        nearest' style. B/2 is half a unit in the last place of the
    two */
    /*        numbers A and SAVEC. Furthermore, A is even, i.e. has last
    bit */
    /*        zero, and SAVEC is odd. Thus adding B/2 to A should not  cha
    nge */
    /*        A, but adding B/2 to SAVEC should change SAVEC. */

    d__1 = b / 2;
    t1 = dlamc3_(&d__1, &a);
    d__1 = b / 2;
    t2 = dlamc3_(&d__1, &savec);
    lieee1 = t1 == a && t2 > savec && lrnd;

    /*        Now find  the  mantissa, t.  It should  be the  integer part
     of */
    /*        log to the base beta of a,  however it is safer to determine
      t */
    /*        by powering.  So we find t as the smallest positive integer
    for */
    /*        which */

    /*           fl( beta**t + 1.0 ) = 1.0. */

    lt = 0;
    a = 1.;
    c = 1.;

    /* +       WHILE( C.EQ.ONE )LOOP */
L30:
    if (c == one)
    {
      ++lt;
      a *= lbeta;
      c = dlamc3_(&a, &one);
      d__1 = -a;
      c = dlamc3_(&c, &d__1);
      goto L30;
    }
    /* +       END WHILE */

  }

  *beta = lbeta;
  *t = lt;
  *rnd = lrnd;
  *ieee1 = lieee1;
  return 0;

  /*     End of DLAMC1 */

} /* dlamc1_ */


/* *********************************************************************** */

/* Subroutine */ int dlamc2_(beta, t, rnd, eps, emin, rmin, emax, rmax)
integer *beta, *t;
logical *rnd;
doublereal *eps;
integer *emin;
doublereal *rmin;
integer *emax;
doublereal *rmax;
{
  /* Initialized data */

  static logical first = TRUE_;
  static logical iwarn = FALSE_;

  /* Format strings */
  static char fmt_9999[] = "(//\002 WARNING. The value EMIN may be incorre\
ct:-\002,\002  EMIN = \002,i8,/\002 If, after inspection, the value EMIN loo\
ks\002,\002 acceptable please comment out \002,/\002 the IF block as marked \
within the code of routine\002,\002 DLAMC2,\002,/\002 otherwise supply EMIN \
explicitly.\002,/)";

  /* System generated locals */
  integer i__1;
  doublereal d__1, d__2, d__3, d__4, d__5;

  /* Builtin functions */
  double pow_di();
  integer s_wsfe(), do_fio(), e_wsfe();

  /* Local variables */
  static logical ieee;
  static doublereal half;
  static logical lrnd;
  static doublereal leps, zero, a, b, c;
  static integer i, lbeta;
  static doublereal rbase;
  static integer lemin, lemax, gnmin;
  static doublereal small;
  static integer gpmin;
  static doublereal third, lrmin, lrmax, sixth;
  extern /* Subroutine */ int dlamc1_();
  extern doublereal dlamc3_();
  static logical lieee1;
  extern /* Subroutine */ int dlamc4_(), dlamc5_();
  static integer lt, ngnmin, ngpmin;
  static doublereal one, two;

  /* Fortran I/O blocks */
  static cilist io___231 = { 0, 6, 0, fmt_9999, 0 };



  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMC2 determines the machine parameters specified in its argument */
  /*  list. */

  /*  Arguments */
  /*  ========= */

  /*  BETA    (output) INTEGER */
  /*          The base of the machine. */

  /*  T       (output) INTEGER */
  /*          The number of ( BETA ) digits in the mantissa. */

  /*  RND     (output) LOGICAL */
  /*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
  /*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
  */
  /*          be a reliable guide to the way in which the machine performs
  */
  /*          its arithmetic. */

  /*  EPS     (output) DOUBLE PRECISION */
  /*          The smallest positive number such that */

  /*             fl( 1.0 - EPS ) .LT. 1.0, */

  /*          where fl denotes the computed value. */

  /*  EMIN    (output) INTEGER */
  /*          The minimum exponent before (gradual) underflow occurs. */

  /*  RMIN    (output) DOUBLE PRECISION */
  /*          The smallest normalized number for the machine, given by */
  /*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
  */
  /*          of BETA. */

  /*  EMAX    (output) INTEGER */
  /*          The maximum exponent before overflow occurs. */

  /*  RMAX    (output) DOUBLE PRECISION */
  /*          The largest positive number for the machine, given by */
  /*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
  */
  /*          value of BETA. */

  /*  Further Details */
  /*  =============== */

  /*  The computation of  EPS  is based on a routine PARANOIA by */
  /*  W. Kahan of the University of California at Berkeley. */

  /* =====================================================================
  */

  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Save statement .. */
  /*     .. */
  /*     .. Data statements .. */
  /*     .. */
  /*     .. Executable Statements .. */

  if (first)
  {
    first = FALSE_;
    zero = 0.;
    one = 1.;
    two = 2.;

    /*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values
     of */
    /*        BETA, T, RND, EPS, EMIN and RMIN. */

    /*        Throughout this routine  we use the function  DLAMC3  to ens
    ure */
    /*        that relevant values are stored  and not held in registers,
     or */
    /*        are not affected by optimizers. */

    /*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
    */

    dlamc1_(&lbeta, &lt, &lrnd, &lieee1);

    /*        Start to find EPS. */

    b = (doublereal) lbeta;
    i__1 = -lt;
    a = pow_di(&b, &i__1);
    leps = a;

    /*        Try some tricks to see whether or not this is the correct  E
    PS. */

    b = two / 3;
    half = one / 2;
    d__1 = -half;
    sixth = dlamc3_(&b, &d__1);
    third = dlamc3_(&sixth, &sixth);
    d__1 = -half;
    b = dlamc3_(&third, &d__1);
    b = dlamc3_(&b, &sixth);
    b = abs(b);
    if (b < leps)
    {
      b = leps;
    }

    leps = 1.;

    /* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
    if (leps > b && b > zero)
    {
      leps = b;
      d__1 = half * leps;
      /* Computing 5th power */
      d__3 = two, d__4 = d__3, d__3 *= d__3;
      /* Computing 2nd power */
      d__5 = leps;
      d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
      c = dlamc3_(&d__1, &d__2);
      d__1 = -c;
      c = dlamc3_(&half, &d__1);
      b = dlamc3_(&half, &c);
      d__1 = -b;
      c = dlamc3_(&half, &d__1);
      b = dlamc3_(&half, &c);
      goto L10;
    }
    /* +       END WHILE */

    if (a < leps)
    {
      leps = a;
    }

    /*        Computation of EPS complete. */

    /*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3
    )). */
    /*        Keep dividing  A by BETA until (gradual) underflow occurs. T
    his */
    /*        is detected when we cannot recover the previous A. */

    rbase = one / lbeta;
    small = one;
    for (i = 1; i <= 3; ++i)
    {
      d__1 = small * rbase;
      small = dlamc3_(&d__1, &zero);
      /* L20: */
    }
    a = dlamc3_(&one, &small);
    dlamc4_(&ngpmin, &one, &lbeta);
    d__1 = -one;
    dlamc4_(&ngnmin, &d__1, &lbeta);
    dlamc4_(&gpmin, &a, &lbeta);
    d__1 = -a;
    dlamc4_(&gnmin, &d__1, &lbeta);
    ieee = FALSE_;

    if (ngpmin == ngnmin && gpmin == gnmin)
    {
      if (ngpmin == gpmin)
      {
        lemin = ngpmin;
        /*            ( Non twos-complement machines, no gradual under
        flow; */
        /*              e.g.,  VAX ) */
      }
      else if (gpmin - ngpmin == 3)
      {
        lemin = ngpmin - 1 + lt;
        ieee = TRUE_;
        /*            ( Non twos-complement machines, with gradual und
        erflow; */
        /*              e.g., IEEE standard followers ) */
      }
      else
      {
        lemin = min(ngpmin, gpmin);
        /*            ( A guess; no known machine ) */
        iwarn = TRUE_;
      }

    }
    else if (ngpmin == gpmin && ngnmin == gnmin)
    {
      if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1)
      {
        lemin = max(ngpmin, ngnmin);
        /*            ( Twos-complement machines, no gradual underflow
        ; */
        /*              e.g., CYBER 205 ) */
      }
      else
      {
        lemin = min(ngpmin, ngnmin);
        /*            ( A guess; no known machine ) */
        iwarn = TRUE_;
      }

    }
    else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
    {
      if (gpmin - min(ngpmin, ngnmin) == 3)
      {
        lemin = max(ngpmin, ngnmin) - 1 + lt;
        /*            ( Twos-complement machines with gradual underflo
        w; */
        /*              no known machine ) */
      }
      else
      {
        lemin = min(ngpmin, ngnmin);
        /*            ( A guess; no known machine ) */
        iwarn = TRUE_;
      }

    }
    else
    {
      /* Computing MIN */
      i__1 = min(ngpmin, ngnmin), i__1 = min(i__1, gpmin);
      lemin = min(i__1, gnmin);
      /*         ( A guess; no known machine ) */
      iwarn = TRUE_;
    }
    /* ** */
    /* Comment out this if block if EMIN is ok */
    if (iwarn)
    {
      first = TRUE_;
      s_wsfe(&io___231);
      do_fio(&c__1, (char *)&lemin, (ftnlen)sizeof(integer));
      e_wsfe();
    }
    /* ** */

    /*        Assume IEEE arithmetic if we found denormalised  numbers abo
    ve, */
    /*        or if arithmetic seems to round in the  IEEE style,  determi
    ned */
    /*        in routine DLAMC1. A true IEEE machine should have both  thi
    ngs */
    /*        true; however, faulty machines may have one or the other. */

    ieee = ieee || lieee1;

    /*        Compute  RMIN by successive division by  BETA. We could comp
    ute */
    /*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow dur
    ing */
    /*        this computation. */

    lrmin = 1.;
    i__1 = 1 - lemin;
    for (i = 1; i <= i__1; ++i)
    {
      d__1 = lrmin * rbase;
      lrmin = dlamc3_(&d__1, &zero);
      /* L30: */
    }

    /*        Finally, call DLAMC5 to compute EMAX and RMAX. */

    dlamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
  }

  *beta = lbeta;
  *t = lt;
  *rnd = lrnd;
  *eps = leps;
  *emin = lemin;
  *rmin = lrmin;
  *emax = lemax;
  *rmax = lrmax;

  return 0;


  /*     End of DLAMC2 */

} /* dlamc2_ */


/* *********************************************************************** */

doublereal dlamc3_(a, b)
doublereal *a, *b;
{
  /* System generated locals */
  doublereal ret_val;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
  */
  /*  the addition of  A  and  B ,  for use in situations where optimizers
  */
  /*  might hold one of these in a register. */

  /*  Arguments */
  /*  ========= */

  /*  A, B    (input) DOUBLE PRECISION */
  /*          The values A and B. */

  /* =====================================================================
  */

  /*     .. Executable Statements .. */

  ret_val = *a + *b;

  return ret_val;

  /*     End of DLAMC3 */

} /* dlamc3_ */


/* *********************************************************************** */

/* Subroutine */ int dlamc4_(emin, start, base)
integer *emin;
doublereal *start;
integer *base;
{
  /* System generated locals */
  integer i__1;
  doublereal d__1;

  /* Local variables */
  static doublereal zero, a;
  static integer i;
  static doublereal rbase, b1, b2, c1, c2, d1, d2;
  extern doublereal dlamc3_();
  static doublereal one;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMC4 is a service routine for DLAMC2. */

  /*  Arguments */
  /*  ========= */

  /*  EMIN    (output) EMIN */
  /*          The minimum exponent before (gradual) underflow, computed by
  */
  /*          setting A = START and dividing by BASE until the previous A */
  /*          can not be recovered. */

  /*  START   (input) DOUBLE PRECISION */
  /*          The starting point for determining EMIN. */

  /*  BASE    (input) INTEGER */
  /*          The base of the machine. */

  /* =====================================================================
  */

  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  a = *start;
  one = 1.;
  rbase = one / *base;
  zero = 0.;
  *emin = 1;
  d__1 = a * rbase;
  b1 = dlamc3_(&d__1, &zero);
  c1 = a;
  c2 = a;
  d1 = a;
  d2 = a;
  /* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND. */
  /*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
  if (c1 == a && c2 == a && d1 == a && d2 == a)
  {
    --(*emin);
    a = b1;
    d__1 = a / *base;
    b1 = dlamc3_(&d__1, &zero);
    d__1 = b1 * *base;
    c1 = dlamc3_(&d__1, &zero);
    d1 = zero;
    i__1 = *base;
    for (i = 1; i <= i__1; ++i)
    {
      d1 += b1;
      /* L20: */
    }
    d__1 = a * rbase;
    b2 = dlamc3_(&d__1, &zero);
    d__1 = b2 / rbase;
    c2 = dlamc3_(&d__1, &zero);
    d2 = zero;
    i__1 = *base;
    for (i = 1; i <= i__1; ++i)
    {
      d2 += b2;
      /* L30: */
    }
    goto L10;
  }
  /* +    END WHILE */

  return 0;

  /*     End of DLAMC4 */

} /* dlamc4_ */


/* *********************************************************************** */

/* Subroutine */ int dlamc5_(beta, p, emin, ieee, emax, rmax)
integer *beta, *p, *emin;
logical *ieee;
integer *emax;
doublereal *rmax;
{
  /* System generated locals */
  integer i__1;
  doublereal d__1;

  /* Local variables */
  static integer lexp;
  static doublereal oldy;
  static integer uexp, i;
  static doublereal y, z;
  static integer nbits;
  extern doublereal dlamc3_();
  static doublereal recbas;
  static integer exbits, expsum, try_;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMC5 attempts to compute RMAX, the largest machine floating-point */
  /*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum */
  /*  approximately to a power of 2.  It will fail on machines where this */
  /*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
  */
  /*  EMAX = 28718).  It will also fail if the value supplied for EMIN is */
  /*  too large (i.e. too close to zero), probably with overflow. */

  /*  Arguments */
  /*  ========= */

  /*  BETA    (input) INTEGER */
  /*          The base of floating-point arithmetic. */

  /*  P       (input) INTEGER */
  /*          The number of base BETA digits in the mantissa of a */
  /*          floating-point value. */

  /*  EMIN    (input) INTEGER */
  /*          The minimum exponent before (gradual) underflow. */

  /*  IEEE    (input) LOGICAL */
  /*          A logical flag specifying whether or not the arithmetic */
  /*          system is thought to comply with the IEEE standard. */

  /*  EMAX    (output) INTEGER */
  /*          The largest exponent before overflow */

  /*  RMAX    (output) DOUBLE PRECISION */
  /*          The largest machine floating-point number. */

  /* =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     First compute LEXP and UEXP, two powers of 2 that bound */
  /*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum */
  /*     approximately to the bound that is closest to abs(EMIN). */
  /*     (EMAX is the exponent of the required number RMAX). */

  lexp = 1;
  exbits = 1;
L10:
  try_ = lexp << 1;
  if (try_ <= -(*emin))
  {
    lexp = try_;
    ++exbits;
    goto L10;
  }
  if (lexp == -(*emin))
  {
    uexp = lexp;
  }
  else
  {
    uexp = try_;
    ++exbits;
  }

  /*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater */
  /*     than or equal to EMIN. EXBITS is the number of bits needed to */
  /*     store the exponent. */

  if (uexp + *emin > -lexp - *emin)
  {
    expsum = lexp << 1;
  }
  else
  {
    expsum = uexp << 1;
  }

  /*     EXPSUM is the exponent range, approximately equal to */
  /*     EMAX - EMIN + 1 . */

  *emax = expsum + *emin - 1;
  nbits = exbits + 1 + *p;

  /*     NBITS is the total number of bits needed to store a */
  /*     floating-point number. */

  if (nbits % 2 == 1 && *beta == 2)
  {

    /*        Either there are an odd number of bits used to store a */
    /*        floating-point number, which is unlikely, or some bits are
    */
    /*        not used in the representation of numbers, which is possible
    , */
    /*        (e.g. Cray machines) or the mantissa has an implicit bit, */
    /*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
     */
    /*        most likely. We have to assume the last alternative. */
    /*        If this is true, then we need to reduce EMAX by one because
    */
    /*        there must be some way of representing zero in an implicit-b
    it */
    /*        system. On machines like Cray, we are reducing EMAX by one
    */
    /*        unnecessarily. */

    --(*emax);
  }

  if (*ieee)
  {

    /*        Assume we are on an IEEE machine which reserves one exponent
     */
    /*        for infinity and NaN. */

    --(*emax);
  }

  /*     Now create RMAX, the largest machine number, which should */
  /*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX . */

  /*     First compute 1.0 - BETA**(-P), being careful that the */
  /*     result is less than 1.0 . */

  recbas = 1. / *beta;
  z = *beta - 1.;
  y = 0.;
  i__1 = *p;
  for (i = 1; i <= i__1; ++i)
  {
    z *= recbas;
    if (y < 1.)
    {
      oldy = y;
    }
    y = dlamc3_(&y, &z);
    /* L20: */
  }
  if (y >= 1.)
  {
    y = oldy;
  }

  /*     Now multiply by BETA**EMAX to get RMAX. */

  i__1 = *emax;
  for (i = 1; i <= i__1; ++i)
  {
    d__1 = y * *beta;
    y = dlamc3_(&d__1, &c_b21);
    /* L30: */
  }

  *rmax = y;
  return 0;

  /*     End of DLAMC5 */

} /* dlamc5_ */

doublereal dlapy2_(x, y)
doublereal *x, *y;
{
  /* System generated locals */
  doublereal ret_val, d__1;

  /* Builtin functions */
  double sqrt();

  /* Local variables */
  static doublereal xabs, yabs, w, z;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
  */
  /*  overflow. */

  /*  Arguments */
  /*  ========= */

  /*  X       (input) DOUBLE PRECISION */
  /*  Y       (input) DOUBLE PRECISION */
  /*          X and Y specify the values x and y. */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  xabs = abs(*x);
  yabs = abs(*y);
  w = max(xabs, yabs);
  z = min(xabs, yabs);
  if (z == 0.)
  {
    ret_val = w;
  }
  else
  {
    /* Computing 2nd power */
    d__1 = z / w;
    ret_val = w * sqrt(d__1 * d__1 + 1.);
  }
  return ret_val;

  /*     End of DLAPY2 */

} /* dlapy2_ */

integer ilaenv_(ispec, name, opts, n1, n2, n3, n4, name_len, opts_len)
integer *ispec;
char *name, *opts;
integer *n1, *n2, *n3, *n4;
ftnlen name_len;
ftnlen opts_len;
{
  /* System generated locals */
  integer ret_val;

  /* Builtin functions */
  /* Subroutine */
  int s_copy();
  integer s_cmp();

  /* Local variables */
  static integer i;
  static logical cname, sname;
  static integer nbmin;
  static char c1[1], c2[2], c3[3], c4[2];
  static integer ic, nb, iz, nx;
  static char subnam[6];


  /*  -- LAPACK auxiliary routine (preliminary version) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     February 20, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  ILAENV is called from the LAPACK routines to choose problem-dependent
  */
  /*  parameters for the local environment.  See ISPEC for a description of
  */
  /*  the parameters. */

  /*  This version provides a set of parameters which should give good, */
  /*  but not optimal, performance on many of the currently available */
  /*  computers.  Users are encouraged to modify this subroutine to set */
  /*  the tuning parameters for their particular machine using the option */
  /*  and problem size information in the arguments. */

  /*  This routine will not function correctly if it is converted to all */
  /*  lower case.  Converting it to all upper case is allowed. */

  /*  Arguments */
  /*  ========= */

  /*  ISPEC   (input) INTEGER */
  /*          Specifies the parameter to be returned as the value of */
  /*          ILAENV. */
  /*          = 1: the optimal blocksize; if this value is 1, an unblocked
  */
  /*               algorithm will give the best performance. */
  /*          = 2: the minimum block size for which the block routine */
  /*               should be used; if the usable block size is less than */
  /*               this value, an unblocked routine should be used. */
  /*          = 3: the crossover point (in a block routine, for N less */
  /*               than this value, an unblocked routine should be used) */
  /*          = 4: the number of shifts, used in the nonsymmetric */
  /*               eigenvalue routines */
  /*          = 5: the minimum column dimension for blocking to be used; */
  /*               rectangular blocks must have dimension at least k by m,
  */
  /*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
  */
  /*          = 6: the crossover point for the SVD (when reducing an m by n
  */
  /*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
  */
  /*               this value, a QR factorization is used first to reduce */
  /*               the matrix to a triangular form.) */
  /*          = 7: the number of processors */
  /*          = 8: the crossover point for the multishift QR and QZ methods
  */
  /*               for nonsymmetric eigenvalue problems. */

  /*  NAME    (input) CHARACTER*(*) */
  /*          The name of the calling subroutine, in either upper case or */
  /*          lower case. */

  /*  OPTS    (input) CHARACTER*(*) */
  /*          The character options to the subroutine NAME, concatenated */
  /*          into a single character string.  For example, UPLO = 'U', */
  /*          TRANS = 'T', and DIAG = 'N' for a triangular routine would */
  /*          be specified as OPTS = 'UTN'. */

  /*  N1      (input) INTEGER */
  /*  N2      (input) INTEGER */
  /*  N3      (input) INTEGER */
  /*  N4      (input) INTEGER */
  /*          Problem dimensions for the subroutine NAME; these may not all
  */
  /*          be required. */

  /* (ILAENV) (output) INTEGER */
  /*          >= 0: the value of the parameter specified by ISPEC */
  /*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
  */

  /*  Further Details */
  /*  =============== */

  /*  The following conventions have been used when calling ILAENV from the
  */
  /*  LAPACK routines: */
  /*  1)  OPTS is a concatenation of all of the character options to */
  /*      subroutine NAME, in the same order that they appear in the */
  /*      argument list for NAME, even if they are not used in determining
  */
  /*      the value of the parameter specified by ISPEC. */
  /*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
  */
  /*      that they appear in the argument list for NAME.  N1 is used */
  /*      first, N2 second, and so on, and unused problem dimensions are */
  /*      passed a value of -1. */
  /*  3)  The parameter value returned by ILAENV is checked for validity in
  */
  /*      the calling subroutine.  For example, ILAENV is used to retrieve
  */
  /*      the optimal blocksize for STRTRI as follows: */

  /*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 ) */
  /*      IF( NB.LE.1 ) NB = MAX( 1, N ) */

  /*  =====================================================================
  */

  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  switch ((int)*ispec)
  {
  case 1:
    goto L100;
  case 2:
    goto L100;
  case 3:
    goto L100;
  case 4:
    goto L400;
  case 5:
    goto L500;
  case 6:
    goto L600;
  case 7:
    goto L700;
  case 8:
    goto L800;
  }

  /*     Invalid value for ISPEC */

  ret_val = -1;
  return ret_val;

L100:

  /*     Convert NAME to upper case if the first character is lower case. */

  ret_val = 1;
  s_copy(subnam, name, 6L, name_len);
  ic = *subnam;
  iz = 'Z';
  if (iz == 90 || iz == 122)
  {

    /*        ASCII character set */

    if (ic >= 97 && ic <= 122)
    {
      *subnam = (char)(ic - 32);
      for (i = 2; i <= 6; ++i)
      {
        ic = (integer) subnam[i - 1];
        if (ic >= 97 && ic <= 122)
        {
          subnam[i - 1] = (char)(ic - 32);
        }
        /* L10: */
      }
    }

  }
  else if (iz == 233 || iz == 169)
  {

    /*        EBCDIC character set */

    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 &&
        ic <= 169)
    {
      *subnam = (char)(ic + 64);
      for (i = 2; i <= 6; ++i)
      {
        ic = (integer) subnam[i - 1];
        if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >=
            162 && ic <= 169)
        {
          subnam[i - 1] = (char)(ic + 64);
        }
        /* L20: */
      }
    }

  }
  else if (iz == 218 || iz == 250)
  {

    /*        Prime machines:  ASCII+128 */

    if (ic >= 225 && ic <= 250)
    {
      *subnam = (char)(ic - 32);
      for (i = 2; i <= 6; ++i)
      {
        ic = (integer) subnam[i - 1];
        if (ic >= 225 && ic <= 250)
        {
          subnam[i - 1] = (char)(ic - 32);
        }
        /* L30: */
      }
    }
  }

  *c1 = *subnam;
  sname = *c1 == 'S' || *c1 == 'D';
  cname = *c1 == 'C' || *c1 == 'Z';
  if (!(cname || sname))
  {
    return ret_val;
  }
  s_copy(c2, subnam + 1, 2L, 2L);
  s_copy(c3, subnam + 3, 3L, 3L);
  s_copy(c4, c3 + 1, 2L, 2L);

  switch ((int)*ispec)
  {
  case 1:
    goto L110;
  case 2:
    goto L200;
  case 3:
    goto L300;
  }

L110:

  /*     ISPEC = 1:  block size */

  /*     In these examples, separate code is provided for setting NB for */
  /*     real and complex.  We assume that NB will take the same value in */
  /*     single or double precision. */

  nb = 1;

  if (s_cmp(c2, "GE", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "TRF", 3L, 3L) == 0)
    {
      if (sname)
      {
        nb = 64;
      }
      else
      {
        nb = 64;
      }
    }
    else if (s_cmp(c3, "QRF", 3L, 3L) == 0 || s_cmp(c3, "RQF", 3L, 3L)
             == 0 || s_cmp(c3, "LQF", 3L, 3L) == 0 || s_cmp(c3, "QLF", 3L,
                 3L) == 0)
    {
      if (sname)
      {
        nb = 32;
      }
      else
      {
        nb = 32;
      }
    }
    else if (s_cmp(c3, "HRD", 3L, 3L) == 0)
    {
      if (sname)
      {
        nb = 32;
      }
      else
      {
        nb = 32;
      }
    }
    else if (s_cmp(c3, "BRD", 3L, 3L) == 0)
    {
      if (sname)
      {
        nb = 32;
      }
      else
      {
        nb = 32;
      }
    }
    else if (s_cmp(c3, "TRI", 3L, 3L) == 0)
    {
      if (sname)
      {
        nb = 64;
      }
      else
      {
        nb = 64;
      }
    }
  }
  else if (s_cmp(c2, "PO", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "TRF", 3L, 3L) == 0)
    {
      if (sname)
      {
        nb = 64;
      }
      else
      {
        nb = 64;
      }
    }
  }
  else if (s_cmp(c2, "SY", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "TRF", 3L, 3L) == 0)
    {
      if (sname)
      {
        nb = 64;
      }
      else
      {
        nb = 64;
      }
    }
    else if (sname && s_cmp(c3, "TRD", 3L, 3L) == 0)
    {
      nb = 1;
    }
    else if (sname && s_cmp(c3, "GST", 3L, 3L) == 0)
    {
      nb = 64;
    }
  }
  else if (cname && s_cmp(c2, "HE", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "TRF", 3L, 3L) == 0)
    {
      nb = 64;
    }
    else if (s_cmp(c3, "TRD", 3L, 3L) == 0)
    {
      nb = 1;
    }
    else if (s_cmp(c3, "GST", 3L, 3L) == 0)
    {
      nb = 64;
    }
  }
  else if (sname && s_cmp(c2, "OR", 2L, 2L) == 0)
  {
    if (*c3 == 'G')
    {
      if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0
          || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
          == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR",
              2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0)
      {
        nb = 32;
      }
    }
    else if (*c3 == 'M')
    {
      if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0
          || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
          == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR",
              2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0)
      {
        nb = 32;
      }
    }
  }
  else if (cname && s_cmp(c2, "UN", 2L, 2L) == 0)
  {
    if (*c3 == 'G')
    {
      if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0
          || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
          == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR",
              2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0)
      {
        nb = 32;
      }
    }
    else if (*c3 == 'M')
    {
      if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0
          || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
          == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR",
              2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0)
      {
        nb = 32;
      }
    }
  }
  else if (s_cmp(c2, "GB", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "TRF", 3L, 3L) == 0)
    {
      if (sname)
      {
        if (*n4 <= 64)
        {
          nb = 1;
        }
        else
        {
          nb = 32;
        }
      }
      else
      {
        if (*n4 <= 64)
        {
          nb = 1;
        }
        else
        {
          nb = 32;
        }
      }
    }
  }
  else if (s_cmp(c2, "PB", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "TRF", 3L, 3L) == 0)
    {
      if (sname)
      {
        if (*n2 <= 64)
        {
          nb = 1;
        }
        else
        {
          nb = 32;
        }
      }
      else
      {
        if (*n2 <= 64)
        {
          nb = 1;
        }
        else
        {
          nb = 32;
        }
      }
    }
  }
  else if (s_cmp(c2, "TR", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "TRI", 3L, 3L) == 0)
    {
      if (sname)
      {
        nb = 64;
      }
      else
      {
        nb = 64;
      }
    }
  }
  else if (s_cmp(c2, "LA", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "UUM", 3L, 3L) == 0)
    {
      if (sname)
      {
        nb = 64;
      }
      else
      {
        nb = 64;
      }
    }
  }
  else if (sname && s_cmp(c2, "ST", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "EBZ", 3L, 3L) == 0)
    {
      nb = 1;
    }
  }
  ret_val = nb;
  return ret_val;

L200:

  /*     ISPEC = 2:  minimum block size */

  nbmin = 2;
  if (s_cmp(c2, "GE", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "QRF", 3L, 3L) == 0 || s_cmp(c3, "RQF", 3L, 3L) == 0 ||
        s_cmp(c3, "LQF", 3L, 3L) == 0 || s_cmp(c3, "QLF", 3L, 3L) ==
        0)
    {
      if (sname)
      {
        nbmin = 2;
      }
      else
      {
        nbmin = 2;
      }
    }
    else if (s_cmp(c3, "HRD", 3L, 3L) == 0)
    {
      if (sname)
      {
        nbmin = 2;
      }
      else
      {
        nbmin = 2;
      }
    }
    else if (s_cmp(c3, "BRD", 3L, 3L) == 0)
    {
      if (sname)
      {
        nbmin = 2;
      }
      else
      {
        nbmin = 2;
      }
    }
    else if (s_cmp(c3, "TRI", 3L, 3L) == 0)
    {
      if (sname)
      {
        nbmin = 2;
      }
      else
      {
        nbmin = 2;
      }
    }
  }
  else if (s_cmp(c2, "SY", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "TRF", 3L, 3L) == 0)
    {
      if (sname)
      {
        nbmin = 2;
      }
      else
      {
        nbmin = 2;
      }
    }
    else if (sname && s_cmp(c3, "TRD", 3L, 3L) == 0)
    {
      nbmin = 2;
    }
  }
  else if (cname && s_cmp(c2, "HE", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "TRD", 3L, 3L) == 0)
    {
      nbmin = 2;
    }
  }
  else if (sname && s_cmp(c2, "OR", 2L, 2L) == 0)
  {
    if (*c3 == 'G')
    {
      if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0
          || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
          == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR",
              2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0)
      {
        nbmin = 2;
      }
    }
    else if (*c3 == 'M')
    {
      if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0
          || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
          == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR",
              2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0)
      {
        nbmin = 2;
      }
    }
  }
  else if (cname && s_cmp(c2, "UN", 2L, 2L) == 0)
  {
    if (*c3 == 'G')
    {
      if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0
          || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
          == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR",
              2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0)
      {
        nbmin = 2;
      }
    }
    else if (*c3 == 'M')
    {
      if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0
          || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
          == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR",
              2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0)
      {
        nbmin = 2;
      }
    }
  }
  ret_val = nbmin;
  return ret_val;

L300:

  /*     ISPEC = 3:  crossover point */

  nx = 0;
  if (s_cmp(c2, "GE", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "QRF", 3L, 3L) == 0 || s_cmp(c3, "RQF", 3L, 3L) == 0 ||
        s_cmp(c3, "LQF", 3L, 3L) == 0 || s_cmp(c3, "QLF", 3L, 3L) ==
        0)
    {
      if (sname)
      {
        nx = 128;
      }
      else
      {
        nx = 128;
      }
    }
    else if (s_cmp(c3, "HRD", 3L, 3L) == 0)
    {
      if (sname)
      {
        nx = 128;
      }
      else
      {
        nx = 128;
      }
    }
    else if (s_cmp(c3, "BRD", 3L, 3L) == 0)
    {
      if (sname)
      {
        nx = 128;
      }
      else
      {
        nx = 128;
      }
    }
  }
  else if (s_cmp(c2, "SY", 2L, 2L) == 0)
  {
    if (sname && s_cmp(c3, "TRD", 3L, 3L) == 0)
    {
      nx = 1;
    }
  }
  else if (cname && s_cmp(c2, "HE", 2L, 2L) == 0)
  {
    if (s_cmp(c3, "TRD", 3L, 3L) == 0)
    {
      nx = 1;
    }
  }
  else if (sname && s_cmp(c2, "OR", 2L, 2L) == 0)
  {
    if (*c3 == 'G')
    {
      if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0
          || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
          == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR",
              2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0)
      {
        nx = 128;
      }
    }
  }
  else if (cname && s_cmp(c2, "UN", 2L, 2L) == 0)
  {
    if (*c3 == 'G')
    {
      if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0
          || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
          == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR",
              2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0)
      {
        nx = 128;
      }
    }
  }
  ret_val = nx;
  return ret_val;

L400:

  /*     ISPEC = 4:  number of shifts (used by xHSEQR) */

  ret_val = 6;
  return ret_val;

L500:

  /*     ISPEC = 5:  minimum column dimension (not used) */

  ret_val = 2;
  return ret_val;

L600:

  /*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

  ret_val = (integer)((real) min(*n1, *n2) * (float)1.6);
  return ret_val;

L700:

  /*     ISPEC = 7:  number of processors (not used) */

  ret_val = 1;
  return ret_val;

L800:

  /*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

  ret_val = 50;
  return ret_val;

  /*     End of ILAENV */

} /* ilaenv_ */

doublereal dlansy_(norm, uplo, n, a, lda, work, norm_len, uplo_len)
char *norm, *uplo;
integer *n;
doublereal *a;
integer *lda;
doublereal *work;
ftnlen norm_len;
ftnlen uplo_len;
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2;
  doublereal ret_val, d__1, d__2, d__3;

  /* Builtin functions */
  double sqrt();

  /* Local variables */
  static doublereal absa;
  static integer i, j;
  static doublereal scale;
  extern logical lsame_();
  static doublereal value;
  extern /* Subroutine */ int dlassq_();
  static doublereal sum;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLANSY  returns the value of the one norm,  or the Frobenius norm, or
  */
  /*  the  infinity norm,  or the  element of  largest absolute value  of a
  */
  /*  real symmetric matrix A. */

  /*  Description */
  /*  =========== */

  /*  DLANSY returns the value */

  /*     DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
  /*              ( */
  /*              ( norm1(A),         NORM = '1', 'O' or 'o' */
  /*              ( */
  /*              ( normI(A),         NORM = 'I' or 'i' */
  /*              ( */
  /*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */

  /*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
  */
  /*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
  */
  /*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
  */
  /*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm. */

  /*  Arguments */
  /*  ========= */

  /*  NORM    (input) CHARACTER*1 */
  /*          Specifies the value to be returned in DLANSY as described */
  /*          above. */

  /*  UPLO    (input) CHARACTER*1 */
  /*          Specifies whether the upper or lower triangular part of the */
  /*          symmetric matrix A is to be referenced. */
  /*          = 'U':  Upper triangular part of A is referenced */
  /*          = 'L':  Lower triangular part of A is referenced */

  /*  N       (input) INTEGER */
  /*          The order of the matrix A.  N >= 0.  When N = 0, DLANSY is */
  /*          set to zero. */

  /*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          The symmetric matrix A.  If UPLO = 'U', the leading n by n */
  /*          upper triangular part of A contains the upper triangular part
  */
  /*          of the matrix A, and the strictly lower triangular part of A
  */
  /*          is not referenced.  If UPLO = 'L', the leading n by n lower */
  /*          triangular part of A contains the lower triangular part of */
  /*          the matrix A, and the strictly upper triangular part of A is
  */
  /*          not referenced. */

  /*  LDA     (input) INTEGER */
  /*          The leading dimension of the array A.  LDA >= max(N,1). */

  /*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK), */
  /*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise, */
  /*          WORK is not referenced. */

  /* =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --work;
  a_dim1 = *lda;
  a_offset = a_dim1 + 1;
  a -= a_offset;

  /* Function Body */
  if (*n == 0)
  {
    value = 0.;
  }
  else if (lsame_(norm, "M", 1L, 1L))
  {

    /*        Find max(abs(A(i,j))). */

    value = 0.;
    if (lsame_(uplo, "U", 1L, 1L))
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        i__2 = j;
        for (i = 1; i <= i__2; ++i)
        {
          /* Computing MAX */
          d__2 = value, d__3 = (d__1 = a[i + j * a_dim1], abs(d__1))
                               ;
          value = max(d__2, d__3);
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
        i__2 = *n;
        for (i = j; i <= i__2; ++i)
        {
          /* Computing MAX */
          d__2 = value, d__3 = (d__1 = a[i + j * a_dim1], abs(d__1))
                               ;
          value = max(d__2, d__3);
          /* L30: */
        }
        /* L40: */
      }
    }
  }
  else if (lsame_(norm, "I", 1L, 1L) || lsame_(norm, "O", 1L, 1L) || *
           norm == '1')
  {

    /*        Find normI(A) ( = norm1(A), since A is symmetric). */

    value = 0.;
    if (lsame_(uplo, "U", 1L, 1L))
    {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        sum = 0.;
        i__2 = j - 1;
        for (i = 1; i <= i__2; ++i)
        {
          absa = (d__1 = a[i + j * a_dim1], abs(d__1));
          sum += absa;
          work[i] += absa;
          /* L50: */
        }
        work[j] = sum + (d__1 = a[j + j * a_dim1], abs(d__1));
        /* L60: */
      }
      i__1 = *n;
      for (i = 1; i <= i__1; ++i)
      {
        /* Computing MAX */
        d__1 = value, d__2 = work[i];
        value = max(d__1, d__2);
        /* L70: */
      }
    }
    else
    {
      i__1 = *n;
      for (i = 1; i <= i__1; ++i)
      {
        work[i] = 0.;
        /* L80: */
      }
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
      {
        sum = work[j] + (d__1 = a[j + j * a_dim1], abs(d__1));
        i__2 = *n;
        for (i = j + 1; i <= i__2; ++i)
        {
          absa = (d__1 = a[i + j * a_dim1], abs(d__1));
          sum += absa;
          work[i] += absa;
          /* L90: */
        }
        value = max(value, sum);
        /* L100: */
      }
    }
  }
  else if (lsame_(norm, "F", 1L, 1L) || lsame_(norm, "E", 1L, 1L))
  {

    /*        Find normF(A). */

    scale = 0.;
    sum = 1.;
    if (lsame_(uplo, "U", 1L, 1L))
    {
      i__1 = *n;
      for (j = 2; j <= i__1; ++j)
      {
        i__2 = j - 1;
        dlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
        /* L110: */
      }
    }
    else
    {
      i__1 = *n - 1;
      for (j = 1; j <= i__1; ++j)
      {
        i__2 = *n - j;
        dlassq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
        /* L120: */
      }
    }
    sum *= 2;
    i__1 = *lda + 1;
    dlassq_(n, &a[a_offset], &i__1, &scale, &sum);
    value = scale * sqrt(sum);
  }

  ret_val = value;
  return ret_val;

  /*     End of DLANSY */

} /* dlansy_ */

/* Subroutine */ int dlassq_(n, x, incx, scale, sumsq)
integer *n;
doublereal *x;
integer *incx;
doublereal *scale, *sumsq;
{
  /* System generated locals */
  integer i__1, i__2;
  doublereal d__1;

  /* Local variables */
  static doublereal absxi;
  static integer ix;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLASSQ  returns the values  scl  and  smsq  such that */

  /*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
  */

  /*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is */
  /*  assumed to be non-negative and  scl  returns the value */

  /*     scl = max( scale, abs( x( i ) ) ). */

  /*  scale and sumsq must be supplied in SCALE and SUMSQ and */
  /*  scl and smsq are overwritten on SCALE and SUMSQ respectively. */

  /*  The routine makes only one pass through the vector x. */

  /*  Arguments */
  /*  ========= */

  /*  N       (input) INTEGER */
  /*          The number of elements to be used from the vector X. */

  /*  X       (input) DOUBLE PRECISION */
  /*          The vector for which a scaled sum of squares is computed. */
  /*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n. */

  /*  INCX    (input) INTEGER */
  /*          The increment between successive values of the vector X. */
  /*          INCX > 0. */

  /*  SCALE   (input/output) DOUBLE PRECISION */
  /*          On entry, the value  scale  in the equation above. */
  /*          On exit, SCALE is overwritten with  scl , the scaling factor
  */
  /*          for the sum of squares. */

  /*  SUMSQ   (input/output) DOUBLE PRECISION */
  /*          On entry, the value  sumsq  in the equation above. */
  /*          On exit, SUMSQ is overwritten with  smsq , the basic sum of */
  /*          squares from which  scl  has been factored out. */

  /* =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /* Parameter adjustments */
  --x;

  /* Function Body */
  if (*n > 0)
  {
    i__1 = (*n - 1) * *incx + 1;
    i__2 = *incx;
    for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2)
    {
      if (x[ix] != 0.)
      {
        absxi = (d__1 = x[ix], abs(d__1));
        if (*scale < absxi)
        {
          /* Computing 2nd power */
          d__1 = *scale / absxi;
          *sumsq = *sumsq * (d__1 * d__1) + 1;
          *scale = absxi;
        }
        else
        {
          /* Computing 2nd power */
          d__1 = absxi / *scale;
          *sumsq += d__1 * d__1;
        }
      }
      /* L10: */
    }
  }
  return 0;

  /*     End of DLASSQ */

} /* dlassq_ */

doublereal dlaran_(iseed)
integer *iseed;
{
  /* System generated locals */
  doublereal ret_val;

  /* Local variables */
  static integer it1, it2, it3, it4;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     February 29, 1992 */

  /*     .. Array Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLARAN returns a random DOUBLE PRECISION number from a uniform (0,1)
  */
  /*  distribution. */

  /*  Arguments */
  /*  ========= */

  /*  ISEED   (input/output) INTEGER array, dimension (4) */
  /*          On entry, the seed of the random number generator; the array
  */
  /*          elements must be between 0 and 4095, and ISEED(4) must be */
  /*          odd. */
  /*          On exit, the seed is updated. */

  /*  Further Details */
  /*  =============== */

  /*  This routine uses a multiplicative congruential method with modulus */
  /*  2**48 and multiplier 33952834046453 (see G.S.Fishman, */
  /*  'Multiplicative congruential random number generators with modulus */
  /*  2**b: an exhaustive analysis for b = 32 and a partial analysis for */
  /*  b = 48', Math. Comp. 189, pp 331-344, 1990). */

  /*  48-bit integers are stored in 4 integer array elements with 12 bits */
  /*  per element. Hence the routine is portable across machines with */
  /*  integers of 32 bits or more. */

  /*  =====================================================================
  */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     multiply the seed by the multiplier modulo 2**48 */

  /* Parameter adjustments */
  --iseed;

  /* Function Body */
  it4 = iseed[4] * 2549;
  it3 = it4 / 4096;
  it4 -= it3 << 12;
  it3 = it3 + iseed[3] * 2549 + iseed[4] * 2508;
  it2 = it3 / 4096;
  it3 -= it2 << 12;
  it2 = it2 + iseed[2] * 2549 + iseed[3] * 2508 + iseed[4] * 322;
  it1 = it2 / 4096;
  it2 -= it1 << 12;
  it1 = it1 + iseed[1] * 2549 + iseed[2] * 2508 + iseed[3] * 322 + iseed[4]
        * 494;
  it1 %= 4096;

  /*     return updated seed */

  iseed[1] = it1;
  iseed[2] = it2;
  iseed[3] = it3;
  iseed[4] = it4;

  /*    convert 48-bit integer to a DOUBLE PRECISION number in the interval
  (0,1)*/

  ret_val = ((doublereal) it1 + ((doublereal) it2 + ((doublereal) it3 + (
                                   doublereal) it4 * 2.44140625e-4) * 2.44140625e-4) * 2.44140625e-4)
            * 2.44140625e-4;
  return ret_val;

  /*     End of DLARAN */

} /* dlaran_ */

logical lsame_(ca, cb, ca_len, cb_len)
char *ca, *cb;
ftnlen ca_len;
ftnlen cb_len;
{
  /* System generated locals */
  logical ret_val;

  /* Local variables */
  static integer inta, intb, zcode;


  /*  -- LAPACK auxiliary routine (version 1.1) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     February 29, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  LSAME returns .TRUE. if CA is the same letter as CB regardless of */
  /*  case. */

  /*  Arguments */
  /*  ========= */

  /*  CA      (input) CHARACTER*1 */
  /*  CB      (input) CHARACTER*1 */
  /*          CA and CB specify the single characters to be compared. */

  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     Test if the characters are equal */

  ret_val = *ca == *cb;
  if (ret_val)
  {
    return ret_val;
  }

  /*     Now test for equivalence if both characters are alphabetic. */

  zcode = 'Z';

  /*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
  /*     machines, on which ICHAR returns a value with bit 8 set. */
  /*     ICHAR('A') on Prime machines returns 193 which is the same as */
  /*     ICHAR('A') on an EBCDIC machine. */

  inta = *ca;
  intb = *cb;

  if (zcode == 90 || zcode == 122)
  {

    /*        ASCII is assumed - ZCODE is the ASCII code of either lower o
    r */
    /*        upper case 'Z'. */

    if (inta >= 97 && inta <= 122)
    {
      inta += -32;
    }
    if (intb >= 97 && intb <= 122)
    {
      intb += -32;
    }

  }
  else if (zcode == 233 || zcode == 169)
  {

    /*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
     or */
    /*        upper case 'Z'. */

    if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta
        >= 162 && inta <= 169)
    {
      inta += 64;
    }
    if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb
        >= 162 && intb <= 169)
    {
      intb += 64;
    }

  }
  else if (zcode == 218 || zcode == 250)
  {

    /*        ASCII is assumed, on Prime machines - ZCODE is the ASCII cod
    e */
    /*        plus 128 of either lower or upper case 'Z'. */

    if (inta >= 225 && inta <= 250)
    {
      inta += -32;
    }
    if (intb >= 225 && intb <= 250)
    {
      intb += -32;
    }
  }
  ret_val = inta == intb;

  /*     RETURN */

  /*     End of LSAME */

  return ret_val;
} /* lsame_ */

logical lsamen_(n, ca, cb, ca_len, cb_len)
integer *n;
char *ca, *cb;
ftnlen ca_len;
ftnlen cb_len;
{
  /* System generated locals */
  integer i__1;
  logical ret_val;

  /* Builtin functions */
  integer i_len();

  /* Local variables */
  static integer i;
  extern logical lsame_();


  /*  -- LAPACK auxiliary routine (version 1.0b) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     February 29, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  LSAMEN  tests if the first N letters of CA are the same as the */
  /*  first N letters of CB, regardless of case. */
  /*  LSAMEN returns .TRUE. if CA and CB are equivalent except for case */
  /*  and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA ) */
  /*  or LEN( CB ) is less than N. */

  /*  Arguments */
  /*  ========= */

  /*  N       (input) INTEGER */
  /*          The number of characters in CA and CB to be compared. */

  /*  CA      (input) CHARACTER*(*) */
  /*  CB      (input) CHARACTER*(*) */
  /*          CA and CB specify two character strings of length at least N.
  */
  /*          Only the first N characters of each string will be accessed.
  */

  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  ret_val = FALSE_;
  if (i_len(ca, ca_len) < *n || i_len(cb, cb_len) < *n)
  {
    goto L20;
  }

  /*     Do for each character in the two strings. */

  i__1 = *n;
  for (i = 1; i <= i__1; ++i)
  {

    /*        Test if the characters are equal using LSAME. */

    if (! lsame_(ca + (i - 1), cb + (i - 1), 1L, 1L))
    {
      goto L20;
    }

    /* L10: */
  }
  ret_val = TRUE_;

L20:
  return ret_val;

  /*     End of LSAMEN */

} /* lsamen_ */

