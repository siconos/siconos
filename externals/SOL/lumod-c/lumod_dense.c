
#include "SiconosBlas.h"
#include "lumod_dense.h"
#include <assert.h>     // for assert
#include <math.h>       // for fabs
#include <stdlib.h>     // for abs
#include "commonlib.h"  // for subvec, MACHINEPREC, TINYNUMBER

/* TODO kill this */
#define ADDR_OFFSET 1

#ifdef __cplusplus
#undef restrict
#include <sys/cdefs.h>  // for __restrict
#define restrict __restrict
#endif

/* Modification for siconos, Olivier Huber:
 * - add restrict keyword to pointer to ensure aliasing
 * - replace BLAS Fortran API with CBLAS
 */

/*
                        LUMOD 2.0  (29 Jun 1999)
      The following routines are members of a package for maintaining
      L*C = U,  an LU factorization of a dense square matrix C:

            LUmod    Lprod     Usolve
            elm      elmgen    LUback   LUforw
      Also needed:
            cblas_daxpy    cblas_dcopy     cblas_ddot    (from the BLAS)

      The first three routines may be used as follows:

            call LUmod ( 1, ...         )    to build up the matrix C,
            call LUmod ( 1,2,3 or 4, ...)    to modify   the matrix C,

            call Lprod ( 1, ... b, x )
            call Usolve( 1, ...    x )       to solve  C x = b,

            call Usolve( 2, ... b    )
            call Lprod ( 2, ... b, z )       to solve  C'x = b.

      BEWARE:  LUmod uses a machine-dependent constant eps.

      Michael Saunders
      Systems Optimization Laboratory,
      Department of EESOR, Stanford University.

      05 Apr 1981: Original version used 2-D arrays and plane rotations.
      09 Oct 1986: Plane rotations optionally replaced by eliminations.
      17 Apr 1990: Mode 2 altered to put new column at end.
      26 Apr 1990: L and U stored by rows in 1-D arrays.
      21 Aug 1991: LUmod created from QRmod3.
      29 Jun 1999: Minor clean-up.
                   Comments added for each mode to show what updates
                   should be made to two lists, say listr[*] and listc[*],
                   that the user will probably need to maintain.
                   They contain the row and column indices of some bigger
                   matrix that are now the rows and columns of C.
      29 Jun 1999: Mode 2 changed back to original: put new column
                   in same place as the column being replaced.
                   This generally needs more arithmetic, but it
                   simplifies updating listr[*] and listc[*],
                   and it keeps C symmetric where possible.
      30 Jun 1999: Mode 4 now replaces a row and column by
                   swapping them with the last row and column.
      ------------------------------------------------------------------
      LUmod  modifies the matrix factorization  L*C = U  where
      L and C are square, L is a product of 2x2 stabilized elimination
      matrices, and U is upper triangular.

      The modifications depend upon 'mode' as follows.
      n is not altered by any of the entries to LUmod.

      ==================================================================
      mode=1.  Expand L and U by adding a new row and column to C.
               The dimension expands from n-1 to n (not from n to n+1).
      ==================================================================
         On entry:
               y[*], z[*]  contain the new row and column.
               y[n]        contains the new diagonal element of C.
         On exit:
               y[*]        is altered.
         Not used:
               krow, kcol, w[*].
         Updates (before call to LUmod):
               n        = n + 1
               listr[n] = new row index
               listc[n] = new column index

      ==================================================================
      mode=2.  Replace the kcol-th column of C.
      ==================================================================
         On entry:
               kcol        says which column of C is being replaced.
               z[*]        contains the new column.
         On exit:
               y[*], w[*]  are altered.
               z[*]        is not altered.
         Not used:
               krow
         Updates:
               listc[kcol] = new column index

      ==================================================================
      mode=3.  Replace the krow-th row of C.
      ==================================================================
         On entry:
               krow        says which row of C is being replaced.
               y[*]        contains the new row.
         On exit:
               y[*], z[*], w[*]  are altered.
         Not used:
               kcol
         Updates:
               listr[krow] = new row index

      ==================================================================
      mode=4.  Shrink L and U by deleting the
               krow-th row and the kcol-th column of C
               (swapping them with the n-th row and column).
               n should be decreased to n-1 by the user after the call.
      ==================================================================
         On entry:
               krow, kcol  say which row and col of C are being replaced.
         On exit:
               y[*], z[*], w[*]  are altered.
         Updates (after call to LUmod):
               listr[krow] = listr[n]
               listc[kcol] = listc[n]
               n        = n - 1


      ==================================================================
      The LU factors are well defined even if C is singular
      (in which case U is also singular).  At some stage the factors
      will be used to solve systems of linear equations of the form
             C x = b     or     C' x = b.
      The diagonals of U should be reasonably different from zero.
      ==================================================================


      Other input parameters:

      maxmod   The maximum dimension of C.

      L[*]     An array of length at least maxmod*maxmod.
               When L has maximum dimension maxmod, its rows are stored
               contiguously in L[*].
               For lower dimensions, each row of L starts in the same
               place but fills only the front of its allowable space.
               Row i of L starts in position (i - 1)*maxmod + 1.

      U[*]     An array of length at least maxmod*(maxmod + 1)/2.
               When U has maximum dimension maxmod, the upper-triangular
               part of its rows are stored contiguously in U[*].
               For lower dimensions, each row of U starts in the same
               place but fills only the front of its allowable space.
               Row i of U starts in position (i-1)*maxmod + (3-i)*i/2.

      w[n]     A work vector.

      ------------------------------------------------------------------ */


void LUmod_dense(int mode, int maxmod, int n, int krow, int kcol,
                 double* restrict L, double* restrict U, double* restrict y, double* restrict z, double* restrict w)
{
  int    first, last, n1, i, j, lastu, lastl, ls, ll, lu, incu;
  double zero = 0.0;
  double one  = 1.0;
  double eps  = MACHINEPREC;  // The machine precision -- A value slightly too large is OK.

  n1 = n - 1;
  assert(n1 >= 0);

  if(mode == 1)
  {
    /*   ---------------------------------------------------------------
         mode = 1.  Add a row y and a column z.
         The LU factors will expand in dimension from n-1 to n.
         The new diagonal element of C is in y[n].
         --------------------------------------------------------------- */
    lastu    = n1*maxmod + (3 - n)*n/2;
    lastl    = n1*maxmod + n;
    ls       = n1*maxmod + 1;
    L[lastl] = one;
    if(n == 1)
    {
      U[lastu] = y[n];
      return;
    }

    /*   Compute L*z and temporarily store it in w;
         (changed to w from last row of L by KE). */

    Lprod_dense(1, maxmod, n1, L+ADDR_OFFSET, z+ADDR_OFFSET, w+ADDR_OFFSET);

    /*       Copy L*z into the new last column of U.
             Border L with zeros. */

    ll     = ls;
    lu     = n;
    incu   = maxmod - 1;
    for(j = 1; j<=n1; j++)
    {
      U[lu] = w[j];
      L[ll] = zero;
      ll++;
      lu   += incu;
      incu--;
    }

    ll     = n;
    for(i = 1; i<=n1; i++)
    {
      L[ll] = zero;
      ll   += maxmod;
    }

    /*   Add row y to the factorization
         using a forward sweep of eliminations. */

    last   = n;
    LUforw_dense(1, last, n, n, maxmod, eps, L, U, y);

  }
  else if(mode == 2)
  {
    /*   ---------------------------------------------------------------
         mode=2.  Replace the kcol-th column of C by the vector z.
         ---------------------------------------------------------------*/

    /*   Compute w = L*z. */

    Lprod_dense(1, maxmod, n, L+ADDR_OFFSET, z+ADDR_OFFSET, w+ADDR_OFFSET);

    /*   Copy the top of w into column kcol of U. */

    lu     = kcol;
    incu   = maxmod - 1;
    for(i = 1; i<=kcol; i++)
    {
      U[lu] = w[i];
      lu   += incu;
      incu--;
    }

    if(kcol < n)
    {

      /*      Find w[last], the last nonzero in the bottom part of w.
              Eliminate elements last-1, last-2, ... kcol+1 of w[*]
              using a partial backward sweep of eliminations. */

      first  = kcol + 1;
      last   = n;
      LUback_dense(first, &last, n, n, maxmod, eps, L, U, y, w);
      y[kcol] = w[last];

      /*      Eliminate elements kcol, kcol+1, ... last-1 of y[*]
              using a partial forward sweep of eliminations. */

      LUforw_dense(kcol, last, n, n, maxmod, eps, L, U, y);
    }

  }
  else if(mode == 3)
  {
    /*   ---------------------------------------------------------------
         mode=3.  Replace the krow-th row of C by the vector y.
         --------------------------------------------------------------- */
    if(n == 1)
    {
      L[1] = one;
      U[1] = y[1];
      return;
    }

    /*   Copy the krow-th column of L into w, and zero the column. */

    ll     = krow;
    for(i = 1; i<=n; i++)
    {
      w[i]  = L[ll];
      L[ll] = zero;
      ll   += maxmod;
    }

    /*   Reduce the krow-th column of L to the unit vector e(last).
         where 'last' is determined by LUback.
         This is done by eliminating elements last-1, last-2, ..., 1
         using a backward sweep of eliminations.
         On exit, row 'last' of U is a spike stored in z, whose first
         nonzero entry is in z[first].  However, z will be discarded. */

    first  = 1;
    last   = n;
    LUback_dense(first, &last, n, n, maxmod, eps, L, U, z, w);

    /*   Replace the 'last' row of L by the krow-th unit vector. */

    ll     = (last - 1)*maxmod;
    for(j = 1; j<=n; j++)
      L[ll + j] = zero;
    L[ll + krow] = one;

    /*   Eliminate the elements of the new row y,
         using a forward sweep of eliminations. */

    LUforw_dense(1, last, n, n, maxmod, eps, L, U, y);

  }

  else if(mode == 4)
  {
    /*   ---------------------------------------------------------------
         mode=4.  Delete the krow-th row and the kcol-th column of C.
                  Replace them by the last row and column respectively.
         --------------------------------------------------------------- */

    /*   First, move the last column into position kcol. */

    if(kcol < n)
    {

      /*      Set w = last column of U. */

      lu     = n;
      incu   = maxmod - 1;
      for(i = 1; i<=n; i++)
      {
        w[i] = U[lu];
        lu  += incu;
        incu--;
      }

      /*      Copy the top of w into column kcol of U. */

      lu     = kcol;
      incu   = maxmod - 1;
      for(i = 1; i<=kcol; i++)
      {
        U[lu] = w[i];
        lu   += incu;
        incu--;
      }

      /*      U now has only n-1 columns.
              Find w[last], the last nonzero in the bottom part of w.
              Eliminate elements last-1, last-2, ... kcol+1 of w[*]
              using a partial backward sweep of eliminations. */

      first  = kcol + 1;
      last   = n;
      LUback_dense(first, &last, n, n1, maxmod, eps, L, U, y, w);
      y[kcol] = w[last];

      /*      Eliminate elements kcol, kcol+1, ... last-1 of y[*]
              using a partial forward sweep of eliminations. */

      LUforw_dense(kcol, last, n, n1, maxmod, eps, L, U, y);
    }

    /*   Now, move the last row into position krow. */

    /*   Swap columns krow and n of L, using w = krow-th column of L. */

    cblas_dcopy(n, &L[subvec(krow)] + ADDR_OFFSET, maxmod, w + ADDR_OFFSET, 1);
    if(krow < n)
      cblas_dcopy(n, &L[subvec(n)] + ADDR_OFFSET, maxmod, &L[subvec(krow)] + ADDR_OFFSET, maxmod);

    /*   Reduce the last column of L (in w) to the unit vector e(n).
         This is done by eliminating elements n-1, n-2, ..., 1
         using a backward sweep of eliminations. */

    last   = - n;
    LUback_dense(1, &last, n, n1, maxmod, eps, L, U, z, w);
//     printvec(n, w, 0);
//     printmatUT(maxmod, n, U, 0);
//     printmatSQ(maxmod, n, L, 0);
  }

  /*     End of LUmod  */
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void Lprod_dense(int mode, int maxmod, int n,
                 double* restrict L, double* restrict y, double* restrict z)
{

  /* This function works with C indexing */

  /* ------------------------------------------------------------------
     If  mode = 1,  Lprod  computes  z = L*y.
     If  mode = 2,  Lprod  computes  z = L[transpose]*y.
     L is stored by rows in L[*].  It is equivalent to storing
     L[transpose] by columns in a 2-D array L[maxmod,n].
     y is not altered.
     ------------------------------------------------------------------ */

  if(mode == 1)
  {
    for(int i = 0, ll = 0; i < n; ++i, ll += maxmod)
    {
      z[i]  = cblas_ddot(n, &L[ll], 1, y, 1);
    }
  }
  else
  {
    /*        call dzero ( n, z, 1 ) */
    for(int i = 0; i<n; ++i)
      z[i]  = 0.;

    for(int i = 0, ll = 0; i < n; ++i, ll += maxmod)
    {
      cblas_daxpy(n, y[i], &L[ll], 1, z, 1);
    }
  }

  /*     End of Lprod */
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void Usolve_dense(int mode, int maxmod, int n, double* restrict U, double* restrict y)
{

  double sum, hold;
  int    lu, incu, nu;

  /* ------------------------------------------------------------------
     If  mode = 1,  Usolve solves U * y[new] = y[old].
     If  mode = 2,  Usolve solves U[transpose] * y[new] = y[old].
     U is upper triangular, stored by rows.
     y is overwritten by the solution.
     ------------------------------------------------------------------ */

  if(mode == 1)
  {
    lu     = (n - 1)*maxmod + (3 - n)*n/2;
    hold   = U[lu];
    y[n]  /= hold;
    incu   = maxmod + 1 - n;
    nu     = 0;

    for(int i = n-1; i>=1; --i)
    {
      nu++;
      incu++;
      lu   -= incu;
      sum   = y[i] - cblas_ddot(nu, &U[lu] + ADDR_OFFSET, 1, &y[i] + ADDR_OFFSET, 1);
      hold  = U[lu];
      y[i]  = sum / hold;
    }
  }
  else
  {
    lu     = 1;
    incu   = maxmod;
    nu     = n - 1;

    for(int i = 1; i<n; ++i)
    {
      hold = U[lu];
      y[i] /= hold;
      cblas_daxpy(nu, -y[i], &U[lu] + ADDR_OFFSET, 1, &y[i] + ADDR_OFFSET, 1);
      lu   += incu;
      incu--;
      nu--;
    }
    hold = U[lu];
    y[n]  /= hold;
  }

  /*    End of Usolve */
}


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void elm_dense(int first, int last, double* restrict x, double* restrict y, double cs, double sn)
{
  /* ------------------------------------------------------------------
     elm  computes the elemental transformation  (x y)*E  and returns
     the result in  (x y),  where the 2 by 2 matrix  E  is defined by
     cs and sn as follows:

       E  =  ( 1  sn )  if  cs >= zero,    E  =  (     1 )  otherwise.
             (     1 )                           ( 1  sn )
     ------------------------------------------------------------------ */

  double  zero = 0.0;

  if(cs < zero)
    cblas_dswap(last, x + ADDR_OFFSET, 1, y + ADDR_OFFSET, 1);
  if(sn != zero)
    cblas_daxpy(last, sn, x + ADDR_OFFSET, 1, y + ADDR_OFFSET, 1);

  /*     End of elm */
}


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void elmgen_dense(double* restrict x, double* restrict y, double eps, double* restrict cs, double* restrict sn)
{

  /* ------------------------------------------------------------------
     elmgen  generates an elimination transformation  E  such that
     (x y)*E  =  (x  0)   or   (y  0),  depending on the relative
     sizes of  x  and  y.
     eps  is an input parameter -- the machine precision.

     CAUTION: It is assumed that the data generating x and y
              are in general well-scaled, so that if both x and y
              are smaller than 'tiny', they are both changed to zero.
              This is an attempt to save work and reduce underflow.
     ------------------------------------------------------------------ */

  double zero = 0.0;
  double one = 1.0;
  double tiny;

  tiny   = eps * TINYNUMBER;
  (*cs)     = zero;

  if(fabs((*x)) >= fabs((*y)))
  {
    if(fabs((*x)) <= tiny)
    {
      (*sn)  =   zero;
      (*x)   =   zero;
    }
    else
      (*sn)  = - (*y)/(*x);
  }
  else
  {
    if(fabs((*y)) <= tiny)
    {
      (*sn)  =   zero;
      (*x)   =   zero;
    }
    else
    {
      (*cs)  = - one;
      (*sn)  = - (*x)/(*y);
      (*x)   = (*y);
    }
  }

  (*y)      = zero;

  /*     End of elmgen */
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void LUback_dense(int first, int* restrict last, int n, int nu,
                  int maxmod, double eps,
                  double* restrict L, double* restrict U, double* restrict y, double* restrict z)
{

  /* ------------------------------------------------------------------
     LUback updates the factors LC = U by performing a backward sweep
     to eliminate all but the 'last' nonzero in the column vector z[*],
     stopping at z[first].

     If 'last' is positive, LUback searches backwards for a nonzero
     element in z and possibly alters 'last' accordingly.
     Otherwise, 'last' will be reset to abs(last) and so used.

     L     is n by n.
     U     is n by nu.
     y[*]  will eventually contain a row spike in row 'last' of U.
     The 'spike' row of L begins at L[ls].

     18 Mar 1990: First version with L and U stored row-wise.
     29 Jun 1999: Save w[last] = zlast at end so column replace
                  can do correct forward sweep.
     ------------------------------------------------------------------ */

  int    i, lz, lu, ll,ls, incu, numu;
  double zero = 0.0;
  double zlast, cs, sn;

  if((*last) > 0)
  {

    /*   Find the last significant element in z[*]. */

    for(i = (*last); i>first; i--)
      if(fabs(z[i]) > eps) break;

    (*last)   = i;
  }
  else
    (*last)   = abs((*last));

  /*  Load the 'last' row of U into the end of y
      and do the backward sweep. */

  zlast  = z[(*last)];
  lu     = ((*last)  - 1)*maxmod + (3 - (*last))*(*last)/2;
  ll     = ((*last)  - 1)*maxmod + 1;
  ls     = ll;
  incu   = maxmod + 1 - (*last);
  numu   = nu     + 1 - (*last);
  if(numu > 0)
    cblas_dcopy(numu, &U[subvec(lu)] + ADDR_OFFSET, 1, &y[subvec((*last))] + ADDR_OFFSET, 1);

  for(lz = (*last) - 1; lz>=first; lz--)
  {
    ll    -= maxmod;
    incu++;
    lu    -= incu;
    y[lz]  = zero;

    /*   See if this element of z is worth eliminating.
         We compare z[lz] with the current last nonzero, zlast. */

    if(fabs(z[lz]) <= eps * fabs(zlast)) continue;

    /*   Generate a 2x2 elimination and apply it to U and L. */

    numu   = nu + 1 - lz;
    elmgen_dense(&zlast, &z[lz], eps, &cs, &sn);
    elm_dense(1, numu, &y[subvec(lz)], &U[subvec(lu)], cs, sn);
    elm_dense(1, n, &L[subvec(ls)], &L[subvec(ll)], cs, sn);
  }

  z[(*last)] = zlast;

  /*   End of LUback */
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void LUforw_dense(int first, int last, int n, int nu, int maxmod,
                  double eps, double* restrict L, double* restrict U, double* restrict y)
{

  /* ------------------------------------------------------------------
     LUforw updates the factors LC = U by performing a forward sweep
     to eliminate subdiagonals in a row spike in row 'last' of U.

     L     is n by n.
     U     is n by nu.
     y[*]  contains the row spike.  The first nonzero to be eliminated
           is in y[first] or later.
     The 'spike' row of L begins at L[ls].

     18 Mar 1990: First version with L and U stored row-wise.
     ------------------------------------------------------------------ */

  int    lu, ll, ls, incu, ly, numu;
  double cs, sn;

  lu     = (first - 1)*maxmod + (3 - first)*first/2;
  ll     = (first - 1)*maxmod + 1;
  ls     = (last  - 1)*maxmod + 1;
  incu   = maxmod + 1 - first;

  for(ly = first; ly<last; ly++)
  {

    /*   See if this element of y is worth eliminating.
         We compare y[ly] with the corresponding diagonal of U. */

    if(fabs(y[ly]) > eps * fabs(U[lu]))
    {

      /*   Generate a 2x2 elimination and apply it to U and L. */

      numu   = nu - ly;
      elmgen_dense(&U[lu], &y[ly], eps, &cs, &sn);
      if(numu > 0)
        elm_dense(1, numu, &U[lu], &y[ly], cs, sn);
      elm_dense(1, n, &L[subvec(ll)], &L[subvec(ls)], cs, sn);
    }

    ll    += maxmod;
    lu    += incu;
    incu--;
  }

  /* Copy the remaining part of y into U. */

  numu     = nu - last + 1;
  if(numu > 0)
    cblas_dcopy(numu, &y[subvec(last)] + ADDR_OFFSET, 1, &U[subvec(lu)] + ADDR_OFFSET, 1);

  /*   End of LUforw */
}


