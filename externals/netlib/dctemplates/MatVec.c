/* MatVec.f -- translated by f2c (version of 20 August 1993  13:15:44).
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

/* Table of constant values */

static integer c__1 = 1;

int matvec_(alpha, x, beta, y)

doublereal *alpha, *x, *beta, *y;
{
  extern int dgemv_();


  /*     This MatVec routine assumes the matrix is in dense format, */
  /*     and uses the BLAS DGEMV. */

  /*     .. Common Blocks .. */
  /*     MAXDIM2 = MAXDIM*MAXDIM. */



  /* Parameter adjustments */
  --y;
  --x;

  /* Executable Statements */
  dgemv_("NOTRANSPOSE", &matdim_1.n, &matdim_1.n, alpha, system_1.a, &
         matdim_1.lda, &x[1], &c__1, beta, &y[1], &c__1, 11L);

  return 0;

} /* matvec_ */


/*     ================================================= */
int matvectrans_(alpha, x, beta, y)

doublereal *alpha, *x, *beta, *y;
{
  extern int dgemv_();


  /*     This MatVec routine assumes the matrix is in dense format, */
  /*     and uses the BLAS DGEMV. */

  /*     .. Common Blocks .. */
  /*     MAXDIM2 = MAXDIM*MAXDIM. */


  /* Parameter adjustments */
  --y;
  --x;

  /* Function Body */
  dgemv_("TRANSPOSE", &matdim_1.n, &matdim_1.n, alpha, system_1.a, &
         matdim_1.lda, &x[1], &c__1, beta, &y[1], &c__1, 9L);

  return 0;

}

