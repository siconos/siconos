#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*!\file lemke_lcp.c


   This subroutine allows the resolution of LCP (Linear Complementary Problem).
   Try \f$(z,w)\f$ such that:

\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
0 \le z \perp w \ge 0\\
\end{array}
\right.
\f$

  here M is an n by n  matrix, q an n-dimensional vector, w an n-dimensional  vector and z an n-dimensional vector.
*/

/*!\fn  lemke_lcp(double vec[], double *qq,int *nn, int *itermax, double *z, double *w, int *it_end, double *res, int *info )


   lemke_lcp is a direct solver for LCP.


   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param qq On enter a pointer over doubles containing the components of the double vector.
   \param nn On enter a pointer over integers, the dimension of the second member.
   \param itermax On enter a pointer over integers, the maximum iterations required.
   \param it_end On enter a pointer over integers, the number of iterations carried out.
   \param res On return a pointer over doubles, the error value.
   \param z On return double vector, the solution of the problem.
   \param w On return double vector, the solution of the problem.
   \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).
   \author Mathieu Renouf
*/

void lexicolemke_lcp(double *vec, double *q , int *nn, int *itermax , double *zlem , int *ispeak ,
                     double *wlem, int *it_end, double *res, int *info)
{


  int i, j, k, drive, block, Ifound;
  int ic, jc, iadj, idlaw;
  int icd, jcd, ian, jan;
  int DIM, DIM2, NC, ITER;
  int nobasis;

  double qs, z0, zb, dblock;
  double pivot, tovip;
  double tmp;
  int *basis;
  double** A;

  DIM = *nn;
  DIM2 = 2 * (DIM + 1);

  basis = (int *)malloc(DIM * sizeof(int));
  A = (double **)malloc(DIM * sizeof(double*));

  for (ic = 0 ; ic < DIM; ++ic)
    A[ic] = (double *)malloc(DIM2 * sizeof(double));

  for (ic = 0 ; ic < DIM; ++ic)
    for (jc = 0 ; jc < DIM2; ++jc)
      A[ic][jc] = 0.0;

  /*! construction of A matrix such as
   * A = [ Id | -d | -M | q ] with d = (1,...1)
   */

  for (ic = 0 ; ic < DIM; ++ic)
    for (jc = 0 ; jc < DIM; ++jc)
      A[ic][jc + DIM + 2] = -vec[DIM * jc + ic];

  for (ic = 0 ; ic < DIM; ++ic) A[ic][0] = q[ic];

  for (ic = 0 ; ic < DIM; ++ic) A[ic][ic + 1 ] =  1.0;
  for (ic = 0 ; ic < DIM; ++ic) A[ic][DIM + 1] = -1.0;

  /* End of construction of A */

  /*! STEP 0
   * qs = min{ q[i], i=1,...,NC }
   */


  qs = q[0];

  for (ic = 1 ; ic < DIM ; ++ic)
  {
    if (q[ic] < qs) qs = q[ic];
  }

  Ifound = 0;

  if (qs >= 0)
  {

    /*! TRIVIAL CASE
     * z = 0 and w = q is solution of LCP(q,M)
     */

    for (ic = 0 ; ic < DIM; ++ic)
    {
      zlem[ic] = 0.0;
      wlem[ic] = q[ic];
      z0 = 0.0;
    }

  }
  else
  {

    for (ic = 0 ; ic < DIM  ; ++ic) basis[ic] = ic + 1;

    drive = DIM + 1;
    block = 0;
    z0 = A[block][0];
    ITER = 0;

    /* Start research of argmin lexico */
    /* With this first step the covering vector enter in the basis */

    for (ic = 1 ; ic < DIM ; ++ic)
    {
      zb = A[ic][0];
      if (zb < z0)
      {
        z0    = zb;
        block = ic;
      }
      else if (zb == z0)
      {
        for (jc = 0 ; jc < DIM ; ++jc)
        {
          dblock = A[block][1 + jc] - A[ic][1 + jc];
          if (dblock < 0)
          {
            break;
          }
          else if (dblock > 0)
          {
            block = ic;
            break;
          }
        }
      }
    }

    /* Stop research of argmin lexico */

    pivot = A[block][drive];
    tovip = 1.0 / pivot;

    /* Pivot < block , drive > */

    A[block][drive] = 1;
    for (ic = 0       ; ic < drive ; ++ic) A[block][ic] = A[block][ic] * tovip;
    for (ic = drive + 1 ; ic < DIM2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

    /* */

    for (ic = 0 ; ic < block ; ++ic)
    {
      tmp = A[ic][drive];
      for (jc = 0 ; jc < DIM2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }
    for (ic = block + 1 ; ic < DIM ; ++ic)
    {
      tmp = A[ic][drive];
      for (jc = 0 ; jc < DIM2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }

    nobasis = basis[block];
    basis[block] = drive;

    while (ITER < *itermax && !Ifound)
    {

      ++ITER;

      if (nobasis < DIM + 1)      drive = nobasis + (DIM + 1);
      else if (nobasis > DIM + 1) drive = nobasis - (DIM + 1);

      /* Start research of argmin lexico for minimum ratio test */

      pivot = 1e20;
      block = -1;

      for (ic = 0 ; ic < DIM ; ++ic)
      {
        zb = A[ic][drive];
        if (zb > 0.0)
        {
          z0 = A[ic][0] / zb;
          if (z0 > pivot) continue;
          if (z0 < pivot)
          {
            pivot = z0;
            block = ic;
          }
          else
          {
            for (jc = 1 ; jc < DIM + 1 ; ++jc)
            {
              dblock = A[block][jc] / pivot - A[ic][jc] / zb;
              if (dblock < 0) break;
              else if (dblock > 0)
              {
                block = ic;
                break;
              }
            }
          }
        }
      }
      if (block == -1) break;

      if (basis[block] == DIM + 1) Ifound = 1;

      /* Pivot < block , drive > */

      pivot = A[block][drive];
      tovip = 1.0 / pivot;
      A[block][drive] = 1;

      for (ic = 0       ; ic < drive ; ++ic) A[block][ic] = A[block][ic] * tovip;
      for (ic = drive + 1 ; ic < DIM2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

      /* */

      for (ic = 0 ; ic < block ; ++ic)
      {
        tmp = A[ic][drive];
        for (jc = 0 ; jc < DIM2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
      }
      for (ic = block + 1 ; ic < DIM ; ++ic)
      {
        tmp = A[ic][drive];
        for (jc = 0 ; jc < DIM2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
      }

      nobasis = basis[block];
      basis[block] = drive;

    }

    for (ic = 0 ; ic < DIM; ++ic)
    {
      drive = basis[ic];
      if (drive < DIM + 1)
      {
        zlem[drive - 1] = 0.0;
        wlem[drive - 1] = A[ic][0];
      }
      else if (drive > DIM + 1)
      {
        zlem[drive - DIM - 2] = A[ic][0];
        wlem[drive - DIM - 2] = 0.0;
      }
    }

  }

  *it_end = ITER;
  if (Ifound) *info = 0;
  else *info = 1;

  free(basis);

  for (i = 0 ; i < DIM ; ++i) free(A[i]);
  free(A);

}
