#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*!\file lcp_lexicolemke.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).
 * Try \f$(z,w)\f$ such that:
 *
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *    0 \le z \perp M z + q = w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 *
 * M is an (n x n) matrix, q , w and z an n-vectors.
 *
 *!\fn  lcp_lexicolemke( double *vec, double *q , int *nn, int *itermax , int *ispeak , double *zlem ,
 *             double *wlem, int *it_end, int *info )
 *
 * lcp_lexicolemke is a direct solver for LCP based on pivoting method principle for degenrate problem.
 * Choice of pivot variable is performed via lexicographic ordering
 * Ref: "The Linear Complementary Problem" Cottle, Pang, Stone (1992)
 *
 * \param double* vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
 * \param double* q       Unchanged parameter which contains the components of the right hand side vector.
 * \param int* nn         Unchanged parameter which represents the dimension of the system.
 * \param int* itermax    Unchanged parameter which represents the maximum number of iterations allowed.
 * \param int* ispeak     Unchanged parameter which represents the output log identifiant
 *                        0 - no output
 *                        0 < identifiant
 *
 * \param double* zlem    Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param double* wlem    Modified parameter which returns the solution of the problem.
 * \param int* info       Modified parameter which returns the termination value
 *                        0 - pivot termination
 *                        1 - maximal pivot number reached
 *
 * \author Mathieu Renouf
 */

void lcp_lexicolemke(double *vec, double *q , int *nn, int *itermax , double *zlem , int *ispeak ,
                     double *wlem, int *it_end , int *info)
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
