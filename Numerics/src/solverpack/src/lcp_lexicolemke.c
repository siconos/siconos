#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*!\file lcp_lexicolemke.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
 * Try \f$(z,w)\f$ such that:\n
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *    M z + q= w\\
 *    0 \le z \perp w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 *
 * where M is an (n x n)-matrix, q , w and z n-vectors.
 *
 * \fn  lcp_lexicolemke( int *nn , double *vec , double *q , double *z , int *info ,
 *                       int *iparamLCP , double *dparamLCP )
 *
 * lcp_lexicolemke is a direct solver for LCP based on pivoting method principle for degenrate problem.
 * Choice of pivot variable is performed via lexicographic ordering
 * Ref: "The Linear Complementary Problem" Cottle, Pang, Stone (1992)
 *
 * Generic lcp parameters:\n
 *
 * \param nn      Unchanged parameter which represents the dimension of the system.
 * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
 * \param q       Unchanged parameter which contains the components of the right hand side vector.
 * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param w       Modified parameter which returns the solution of the problem.
 * \param info    Modified parameter which returns the termination value\n
 *                0 - convergence\n
 *                1 - iter = itermax\n
 *                2 - negative diagonal term\n
 *
 * Specific Lexico Lemke parameters:\n
 *
 * \param iparamLCP[0] = itermax Input unchanged parameter which represents the maximum number of pivots allowed.
 * \param iparamLCP[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
 *                       0 - no output\n
 *                       0 < active screen output\n
 * \param iparamLCP[2] = it_end  Output modified parameter which returns the number of pivots performed by the algorithm.
 *
 *
 * \author Mathieu Renouf
 *
 */

void lcp_lexicolemke(int *nn , double *vec , double *q , double *zlem , double *wlem , int *info ,
                     int *iparamLCP , double *dparamLCP)
{

  int i, drive, block, Ifound;
  int ic, jc;
  int dim, dim2, ITER;
  int nobasis;
  int itermax, ispeak;

  double qs, z0, zb, dblock;
  double pivot, tovip;
  double tmp;
  int *basis;
  double** A;

  dim = *nn;
  dim2 = 2 * (dim + 1);

  /*input*/

  itermax = iparamLCP[0];
  ispeak  = iparamLCP[1];

  /*output*/

  iparamLCP[2] = 0;

  /* Allocation */

  basis = (int *)malloc(dim * sizeof(int));
  A = (double **)malloc(dim * sizeof(double*));

  for (ic = 0 ; ic < dim; ++ic)
    A[ic] = (double *)malloc(dim2 * sizeof(double));

  for (ic = 0 ; ic < dim; ++ic)
    for (jc = 0 ; jc < dim2; ++jc)
      A[ic][jc] = 0.0;

  /*! construction of A matrix such as
   * A = [ Id | -d | -M | q ] with d = (1,...1)
   */

  for (ic = 0 ; ic < dim; ++ic)
    for (jc = 0 ; jc < dim; ++jc)
      A[ic][jc + dim + 2] = -vec[dim * jc + ic];

  for (ic = 0 ; ic < dim; ++ic) A[ic][0] = q[ic];

  for (ic = 0 ; ic < dim; ++ic) A[ic][ic + 1 ] =  1.0;
  for (ic = 0 ; ic < dim; ++ic) A[ic][dim + 1] = -1.0;

  /* End of construction of A */

  /*! STEP 0
   * qs = min{ q[i], i=1,...,NC }
   */


  qs = q[0];

  for (ic = 1 ; ic < dim ; ++ic)
  {
    if (q[ic] < qs) qs = q[ic];
  }

  Ifound = 0;

  if (qs >= 0)
  {

    /*! TRIVIAL CASE
     * z = 0 and w = q is solution of LCP(q,M)
     */

    for (ic = 0 ; ic < dim; ++ic)
    {
      zlem[ic] = 0.0;
      wlem[ic] = q[ic];
      z0 = 0.0;
    }

  }
  else
  {

    for (ic = 0 ; ic < dim  ; ++ic) basis[ic] = ic + 1;

    drive = dim + 1;
    block = 0;
    z0 = A[block][0];
    ITER = 0;

    /* Start research of argmin lexico */
    /* With this first step the covering vector enter in the basis */

    for (ic = 1 ; ic < dim ; ++ic)
    {
      zb = A[ic][0];
      if (zb < z0)
      {
        z0    = zb;
        block = ic;
      }
      else if (zb == z0)
      {
        for (jc = 0 ; jc < dim ; ++jc)
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
    for (ic = drive + 1 ; ic < dim2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

    /* */

    for (ic = 0 ; ic < block ; ++ic)
    {
      tmp = A[ic][drive];
      for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }
    for (ic = block + 1 ; ic < dim ; ++ic)
    {
      tmp = A[ic][drive];
      for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }

    nobasis = basis[block];
    basis[block] = drive;

    while (ITER < itermax && !Ifound)
    {

      ++ITER;

      if (nobasis < dim + 1)      drive = nobasis + (dim + 1);
      else if (nobasis > dim + 1) drive = nobasis - (dim + 1);

      /* Start research of argmin lexico for minimum ratio test */

      pivot = 1e20;
      block = -1;

      for (ic = 0 ; ic < dim ; ++ic)
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
            for (jc = 1 ; jc < dim + 1 ; ++jc)
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

      if (basis[block] == dim + 1) Ifound = 1;

      /* Pivot < block , drive > */

      pivot = A[block][drive];
      tovip = 1.0 / pivot;
      A[block][drive] = 1;

      for (ic = 0       ; ic < drive ; ++ic) A[block][ic] = A[block][ic] * tovip;
      for (ic = drive + 1 ; ic < dim2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

      /* */

      for (ic = 0 ; ic < block ; ++ic)
      {
        tmp = A[ic][drive];
        for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
      }
      for (ic = block + 1 ; ic < dim ; ++ic)
      {
        tmp = A[ic][drive];
        for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
      }

      nobasis = basis[block];
      basis[block] = drive;

    }

    for (ic = 0 ; ic < dim; ++ic)
    {
      drive = basis[ic];
      if (drive < dim + 1)
      {
        zlem[drive - 1] = 0.0;
        wlem[drive - 1] = A[ic][0];
      }
      else if (drive > dim + 1)
      {
        zlem[drive - dim - 2] = A[ic][0];
        wlem[drive - dim - 2] = 0.0;
      }
    }

  }

  iparamLCP[2] = ITER;

  if (Ifound) *info = 0;
  else *info = 1;

  free(basis);

  for (i = 0 ; i < dim ; ++i) free(A[i]);
  free(A);

}
