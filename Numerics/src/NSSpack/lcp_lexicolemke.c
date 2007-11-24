/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/

/*!\file lcp_lexicolemke.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
 * Try \f$(z,w)\f$ such that:\n
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *    w - M z = q\\
 *    0 \le z \perp w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 *
 * where M is an (\f$nn \times nn\f$)-matrix, q , w and z nn-vectors.
 */


/*!\fn  void lcp_lexicolemke( int *nn , double *vec , double *q , double *zlem , double *wlem , int *info , int *iparamLCP , double *dparamLCP )

  lcp_lexicolemke is a direct solver for LCP based on pivoting method principle for degenerate problem.\n
  Choice of pivot variable is performed via lexicographic ordering
  Ref: "The Linear Complementary Problem" Cottle, Pang, Stone (1992)\n


  \param nn      On enter, an integer which represents the dimension of the system.
  \param vec     On enter, a (\f$nn\times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
  \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
  \param zlem    On return, a nn-vector of doubles which contains the solution of the problem.
  \param wlem    On return, a nn-vector of doubles which contains the solution of the problem.
  \param info    On return, an integer which returns the termination value:\n
                 0 : convergence\n
                 1 : iter = itermax\n
                 2 : negative diagonal term\n

  \param iparamLCP  On enter/return, a vetor of integers:\n
                 - iparamLCP[0] = itermax On enter, the maximum number of pivots allowed.
                 - iparamLCP[1] = ispeak  On enter, the output log identifiant:\n
                        0 : no output\n
                        >0: active screen output\n
                 - iparamLCP[2] = it_end  On return, the number of pivots performed by the algorithm.

  \param dparamLCP  On enter/return, a vetor of doubles (not used).\n



  \author Mathieu Renouf

 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

void lcp_lexicolemke(int *nn , double *vec , double *q , double *zlem , double *wlem , int *info , int *iparamLCP , double *dparamLCP)
{

  int i, drive, block, Ifound;
  int ic, jc;
  int dim, dim2, ITER;
  int nobasis;
  int itermax, ispeak;

  double z0, zb, dblock;
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

  /* construction of A matrix such as
   * A = [ q | Id | -d | -M ] with d = (1,...1)
   */

  for (ic = 0 ; ic < dim; ++ic)
    for (jc = 0 ; jc < dim; ++jc)
      A[ic][jc + dim + 2] = -vec[dim * jc + ic];

  for (ic = 0 ; ic < dim; ++ic) A[ic][0] = q[ic];

  for (ic = 0 ; ic < dim; ++ic) A[ic][ic + 1 ] =  1.0;
  for (ic = 0 ; ic < dim; ++ic) A[ic][dim + 1] = -1.0;

  /* End of construction of A */

  Ifound = 0;


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

  } /* end while*/

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


  iparamLCP[2] = ITER;

  if (Ifound) *info = 0;
  else *info = 1;

  free(basis);

  for (i = 0 ; i < dim ; ++i) free(A[i]);
  free(A);

}
