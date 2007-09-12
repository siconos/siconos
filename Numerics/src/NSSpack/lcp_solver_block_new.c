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
/*!\file lcp_solver_block_new.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
 * Try \f$(z,w)\f$ such that:\n
 * \f$
 *    0 \le z \perp Mz + q = w \ge 0
 * \f$
 *
 * Here M is an (\f$dn \times dn\f$) matrix , q , w and z dn-vector.\n
 *
 * This system of equalities and inequalities is solved thanks to @ref block_lcp solvers. */

/*!\fn int lcp_solver_block(SparseBlockStructuredMatrix *blmat, double *q, method *pt , double *z , double *w , int *it_end ,               int *itt_end ,double *res )
*  lcp_solver_block is a generic interface allowing the call of one of the block LCP solvers.\n
*
* \param blmat    On enter, the sparse block matrix M of the LCP
* \param q        On enter, the vector of doubles of the LCP
* \param pt       On enter, a union which contains a LCP structure
*
* \param z        On enter/return, a vector of doubles which contains the initial iterate for the LCP(q,M) and returns the solution of the problem.
* \param w        On return, a vector of doubles which returns the solution of the problem.
* \param it_end   On return, an integer which returns the final number of block iterations
* \param itt_end  On return, an integer which returns the total number of local LCP iterations or pivots
* \param res      On return, a double which returns the final value of error criteria.
*
* \return info    Integer identifiant for the solver result\n
*                 0 : convergence\n
*                 >0 : no convergence (see solver for specific info value)\n
*
* lcp_solver_block is a generic interface which consider LCP with a sparse block structure. The global LCP is solved
* as a succession of local LCP solved via lcp_solver.\n
*
* list Keywords to call solvers:
*
*   - Lemke    for lcp_lexicolemke
*   - PGS      for lcp_pgs
*   - RPGS     for lcp_rpgs
*   - CPG      for lcp_cpg
*   - QP       for lcp_qp
*   - NSQP     for lcp_nsqp
*   - Latin    for lcp_latin
*   - Newton   for lcp_newton_min
*
* Data file example:\n
* If we consider the matrix M and the right-hand-side q defined as
*
* \f$
* M=\left[\begin{array}{cccc|cc|cc}
*          1 & 2 & 0 & 4   & 3 &-1   & 0 & 0\\
*          2 & 1 & 0 & 0   & 4 & 1   & 0 & 0\\
*          0 & 0 & 1 &-1   & 0 & 0   & 0 & 0\\
*          5 & 0 &-1 & 6   & 0 & 6   & 0 & 0\\
*          \hline
*          0 & 0 & 0 & 0   & 1 & 0   & 0 & 5\\
*          0 & 0 & 0 & 0   & 0 & 2   & 0 & 2\\
*          \hline
*          0 & 0 & 2 & 1   & 0 & 0   & 2 & 2\\
*          0 & 0 & 2 & 2   & 0 & 0   & -1 & 2\\
*        \end{array}\right] \quad, q=\left[\begin{array}{c}-1\\-1\\0\\-1\\\hline 1\\0\\\hline -1\\2\end{array}\right].
* \f$
*
* then
* - the number of non null blocks is 6 (blmat.nbblocks) and the number of diagonal blocks is 3 (blmat.size)
* - the vector blmat.blocksize of diagonal blocks sizes is equal to [4,2,2]
* - the vectors blmat.ColumnIndex and blmat.RowIndex of non null blocks indices are equal to [0,1,1,2,0,2] and [0,0,1,1,2,2]
* - the blmat.block contains all non null block matrices stored in Fortran order (column by column) as\n
*   blmat.block[0] = {1,2,0,5,2,1,0,0,0,0,1,-1,4,0,-1,6}\n
*   blmat.block[1] = {3,4,0,0,-1,1,0,6}\n
*   ...\n
*   blmat.block[5] = {2,-1,2,2}
* \author Mathieu Renouf & Pascal Denoyelle
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "NSSpack.h"
#endif
#include "blaslapack.h"

int lcp_solver_block(SparseBlockStructuredMatrix *blmat, double *q, method *pt , double *z , double *w , int *it_end ,
                     int *itt_end , double *res)
{

  int info;
  int n, nbbl, nbblrow, blsizemax;
  int rowprecbl, rowcurbl, indicrow, rowsize;
  int colprecbl, colcurbl, indiccol, colsize;
  int i, j, k;
  int iter, itermax, totaliter;
  int info1;
  int incx, incy;

  double a1, b1, tol;
  double qs, err,/* num, */den;

  double *adrcurbl, *adrbldiag;
  double *rhs;
  /*  double *ww;*/

  char NOTRANS = 'N';

  /*   char savename[40]; */

  *it_end   = 0;
  *itt_end   = 0;
  *res      = 0.0;
  info      = 1;
  info1     = 1;
  itermax   = pt->lcp.itermax;
  tol       = pt->lcp.tol;
  totaliter = 0;

  nbbl = blmat->nbblocks;
  if (nbbl < 1)
  {
    printf(" Problem : Null LCP block matrix !!! \n");
    return 1;
  }

  nbblrow = blmat->size;
  n = 0;
  blsizemax = 0;
  for (i = 0 ; i < nbblrow ; i++)
  {
    k = blmat->blocksize[i];
    n += k;
    if (k > blsizemax) blsizemax = k;
  }

  i = 0;
  while ((i < (n - 1)) && (q[i] >= 0.)) i++;
  if ((i == (n - 1)) && (q[n - 1] >= 0.))
  {
    /* TRIVIAL CASE : q >= 0
     * z = 0 and w = q is solution of LCP(q,M)
     */
    for (j = 0 ; j < n; j++)
    {
      z[j] = 0.0;
      w[j] = q[j];
    }
    if (pt->lcp.chat > 0) printf("Trivial case of block LCP : positive vector q \n");
    return 0;
  }

  /*  ww   = ( double* )malloc(         n * sizeof( double ) );*/
  rhs  = (double*)malloc(blsizemax * sizeof(double));

  incx = 1;
  qs = dnrm2_((integer *)&n , q , (integer *)&incx);
  den = 1.0 / qs;

  /* Initialization of z and w */

  //   for( i = 0 ; i < n ; i++ ) z[i] = 0.;

  incx = 1;
  incy = 1;
  //   dcopy_( (integer *)&n , q , (integer *)&incx , w , (integer *)&incy );


  /*   test on matrix structure : is there null blocks rows ? */
  rowprecbl = -1;
  for (i = 0 ; i < nbbl ; i++)
  {
    rowcurbl = blmat->RowIndex[i];
    if (rowcurbl > rowprecbl + 1)
    {
      printf(" Null blocks row in LCP matrix !!!\n");
      free(rhs);
      return 1;
    }
    rowprecbl = rowcurbl;
  }
  if (rowcurbl != nbblrow - 1)
  {
    printf(" Null blocks row in LCP matrix !!!\n");
    free(rhs);
    return 1;
  }

  iter = 0;
  err  = 1.;

  /*  while( ( iter < itermax ) && ( err > tol ) ){*/
  while ((iter < itermax) && (info))
  {

    iter++;

    /*    incx = 1;
        incy = 1;
        dcopy_( (integer *)&n , w , (integer *)&incx , ww , (integer *)&incy );*/

    rowprecbl = -1;
    rowsize = 0;
    indicrow = 0;

    for (i = 0 ; i < nbbl ; i++)
    {
      rowcurbl = blmat->RowIndex[i];
      colcurbl = blmat->ColumnIndex[i];
      if (rowcurbl != rowprecbl)   /*  new row  */
      {
        if (rowprecbl != -1)
        {
          if (adrbldiag != NULL)
          {
            /* Local LCP resolution  */
            pt->lcp.iter = 0;
            info1 = lcp_solver(adrbldiag , rhs , &rowsize , pt , &z[indicrow] , &w[indicrow]);
            totaliter += pt->lcp.iter;
            if (info1 > 0)
            {
              if (pt->lcp.chat > 0)
              {
                printf(" Sub LCP solver failed at global iteration %d in lcp_solver_block", iter);
                printf(" for block(%d,%d)\n", rowprecbl, rowprecbl);
              }
              free(rhs);
              return 1;
            }
          }
          else
          {
            printf("NULL diagonal block in LCP !!!\n");
            free(rhs);
            return 1;
          }
        }
        adrbldiag = NULL;

        indiccol = 0;
        for (j = 0 ; j < colcurbl ; j++) indiccol += blmat->blocksize[j];
        colprecbl = colcurbl;

        indicrow += rowsize;
        for (j = rowprecbl + 1 ; j < rowcurbl ; j++) indicrow += blmat->blocksize[j];

        rowprecbl = rowcurbl;
        rowsize = blmat->blocksize[rowcurbl];
        for (j = 0 ; j < rowsize ; j++) rhs[j] = q[indicrow + j];
      }
      if (rowcurbl == colcurbl)   /* diagonal block  */
      {
        adrbldiag = blmat->block[i];
      }
      else                        /* extra diagonal block  */
      {
        for (j = colprecbl ; j < colcurbl ; j++) indiccol += blmat->blocksize[j];
        colprecbl = colcurbl;

        rowsize = blmat->blocksize[rowcurbl];
        colsize = blmat->blocksize[colcurbl];
        adrcurbl = blmat->block[i];

        a1 = 1.;
        b1 = 1.;
        incx = 1;
        incy = 1;
        dgemv_(&NOTRANS , (integer *)&rowsize , (integer *)&colsize , &a1 , adrcurbl , (integer *)&rowsize , &z[indiccol] ,
               (integer *)&incx , &b1 , rhs , (integer *)&incy);
      }
    }

    if (adrbldiag != NULL)
    {
      /* Local LCP resolution  */


      /* strcpy(savename,pt->lcp.name); */
      /*strcpy(pt->lcp.name,"NLGS"); */

      pt->lcp.iter = 0;
      info1 = lcp_solver(adrbldiag , rhs , &rowsize , pt , &z[indicrow] , &w[indicrow]);

      /* strcpy(pt->lcp.name,savename); */

      totaliter += pt->lcp.iter;
      if (info1 > 0)
      {
        if (pt->lcp.chat > 0)
        {
          printf(" Sub LCP solver failed at global iteration %d in lcp_solver_block", iter);
          printf(" for block(%d,%d)\n", rowprecbl, rowprecbl);
        }
        free(rhs);
        return 1;
      }
    }
    else
    {
      printf("NULL diagonal block in LCP !!!\n");
      free(rhs);
      return 1;
    }

    /* **** Criterium convergence **** */

    info = filter_result_LCP_block(blmat, q, z, tol, pt->lcp.chat, w);
    /*    qs   = -1;
        incx =  1;
        incy =  1;

        daxpy_( (integer *)&n , &qs , w , (integer *)&incx , ww , (integer *)&incy );

        num = dnrm2_( (integer *)&n, ww , (integer *)&incx );
        err = num*den;*/
    /* **** ********************* **** */
  }

  *it_end  = iter;
  *itt_end = totaliter;
  *res     = err;

  /*  free(ww);*/
  free(rhs);

  if (pt->lcp.chat > 0)
  {
    /*    if( err > tol ){*/
    /*      printf( " The residue is : %g \n", err );*/
    /*      return 1;*/
    if (info)
    {
      printf(" No convergence of lcp_solver_block after %d iterations\n" , iter);
    }
    else
    {
      printf(" Convergence of lcp_solver_block after %d iterations\n" , iter);
      /*      printf( " The residue is : %g \n", err );*/
    }
  }
  /*  else{
      if( err > tol ) return 1;
    }*/

  /*  info = filter_result_LCP_block(blmat,q,z,tol,pt->lcp.chat,w);*/

  return info;

}
