/* Siconos-Numerics, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
/*!\file lcp_solver_block_pred_vec.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
 * Try \f$(z,w)\f$ such that:\n
 * \f$
 *    0 \le z \perp Mz + q = w \ge 0
 * \f$
 *
 * Here M is an (\f$dn \times dn\f$) matrix , q , w and z dn-vector.\n
 *
 * This system of equalities and inequalities is solved thanks to  block_lcp solvers.
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include "LA.h"

int lcp_solver_block_pred_vec(SparseBlockStructuredMatrix *blmat,
                              SparseBlockStructuredMatrixPred *blmatpred,
                              int nbmethod,
                              int maxiterglob, double tolglob,
                              double *q, method_lcp **ptvec ,
                              double *z , double *w , int *it_end , int *itt_end , double *res)
{

  static int firsttime = 1;
  int info;
  int n, nbbl, nbblrow, blsizemax;
  int rowprecbl, rowcurbl, indicrow, rowsize;
  int colprecbl, colcurbl, indiccol, colsize;
  int i, j, k;
  int iter, totaliter;
  int info1;
  int incx, incy;

  double a1, b1;
  double qs, err, den;

  double *adrcurbl, *adrbldiag;
  double *rhs;
  /*  double *ww;*/

  int iterrow0;
  method_lcp *pt;
  int numsolver;
  int *soltype;
  int *indic;
  int *indicop;
  double *submatlcp;
  double *submatlcpop;
  int *ipiv;
  int *sizesublcp;
  int *sizesublcpop;
  double *subq;
  double *bufz;
  double *newz;
  double *workspace;

  int sizelcp;

  *it_end   = 0;
  *itt_end   = 0;
  *res      = 0.0;
  info      = 1;
  info1     = 1;
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

  sizelcp = n;

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
  qs = DNRM2(n , q , incx);
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
  iterrow0 = 0;

  while ((iter < maxiterglob) && (info))
  {

    iter++;

    /*    incx = 1;
        incy = 1;
        dcopy_( (integer *)&n , w , (integer *)&incx , ww , (integer *)&incy );*/

    rowprecbl = -1;
    rowsize = 0;
    indicrow = 0;
    numsolver = 0;
    pt = ptvec[numsolver];

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

            indic           = blmatpred->indic[rowprecbl];
            indicop         = blmatpred->indicop[rowprecbl];
            submatlcp       = blmatpred->submatlcp[rowprecbl];
            submatlcpop     = blmatpred->submatlcpop[rowprecbl];
            ipiv            = blmatpred->ipiv[rowprecbl];
            sizesublcp      = &(blmatpred->sizesublcp[rowprecbl]);
            sizesublcpop    = &(blmatpred->sizesublcpop[rowprecbl]);
            subq            = blmatpred->subq[rowprecbl];
            bufz            = blmatpred->bufz[rowprecbl];
            newz            = blmatpred->newz[rowprecbl];
            workspace       = blmatpred->workspace[rowprecbl];

            info1 = lcp_solver_pred(adrbldiag , rhs , &rowsize , pt , &z[indicrow] , &w[indicrow] ,
                                    firsttime, soltype , indic , indicop , submatlcp , submatlcpop , ipiv , sizesublcp , sizesublcpop ,
                                    subq , bufz , newz , workspace);

            totaliter += pt->lcp.iter;
            if (rowprecbl == 0) iterrow0 += pt->lcp.iter;

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

            if (numsolver + 1 < nbmethod)
            {
              numsolver++;
              pt = ptvec[numsolver];
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
      else                      /* extra diagonal block  */
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
        DGEMV(LA_NOTRANS, rowsize , colsize , a1 , adrcurbl , rowsize , &z[indiccol] , incx , b1 , rhs , incy);
      }
    }

    if (adrbldiag != NULL)
    {

      /*        strcpy(pt->lcp.name,"PGS"); */
      /*        printf("row : %d , solver %s\n",rowprecbl,pt->lcp.name);*/
      /* Local LCP resolution  */
      pt->lcp.iter = 0;
      /*
              info1 = lcp_solver(adrbldiag , rhs , &rowsize , pt , &z[indicrow] , &w[indicrow] );

              strcpy(pt->lcp.name,savename);
      */
      indic           = blmatpred->indic[rowprecbl];
      indicop         = blmatpred->indicop[rowprecbl];
      submatlcp       = blmatpred->submatlcp[rowprecbl];
      submatlcpop     = blmatpred->submatlcpop[rowprecbl];
      ipiv            = blmatpred->ipiv[rowprecbl];
      sizesublcp      = &(blmatpred->sizesublcp[rowprecbl]);
      sizesublcpop    = &(blmatpred->sizesublcpop[rowprecbl]);
      subq            = blmatpred->subq[rowprecbl];
      bufz            = blmatpred->bufz[rowprecbl];
      newz            = blmatpred->newz[rowprecbl];
      workspace       = blmatpred->workspace[rowprecbl];

      info1 = lcp_solver_pred(adrbldiag , rhs , &rowsize , pt , &z[indicrow] , &w[indicrow] ,
                              firsttime, soltype , indic , indicop , submatlcp , submatlcpop , ipiv , sizesublcp , sizesublcpop ,
                              subq , bufz , newz , workspace);

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

    firsttime = 0;
    /* **** Criterium convergence **** */
    /*    info = filter_result_LCP_block(blmat,q,z,tolglob,pt->lcp.chat,w);*/
    info = filter_result_LCP(sizelcp, blmat->vec, q, z, tolglob, pt->lcp.chat, w);

    /*    qs   = -1;
        incx =  1;
        incy =  1;

        daxpy_( (integer *)&n , &qs , w , (integer *)&incx , ww , (integer *)&incy );

        num = dnrm2_( (integer *)&n, ww , (integer *)&incx );
        err = num*den;*/
    /* **** ********************* **** */
  }
  /*  printf("fin iteration solver bloc : %d iterations globales , %d iterations totales , %d iterations 1�e rang�\n",iter,totaliter,iterrow0);*/

  *it_end  = iter;
  *itt_end = totaliter;
  *res     = err;

  /*  free(ww);*/
  free(rhs);

  if (pt->lcp.chat > 0)
  {
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

  /*  info = filter_result_LCP_block(blmat,q,z,tolglob,pt->lcp.chat,w);*/

  return info;

}
