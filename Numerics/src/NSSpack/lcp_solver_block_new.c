/* Siconos-Numerics version 1.2.0, Copyright INRIA 2005-2006.
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
/*!\file lcp_solver_block.c
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

/*!\fn int lcp_solver_block( int *inb , int *iid , double *vec , double *q , int *dn , int *db , method *pt , double *z , double *w , int *it_end , int *itt_end ,double *res )
*  lcp_solver_block is a generic interface allowing the call of one of the block LCP solvers.\n
*
* \param inb      On enter, an integer which contains the number of non nul block element on the block row.
* \param iid      On enter, a vector of integers which contains the list of active block on each row.
* \param vec      On enter, a vector of doubles which contains the components of the block matrices. Each block
*                 matrix is stored as a fortran matrix and all matrices are stored also as a fortran allocation.
* \param q        On enter, a vector of doubles which contains the components of the right hand side.

* \param dn       On enter, an integer which contains the dimension of the problem.
* \param db       On enter, an integer which contains the dimension of the block matrices.
* \param pt       On enter, a union which contains a LCP structure
*
* \param z        On enter/return, a vector of doubles which contains the initial iterate for the LCP(q,M) and returns the solution of the problem.
* \param w        On return, a vector of doubles which returns the solution of the problem.
* \param it_end   On return, an integer which returns the final number of iterations or pivots
* \param itt_end  On return, an interger which returns the total number of iterations or pivots.
* \param res      On return, a doubles which returns the final value of error criteria.
*
* \return info    Integer identifiant for the solver result\n
*                 0 : convergence\n
*                 >0 : no convergence (see solver for specific info value)\n
*
* lcp_solver_block is a generic interface which consider LCP with a block structure. The global LCP is solved
* as a succession of local LCP solved via lcp_solver.\n
*
* list Keywords to call solvers:
*
*   - Lemke    for lcp_lexicolemke
*   - NLGS     for lcp_nlgs
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
* M=\left[\begin{array}{cc|cc|cc|cc}
*          1 & 2 & 0 & 0 & 3 &-1 & 0 & 0\\
*          2 & 1 & 0 & 0 & 4 & 1 & 0 & 0\\
*          \hline
*          0 & 0 & 1 &-1 & 0 & 0 & 2 & 2\\
*          0 & 0 &-1 & 6 & 0 & 0 & 1 & 2\\
*          \hline
*          3 & 4 & 0 & 0 & 1 & 0 & 0 & 0\\
*         -1 & 3 & 0 & 0 & 0 & 2 & 0 & 0\\
*          \hline
*          0 & 0 & 2 & 1 & 0 & 0 & 2 & 2\\
*          0 & 0 & 2 & 2 & 0 & 0 & 2 & 2\\
*        \end{array}\right] \quad, q=\left[\begin{array}{c}-1\\-1\\\hline 0\\-1\\\hline 1\\0\\\hline -1\\2\end{array}\right].
* \f$
*
* then
* - the number of block is 4 (=dn) and the dimension of block is 2 (=db)
* - the vector inb[dn] is equal to [2,2,2,2]
* - the vector iid[dn*db] is equal to [1,3,2,4,1,3,2,4]
* - the vector vec contains all block matrices stored as\n
*   vec = { M11[1][1],M11[2][1],M11[1][2],M11[2][2],M12[1][1],M12[2][1],..., \n
*           M43[1][2],M43[2][2],M44[1][1],M44[2][1],M44[1][2],M44[2][2]}\n
* \author Mathieu Renouf
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "NSSpack.h"
#endif
#include "blaslapack.h"

/*
 * Pointer function
 */

void (*local_solver)(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                     int *iparamLCP , double *dparamLCP) = NULL;

/* int lcp_solver_block( int *inb , int *iid , double *vec , double *q , int *dn , int *db , method *pt , double *z , double *w , int *it_end ,
int *itt_end ,double *res ) */
int lcp_solver_block(SparseBlockStructuredMatrix *blmat, double *q, method *pt , double *z , double *w , int *it_end ,
                     int *itt_end , double *res)

{
  const char mot1[15] = "LexicoLemke", mot2[10] = "NLGS", mot3[10] = "CPG";
  const char mot4[10] = "QP"         , mot5[10] = "NSQP", mot6[10] = "NewtonMin";
  const char mot7[10] = "Latin";

  int info;
  int n, na, nbbl, nbblrow, blsizemax;
  int rowprecbl, rowcurbl, indicrow, rowsize;
  int colprecbl, colcurbl, indiccol, colsize;
  int i, j, k;
  int iter, itermax, totaliter;
  int info1;
  int incx, incy;

  double a1, b1, tol;
  double qs, err, num, den;

  double *adrcurbl, *adrbldiag;
  double *ww, *rhs;

  char NOTRANS = 'N';

  int     iparamLCP[5];
  double  dparamLCP[5];

  *it_end   = 0;
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

  for (i = 0 ; i < 5 ; ++i) iparamLCP[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparamLCP[i] = 0.0;

  if (strcmp(pt->lcp.name , mot1) == 0)
  {
    /* Lexico Lemke */
    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    local_solver = &lcp_lexicolemke;
  }
  else if (strcmp(pt->lcp.name , mot2) == 0)
  {
    /* NLGS */
    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.relax;
    local_solver = &lcp_nlgs;
  }
  else if (strcmp(pt->lcp.name , mot3) == 0)
  {
    /* CPG */
    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    local_solver = &lcp_cpg;
  }
  else if (strcmp(pt->lcp.name , mot4) == 0)
  {
    /* QP */
    dparamLCP[0] = pt->lcp.tol;
    local_solver = &lcp_qp;
  }
  else if (strcmp(pt->lcp.name , mot5) == 0)
  {
    /* NSQP */
    dparamLCP[0] = pt->lcp.tol;
    local_solver = &lcp_nsqp;
  }
  else if (strcmp(pt->lcp.name , mot6) == 0)
  {
    /* Newton Min */
    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    local_solver = &lcp_newton_min;
  }
  else if (strcmp(pt->lcp.name , mot7) == 0)
  {
    /* Latin */
    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.k_latin;
    local_solver = &lcp_latin;
  }
  else printf("Warning : Unknown solver : %s\n", pt->lcp.name);

  ww   = (double*)malloc(n * sizeof(double));
  rhs  = (double*)malloc(blsizemax * sizeof(double));

  /* Check for non trivial case */

  incx = 1;
  qs = dnrm2_((integer *)&n , q , (integer *)&incx);

  if (qs > 1e-16) den = 1.0 / qs;
  else
  {
    for (i = 0 ; i < n ; ++i)
    {
      w[i] = 0.;
      z[i] = 0.;
    }
    info = 0;
    return info;
  }


  /* Initialization of z and w */

  for (i = 0 ; i < n ; i++) z[i] = 0.;

  incx = 1;
  incy = 1;
  dcopy_((integer *)&n , q , (integer *)&incx , w , (integer *)&incy);


  /*   test on matrix structure : is there null blocks rows ? */
  rowprecbl = -1;
  for (i = 0 ; i < nbbl ; i++)
  {
    rowcurbl = blmat->RowIndex[i];
    if (rowcurbl > rowprecbl + 1)
    {
      printf(" Null blocks row in LCP matrix !!!\n");
      return 1;
    }
    rowprecbl = rowcurbl;
  }
  if (rowcurbl != nbblrow - 1)
  {
    printf(" Null blocks row in LCP matrix !!!\n");
    return 1;
  }

  iter = 0;
  err  = 1.;

  while ((iter < itermax) && (err > tol))
  {

    iter++;

    incx = 1;
    incy = 1;

    dcopy_((integer *)&n , w , (integer *)&incx , ww , (integer *)&incy);

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
            /* Local LCP resolution with NLGS algorithm */
            (*local_solver)(&rowsize , adrbldiag , rhs , &z[indicrow] , &w[indicrow] , &info1 , iparamLCP , dparamLCP);
          }
          else
          {
            printf("NULL diagonal block in LCP !!!\n");
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
      /* Local LCP resolution with NLGS algorithm */
      (*local_solver)(&rowsize , adrbldiag , rhs , &z[indicrow] , &w[indicrow] , &info1 , iparamLCP , dparamLCP);
    }
    else
    {
      printf("NULL diagonal block in LCP !!!\n");
      return 1;
    }

    totaliter += iparamLCP[2];

    /* **** Criterium convergence **** */

    qs   = -1;
    incx =  1;
    incy =  1;

    daxpy_((integer *)&n , &qs , w , (integer *)&incx , ww , (integer *)&incy);

    num = dnrm2_((integer *)&n, ww , (integer *)&incx);
    err = num * den;
    /* **** ********************* **** */
  }

  *it_end  = iter;
  *itt_end = totaliter;
  *res     = err;

  free(ww);
  free(rhs);
  local_solver = NULL;

  if (pt->lcp.chat > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of NLGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      info = 1;
    }
    else
    {
      printf(" Convergence of NLGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      info = 0;
    }
  }
  else
  {
    if (err > tol) info = 1;
    else info = 0;
  }

  return info;

}
