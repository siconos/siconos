/*!\file lcp_solver_block.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
 * Try \f$(z,w)\f$ such that:\n
 * \f$
 *    0 \le z \perp Mz + q = w \ge 0
 * \f$
 *
 * Here M is an (n x n) matrix , q , w and z n-vector.\n
 * This system of equalities and inequalities is solved thanks to @ref block_lcp solvers.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "solverpack.h"
#endif // MEXFLAG
#include "blaslapack.h"

/*!\fn int lcp_solver_block( int *inb , int *iid , double *vec, double *q , int *n , int *nb ,
 *                           methode *pt , double *z , double *w , int *it_end , double *res ){
 *
 * \brief solve_block_lcp is a generic interface allowing the call of one of the block LCP solvers.
 *
 * \param int* inb      Unchanged parameter which contains the number of non nul block element on the block row.
 * \param int* iid      Unchanged parameter which contains the list of active block on each row.
 * \param double* vec   Unchanged parameter which contains the components of the block matrices. Each block
 *                      matrix is stored as a fortran matrix and all matrices are stored also as a fortran allocation.
 * \param double* q     Unchanged parameter which contains the components of the right hand side.

 * \param int* n        Unchanged parameter which contains the dimension of the problem.
 * \param int* nb       Unchanged parameter which contains the dimension of the block matrices.
 * \param methode* pt   Unchanged parameter which contains a LCP structure
 *
 * \param double* z     Modified parameter which contains the initial iterate for the LCP(q,M)
 *                      and returns the solution of the problem.
 * \param double* w     Modified parameter which returns the solution of the problem.
 * \param int* it_end   Modified parameter which returns the final number of iterations or pivots
 * \param int* res      Modified parameter which returns the final value of error criteria.
 *
 * \return info         Integer identifiant for the solver result\n
 *                      0 - convergence\n
 *                      0 < no convergence (see solver for specific info value)\n
 *
 * \author Mathieu Renouf
 *
 * solve_block_lcp is a generic interface which consider LCP with a block structure. The global LCP is solved
 * as a succession of local LCP solved via lcp_solver.\n
 *
 * - list Keywords to call solvers:
 *
 *   - sub Lemke for lcp_lexicolemke
 *   - sub NLGS  for lcp_nlgs
 *   - sub CPG   for lcp_cpg
 *   - sub QP    for lcp_qp
 *   - sub NSQP  for lcp_nsqp
 *
 */

/*
 * Pointer function
 */

void (*local_solver)(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                     int *iparamLCP , double *dparamLCP) = NULL;
/*
 */

int lcp_solver_block(int *inb , int *iid , double *vec , double *q , int *dn , int *db , method *pt , double *z ,  /* in  */
                     double *w , int *it_end , int *itt_end , double *res)                                                      /* out */
{

  const char mot1[15] = "LexicoLemke", mot2[10] = "NLGS", mot3[10] = "CPG";
  const char mot4[10] = "QP"         , mot5[10] = "NSQP", mot6[10] = "NewtonMin";
  const char mot7[10] = "Latin";

  int info;
  int n, na, db2, db10;
  int i, j, k, il, iblock;
  int iter, itermax, totaliter;
  int info1;
  int incx, incy;

  double a1, b1, tol;
  double qs, err, num, den;

  double *ww, *rhs, *zl;

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

  n   = (*db) * (*dn);
  db2 = (*db) * 2;
  db10 = (*db) * 10;

  for (i = 0 ; i < 5 ; ++i) iparamLCP[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparamLCP[i] = 0.0;

  if (strcmp(pt->lcp.name , mot1) == 0)
  {
    /* Lexico Lemke */
    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;
    local_solver = &lcp_lexicolemke;
  }
  else if (strcmp(pt->lcp.name , mot2) == 0)
  {
    /* NLGS */
    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.relax;
    local_solver = &lcp_nlgs;
  }
  else if (strcmp(pt->lcp.name , mot3) == 0)
  {
    /* CPG */
    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;
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
    iparamLCP[1] = pt->lcp.iout;
    dparamLCP[0] = pt->lcp.tol;
    local_solver = &lcp_newton_min;
  }
  else if (strcmp(pt->lcp.name , mot7) == 0)
  {
    /* Latin */
    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.k_latin;
    local_solver = &lcp_latin;
  }
  else printf("Warning : Unknown solver : %s\n", pt->lcp.name);

  ww   = (double*)malloc(n * sizeof(double));
  rhs  = (double*)malloc((*db) * sizeof(double));
  zl   = (double*)malloc((*db) * sizeof(double));

  /* Check for non trivial case */

  incx = 1;
  qs = dnrm2_(&n , q , &incx);

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

  for (i = 0 ; i < n ; ++i)
  {
    ww[i] = 0.;
    w[i]  = 0.;
  }

  /* Intialization of w */

  incx = 1;
  incy = 1;
  dcopy_(&n , q , &incx , w , &incy);

  iter = 0;
  err  = 1.;

  while ((iter < itermax) && (err > tol))
  {

    ++iter;
    iblock = 0;

    incx = 1;
    incy = 1;

    dcopy_(&n , w , &incx , ww , &incy);

    for (i = 0 ; i < *dn ; ++i)
    {

      na = inb[i];

      incx = 1;
      incy = 1;

      for (j = 0 ; j < *db ; ++j) rhs[j] = q[(*db) * i + j];

      for (j = 0 ; j < na ; ++j)
      {

        k = iid[iblock] - 1;
        if (i == k)
        {
          il = iblock;
          ++iblock;
          continue;
        }
        a1 = 1.;
        b1 = 1.;
        incx = 1;
        incy = 1;

        dgemv_(&NOTRANS , db , db , &a1 , &vec[iblock * db2] , db , &z[(*db)*k] , &incx , &b1 , rhs , &incy);
        ++iblock;

      }

      /* Local LCP resolution with NLGS algorithm */

      (*local_solver)(db , &vec[il * db2] , rhs , &z[(*db)*i] , &w[(*db)*i] , &info1 , iparamLCP , dparamLCP);

      totaliter += iparamLCP[2];
    }

    /* **** Criterium convergence **** */

    qs   = -1;
    incx =  1;
    incy =  1;

    daxpy_(&n , &qs , w , &incx , ww , &incy);

    num = dnrm2_(&n, ww , &incx);
    err = num * den;

    /* **** ********************* **** */

  }

  *it_end  = iter;
  *itt_end = totaliter;
  *res     = err;

  free(zl);
  free(ww);
  free(rhs);
  local_solver = NULL;

  if (pt->lcp.iout > 0)
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
