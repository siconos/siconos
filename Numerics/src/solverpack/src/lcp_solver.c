/*!\file lcp_solver.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).
 * Try \f$(z,w)\f$ such that:
 *
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *    0 \le z \perp Mz + q = w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 *
 * M is an ( n x n ) matrix, q , w and z n-vector. This system of equalities and inequalities
 * is solved thanks to @ref lcp solvers. The routine's call is due to the function lcp_solver.c.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "solverpack.h"
#endif

/*!\fn int lcp_solver( double *vec , double *q , int *nn , method *pt , double *z , double *w , int *it_end , double *res )
 *
 * lcp_solver is a generic interface allowing the call of one of the LCP solvers.
 *
 * \param double* vec  Unchanged parameter which contains the components of the LCP matrix with a Fortran storage.
 * \param double* q    Unchanged parameter which contains the components of the constant right hand side vector.
 * \param int* nn      Unchanged parameter which represents the dimension of the LCP problem.
 * \param method* pt   Unchanged parameter which represents the LCP structure.
 * \param double* z    Modified parameter which contains the initial value of the LCP and returns the solution of the problem.
 * \param double* w    Modified parameter which returns the complementary solution of the problem.
 * \param it_end       Modified parameter which returns the number of iterations performed by the algorithm.
 * \param res          Modified parameter which returns the final error value.
 *
 * \return integer     0 - successful\n
 *                     0 >  - otherwise (see specific solvers for more information about the log info)
 *
 * \author Nineb Sheherazade & Mathieu Renouf
 */

int lcp_solver(double *vec, double *q , int *n , method *pt , double *z , double *w , int *it_end , double *res)
{

  const char mot1[10] = "Lemke", mot2[10] = "NLGS", mot3[10] = "CPG";
  const char mot4[10] = "Latin", mot5[10] = "QP", mot6[10] = "NSQP";
  const char mot7[15] = "LexicoLemke";

  int i, info = 1;

  int     iparamLCP[5];
  double  dparamLCP[5];

  for (i = 0 ; i < 5 ; ++i) iparamLCP[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparamLCP[i] = 0.0;

  *it_end = 0;
  *res    = 0.0;

  if (strcmp(pt->lcp.name , mot1) == 0)

    lcp_lemke(vec , q , n , &pt->lcp.itermax , z ,   /* in  */
              w , it_end , res , &info);            /* out */

  /* *** LCP signature *** */

  /* **** Latin Solver **** */

  else if (strcmp(pt->lcp.name , mot4) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.k_latin;

    lcp_latin(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    *it_end = iparamLCP[2];
    *res    = dparamLCP[2];

  }
  /* **** NLGS Solver **** */

  else if (strcmp(pt->lcp.name , mot2) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.relax;

    lcp_nlgs(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    *it_end = iparamLCP[2];
    *res    = dparamLCP[2];

  }

  /* **** CPG Solver **** */

  else if (strcmp(pt->lcp.name , mot3) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;
    dparamLCP[0] = pt->lcp.tol;

    lcp_cpg(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    *it_end = iparamLCP[2];
    *res    = dparamLCP[1];

  }

  /* ***** QP Solver ***** */

  else if (strcmp(pt->lcp.name , mot5) == 0)
  {

    // We assume that the LCP matrix M is symmetric

    dparamLCP[0] = pt->lcp.tol;

    lcp_qp(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

  }

  /* **** NSQP Solver **** */

  else if (strcmp(pt->lcp.name , mot6) == 0)
  {

    // We assume that the LCP matrix M is not symmetric

    dparamLCP[0] = pt->lcp.tol;

    lcp_nsqp(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

  }
  else if (strcmp(pt->lcp.name , mot7) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;

    lcp_lexicolemke(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    *it_end = iparamLCP[2];
  }
  else printf("Warning : Unknown solver : %s\n", pt->lcp.name);

  return info;

}
