/*!\file lcp_solver.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
 * Try \f$(z,w)\f$ such that:\n
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
 *!\fn int lcp_solver( double *vec , double *q , int *nn , method *pt , double *z , double *w )
 *
 * lcp_solver is a generic interface allowing the call of one of the LCP solvers.
 *
 * \param vec          Unchanged parameter which contains the components of the LCP matrix with a Fortran storage.
 * \param q            Unchanged parameter which contains the components of the constant right hand side vector.
 * \param nn           Unchanged parameter which represents the dimension of the LCP problem.
 * \param pt           Unchanged parameter which represents the LCP structure.
 * \param z            Modified parameter which contains the initial value of the LCP and returns the solution of the problem.
 * \param w            Modified parameter which returns the complementary solution of the problem.
 *
 * \return integer     0 - successful\n
 *                     0 >  - otherwise (see specific solvers for more information about the log info)
 *
 * \author Nineb Sheherazade & Mathieu Renouf
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "solverpack.h"
#endif

int lcp_solver(double *vec, double *q , int *n , method *pt , double *z , double *w)
{

  const char lcpkey1[10] = "Lemke", lcpkey2[10] = "NLGS", lcpkey3[10] = "CPG";
  const char lcpkey4[10] = "Latin", lcpkey5[10] = "QP", lcpkey6[10] = "NSQP";
  const char lcpkey7[15] = "LexicoLemke", lcpkey8[15] = "NewtonMin";

  int i, info = 1;

  int     iparamLCP[5];
  double  dparamLCP[5];

  for (i = 0 ; i < 5 ; ++i) iparamLCP[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparamLCP[i] = 0.0;

  if (strcmp(pt->lcp.name , lcpkey1) == 0)

    lcp_lemke(vec , q , n , &pt->lcp.itermax , z ,                   /* in  */
              w , &pt->lcp.iter , &pt->lcp.err , &info);            /* out */

  /* *** LCP signature *** */

  /* **** Latin Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey4) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.k_latin;

    lcp_latin(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }
  /* **** NLGS Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey2) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.relax;

    lcp_nlgs(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }

  /* **** CPG Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey3) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;
    dparamLCP[0] = pt->lcp.tol;

    lcp_cpg(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[1];

  }

  /* ***** QP Solver ***** */

  else if (strcmp(pt->lcp.name , lcpkey5) == 0)
  {

    // We assume that the LCP matrix M is symmetric

    dparamLCP[0] = pt->lcp.tol;

    lcp_qp(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

  }

  /* **** NSQP Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey6) == 0)
  {

    // We assume that the LCP matrix M is not symmetric

    dparamLCP[0] = pt->lcp.tol;

    lcp_nsqp(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

  }
  else if (strcmp(pt->lcp.name , lcpkey7) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;

    lcp_lexicolemke(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];

  }
  else if (strcmp(pt->lcp.name , lcpkey8) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.iout;
    dparamLCP[0] = pt->lcp.tol;

    lcp_newton_min(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[1];

  }
  else printf("Warning : Unknown solver : %s\n", pt->lcp.name);

  return info;

}
