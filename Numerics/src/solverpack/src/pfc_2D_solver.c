/*!\file pfc_2D_solver.c
 *
 * This subroutine allows the primal resolution of contact problems with friction.
 *
 * Try \f$(z,w)\f$ such that:\n
 *
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *    M z + q = w \\
 *    0 \le z_n \perp w_n \ge 0\\
 *    -w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
 *    \end{array}
 *   \right.
 *  \f$
 *
 * here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.
 *
 * This system of equations and inequalities is solved thanks to @ref pfc solvers.
 * The routine's call is due to the function pfc_2D_solver.c.
 *
 * \fn int pfc_2D_solver( double *vec , double *q ,int *n , method *pt , double *z , double *w )
 *
 *  pfc_2D_solver is a generic interface allowing the call of one of the PFC solvers.
 *
 *  \param vec  components of the double matrix with a fortran allocation.
 *  \param q    the components of the second member of the system.
 *  \param nn   the dimension of the second member.
 *  \param pt   structure
 *  \param z    the solution of the problem.
 *  \param w    the complementarity solution of the problem.
 *
 *  \return     result (0 is successful otherwise 1).
 *
 * \author Nineb Sheherazade & Mathieu Renouf.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "solverpack.h"
#endif

int pfc_2D_solver(double *vec , double *q , int *n , method *pt , double *z , double *w)
{

  char pfckey1[10] = "NLGS", pfckey2[10] = "CPG", pfckey3[10] = "Latin";

  int i, info;

  int     iparamLCP[5];
  double  dparamLCP[5];

  clock_t t1, t2;

  for (i = 0 ; i < 5 ; ++i) iparamLCP[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparamLCP[i] = 0.0;

  info    = -1;

  t1 = clock();

  if (strcmp(pt->pfc.name , pfckey1) == 0)
  {

    iparamLCP[0] = pt->pfc.itermax;
    iparamLCP[1] = pt->pfc.iout;
    dparamLCP[0] = pt->pfc.mu;
    dparamLCP[1] = pt->pfc.tol;

    pfc_2D_nlgs(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->pfc.iter = iparamLCP[2];
    pt->pfc.err  = dparamLCP[2];

  }
  else if (strcmp(pt->pfc.name , pfckey2) == 0)
  {

    iparamLCP[0] = pt->pfc.itermax;
    iparamLCP[1] = pt->pfc.iout;
    dparamLCP[0] = pt->pfc.mu;
    dparamLCP[1] = pt->pfc.tol;

    pfc_2D_cpg(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->pfc.iter = iparamLCP[2];
    pt->pfc.err  = dparamLCP[2];

  }
  else if (strcmp(pt->pfc.name , pfckey3) == 0)
  {

    iparamLCP[0] = pt->pfc.itermax;
    iparamLCP[1] = pt->pfc.iout;
    dparamLCP[0] = pt->pfc.mu;
    dparamLCP[1] = pt->pfc.tol;
    dparamLCP[2] = pt->pfc.k_latin;

    pfc_2D_latin(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->pfc.iter = iparamLCP[2];
    pt->pfc.err  = dparamLCP[3];

  }
  else printf("Warning : Unknown solving method : %s\n", pt->pfc.name);

  t2 = clock();

  printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

  return info;

}
