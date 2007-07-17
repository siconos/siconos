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
/*!\file lcp_solver.c

  This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
  Try \f$(z,w)\f$ such that:\n

  \f$
   \left\lbrace
    \begin{array}{l}
     w - M z = q\\
     0 \le z \perp  w \ge 0\\
    \end{array}
   \right.
  \f$

  M is an (\f$ n \times n\f$ ) matrix, q , w and z n-vector. This system of equalities and inequalities
  is solved thanks to @ref lcp solvers. The routine's call is due to the function lcp_solver.c.

 !\fn int lcp_solver( double *vec , double *q , int *n , method *pt , double *z , double *w )

  lcp_solver is a generic interface allowing the call of one of the LCP solvers.

  \param vec          On enter, a (\f$n \times n\f$)-vector of doubles which contains the components of the LCP matrix with a Fortran storage.
  \param q            On enter, a n-vector of doubles which contains the components of the constant right hand side vector.
  \param n            On enter, an integer which represents the dimension of the LCP problem.
  \param pt           On enter, a union containing the LCP structure.
  \n \n
  \param z            On return, a n-vector of doubles which contains the solution of the problem.
  \param w            On return, a n-vector of doubles which contains the complementary solution of the problem.



  \return integer
                   - 0 : successful\n
                   - >0 : otherwise (see specific solvers for more information about the log info)

  \author Nineb Sheherazade & Mathieu Renouf
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "blaslapack.h"

#ifndef MEXFLAG
#include "NSSpack.h"
#endif

int lcp_solver(double *vec, double *q , int *n , method *pt , double *z , double *w)
{

  const char lcpkey1[10] = "Lemke", lcpkey2[10] = "PGS", lcpkey3[10] = "CPG";
  const char lcpkey4[10] = "Latin", lcpkey5[10] = "QP", lcpkey6[10] = "NSQP";
  const char lcpkey7[15] = "LexicoLemke", lcpkey8[15] = "NewtonMin";
  const char lcpkey9[15] = "Latin_w", lcpkey10[15] = "NewtonFB", lcpkey11[15] = "PSOR";
  const char lcpkey12[10] = "NLGS";
  const char lcpkey13[10] = "RPGS";

  // Remark: Lemke = LexicoLemke. Only one solver is called: lexicoLemke.

  int i, j, info = 1;

  int     iparamLCP[5];
  double  dparamLCP[5];

  for (i = 0 ; i < 5 ; ++i) iparamLCP[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparamLCP[i] = 0.0;

  /*  limqpos = -1e-16 / sqrt((double) *n); */
  i = 0;
  while ((i < (*n - 1)) && (q[i] >= 0.)) i++;
  if ((i == (*n - 1)) && (q[*n - 1] >= 0.))
  {
    /* TRIVIAL CASE : q >= 0
     * z = 0 and w = q is solution of LCP(q,M)
     */
    for (j = 0 ; j < *n; j++)
    {
      z[j] = 0.0;
      w[j] = q[j];
    }
    info = 0;
    pt->lcp.iter = 0;
    pt->lcp.err  = 0.;
    if (pt->lcp.chat > 0) printf("Trivial case of LCP : positive vector q \n");
    return info;
  }

  if (strcmp(pt->lcp.name , lcpkey1) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;

    lcp_lexicolemke(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
  }
  /* **** Latin Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey4) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.k_latin;

    lcp_latin(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }

  /* **** Latin_w Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey9) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.k_latin;
    dparamLCP[3] = pt->lcp.relax;


    lcp_latin_w(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }

  /* **** PGS Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey2) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    /* dparamLCP[1] = pt->lcp.relax;*/

    lcp_pgs(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }
  /* **** NLGS Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey12) == 0)
  {

    printf("Warning: NLGS method is obsolete. Use PGS instead.\n");

  }
  /* **** SOR Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey11) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.relax;

    lcp_psor(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }
  /* **** CPG Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey3) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;

    lcp_cpg(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[1];

  }

  /* ***** QP Solver ***** */

  else if (strcmp(pt->lcp.name , lcpkey5) == 0)
  {

    /* We assume that the LCP matrix M is symmetric*/

    dparamLCP[0] = pt->lcp.tol;

    lcp_qp(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

  }

  /* **** NSQP Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey6) == 0)
  {

    /* We assume that the LCP matrix M is not symmetric*/

    dparamLCP[0] = pt->lcp.tol;

    lcp_nsqp(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

  }
  else if (strcmp(pt->lcp.name , lcpkey7) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;

    lcp_lexicolemke(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];

  }
  else if (strcmp(pt->lcp.name , lcpkey8) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;

    lcp_newton_min(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[1];

  }
  else if (strcmp(pt->lcp.name , lcpkey10) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;

    lcp_newton_FB(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[1];

  }
  /* **** RPGS (Regularized Projected Gauss-Seidel) Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey13) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.rho;
    /* dparamLCP[2] = pt->lcp.relax;*/

    lcp_rpgs(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[3];

  }

  else printf("Warning : Unknown solver : %s\n", pt->lcp.name);

  /* Checking validity of z found  */
  /*  if (info == 0) info = filter_result_LCP(*n,vec,q,z,pt->lcp.tol,pt->lcp.chat,w);*/
  if ((strcmp(pt->lcp.name , lcpkey2) != 0) && (strcmp(pt->lcp.name , lcpkey13) != 0))
    info = filter_result_LCP(*n, vec, q, z, pt->lcp.tol, pt->lcp.chat, w);

  return info;

}
