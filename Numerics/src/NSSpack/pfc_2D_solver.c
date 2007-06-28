/* Siconos-Numerics version 2.1.0, Copyright INRIA 2005-2006.
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
/*!\file pfc_2D_solver.c
 *
 * This subroutine allows the primal resolution of contact problems with friction in the 2D case (PFC_2D).
 *
 * Try \f$(z,w)\f$ such that:\n
 *
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *    w - M z = q \\
 *    0 \le z_n \perp w_n \ge 0\\
 *    -w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
 *    \end{array}
 *   \right.
 *  \f$
 *
 * here M is an (n \f$\times\f$n)-matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.
 *
 * This system of equations and inequalities is solved thanks to @ref pfc_2D solvers.
 * The routine's call is due to the function pfc_2D_solver.c.\n\n
 *
 * !\fn int pfc_2D_solver( double *vec , double *q ,int *n , method *pt , double *z , double *w )
 *
 *
 *
 *  \param vec  On enter, a (n \f$\times\f$n)-vector of doubles which contains the components of the double matrix with a fortran allocation.
 *  \param q    On enter, a n-vector of doubles containing the components of the second member of the system.
 *  \param n    On enter, an integer, the dimension of the second member.
 *  \param pt   On enter, a union (::method) containing the PFC_2D structure.
 *   \n \n
 *  \param z    On return, a n-vector of doubles containing the solution of the problem.
 *  \param w    On return, a n-vector of doubles containing the solution of the problem.
 *
 *  \return     integer
 *                       - 0: successful,
 *                       - otherwise (see specific solvers for more information about the
 *                           termination reason).
 *
 * \author Nineb Sheherazade.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "NSSpack.h"
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

  if (strcmp(pt->pfc_2D.name , pfckey1) == 0)
  {

    iparamLCP[0] = pt->pfc_2D.itermax;
    iparamLCP[1] = pt->pfc_2D.chat;
    dparamLCP[0] = pt->pfc_2D.mu;
    dparamLCP[1] = pt->pfc_2D.tol;

    pfc_2D_nlgs(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->pfc_2D.iter = iparamLCP[2];
    pt->pfc_2D.err  = dparamLCP[2];

  }
  else if (strcmp(pt->pfc_2D.name , pfckey2) == 0)
  {

    iparamLCP[0] = pt->pfc_2D.itermax;
    iparamLCP[1] = pt->pfc_2D.chat;
    dparamLCP[0] = pt->pfc_2D.mu;
    dparamLCP[1] = pt->pfc_2D.tol;

    pfc_2D_cpg(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->pfc_2D.iter = iparamLCP[2];
    pt->pfc_2D.err  = dparamLCP[2];

  }
  else if (strcmp(pt->pfc_2D.name , pfckey3) == 0)
  {

    iparamLCP[0] = pt->pfc_2D.itermax;
    iparamLCP[1] = pt->pfc_2D.chat;
    dparamLCP[0] = pt->pfc_2D.mu;
    dparamLCP[1] = pt->pfc_2D.tol;
    dparamLCP[2] = pt->pfc_2D.k_latin;

    pfc_2D_latin(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->pfc_2D.iter = iparamLCP[2];
    pt->pfc_2D.err  = dparamLCP[3];

  }
  else printf("Warning : Unknown solving method : %s\n", pt->pfc_2D.name);

  t2 = clock();

  return info;

}
