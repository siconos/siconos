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

/*!\file pr_solver.c
 *
 * This subroutine allows the primal resolution of relay problems (PR).\n
 *
 * Try \f$(z,w)\f$ such that:\n
 *
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *     w - M z = q\\
 *     -w \in \partial\psi_{[-b, a]}(z)\\
 *   \end{array}
 *  \right.
 * \f$

 * here M is an (\f$ nn\times nn\f$)-matrix, q an nn-dimensional vector, z an nn-dimensional  vector
 *and w an nn-dimensional vector.\n
 *This system of equations and inequalities is solved thanks to @ref pr solvers.
 *The routine's call is due to the function pr_solver.c.

 * \fn int pr_solver(double *vec,double *q,int *nn, method *pt,double *z,double *w)
 *
 *
 * pr_solver is a generic interface allowing the call of one of the PR solvers.
 *
 *
 *  \param vec       On enter, a (\f$nn \times nn\f$)-vector of doubles containing the components of the  matrix with a fortran90 allocation.
 *  \param q         On enter, a nn-vector of doubles containing the components of the second member of the system.
 *  \param nn        On enter, an integer which represents the dimension of the second member.
 *  \param pt        On enter, a union (::method) containing the PR structure.
 *  \n \n
 *  \param z         On return, a nn-vector of doubles which contains the solution of the problem.
 *  \param w         On return, a nn-vector of doubles which contains the solution of the problem.
 *
 *   \return integer : the termination reason\n
 *                    - 0 : successful\n
 *                    - otherwise : see specific solvers for more information about the log info.
 *
 *   \author Nineb Sheherazade.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MEXFLAG
#include "NSSpack.h"
#endif
#include <time.h>


int pr_solver(double *vec, double *q, int *nn, method *pt, double *z, double *w)
{


  int info = -1, it_end;

  char prkey1[10] = "NLGS", prkey2[10] = "Latin";

  double res;

  clock_t t1, t2;


  t1 = clock();


  if (strcmp(pt->pr.name , prkey1) == 0)
  {
    pr_nlgs(vec, q, nn, pt->pr.a, pt->pr.b, & pt->pr.itermax, & pt->pr.tol, &pt->pr.chat, z, w, &it_end, &res, &info);

    pt->pr.err = res;
    pt->pr.iter = it_end;

  }
  else if (strcmp(pt->pr.name, prkey2) == 0)
  {
    pr_latin(vec, q, nn, &pt->pr.k_latin, pt->pr.a, pt->pr.b, &pt->pr.itermax, &pt->pr.tol, &pt->pr.chat, z, w, &it_end, &res, &info);

    pt->pr.err = res;
    pt->pr.iter = it_end;
  }

  else printf("Warning : Unknown solving method : %s\n", pt->pr.name);

  t2 = clock();


  printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

  return info;
}
