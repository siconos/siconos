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
/*!\file dr_solver.c
 *
 * This subroutine allows the dual resolution of relay problems (DR).\n
 *
 * Try \f$(z,w)\f$ such that:\n
 *  \f$
 *   \left\lbrace
 *    \begin{array}{l}
 *      w - M z = q\\
 *     -z \in \partial\psi_{[-b, a]}(w)\\
 *    \end{array}
 *   \right.
 *  \f$
 *
 * where M is an (nn \f$\times\f$nn)-matrix, q , w and z nn-vectors.\n
 *
 * This system of equations and inequalities is solved thanks to @ref dr solvers.
 * The routine's call is due to the function pr_solver.c.\n
 *
 *
 * \fn int dr_solver( double *vec , double *q ,int *nn , method *pt , double *z , double *w )
 *
 *\n
 * dr_solver is a generic interface allowing the call of one of the DR solvers.\n

 * \param vec      On enter, (nn \f$\times\f$nn)-vector of doubles containing the components of the double matrix with a fortran90 allocation.
 * \param q        On enter, a nn-vector of doubles containing the components of the second member of the system.
 * \param nn       On enter, an integer, the dimension of the second member.
 * \param pt       On enter, a union (::method) containing the DR structure.\n \n
 * \param z        On return, a nn-vector of doubles, the solution of the problem.
 * \param w        On return, a nn-vector of doubles, the solution of the problem.
 *
 *
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

int dr_solver(double *vec , double *q , int *nn , method *pt , double *z , double *w)
{

  int info = -1, it_end;

  double res;

  char drkey1[10] = "NLGS" , drkey2[10] = "Latin"/* ,drkey3[10] = "CPG" */;




  clock_t t1, t2;

  t1 = clock();

  if (strcmp(pt->dr.name , drkey2) == 0)
  {

    dr_latin(vec , q , nn , &pt->dr.k_latin , pt->dr.a , pt->dr.b , &pt->dr.itermax , &pt->dr.tol , &pt->dr.chat, z , w , &it_end , &res , &info);

    pt->dr.err = res;
    pt->dr.iter = it_end;


  }
  else if (strcmp(pt->dr.name , drkey1) == 0)
  {

    dr_nlgs(vec , q , nn , pt->dr.a , pt->dr.b , &pt->dr.itermax , &pt->dr.tol , &pt->dr.chat, z , w , &it_end , &res , &info);

    pt->dr.err = res;
    pt->dr.iter = it_end;
  }
  else printf("Warning : Unknown solving method : %s\n", pt->dr.name);

  t2 = clock();
  printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

  return info;

}
