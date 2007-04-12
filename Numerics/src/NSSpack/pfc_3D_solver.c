/* Siconos-Numerics version 2.0.1, Copyright INRIA 2005-2006.
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
/*!\file pfc_3D_solver.c
 *
 *  This subroutine allows the primal resolution of contact problems with friction.\n
 *
 *  Try \f$(z,w)\f$ such that:\n
 *
 *  \f$
 *   \left\lbrace
 *    \begin{array}{l}
 *     M z + q = w \\
 *     0 \le z_n \perp w_n \ge 0\\
 *     -w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
 *     \end{array}
 *    \right.
 *   \f$
 *
 *  here M is an (3n x 3n) matrix, q, z and w 3n-vectors.\n
 *
 *  This system of equations and inequalities is solved thanks to @ref pfc_3D solvers.
 *  The routine's call is due to the function pfc_3D_solver.c.
 *
 *  \fn int pfc_3D_solver( double *vec , double *q ,int *n , method *pt , double *z , double *w )
 *
 *  \brief pfc_3D_solver is a generic interface allowing the call of one of the PFC solvers.
 *
 *  \param vec  components of the double matrix with a fortran allocation.
 *  \param q    the components of the second member of the system.
 *  \param 3*n  the dimension of the second member.
 *  \param pt   structure
 *  \param z    the solution of the problem.
 *  \param w    the complementarity solution of the problem.
 *
 *  \return     result (0 is successful otherwise 1).
 *
 *  \author Mathieu Renouf.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "NSSpack.h"
#endif

int pfc_3D_solver(double *vec , double *q , int *n , method *pt , double *z , double *w)
{

  char pfckey1[10] = "NLGS", pfckey2[10] = "CPG", pfckey3[15] = "NLGSNEWTON";

  int i, info;

  int     iparamLCP[5];
  double  dparamLCP[5];

  clock_t t1, t2;

  for (i = 0 ; i < 5 ; ++i) iparamLCP[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparamLCP[i] = 0.0;

  info    = -1;

  t1 = clock();

  if (strcmp(pt->pfc_3D.name , pfckey1) == 0)
  {

    iparamLCP[0] = pt->pfc_3D.itermax;
    iparamLCP[1] = pt->pfc_3D.chat;
    dparamLCP[0] = pt->pfc_3D.mu;
    dparamLCP[1] = pt->pfc_3D.tol;

    pfc_3D_nlgs(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->pfc_3D.iter = iparamLCP[2];
    pt->pfc_3D.err  = dparamLCP[2];

  }
  /*
  else if( strcmp( pt->pfc_3D.name , pfckey2 ) == 0 ){

    iparamLCP[0] = pt->pfc_3D.itermax;
    iparamLCP[1] = pt->pfc_3D.chat;
    dparamLCP[0] = pt->pfc_3D.mu;
    dparamLCP[1] = pt->pfc_3D.tol;

    pfc_3D_cpg( n , vec , q , z , w , &info , iparamLCP , dparamLCP );

    pt->pfc_3D.iter = iparamLCP[2];
    pt->pfc_3D.err  = dparamLCP[2];

  }
  */
  else if (strcmp(pt->pfc_3D.name , pfckey3) == 0)
  {

    iparamLCP[0] = pt->pfc_3D.itermax;
    iparamLCP[1] = pt->pfc_3D.chat;
    dparamLCP[0] = pt->pfc_3D.mu;
    dparamLCP[1] = pt->pfc_3D.tol;

    pfc_3D_nlgsnewton(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->pfc_3D.iter = iparamLCP[2];
    pt->pfc_3D.err  = dparamLCP[2];

  }

  else printf("Warning : Unknown solving method : %s\n", pt->pfc_3D.name);

  t2 = clock();

  //printf("%.4lf seconds of processing\n", (t2-t1)/(double)CLOCKS_PER_SEC);

  return info;

}
