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
 *  here M is an (n x n) matrix, q, z and w n-vectors.\n
 *
 *  This system of equations and inequalities is solved thanks to @ref pfc_3D solvers.
 *  The routine's call is due to the function pfc_3D_solver.c.
 *
 *  \fn int pfc_3D_solver( int nc, double *vec, double *q ,method *pt , double *z , double *w, double *mu )
 *
 *  \brief pfc_3D_solver is a generic interface allowing the call of one of the PFC solvers.
 *
 *  \param nc   the number of contacts. The dimension of the system is 3*nc.
 *  \param vec  components of the double matrix with a fortran allocation.
 *  \param q    the components of the second member of the system.
 *  \param pt   structure
 *  \param z    the solution of the problem.
 *  \param w    the complementarity solution of the problem.
 *  \param mu   the list of friction coefficients. mu[i] corresponds to contact number i.
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
#include "NonSmoothDrivers.h"
#endif

int pfc_3D_driver(int nc, double *vec , double *q , method *pt , double *z , double *w, double *mu)
{
  /* Solver name */
  char pfckey1[10] = "NLGS", /*pfckey2[10]="CPG",*/ pfckey3[15] = "NSGS";

  int info;
  /*
   *  iparam[0] = itermax, the maximum number of iterations allowed.
   *  iparam[1] = ispeak, the output log identifiant\n
   *                       0 - no output\n
   *                       1 - active screen output\n
   *  iparam[2] = the number of iterations performed by the algorithm.
   *  iparam[3] = local value for itermax
   *  iparam[4] = local number of iterations perfomed
   *  iparam[5] = local formulation
   *  iparam[6] = local solver
   */
  int iparam[7];
  /*
   *  dparam[0] = tol     Input unchanged parameter which represents the tolerance required.
   *  dparam[1] = error   Output modified parameter which returns the final error value.
   *  dparam[2] = local tolerance
   *  dparam[3] = local error
   */
  double  dparam[4];

  clock_t t1;

  info    = -1;
  /*   printf("---------------- contact points %i -----------------\n",*n);  */
  t1 = clock();

  if (strcmp(pt->pfc_3D.name , pfckey1) == 0)
  {

    iparam[0] = pt->pfc_3D.itermax;
    iparam[1] = pt->pfc_3D.chat;
    dparam[0] = pt->pfc_3D.tol;

    pfc_3D_nlgs(nc, vec , q , z , w , mu, &info , iparam , dparam);

    pt->pfc_3D.iter = iparam[2];
    pt->pfc_3D.err  = dparam[1];

  }
  /*
  else if( strcmp( pt->pfc_3D.name , pfckey2 ) == 0 ){

    iparam[0] = pt->pfc_3D.itermax;
    iparam[1] = pt->pfc_3D.chat;
    dparam[0] = pt->pfc_3D.tol;

    pfc_3D_cpg(nc, vec , q , z , w , mu, &info , iparam , dparam );

    pt->pfc_3D.iter = iparam[2];
    pt->pfc_3D.err  = dparam[1];

  }
  */
  else if (strcmp(pt->pfc_3D.name , pfckey3) == 0)
  {

    iparam[0] = pt->pfc_3D.itermax;
    iparam[1] = pt->pfc_3D.chat;
    dparam[0] = pt->pfc_3D.tol;

    /* Local tolerance and max number of iterations are set to global ones.*/
    iparam[3] = iparam[0];
    dparam[2] = dparam[0];
    /* \todo: set them in Kernel */
    /*     iparam[3] = pt->pfc_3D.local_itermax; */
    /*     dparam[2] = pt->pfc_3D.local_tol; */

    iparam[5] = pt->pfc_3D.local_formulation;
    iparam[6] = pt->pfc_3D.local_solver;
    /* \todo: set [5] and [6] in Kernel */

    pfc_3D_nsgs(nc, vec , q , z , w , mu, &info , iparam , dparam);

    /* Get output informations: local/global number of iterations and errors. */
    pt->pfc_3D.iter = iparam[2];
    pt->pfc_3D.err  = dparam[1];
    pt->pfc_3D.local_err = dparam[3] ;
    pt->pfc_3D.local_iter = iparam[4] ;
    /** Note: at the time the last computed value of local number of iterations and local error are saved in i/dparam[4].*/
  }

  else printf("Warning : Unknown solving method : %s\n", pt->pfc_3D.name);

  /*  t2 = clock(); */
  /* printf("%.4lf seconds of processing\n", (t2-t1)/(double)CLOCKS_PER_SEC); */
  /* printf("%.4lf \n", (t2-t1)/(double)CLOCKS_PER_SEC); */

  return info;

}
