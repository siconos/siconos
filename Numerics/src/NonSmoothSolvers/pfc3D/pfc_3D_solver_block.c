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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif

int pfc_3D_driver_block(int n, SparseBlockStructuredMatrix *M , double *q , method *pt , double *z , double *w, double *mu)
{
  /* Solver name */
  char pfckey1[10] = "NSGS";

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

  info    = -1;

  /* NSGS Solver */
  if (strcmp(pt->pfc_3D.name , pfckey1) == 0)
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

    //    pfc_3D_nsgs_block2(n, M, q , z , w , mu, &info , iparam , dparam );

    /* Get output informations: local/global number of iterations and errors. */
    pt->pfc_3D.iter = iparam[2];
    pt->pfc_3D.err  = dparam[1];
    pt->pfc_3D.local_err = dparam[3] ;
    pt->pfc_3D.local_iter = iparam[4] ;
    /** Note: at the time the last computed value of local number of iterations and local error are saved in iparam[4]/dparam[3].*/
  }
  else
    printf("Warning : Unknown (or not yet implemented with block storage) solving method : %s\n", pt->pfc_3D.name);

  return info;

}
