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
#ifndef FRICTIONCONTACT3DSOLVERS_H
#define FRICTIONCONTACT3DSOLVERS_H

#include "FrictionContact3D_AlartCurnier.h"
#include "FrictionContact3D_projection.h"
#include "FrictionContact3D_Newton.h"
#include "FrictionContact3D2NCP_Glocker.h"

/*!\file FrictionContact3D_Solvers.h
  Subroutines for the resolution of contact problems with friction (3-dimensional case).\n

  \author Mathieu Renouf, Franck Perignon.

*/

/** pointer to function used to call local solver */
typedef void (*SolverPtr)(int, int, double*, int*, double*);

/** pointer to function used to post-processed results after a call to the (local) solver */
typedef void (*PostSolverPtr)(int, double*);

/** pointer to function used to update velocity and compute error */
typedef void (*ComputeErrorPtr)(int, double*, double*, double*);

/** pointer to function used to free memory for objects used in solvers */
typedef void (*FreeSolverPtr)();

/** Structure used to send/get parameters for Friction Contact 3D problems
    \param chat       verbose mode ( 0: off, >0: on)
    \param name       name of the solver.
    \param itermax    maximum number of iterations.
    \param tol        convergence criteria value.
    \param iter       final number of iterations.
    \param err      final value of error criteria
    \param local_name int that characterize the local solver or formulation
    \param local_itermax local maximum number of iterations.
    \param local_tol   local convergence criteria value.
    \param local_iter       final number of iterations.
    \param local_err        final value of error criteria
*/
typedef struct
{
  int    chat;
  char   name[64];
  int    itermax;
  double tol;
  int    iter;
  double err;
  int    local_name;
  int    local_itermax;
  double local_tol;
  int    local_iter;
  double local_err;
} method_frictionContact_3D;


#include "FrictionContact3D_Problem.h"
#include "Numerics_Options.h"
#include "Solver_Options.h"


#ifdef __cplusplus
extern "C" {
#endif

  /** Non-Smooth Gauss Seidel solver for friction-contact 3D problem
      \param nc, number of contacts (dim of the problem n = 3*nc)
      \param M global matrix (n*n)
      \param q global vector (n)
      \param reaction global vector (n), in-out parameter
      \param velocity global vector (n), in-out parameter
      \param mu, vector of the friction coefficients (size nc)
      \param information about solver result
      \param int vector of parameters (max. iteration number ...)
      \param double vector of parameters (tolerance ...)
  */
  //  void frictionContact3D_nsgs(int, double*, double*, double*, double*, double*, int *, int *, double*);
  void frictionContact3D_nsgs(FrictionContact3D_Problem* problem, double *reaction, double *velocity, int* info, Solver_Options* options);

  /** Non-Smooth Gauss Seidel solver for friction-contact 3D problem, with sparse-block storage for M
      \param nc, number of contacts (dim of the problem n = 3*nc)
      \param M global matrix (n*n)
      \param q global vector (n)
      \param reaction global vector (n), in-out parameter
      \param velocity global vector (n), in-out parameter
      \param mu, vector of the friction coefficients (size nc)
      \param information about solver result:\n
      0 - convergence\n
      1 - iter = itermax, ie the simulation reached the maximum number of iterations allowed\n
      2 - negative diagonal term(s) in M.\n
      \param int vector of parameters (max. iteration number ...)\n
      iparam[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.\n
      iparam[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
      0 - no output\n
      1 - active screen output\n
      iparam[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.\n
      iparam[3] = local iter_max\n
      iparam[4] = iter local (output)\n
      iparam[5] = local formulation (0: Alart-Curnier, 1: Fischer-Burmeister)\n
      iparam[6] = local solver (0: projection, 1: Newton). Projection only for AC case.\n
      \param double vector of parameters (tolerance ...)\n
      dparam[0] = tol     Input unchanged parameter which represents the tolerance required.\n
      dparam[1] = error   Output modified parameter which returns the final error value.\n
      dparam[2] = local tolerance\n
      dparam[3] = Output modified parameter which returns the local error\n
  */
  void frictionContact3D_nsgs_SBS(int, SparseBlockStructuredMatrix*, double*, double*, double*, double*, int *, int *, double*);

  /** Check for trivial solution in the friction-contact 3D problem
      \param dim of the problem
      \param q global vector (n)
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameter
      \param int vector of parameters (max. iteration number ...)
      \param double vector of parameters (tolerance ...)
      \return int =0 if a trivial solution has been found, else = -1
  */
  int checkTrivialCase(int, double*, double*, double*, int*, double*);

#ifdef __cplusplus
}
#endif

#endif
