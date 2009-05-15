/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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

/*!\file FrictionContact3D_Solvers.h
  \brief Subroutines for the resolution of contact problems with friction (3-dimensional case).\n

  \author Mathieu Renouf, Franck Perignon.

*/

/*! \page FC3DSolvers Friction-Contact 3D problems Solvers

This page gives an overview of the available solvers for friction-contact (3D) problems and their required parameters.

For each solver, the input argument are:
- a FrictionContact_Problem
- the unknowns (reaction,velocity)
- info, the termination value (0: convergence, >0 problem which depends on the solver)
- a Solver_Options structure, which handles iparam and dparam

\section fc3Dnsgs Non-Smooth Gauss Seidel Solver

\bf function: frictionContact3D_nsgs()
\bf parameters:


*/

#include "FrictionContact_Problem.h"
#include "Numerics_Options.h"
#include "Solver_Options.h"
#include "FrictionContact3D_AlartCurnier.h"
#include "FrictionContact3D_projection.h"
#include "FrictionContact3D_Newton.h"
#include "FrictionContact3D_FixedP.h"
#include "FrictionContact3D2NCP_Glocker.h"

/** pointer to function used to call local solver */
typedef void (*SolverPtr)(int, int, double*, int*, double*);

/** pointer to function used to post-processed results after a call to the (local) solver */
typedef void (*PostSolverPtr)(int, double*);

/** pointer to function used to update velocity and compute error */
typedef void (*ComputeErrorPtr)(FrictionContact_Problem*, double*, double*, double, double*);

/** pointer to function used to free memory for objects used in solvers */
typedef void (*FreeSolverPtr)();

/** pointer to function used to call internal solver for proximal point solver */
typedef void (*internalSolverPtr)(FrictionContact_Problem*, double*, double*, int *, Solver_Options *);

#ifdef __cplusplus
extern "C" {
#endif

  /** Non-Smooth Gauss Seidel solver for friction-contact 3D problem
      \param problem, the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      iparam[4] : localsolver choice 0: projection on Cone, 1: Newton/AlartCurnier,  2: projection on Cone with local iteration, 2: projection on Disk  with diagonalization,
  */

  void frictionContact3D_nsgs(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, Solver_Options* options);

  /** Non-Smooth Gauss Seidel in velocity solver for friction-contact 3D problem
     \param problem, the friction-contact 3D problem to solve
     \param velocity global vector (n), in-out parameter
     \param reaction global vector (n), in-out parameters
     \param info return 0 if the solution is found
     \param options the solver options :
     iparam[0] : Maximum iteration number
     iparam[4] : localsolver choice 0: projection on Cone, 1: Newton/AlartCurnier,  2: projection on Cone with local iteration, 2: projection on Disk  with diagonalization,
  */

  void frictionContact3D_nsgs_velocity(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, Solver_Options* options);

  /** Proximal point solver for friction-contact 3D problem
      \param problem, the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      iparam[4] : localsolver choice 0: projection on Cone, 1: Newton/AlartCurnier,  2: projection on Cone with local iteration, 2: projection on Disk  with diagonalization,
  */
  void frictionContact3D_proximal(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, Solver_Options* options);


  /** Projected Gradient solver for friction-contact 3D problem
      \param problem, the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
  */
  void frictionContact3D_projectedgradient(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, Solver_Options* options);

  /** Check for trivial solution in the friction-contact 3D problem
      \param dim of the problem
      \param q global vector (n)
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param int vector of parameters (max. iteration number ...)
      \param double vector of parameters (tolerance ...)
      \return int =0 if a trivial solution has been found, else = -1
  */
  int checkTrivialCase(int, double*, double*, double*, int*, double*);

#ifdef __cplusplus
}
#endif

#endif
