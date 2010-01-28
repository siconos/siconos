/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef FRICTIONCONTACT3DSOLVERS_H
#define FRICTIONCONTACT3DSOLVERS_H

/*!\file FrictionContact3D_Solvers.h
  \brief Subroutines for the resolution of contact problems with friction (3-dimensional case).\n

*/

/*! \page FC3DSolvers Friction-Contact 3D problems Solvers

This page gives an overview of the available solvers for friction-contact (3D) problems and their required parameters.

For each solver, the input argument are:
- a FrictionContact_Problem
- the unknowns (reaction,velocity)
- info, the termination value (0: convergence, >0 problem which depends on the solver)
- a SolverOptions structure, which handles iparam and dparam

\section fc3Dnsgs Non-Smooth Gauss Seidel Solver

\bf function: frictionContact3D_nsgs()
\bf parameters:


*/

#include "FrictionContact_Problem.h"
#include "Numerics_Options.h"
#include "SolverOptions.h"
#include "FrictionContact3D_AlartCurnier.h"
#include "FrictionContact3D_projection.h"
#include "FrictionContact3D_Newton.h"
#include "FrictionContact3D_NCPGlockerFixedPoint.h"
#include "FrictionContact3D2NCP_Glocker.h"

/** pointer to function used to call local solver */
typedef void (*SolverPtr)(int, int, double*, SolverOptions *);

/** pointer to function used to post-processed results after a call to the (local) solver */
typedef void (*PostSolverPtr)(int, double*);

/** pointer to function used to update velocity and compute error */
typedef void (*ComputeErrorPtr)(FrictionContact_Problem*, double*, double*, double, double*);

/** pointer to function used to free memory for objects used in solvers */
typedef void (*FreeSolverPtr)();

/** pointer to function used to call internal solver for proximal point solver */
typedef void (*internalSolverPtr)(FrictionContact_Problem*, double*, double*, int *, SolverOptions *);

#ifdef __cplusplus
extern "C" {
#endif

  /** General interface to solvers for friction-contact 3D problem
  \param[in] , the structure which handles the Friction-Contact problem
  \param[in-out] , reaction global vector (n)
  \param[in-out] , velocity global vector (n)
  \param[in,out] options structure used to define the solver(s) and their parameters
  \param[in] general options for Numerics (verbose mode ...)
  \return result (0 if successful otherwise 1).
  */
  int frictionContact3D_driver(FrictionContact_Problem* problem, double *reaction , double *velocity, SolverOptions* options, Numerics_Options* global_options);

  /** set the default solver parameters and perform memory allocation for FrictionContact3D
      \param SolverOptions * the pointer to the options to set
      \param char * the string which identify the solver
  */
  int frictionContact3D_setDefaultSolverOptions(SolverOptions* options, char *);


  /** Non-Smooth Gauss Seidel solver for friction-contact 3D problem
      \param problem, the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      iparam[1] : error computation. 0 : Complete error computation with velocity computation 1: Light error computation with incremental values on reaction verification of absolute error at the end 2: only light error computation (velocity not computed)
      The internal (local) solver must set by the SolverOptions options[1]
  */

  void frictionContact3D_nsgs(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for NSGS
      \param SolverOptions ** the pointer to the array of options to set
  */
  int frictionContact3D_nsgs_setDefaultSolverOptions(SolverOptions* options);




  /** Non-Smooth Gauss Seidel in velocity solver for friction-contact 3D problem
     \param problem, the friction-contact 3D problem to solve
     \param velocity global vector (n), in-out parameter
     \param reaction global vector (n), in-out parameters
     \param info return 0 if the solution is found
     \param options the solver options :
     iparam[0] : Maximum iteration number
     The internal (local) solver must set by the SolverOptions options[1]
  */

  void frictionContact3D_nsgs_velocity(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for NSGSV
      \param SolverOptions ** the pointer to the array of options to set
  */
  int frictionContact3D_nsgs_velocity_setDefaultSolverOptions(SolverOptions* options);

  /** Proximal point solver for friction-contact 3D problem
      \param problem, the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      The internal (local) solver must set by the SolverOptions options[1]
  */
  void frictionContact3D_proximal(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for PROX
    \param SolverOptions ** the pointer to the array of options to set
  */
  int frictionContact3D_proximal_setDefaultSolverOptions(SolverOptions* options);

  /** Fixed point solver for friction-contact 3D problem based on the Tresca
  problem with fixed friction threshold
    \param problem, the friction-contact 3D problem to solve
    \param velocity global vector (n), in-out parameter
    \param reaction global vector (n), in-out parameters
    \param info return 0 if the solution is found
    \param options the solver options :
    iparam[0] : Maximum iteration number
    The internal (local) solver must set by the SolverOptions options[1]
  */
  void frictionContact3D_TrescaFixedPoint(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);


  /** set the default solver parameters and perform memory allocation for TFP
    \param SolverOptions ** the pointer to the array of options to set
  */
  int frictionContact3D_TrescaFixedPoint_setDefaultSolverOptions(SolverOptions* options);

  /**Fixed Point solver for friction-contact 3D problem based on the De Saxce Formulation
      \param problem, the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
  */
  void frictionContact3D_DeSaxceFixedPoint(FrictionContact_Problem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for DSFP
    \param SolverOptions ** the pointer to the array of options to set
  */
  int frictionContact3D_DeSaxceFixedPoint_setDefaultSolverOptions(SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for AlartCurnierNewton
    \param SolverOptions ** the pointer to the array of options to set
  */
  int frictionContact3D_AlartCurnierNewton_setDefaultSolverOptions(SolverOptions* options);

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
