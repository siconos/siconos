/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
- a FrictionContactProblem
- the unknowns (reaction,velocity)
- info, the termination value (0: convergence, >0 problem which depends on the solver)
- a SolverOptions structure, which handles iparam and dparam

\section fc3Dnsgs Non-Smooth Gauss Seidel Solver

 function: frictionContact3D_nsgs()
 parameters:


*/

#include "FrictionContactProblem.h"
#include "NumericsOptions.h"
#include "SolverOptions.h"
#include "FrictionContact3D_AlartCurnier.h"
#include "FrictionContact3D_projection.h"
#include "FrictionContact3D_Newton.h"
#include "FrictionContact3D_NCPGlockerFixedPoint.h"
#include "FrictionContact3D2NCP_Glocker.h"
#include "FrictionContact3D_globalAlartCurnier.h"
#include "FrictionContact3D_unitary_enumerative.h"
#include "Friction_cst.h"
#include "SiconosCompat.h"


/** pointer to function used to call local solver */
typedef int (*SolverPtr)(FrictionContactProblem*, double*, SolverOptions *);

/** pointer to function used to update local problem */
typedef void (*UpdatePtr)(int, FrictionContactProblem*, FrictionContactProblem*, double*, SolverOptions *);

/** pointer to function used to post-processed results after a call to the (local) solver */
typedef void (*PostSolverPtr)(int, double*);

/** pointer to function used to update velocity and compute error */
typedef void (*ComputeErrorPtr)(FrictionContactProblem*, double*, double*, double, SolverOptions*,  double*);

/** pointer to function used to free memory for objects used in solvers */
typedef void (*FreeSolverPtr)();

/** pointer to function used to free memory for objects used in nsgs solvers */
typedef void (*FreeSolverNSGSPtr)(FrictionContactProblem*);

/** pointer to function used to call internal solver for proximal point solver */
typedef void (*internalSolverPtr)(FrictionContactProblem*, double*, double*, int *, SolverOptions *);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** General interface to solvers for friction-contact 3D problem
  \param[in] problem the structure which handles the Friction-Contact problem
  \param[in,out] reaction global vector (n)
  \param[in,out] velocity global vector (n)
  \param[in,out] options structure used to define the solver(s) and their parameters
  \param[in] global_options for Numerics (verbose mode ...)
  \return result (0 if successful otherwise 1).
  */
  int frictionContact3D_driver(FrictionContactProblem* problem, double *reaction , double *velocity, SolverOptions* options, NumericsOptions* global_options);

  /** set the default solver parameters and perform memory allocation for FrictionContact3D
      \param options the pointer to the options to set
      \param solverId the identifier of the solver
  */
  int frictionContact3D_setDefaultSolverOptions(SolverOptions* options, int solverId);


  /** Non-Smooth Gauss Seidel solver for friction-contact 3D problem
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      iparam[1] : error computation. 0 : Complete error computation with velocity computation 1: Light error computation with incremental values on reaction verification of absolute error at the end 2: only light error computation (velocity not computed)
      The internal (local) solver must set by the SolverOptions options[1]
  */

  void frictionContact3D_nsgs(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  void frictionContact3D_nsgs_fillMLocal(FrictionContactProblem * problem, FrictionContactProblem * localproblem, int contact);

  void frictionContact3D_nsgs_computeqLocal(FrictionContactProblem * problem, FrictionContactProblem * localproblem, double * reaction, int contact);


  /** set the default solver parameters and perform memory allocation for NSGS
      \param options the pointer to the array of options to set
  */
  int frictionContact3D_nsgs_setDefaultSolverOptions(SolverOptions* options);




  /** Non-Smooth Gauss Seidel in velocity solver for friction-contact 3D problem
     \param problem the friction-contact 3D problem to solve
     \param velocity global vector (n), in-out parameter
     \param reaction global vector (n), in-out parameters
     \param info return 0 if the solution is found
     \param options the solver options :
     iparam[0] : Maximum iteration number
     The internal (local) solver must set by the SolverOptions options[1]
  */

  void frictionContact3D_nsgs_velocity(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for NSGSV
      \param options the pointer to the array of options to set
  */
  int frictionContact3D_nsgs_velocity_setDefaultSolverOptions(SolverOptions* options);

  /** Proximal point solver for friction-contact 3D problem
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      The internal (local) solver must set by the SolverOptions options[1]
  */
  void frictionContact3D_proximal(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for PROX
    \param  options the pointer to the array of options to set
  */
  int frictionContact3D_proximal_setDefaultSolverOptions(SolverOptions* options);

  /** Fixed point solver for friction-contact 3D problem based on the Tresca
  problem with fixed friction threshold
    \param problem the friction-contact 3D problem to solve
    \param velocity global vector (n), in-out parameter
    \param reaction global vector (n), in-out parameters
    \param info return 0 if the solution is found
    \param options the solver options :
    iparam[0] : Maximum iteration number
    The internal (local) solver must set by the SolverOptions options[1] : possible internal solvers is NSGS.
  */
  void frictionContact3D_TrescaFixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);


  /** set the default solver parameters and perform memory allocation for TFP
   *  \param options the pointer to the array of options to set
   */
  int frictionContact3D_TrescaFixedPoint_setDefaultSolverOptions(SolverOptions* options);

  /** Projected Gradient on Cylinder solver for  Friction-contact 3D problem
   * \param problem the friction-contact 3D problem to solve
   *  \param velocity global vector (n), in-out parameter
   *   \param reaction global vector (n), in-out parameters
   *   \param info return 0 if the solution is found
   *   \param options the solver options :
   *   iparam[0] : Maximum iteration number
   *   if dparam[3] >0 = rho
   *   if dparam[3] <= 0 then  a line-search is performed. iparam[2] is the maximum number of iteration is the line--search.
   *   The internal (local) solver must set by the SolverOptions options->internalsolvers.
  */
  void frictionContact3D_ProjectedGradientOnCylinder(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for NSGS
      \param options the pointer to the array of options to set
  */
  int frictionContact3D_ProjectedGradientOnCylinder_setDefaultSolverOptions(SolverOptions* options);

  /**Fixed Point solver for friction-contact 3D problem based on the De Saxce Formulation
      \param problem : the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      dparam[3] : rho . if dparam[3] >0 then rho=dparam[3] otherwise a computataion of rho is assumed.
  */
  void frictionContact3D_DeSaxceFixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for DSFP
    \param options the pointer to the array of options to set
  */
  int frictionContact3D_DeSaxceFixedPoint_setDefaultSolverOptions(SolverOptions* options);

  /**Extra Gradient solver for friction-contact 3D problem based on the De Saxce Formulation
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      dparam[3] : rho >0
  */
  void frictionContact3D_ExtraGradient(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for EG
    \param options the pointer to the array of options to set
  */
  int frictionContact3D_ExtraGradient_setDefaultSolverOptions(SolverOptions* options);

  /** Hyperplane Projection solver for friction-contact 3D problem based on the De Saxce Formulation
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      dparam[3] : rho >0
  */
  void frictionContact3D_HyperplaneProjection(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for EG
    \param options the pointer to the array of options to set
  */
  int frictionContact3D_HyperplaneProjection_setDefaultSolverOptions(SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for AlartCurnierNewton
    \param options the pointer to the array of options to set
  */
  int frictionContact3D_AlartCurnierNewton_setDefaultSolverOptions(SolverOptions* options);

  /** Check for trivial solution in the friction-contact 3D problem
      \param problem FrictionContactProblem*  the problem
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param options the pointer to the array of options to set
      \return info  =0 if a trivial solution has been found, else = -1
  */
  int checkTrivialCase(FrictionContactProblem* problem , double* velocity, double* reaction, SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
