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
#ifndef PRIMALFRICTIONCONTACT3DSOLVERS_H
#define PRIMALFRICTIONCONTACT3DSOLVERS_H

/*!\file PrimalFrictionContact3D_Solvers.h
  Subroutines for the resolution of contact problems with friction (3-dimensional case).\n

  \author Vincent Acary

*/

/*! \page PrimalFC3DSolvers Primal Friction-Contact 3D problems Solvers

This page gives an overview of the available solvers for friction-contact (3D) problems and their required parameters.

For each solver, the input argument are:
- a FrictionContactProblem
- the unknowns (reaction,velocity)
- info, the termination value (0: convergence, >0 problem which depends on the solver)
- a SolverOptions structure, which handles iparam and dparam

\section pfc3Dnsgs Non-Smooth Gauss Seidel Solver

\bf function: frictionContact3D_nsgs()
\bf parameters:


*/
#include "PrimalFrictionContactProblem.h"
#include "NumericsOptions.h"
#include "SolverOptions.h"
#include "Friction_cst.h"

typedef void (*SolverPrimalPtr)(int, int, double*, int*, double*);
typedef void (*PostSolverPrimalPtr)(int, double*);
typedef void (*ComputeErrorPrimalPtr)(PrimalFrictionContactProblem*, double*, double*, double *, double, double*);
typedef void (*FreeSolverPrimalPtr)(PrimalFrictionContactProblem*);




#ifdef __cplusplus
extern "C"
{
#endif
  /** General interface to solvers for primal friction-contact 3D problem
    \param[in] , the structure which handles the Friction-Contact problem
    \param[in-out] , reaction global vector (n)
    \param[in-out] , velocity global vector (n)
    \param[in,out] options structure used to define the solver(s) and their parameters
    \param[in] general options for Numerics (verbose mode ...)
    \return result (0 if successful otherwise 1).
  */
  int primalFrictionContact3D_driver(PrimalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, SolverOptions* options, NumericsOptions* global_options);

  /** set the default solver parameters and perform memory allocation for PrimalFrictionContact3D
      \param SolverOptions ** the pointer to the array of options to set
      \param int identifier of the solver
  */
  int primalFrictionContact3D_setDefaultSolverOptions(SolverOptions* options, int solverId);

  /** Check for trivial solution in the friction-contact 3D problem
       \param dim of the problem
       \param q global vector (n)
       \param velocity global vector (n), in-out parameter
       \param reaction global vector (n), in-out parameters
       \param int vector of parameters (max. iteration number ...)
       \param double vector of parameters (tolerance ...)
       \return int =0 if a trivial solution has been found, else = -1
   */
  int checkTrivialCasePrimal(int, double*, double*, double*, double*, int*, double*);

  /** Non-Smooth Gauss Seidel solver with reformulation for friction-contact 3D problem
      \param problem, the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param globalVelocity global vector (m), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      iparam[4] : localsolver choice 0: projection on Cone, 1: Newton/AlartCurnier,  2: projection on Cone with local iteration, 2: projection on Disk  with diagonalization,
      dparam[0] : tolerance
      dparam[2] : localtolerance
      dparam[1] : (out) error
  */
  void primalFrictionContact3D_nsgs_wr(PrimalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int* info,  SolverOptions* options);

  int primalFrictionContact3D_nsgs_wr_setDefaultSolverOptions(SolverOptions* options);

  int primalFrictionContact3D_globalAlartCurnier_wr_setDefaultSolverOptions(SolverOptions* options);

  void  primalFrictionContact3D_globalAlartCurnier_wr(PrimalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int *info, SolverOptions* options);

  /** Proximal point solver with reformulation for friction-contact 3D problem
    \param problem, the friction-contact 3D problem to solve
    \param velocity global vector (n), in-out parameter
    \param reaction global vector (n), in-out parameters
    \param globalVelocity global vector (m), in-out parameters
    \param info return 0 if the solution is found
    \param options the solver options :
    iparam[0] : Maximum iteration number
    iparam[4] : localsolver choice 0: projection on Cone, 1: Newton/AlartCurnier,  2: projection on Cone with local iteration, 2: projection on Disk  with diagonalization,
    dparam[0] : tolerance
    dparam[2] : localtolerance
    dparam[1] : (out) error
  */
  void primalFrictionContact3D_proximal_wr(PrimalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int* info,  SolverOptions* options);

  int primalFrictionContact3D_proximal_wr_setDefaultSolverOptions(SolverOptions* options);

  /** Fixed Point iteration on De Saxe formulation solver with reformulation for friction-contact 3D problem
     \param problem, the friction-contact 3D problem to solve
     \param velocity global vector (n), in-out parameter
     \param reaction global vector (n), in-out parameters
     \param globalVelocity global vector (m), in-out parameters
     \param info return 0 if the solution is found
     \param options the solver options :
     iparam[0] : Maximum iteration number
     dparam[0] : tolerance
     dparam[2] : localtolerance
     dparam[1] : (out) error
  */
  void primalFrictionContact3D_DeSaxceFixedPoint_wr(PrimalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int* info,  SolverOptions* options);

  int primalFrictionContact3D_DeSaxceFixedPoint_setDefaultSolverOptions(SolverOptions* options);

  /** Fied Point iteration on Tresca Friction Cylinder with reformulation for friction-contact 3D problem
     \param problem, the friction-contact 3D problem to solve
     \param velocity global vector (n), in-out parameter
     \param reaction global vector (n), in-out parameters
     \param globalVelocity global vector (m), in-out parameters
     \param info return 0 if the solution is found
     \param options the solver options :
     iparam[0] : Maximum iteration number
     dparam[0] : tolerance
     dparam[2] : localtolerance
     dparam[1] : (out) error
  */
  void primalFrictionContact3D_TrescaFixedPoint_wr(PrimalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int* info,  SolverOptions* options);

  int primalFrictionContact3D_TrescaFixedPoint_setDefaultSolverOptions(SolverOptions* options);

  /**  Non-Smooth Gauss Seidel solver  for friction-contact 3D problem with iteration on velocities
        \param problem, the friction-contact 3D problem to solve
        \param velocity global vector (n), in-out parameter
        \param reaction global vector (n), in-out parameters
        \param globalVelocity global vector (m), in-out parameters
        \param info return 0 if the solution is found
        \param options the solver options :
        iparam[0] : Maximum iteration number
        dparam[0] : tolerance
        dparam[2] : localtolerance
        dparam[1] : (out) error
    */
  void  primalFrictionContact3D_nsgs_velocity_wr(PrimalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int *info, SolverOptions* options);

  int primalFrictionContact3D_nsgs_velocity_wr_setDefaultSolverOptions(SolverOptions* options);
  /** Non-Smooth Gauss Seidel solver  for friction-contact 3D problem
        \param problem, the friction-contact 3D problem to solve
        \param velocity global vector (n), in-out parameter
        \param reaction global vector (n), in-out parameters
        \param globalVelocity global vector (m), in-out parameters
        \param info return 0 if the solution is found
        \param options the solver options :
        iparam[0] : Maximum iteration number
        iparam[4] ; local strategy
        dparam[0] : tolerance
        dparam[2] : localtolerance
        dparam[1] : (out) error
        \todo Implement ProdTransSBM
        \todo Improve the splitting Algorithm with a smaller granularity
        \todo Use a global projection perhaps
    */
  void primalFrictionContact3D_nsgs(PrimalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int* info, SolverOptions* options);



#ifdef __cplusplus
}
#endif

#endif
