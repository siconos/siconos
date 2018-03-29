/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#ifndef GLOBALFRICTIONCONTACT3DSOLVERS_H
#define GLOBALFRICTIONCONTACT3DSOLVERS_H

/*!\file gfc3d_Solvers.h
  Subroutines for the resolution of contact problems with friction (3-dimensional case).\n

  \author Vincent Acary

*/

#include "GlobalFrictionContactProblem.h"
#include "SolverOptions.h"
#include "Friction_cst.h"
#include "gfc3d_nonsmooth_Newton_AlartCurnier.h"

typedef void (*SolverGlobalPtr)(int, int, double*, int*, double*);
typedef void (*PostSolverGlobalPtr)(int, double*);
typedef void (*ComputeErrorGlobalPtr)(GlobalFrictionContactProblem*, double*, double*, double *, double, SolverOptions*, double, double*);
typedef void (*FreeSolverGlobalPtr)(GlobalFrictionContactProblem*);




#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** set the default solver parameters and perform memory allocation for gfc3d
      \param options the pointer to the array of options to set
      \param solverId int identifier of the solver
  */
  int gfc3d_setDefaultSolverOptions(SolverOptions* options, int solverId);

  void gfc3d_set_internalsolver_tolerance(GlobalFrictionContactProblem* problem,
                                          SolverOptions* options,
                                          SolverOptions* internalsolver_options,
                                          double error);

  
  /** Check for trivial solution in the friction-contact 3D problem
       \param dim of the problem
       \param q global vector (n)
       \param velocity global vector (n), in-out parameter
       \param reaction global vector (n), in-out parameters
       \param globalVelocity the velocity in global coordinates
       \param options the pointer to the array of options to set
       \return int =0 if a trivial solution has been found, else = -1
   */
  int gfc3d_checkTrivialCaseGlobal(int dim, double* q, double* velocity, double*reaction, double* globalVelocity, SolverOptions* options);

  /** Non-Smooth Gauss Seidel solver with reformulation for friction-contact 3D problem
      \param problem the friction-contact 3D problem to solve
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
  void gfc3d_nsgs_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int* info,  SolverOptions* options);

  int gfc3d_nsgs_wr_setDefaultSolverOptions(SolverOptions* options);

  int gfc3d_nonsmooth_Newton_AlartCurnier_wr_setDefaultSolverOptions(SolverOptions* options);

  void  gfc3d_nonsmooth_Newton_AlartCurnier_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int *info, SolverOptions* options);

  /** Proximal point solver with reformulation for friction-contact 3D problem
    \param problem the friction-contact 3D problem to solve
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
  void gfc3d_proximal_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int* info,  SolverOptions* options);

  int gfc3d_proximal_wr_setDefaultSolverOptions(SolverOptions* options);

  /** Fixed Point iteration on De Saxe formulation solver with reformulation for friction-contact 3D problem
     \param problem the friction-contact 3D problem to solve
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
  void gfc3d_DeSaxceFixedPoint_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int* info,  SolverOptions* options);

  int gfc3d_DeSaxceFixedPoint_setDefaultSolverOptions(SolverOptions* options);

  /** Fied Point iteration on Tresca Friction Cylinder with reformulation for friction-contact 3D problem
     \param problem the friction-contact 3D problem to solve
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
  void gfc3d_TrescaFixedPoint_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int* info,  SolverOptions* options);

  int gfc3d_TrescaFixedPoint_setDefaultSolverOptions(SolverOptions* options);

  /**  Non-Smooth Gauss Seidel solver  for friction-contact 3D problem with iteration on velocities
        \param problem the friction-contact 3D problem to solve
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
  void  gfc3d_nsgs_velocity_wr(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int *info, SolverOptions* options);

  int gfc3d_nsgs_velocity_wr_setDefaultSolverOptions(SolverOptions* options);

  /** Non-Smooth Gauss Seidel solver  for friction-contact 3D problem
        \param problem the friction-contact 3D problem to solve
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
  void gfc3d_nsgs(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int* info, SolverOptions* options);

   /** Solver based on the fixed-point iteration proposed by Cadoux for friction-contact 3D problem
        \param problem the friction-contact 3D problem to solve
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
    */
  void gfc3d_ACLMFixedPoint(GlobalFrictionContactProblem*  problem, double*  reaction, double*  velocity,
                            double*  globalVelocity, int*  info, SolverOptions* options);
  
  int gfc3d_ACLMFixedPoint_setDefaultSolverOptions(SolverOptions* options);

  /** solver using PATH (via GAMS) for friction-contact 3D problem based on an AVI reformulation
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options
  */
  void gfc3d_AVI_gams_path(GlobalFrictionContactProblem* problem, double *reaction,
                           double *velocity, int* info, SolverOptions* options);

  /** solver using PATHVI (via GAMS) for friction-contact 3D problem based on an AVI reformulation
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options
  */
  void gfc3d_AVI_gams_pathvi(GlobalFrictionContactProblem* problem, double *reaction,
                             double *velocity, int* info, SolverOptions* options);


  void gfc3d_nonsmooth_Newton_AlartCurnier(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double *globalVelocity, int *info, SolverOptions* options);

  void gfc3d_VI_ExtraGradient(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int* info, SolverOptions* options);
  
  int gfc3d_VI_ExtraGradient_setDefaultSolverOptions(SolverOptions* options);
  void gfc3d_VI_FixedPointProjection(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, double* globalVelocity, int* info, SolverOptions* options);
  
  int gfc3d_VI_FixedPointProjection_setDefaultSolverOptions(SolverOptions* options);

  


  void gfc3d_ADMM(GlobalFrictionContactProblem*  problem, double*  reaction,
                  double*  velocity, double*  globalVelocity,
                  int*  info, SolverOptions*  options);

  void gfc3d_ADMM_init(GlobalFrictionContactProblem* problem, SolverOptions* options);
  
  void gfc3d_ADMM_free(GlobalFrictionContactProblem* problem, SolverOptions* options);

  int gfc3d_ADMM_setDefaultSolverOptions(SolverOptions* options);
  
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
