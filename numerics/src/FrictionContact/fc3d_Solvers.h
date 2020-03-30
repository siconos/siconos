/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#ifndef FRICTIONCONTACT3DSOLVERS_H
#define FRICTIONCONTACT3DSOLVERS_H

/*!\file fc3d_Solvers.h
  \brief Subroutines for the resolution of contact problems with friction (3-dimensional case).

*/

#include "FrictionContactProblem.h"
#include "SolverOptions.h"
#include "fc3d_AlartCurnier_functions.h"
#include "fc3d_projection.h"
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"
#include "fc3d_NCPGlockerFixedPoint.h"
#include "fc3d_2NCP_Glocker.h"
#include "fc3d_nonsmooth_Newton_AlartCurnier.h"
#include "fc3d_nonsmooth_Newton_FischerBurmeister.h"
#include "fc3d_unitary_enumerative.h"
#include "Friction_cst.h"
#include "fc3d_nonsmooth_Newton_natural_map.h"
#include "fc3d_local_problem_tools.h"

/** pointer to function used to call local solver */
typedef int (*SolverPtr)(FrictionContactProblem*, double*, SolverOptions *);

/** pointer to function used to update local problem */
typedef void (*UpdatePtr)(int, FrictionContactProblem*, FrictionContactProblem*, double*, SolverOptions *);

/** pointer to function used to post-processed results after a call to the (local) solver */
typedef void (*PostSolverPtr)(int, double*);

/** pointer to function used to update velocity and compute error */
typedef void (*ComputeErrorPtr)(FrictionContactProblem*, double*, double*, double, SolverOptions*, double,  double*);

/** pointer to function used to free memory for objects used in solvers */
typedef void (*FreeSolverPtr)();

/** pointer to function used to free memory for objects used in nsgs solvers */
typedef void (*FreeSolverNSGSPtr)(FrictionContactProblem*, FrictionContactProblem*, SolverOptions*  );

/** pointer to function used to call internal solver for proximal point solver */
typedef void (*internalSolverPtr)(FrictionContactProblem*, double*, double*, int *, SolverOptions *);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Non-Smooth Gauss Seidel solver for friction-contact 3D problem
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      [in] iparam[0] : Maximum iteration number

      [in] iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION (7)] : error computation method :
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL (0) : Full error computation with velocity computation
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (1) : Light error computation with incremental values on reaction verification of absolute error at the end
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT (2) : only light error computation (velocity not computed)
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE (3) :  we adapt the frequency of the full erro evaluation.

      [in] iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION(14)] : filter local solution if the local error is greater than 1.0
          SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_FALSE (0) the filter is not applied
          SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE  (1) the filter is applied

      [in] iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION(4)] : method uses overrelaxation
          SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE (0) relaxation is not used,
          SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE  (1) relaxation is used with parameter dparam[8],

      [in] iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE(5)] : shuffle the contact indices in the loop
          SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE (0) : no shuffle
          SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE (1) : shuffle only at the beginning
          SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP (2) : shuffle in each iteration

      [in] iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED(6)] : seed for the random generator in shuffling  contacts

      [out] iparam[SICONOS_IPARAM_ITER_DONE(1)] = iter number of performed iterations

      [in]  iparam[8] = error computation frequency

      [in]  dparam[SICONOS_DPARAM_TOL(0)] user tolerance on the loop
      [in]  dparam[8]  the relaxation parameter omega
      [out] dparam[SICONOS_DPARAM_RESIDU(1)]  reached error

      The internal (local) solver must set by the SolverOptions options[1]

  */
  void fc3d_nsgs(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  void fc3d_nsgs_initialize_local_solver(SolverPtr* solve, UpdatePtr* update, FreeSolverNSGSPtr* freeSolver, ComputeErrorPtr* computeError,
                                         FrictionContactProblem* problem, FrictionContactProblem* localproblem,
                                         SolverOptions * options);

  void fc3d_admm(FrictionContactProblem*  problem, double*  reaction,
                 double*  velocity,
                 int*  info, SolverOptions*  options);

  void fc3d_admm_init(FrictionContactProblem* problem, SolverOptions* options);
  void fc3d_admm_free(FrictionContactProblem* problem, SolverOptions* options);
  /** Non-Smooth Gauss Seidel in velocity solver for friction-contact 3D problem
     \param problem the friction-contact 3D problem to solve
     \param velocity global vector (n), in-out parameter
     \param reaction global vector (n), in-out parameters
     \param info return 0 if the solution is found
     \param options the solver options
  */

  void fc3d_nsgs_velocity(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** Proximal point solver for friction-contact 3D problem
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
  */
  void fc3d_proximal(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  
  /** Fixed point solver for friction-contact 3D problem based on the Tresca
  problem with fixed friction threshold
    \param problem the friction-contact 3D problem to solve
    \param velocity global vector (n), in-out parameter
    \param reaction global vector (n), in-out parameters
    \param info return 0 if the solution is found
    \param options the solver options :
  */
  void fc3d_TrescaFixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);


  /** Fixed point solver for friction-contact 3D problem based on the Panagiotopoulos
  method based on an alternative technique between the normal problem and the tangential one.
    \param problem the friction-contact 3D problem to solve
    \param velocity global vector (n), in-out parameter
    \param reaction global vector (n), in-out parameters
    \param info return 0 if the solution is found
    \param options the solver options
  */
  void fc3d_Panagiotopoulos_FixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);


  void fc3d_SOCLCP(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);
  
  /** Fixed point solver for friction-contact 3D problem based on the ACLM
    \param problem the friction-contact 3D problem to solve
    \param velocity global vector (n), in-out parameter
    \param reaction global vector (n), in-out parameters
    \param info return 0 if the solution is found
    \param options the solver options :
  */
  void fc3d_ACLMFixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);


  /** Projected Gradient on Cylinder solver for  Friction-contact 3D problem
   * \param problem the friction-contact 3D problem to solve
   *  \param velocity global vector (n), in-out parameter
   *   \param reaction global vector (n), in-out parameters
   *   \param info return 0 if the solution is found
   *   \param options the solver options :
   *   if dparam[3] >0 = rho
   *   if dparam[3] <= 0 then  a line-search is performed. iparam[2] is the maximum number of iteration is the line--search.
   *   The internal (local) solver must set by the SolverOptions options->internalsolvers.
  */
  void fc3d_ConvexQP_ProjectedGradient_Cylinder(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /**Fixed Point solver for friction-contact 3D problem based on the De Saxce Formulation
      \param problem : the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      dparam[3] : rho . if dparam[3] >0 then rho=dparam[3] otherwise a computataion of rho is assumed.
  */
  void fc3d_DeSaxceFixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /**Fixed Point Projection solver for friction-contact 3D problem based on the De Saxce Formulation
      \param problem : the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      dparam[3] : rho . if dparam[3] >0 then rho=dparam[3] otherwise a computataion of rho is assumed.
  */
  void fc3d_fixedPointProjection(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /**Fixed Point solver for friction-contact 3D problem based on the VI reformulation
      \param problem : the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
  */
  void fc3d_VI_FixedPointProjection(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);
  
  void fc3d_VI_FixedPointProjection_Cylinder(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /**Extra Gradient solver for friction-contact 3D problem based on the De Saxce Formulation
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options
  */
  void fc3d_ExtraGradient(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /**Extra Gradient solver (VI_EG) for friction-contact 3D problem based on a VI reformulation
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options
  */
  void fc3d_VI_ExtraGradient(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** Hyperplane Projection solver for friction-contact 3D problem based on the De Saxce Formulation
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options
  */
  void fc3d_HyperplaneProjection(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** solver using PATH (via GAMS) for friction-contact 3D problem based on an AVI reformulation
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options
  */
  void fc3d_AVI_gams_path(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** solver using PATHVI (via GAMS) for friction-contact 3D problem based on an AVI reformulation
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options
  */
  void fc3d_AVI_gams_pathvi(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** solver using PATH (via GAMS) for friction-contact 3D problem based on an AVI reformulation
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options
  */
  void fc3d_lcp_gams_path(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  /** solver using PATHVI (via GAMS) for friction-contact 3D problem based on an AVI reformulation
      \param problem the friction-contact 3D problem to solve
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options
  */
  void fc3d_lcp_gams_pathvi(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);


  /** Check for trivial solution in the friction-contact 3D problem
      \param problem FrictionContactProblem*  the problem
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param options the pointer to the array of options to set
      \return info  =0 if a trivial solution has been found, else = -1
  */
  int fc3d_checkTrivialCase(FrictionContactProblem* problem , double* velocity, double* reaction, SolverOptions* options);


  void fc3d_nonsmooth_Newton_AlartCurnier2(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);


  void fc3d_set_internalsolver_tolerance(FrictionContactProblem* problem,
                                         SolverOptions* options,
                                         SolverOptions* internalsolver_options,
                                         double error);


  /** \addtogroup SetSolverOptions @{
   */
  void fc3d_nsgs_set_default(SolverOptions* options);
  void fc3d_nsgs_velocity_set_default(SolverOptions* options);
  void fc3d_proximal_set_default(SolverOptions* options);
  void fc3d_tfp_set_default(SolverOptions* options);
  void fc3d_nsn_ac_set_default(SolverOptions* options);
  void fc3d_dsfp_set_default(SolverOptions* options);
  void fc3d_hp_set_default(SolverOptions* options);
  void fc3d_fpp_set_default(SolverOptions* options);
  void fc3d_eg_set_default(SolverOptions* options);
  void fc3d_nsn_fb_set_default(SolverOptions* options);
  void fc3d_aclmfp_set_default(SolverOptions* options);
  void fc3d_nsn_nm_set_default(SolverOptions* options);
  void fc3d_pfp_set_default(SolverOptions* options);
  void fc3d_admm_set_default(SolverOptions* options);
  void fc3d_onecontact_nsn_set_default(SolverOptions* options);
  void fc3d_onecontact_nsn_gp_set_default(SolverOptions* options);
  void fc3d_poc_set_default(SolverOptions* options);



  /** @} */

  
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
