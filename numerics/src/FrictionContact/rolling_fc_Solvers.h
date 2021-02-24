/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#ifndef ROLLINGFRICTIONCONTACT3DSOLVERS_H
#define ROLLINGFRICTIONCONTACT3DSOLVERS_H

/*!\file fc3d_Solvers.h
  \brief Subroutines for the resolution of contact problems with rolling friction (3-dimensional and 2-dimensional case).

*/

#include "RollingFrictionContactProblem.h"
#include "SolverOptions.h"
#include "rolling_fc3d_projection.h"
#include "rolling_fc3d_local_problem_tools.h"
#include "Friction_cst.h"

/** pointer to function used to call local solver */
typedef int (*RollingSolverPtr)(RollingFrictionContactProblem*, double*, SolverOptions *);

/** pointer to function used to update local problem */
typedef void (*RollingUpdatePtr)(int, RollingFrictionContactProblem*, RollingFrictionContactProblem*, double*, SolverOptions *);

/** pointer to function used to post-processed results after a call to the (local) solver */
typedef void (*RollingPostSolverPtr)(int, double*);

/** pointer to function used to update velocity and compute error */
typedef void (*RollingComputeErrorPtr)(RollingFrictionContactProblem*, double*, double*, double, SolverOptions*, double,  double*);

/** pointer to function used to free memory for objects used in solvers */
typedef void (*RollingFreeSolverPtr)();

/** pointer to function used to free memory for objects used in nsgs solvers */
typedef void (*RollingFreeSolverNSGSPtr)(RollingFrictionContactProblem*, RollingFrictionContactProblem*, SolverOptions*  );

/** pointer to function used to call internal solver for proximal point solver */
typedef void (*internalRollingSolverPtr)(RollingFrictionContactProblem*, double*, double*, int *, SolverOptions *);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif


  /** Non-Smooth Gauss Seidel solver for Rolling friction-contact 3D problem
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

  void rolling_fc3d_nsgs(RollingFrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  void rolling_fc3d_nsgs_initialize_local_solver(RollingSolverPtr* solve,
                                                 RollingUpdatePtr* update,
                                                 RollingFreeSolverNSGSPtr* freeSolver,
                                                 RollingComputeErrorPtr* computeError,
                                                 RollingFrictionContactProblem* problem,
                                                 RollingFrictionContactProblem* localproblem,
                                                 SolverOptions * options);

  /** Check for trivial solution in the rolling friction-contact 3D problem
      \param problem FrictionContactProblem*  the problem
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param options the pointer to the array of options to set
      \return info  =0 if a trivial solution has been found, else = -1
  */

  int rolling_fc3d_checkTrivialCase(RollingFrictionContactProblem* problem , double* velocity, double* reaction, SolverOptions* options);

  void rolling_fc3d_set_internalsolver_tolerance(RollingFrictionContactProblem* problem,
                                                 SolverOptions* options,
                                                 SolverOptions* internalsolver_options,
                                                 double error);


  /** \addtogroup SetSolverOptions @{
   */
  void rfc3d_nsgs_set_default(SolverOptions* options);
  void rfc3d_poc_withLocalIteration_set_default(SolverOptions* options);
  void rfc3d_poc_set_default(SolverOptions* options);


  void rolling_fc3d_admm(RollingFrictionContactProblem*  problem, double*  reaction,
                         double*  velocity,
                         int*  info, SolverOptions*  options);

  void rolling_fc3d_admm_init(RollingFrictionContactProblem* problem, SolverOptions* options);
  void rolling_fc3d_admm_free(RollingFrictionContactProblem* problem, SolverOptions* options);
  void rolling_fc3d_admm_set_default(SolverOptions* options);

  /** Non-Smooth Gauss Seidel solver for Rolling friction-contact 2D problem
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

  void rolling_fc2d_nsgs(RollingFrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options);

  void  rolling_fc2d_nsgs_initialize_local_solver(RollingSolverPtr* solve,
                                                  RollingUpdatePtr* update,
                                                  RollingFreeSolverNSGSPtr* freeSolver,
                                                  RollingComputeErrorPtr* computeError,
                                                  RollingFrictionContactProblem* problem,
                                                  RollingFrictionContactProblem* localproblem,
                                                  SolverOptions * options);

  /** Check for trivial solution in the friction-contact 3D problem
      \param problem FrictionContactProblem*  the problem
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameters
      \param options the pointer to the array of options to set
      \return info  =0 if a trivial solution has been found, else = -1
  */

  int rolling_fc2d_checkTrivialCase(RollingFrictionContactProblem* problem , double* velocity, double* reaction, SolverOptions* options);

  void rolling_fc2d_set_internalsolver_tolerance(RollingFrictionContactProblem* problem,
                                                 SolverOptions* options,
                                                 SolverOptions* internalsolver_options,
                                                 double error);


  /** \addtogroup SetSolverOptions @{
   */
  void rfc2d_nsgs_set_default(SolverOptions* options);
  void rfc2d_poc_withLocalIteration_set_default(SolverOptions* options);
  void rfc2d_poc_set_default(SolverOptions* options);

  /** @} */


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
