/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#ifndef GLOBALROLLINGFRICTIONCONTACT3DSOLVERS_H
#define GLOBALROLLINGFRICTIONCONTACT3DSOLVERS_H

/*!\file grfc3d_Solvers.h
  \brief Subroutines for the resolution of contact problems with rolling friction (3-dimensional and 2-dimensional case).

*/

#include "GlobalRollingFrictionContactProblem.h"
#include "SolverOptions.h"
#include "Friction_cst.h"

/** pointer to function used to update velocity and compute error */
typedef void (* ComputeErrorGlobalRollingPtr)(GlobalRollingFrictionContactProblem* ,
                                      double * , double * , double* ,
                                      double , double * , int );

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Non-Smooth Gauss Seidel solver with reformulation for rolling friction-contact 3D problem
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
  void grfc3d_nsgs_wr(GlobalRollingFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, int* info,  SolverOptions* options);


  int grfc3d_checkTrivialCaseGlobal(int n, double* q, double* velocity, double* reaction, double * globalVelocity, SolverOptions* options);


  /* initialize solver (allocate memory) */
  void grfc3d_IPM_init(GlobalRollingFrictionContactProblem* problem, SolverOptions* options);


  /* deallocate memory */
  void grfc3d_IPM_free(GlobalRollingFrictionContactProblem* problem, SolverOptions* options);


  /* setup default solver parameters */
  void grfc3d_IPM_set_default(SolverOptions* options);


  /**
   * solver based on Interior Point Method (IPM) for Rolling friction-contact 3D problem based on an AVI reformulation
   * Vincent Acary, Paul Armand, Hoang Minh NGUYEN. High-accuracy computation of rolling friction contact problems. 2022.
   * https://hal.inria.fr/hal-03741048
   */
  void grfc3d_IPM(GlobalRollingFrictionContactProblem*  problem, double*  reaction,
                  double*  velocity, double*  globalVelocity,
                  int*  info, SolverOptions*  options);
  
  /* /\** \addtogroup SetSolverOptions @{ */
  /*  *\/ */
  /* void grfc3d_nsgs_sr_set_default(SolverOptions* options); */
  /* /\** @} *\/ */


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
