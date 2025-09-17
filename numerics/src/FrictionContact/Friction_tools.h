/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#ifndef FRICTIONCONTACT3DTOOLS_H
#define FRICTIONCONTACT3DTOOLS_H

/*!\file fc3d_tools.h
  Shared subroutines for the resolution of contact problems with friction
  (3-dimensional case).
*/
#include "FrictionContactProblem.h"


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif
  
/** pointer to function used to update velocity and compute error */
typedef void (*ComputeErrorPtr)(FrictionContactProblem *, double *, double *, double,
                                SolverOptions *, double, double *);

/** pointer to function used to call internal solver for proximal point solver
 */
typedef void (*internalSolverPtr)(FrictionContactProblem *, double *, double *, int *,
                                  SolverOptions *);

/** pointer to function used to free memory for objects used in nsgs solvers */
typedef void (*FreeSolverPtr)();


void fc3d_set_internalsolver_tolerance(FrictionContactProblem *problem, SolverOptions *options,
                                       SolverOptions *internalsolver_options, double error);

/**
    Check for trivial solution in the friction-contact 3D problem

    \param problem FrictionContactProblem*  the problem
    \param velocity global vector (n), in-out parameter
    \param reaction global vector (n), in-out parameters
    \param options the pointer to the array of options to set
    \return info  =0 if a trivial solution has been found, else = -1
*/
int fc3d_checkTrivialCase(FrictionContactProblem *problem, double *velocity, double *reaction,
                          SolverOptions *options);

  

/** @} */

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
