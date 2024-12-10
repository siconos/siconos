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
#ifndef MOHRCOULOMB2DSOLVERS_H
#define MOHRCOULOMB2DSOLVERS_H

/*!\file mc2d_Solvers.h
  Subroutines for the resolution of Mohr Coulomb 2D plasticity
*/

#include "MohrCoulomb2DProblem.h"
#include "Plasticity_cst.h"
#include "SolverOptions.h"

/** pointer to function used to update plastic_strain_rate and compute error */
typedef void (*mc2d_ComputeErrorPtr)(MohrCoulomb2DProblem *, double *, double *, double,
                                SolverOptions *, double, double *);

/** pointer to function used to call internal solver for proximal point solver
 */
typedef void (*mc2d_internalSolverPtr)(MohrCoulomb2DProblem *, double *, double *, int *,
                                  SolverOptions *);

/** pointer to function used to free memory for objects used in nsgs solvers */
typedef void (*FreeSolverPtr)();

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/**
    Non-Smooth Gauss Seidel solver for Mohr Coulomb 2D problem

    \param problem the Mohr Coulomb 2D problem to solve
    \param plastic_strain_rate global vector (n), in-out parameter
    \param stress global vector (n), in-out parameters
    \param info return 0 if the solution is found
    \param options the solver options :
    [in] iparam[0] : Maximum iteration number

    [in] iparam[PLASTICITY_IPARAM_ERROR_EVALUATION (7)] : error
    computation method : PLASTICITY_NSGS_ERROR_EVALUATION_FULL (0) :
    Full error computation with plastic_strain_rate computation
    PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (1) : Light
    error computation with incremental values on stress verification of
    absolute error at the end PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT
    (2) : only light error computation (plastic_strain_rate not computed)
    PLASTICITY_NSGS_ERROR_EVALUATION_ADAPTIVE (3) :  we adapt the
    frequency of the full erro evaluation.

    [in] iparam[PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION(14)] : filter
    local solution if the local error is greater than 1.0
    PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION_FALSE (0) the filter is not
    applied PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION_TRUE  (1) the filter
    is applied

    [in] iparam[PLASTICITY_NSGS_RELAXATION(4)] : method uses
    overrelaxation PLASTICITY_NSGS_RELAXATION_FALSE (0) relaxation is
    not used, PLASTICITY_NSGS_RELAXATION_TRUE  (1) relaxation is used
    with parameter dparam[8],

    [in] iparam[PLASTICITY_NSGS_SHUFFLE(5)] : shuffle the contact
    indices in the loop PLASTICITY_NSGS_SHUFFLE_FALSE (0) : no shuffle
    PLASTICITY_NSGS_SHUFFLE_TRUE (1) : shuffle only at the beginning
    PLASTICITY_NSGS_SHUFFLE_TRUE_EACH_LOOP (2) : shuffle in each
    iteration

    [in] iparam[PLASTICITY_NSGS_SHUFFLE_SEED(6)] : seed for the random
    generator in shuffling  contacts

    [out] iparam[SICONOS_IPARAM_ITER_DONE(1)] = iter number of performed
    iterations

    [in]  iparam[8] = error computation frequency

    [in]  dparam[SICONOS_DPARAM_TOL(0)] user tolerance on the loop
    [in]  dparam[8]  the relaxation parameter omega
    [out] dparam[SICONOS_DPARAM_RESIDU(1)]  reached error

    The internal (local) solver must set by the SolverOptions options[1]

*/
void mc2d_nsgs(MohrCoulomb2DProblem *problem, double *stress, double *plastic_strain_rate,
               int *info, SolverOptions *options);

/**
    Check for trivial solution in the Mohr Coulomb 2D problem

    \param problem MohrCoulomb2DProblem*  the problem
    \param plastic_strain_rate global vector (n), in-out parameter
    \param stress global vector (n), in-out parameters
    \param options the pointer to the array of options to set
    \return info  =0 if a trivial solution has been found, else = -1
*/
int mc2d_checkTrivialCase(MohrCoulomb2DProblem *problem, double *plastic_strain_rate,
                          double *stress, SolverOptions *options);

/** \addtogroup SetSolverOptions
 * @{
 */
void mc2d_nsgs_set_default(SolverOptions *options);
void mc2d_onecontact_nsn_set_default(SolverOptions *options);
void mc2d_onecontact_nsn_gp_set_default(SolverOptions *options);
void mc2d_poc_set_default(SolverOptions *options);
void mc2d_set_internalsolver_tolerance(MohrCoulomb2DProblem *problem, SolverOptions *options,
                                       SolverOptions *internalsolver_options, double error);


  
/** @} */

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
