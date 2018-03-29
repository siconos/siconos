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
#ifndef SOCLCPSOLVERS_H
#define SOCLCPSOLVERS_H

/*!\file SOCLCP_Solvers.h
  \brief Subroutines for the resolution of Second Order Cone Linear Complementarity Problem (SOCLCP).\n

*/

#include "SecondOrderConeLinearComplementarityProblem.h"
#include "SolverOptions.h"

/* #include "soclcp_projection.h" */
/* #include "soclcp_Newton.h" */
/* #include "soclcp_AlartCurnier.h" */
/* #include "soclcp_NCPGlockerFixedPoint.h" */
/* #include "soclcp2NCP_Glocker.h" */
/* #include "soclcp_nonsmooth_Newton_AlartCurnier.h" */
/* #include "soclcp_nonsmooth_Newton_FischerBurmeister.h" */
/* #include "soclcp_unitary_enumerative.h" */

#include "SOCLCP_cst.h"
#include "SiconosCompat.h"


/** pointer to function used to call local solver */
typedef int (*Solver_soclcp_Ptr)(SecondOrderConeLinearComplementarityProblem*, double*, SolverOptions *);

/** pointer to function used to update local problem */
typedef void (*Update_soclcp_Ptr)(int, SecondOrderConeLinearComplementarityProblem*, SecondOrderConeLinearComplementarityProblem*, double*, SolverOptions *);

/** pointer to function used to post-processed results after a call to the (local) solver */
typedef void (*PostSolver_soclcp_Ptr)(int, double*);

/** pointer to function used to update v and compute error */
typedef void (*ComputeError_soclcp_Ptr)(SecondOrderConeLinearComplementarityProblem*, double*, double*, double, SolverOptions*,  double*);

/** pointer to function used to free memory for objects used in solvers */
typedef void (*FreeSolver_soclcp_Ptr)();

/** pointer to function used to free memory for objects used in nsgs solvers */
typedef void (*FreeSolverNSGS_soclcp_Ptr)(SecondOrderConeLinearComplementarityProblem*, SecondOrderConeLinearComplementarityProblem*, SolverOptions*);

/** pointer to function used to call internal solver for proximal point solver */
typedef void (*internalSolver_soclcp_Ptr)(SecondOrderConeLinearComplementarityProblem*, double*, double*, int *, SolverOptions *);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

/** set the default solver parameters and perform memory allocation for soclcp
    \param options the pointer to the options to set
    \param solverId the identifier of the solver
*/
int soclcp_setDefaultSolverOptions(SolverOptions* options, int solverId);

/** Non-Smooth Gauss Seidel solver for SOCLCP problem
    \param problem the SOCLCP problem to solve
    \param v global vector (n), in-out parameter
    \param r global vector (n), in-out parameters
    \param info return 0 if the solution is found
    \param options the solver options :
    [in] iparam[0] : Maximum iteration number
    [in] iparam[1] : error computation method :
        0 : Complete error computation with v computation
        1: Light error computation with incremental values on r verification of absolute error at the end
        2: only light error computation (v not computed)
    [out]iparam[7] = iter number of performed iterations
    [in] iparam[8] : method uses overrelaxation
    [in] iparam[9] : shuffle the contact indices in the loop

    [in]  dparam[0]  user tolerance on the loop
    [in]  dparam[8]  the relaxation parameter omega
    [out] dparam[1]  reached error

    The internal (local) solver must set by the SolverOptions options[1]

*/
  void soclcp_nsgs(SecondOrderConeLinearComplementarityProblem* problem, double *r, double *v, int* info, SolverOptions* options);

  void soclcp_nsgs_fillMLocal(SecondOrderConeLinearComplementarityProblem * problem, SecondOrderConeLinearComplementarityProblem * localproblem, int contact);

  void soclcp_nsgs_computeqLocal(SecondOrderConeLinearComplementarityProblem * problem, SecondOrderConeLinearComplementarityProblem * localproblem, double * r, int contact, SolverOptions * options);


/** set the default solver parameters and perform memory allocation for NSGS
    \param options the pointer to the array of options to set
*/
int soclcp_nsgs_setDefaultSolverOptions(SolverOptions* options);




/* /\** Non-Smooth Gauss Seidel in v solver for SOCLCP problem */
/*    \param problem the SOCLCP problem to solve */
/*    \param v global vector (n), in-out parameter */
/*    \param r global vector (n), in-out parameters */
/*    \param info return 0 if the solution is found */
/*    \param options the solver options : */
/*    iparam[0] : Maximum iteration number */
/*    The internal (local) solver must set by the SolverOptions options[1] */
/* *\/ */

/* void soclcp_nsgs_v(SecondOrderConeLinearComplementarityProblem* problem, double *r, double *v, int* info, SolverOptions* options); */

/* /\** set the default solver parameters and perform memory allocation for NSGSV */
/*     \param options the pointer to the array of options to set */
/* *\/ */
/* int soclcp_nsgs_v_setDefaultSolverOptions(SolverOptions* options); */

/* /\** Proximal point solver for SOCLCP problem */
/*     \param problem the SOCLCP problem to solve */
/*     \param v global vector (n), in-out parameter */
/*     \param r global vector (n), in-out parameters */
/*     \param info return 0 if the solution is found */
/*     \param options the solver options : */
/*     iparam[0] : Maximum iteration number */
/*     The internal (local) solver must set by the SolverOptions options[1] */
/* *\/ */
/* void soclcp_proximal(SecondOrderConeLinearComplementarityProblem* problem, double *r, double *v, int* info, SolverOptions* options); */

/* /\** set the default solver parameters and perform memory allocation for PROX */
/*   \param  options the pointer to the array of options to set */
/* *\/ */
/* int soclcp_proximal_setDefaultSolverOptions(SolverOptions* options); */

/* /\** Fixed point solver for SOCLCP problem based on the Tresca */
/* problem with fixed friction threshold */
/*   \param problem the SOCLCP problem to solve */
/*   \param v global vector (n), in-out parameter */
/*   \param r global vector (n), in-out parameters */
/*   \param info return 0 if the solution is found */
/*   \param options the solver options : */
/*   iparam[0] : Maximum iteration number */
/*   The internal (local) solver must set by the SolverOptions options[1] : possible internal solvers is NSGS. */
/* *\/ */
/* void soclcp_TrescaFixedPoint(SecondOrderConeLinearComplementarityProblem* problem, double *r, double *v, int* info, SolverOptions* options); */


/* /\** set the default solver parameters and perform memory allocation for TFP */
/*  *  \param options the pointer to the array of options to set */
/*  *\/ */
/* int soclcp_TrescaFixedPoint_setDefaultSolverOptions(SolverOptions* options); */

/* /\** Projected Gradient on Cylinder solver for  Friction-contact 3D problem */
/*  * \param problem the SOCLCP problem to solve */
/*  *  \param v global vector (n), in-out parameter */
/*  *   \param r global vector (n), in-out parameters */
/*  *   \param info return 0 if the solution is found */
/*  *   \param options the solver options : */
/*  *   iparam[0] : Maximum iteration number */
/*  *   if dparam[3] >0 = rho */
/*  *   if dparam[3] <= 0 then  a line-search is performed. iparam[2] is the maximum number of iteration is the line--search. */
/*  *   The internal (local) solver must set by the SolverOptions options->internalsolvers. */
/* *\/ */
/* void soclcp_ProjectedGradientOnCylinder(SecondOrderConeLinearComplementarityProblem* problem, double *r, double *v, int* info, SolverOptions* options); */

void soclcp_VI_FixedPointProjection(SecondOrderConeLinearComplementarityProblem* problem, double *r, double *v, int* info, SolverOptions* options);

/** set the default solver parameters and perform memory allocation for DSFP
  \param options the pointer to the array of options to set
*/
int soclcp_VI_FixedPointProjection_setDefaultSolverOptions(SolverOptions* options);

/**Extra Gradient solver (VI_EG) for SOCLCP problem based on a VI reformulation
    \param problem the SOCLCP problem to solve
    \param v global vector (n), in-out parameter
    \param r global vector (n), in-out parameters
    \param info return 0 if the solution is found
    \param options the solver options :
    iparam[0] : Maximum iteration number
    dparam[3] : rho >0
*/
void soclcp_VI_ExtraGradient(SecondOrderConeLinearComplementarityProblem* problem, double *r, double *v, int* info, SolverOptions* options);

/** set the default solver parameters and perform memory allocation for VI_EG
  \param options the pointer to the array of options to set
*/
 int soclcp_VI_ExtraGradient_setDefaultSolverOptions(SolverOptions* options);

/* /\** Hyperplane Projection solver for SOCLCP problem based on the De Saxce Formulation */
/*     \param problem the SOCLCP problem to solve */
/*     \param v global vector (n), in-out parameter */
/*     \param r global vector (n), in-out parameters */
/*     \param info return 0 if the solution is found */
/*     \param options the solver options : */
/*     iparam[0] : Maximum iteration number */
/*     dparam[3] : rho >0 */
/* *\/ */
/* void soclcp_HyperplaneProjection(SecondOrderConeLinearComplementarityProblem* problem, double *r, double *v, int* info, SolverOptions* options); */

/* /\** set the default solver parameters and perform memory allocation for EG */
/*   \param options the pointer to the array of options to set */
/* *\/ */
/* int soclcp_HyperplaneProjection_setDefaultSolverOptions(SolverOptions* options); */

/* /\** set the default solver parameters and perform memory allocation for AlartCurnierNewton */
/*   \param options the pointer to the array of options to set */
/* *\/ */
/* int soclcp_AlartCurnierNewton_setDefaultSolverOptions(SolverOptions* options); */

/** Check for trivial solution in the SOCLCP problem
    \param problem SecondOrderConeLinearComplementarityProblem*  the problem
    \param v global vector (n), in-out parameter
    \param r global vector (n), in-out parameters
    \param options the pointer to the array of options to set
    \return info  =0 if a trivial solution has been found, else = -1
*/
int soclcp_checkTrivialCase(SecondOrderConeLinearComplementarityProblem* problem , double* v, double* r, SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
