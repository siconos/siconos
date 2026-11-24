/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef PLASTICITY2D_local_problem_tools_H
#define PLASTICITY2D_local_problem_tools_H

/*!\file

 */
#include "NumericsFwd.h"  // for PlasticityProblem

#if defined(__cplusplus)
extern "C" {
#endif

/** pointer to function used to call local solver */
typedef int (*SolverPtr)(PlasticityProblem*, double*, SolverOptions*);

/** pointer to function used to update local problem */
typedef void (*UpdatePtr)(int, PlasticityProblem*, PlasticityProblem*, double*,
                          SolverOptions*);

/** pointer to function used to post-processed results after a call to the
 * (local) solver */
typedef void (*PostSolverPtr)(int, double*);

/** pointer to function used to free memory for objects used in nsgs solvers */
typedef void (*FreeLocalPlasticitySolverPtr)(PlasticityProblem*, PlasticityProblem*,
                                             SolverOptions*);

typedef void (*CopyLocalReactionPtr)(double*, double*);

typedef void (*PerformRelaxationPtr)(double*, double*, double);

typedef double (*LightErrorSquaredPtr)(double*, double*);

typedef double (*SquaredNormPtr)(double*);

struct LocalPLASTICITY_2DProblemFunctionToolkit {
  SolverPtr local_solver;
  UpdatePtr update_local_problem;
  PostSolverPtr post_processed_local_result;
  FreeLocalPlasticitySolverPtr free_local_solver;
  CopyLocalReactionPtr copy_local_reaction;
  PerformRelaxationPtr perform_relaxation;
  LightErrorSquaredPtr light_error_squared;
  SquaredNormPtr squared_norm;
};

struct LocalPLASTICITY_2DProblemFunctionToolkit* localPLASTICITY_2DProblemFunctionToolkit_new(
    void);

void localPLASTICITY_2DProblemFunctionToolkit_display(
    struct LocalPLASTICITY_2DProblemFunctionToolkit*);

PlasticityProblem* plasticity_2d_local_problem_allocate(PlasticityProblem* problem);

void plasticity_2d_local_problem_free(PlasticityProblem* localproblem,
                                      PlasticityProblem* problem);
void plasticity_2d_local_problem_compute_q(PlasticityProblem* problem,
                                           PlasticityProblem* localproblem, double* reaction,
                                           int contact);
void plasticity_2d_local_problem_fill_M(PlasticityProblem* problem,
                                        PlasticityProblem* localproblem, int contact);

#if defined(__cplusplus)
}
#endif

#endif
