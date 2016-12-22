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

#ifndef NEWTON_METHODS_H
#define NEWTON_METHODS_H

/*!\file Newton_methods.h
 * \brief Data structure and function for using Newton based solvers
 *
 * The reference used for the implementation is "Finite-Dimensional Variational
 * Inequalities and Complementarity Problems" by Facchinei and Pang.
 *
 * More precisely, the function newton_LSA() is algorithm VFBLSA.
 *
 * \author Olivier Huber
 */

#include "SolverOptions.h"
#include "SiconosConfig.h"
#include "NumericsMatrix.h"

#include <stdbool.h>

typedef void (*compute_F_ptr) (void* data_opaque, double* z, double* F);
typedef void (*compute_F_merit_ptr) (void* data_opaque, double* z, double* F, double* F_merit);

/** \struct functions_LSA Newton_methods.h
 * Struct holding the necessary pointers to functions needed by the
 * newton_LSA() procedure.
 */
typedef struct {
  compute_F_ptr compute_F; /**< function to evaluate w = F(z) */
  compute_F_merit_ptr compute_F_merit; /**< function to evaluate F_merit(z) (e.g. F_FB, F_{min}, ...) */
  void (*compute_H)(void* data_opaque, double* z, double* w, double* workV1, double* workV2, NumericsMatrix* H); /**< function to get an element H of T */
  void (*compute_error)(void* data_opaque, double* z, double* w, double* nabla_theta, double tol, double* err); /**< function to compute the error */
  void (*compute_RHS_desc)(void* data_opaque, double* z, double* w, double* F_desc); /**< function to evaluate F_desc(z) (e.g. F_FB, F_{min}, ...), optional */
  void (*compute_H_desc)(void* data_opaque, double* z, double* w, double* workV1, double* workV2, NumericsMatrix* H_desc); /**< function to get an element H_desc of T_desc, optional */
  int (*compute_descent_direction)(void* data_opaque, double* z, double* w, double* descent_dir, SolverOptions* options); /**< function to get the descent direction, used for instance in the Newton-Josephy method */
  void (*compute_JacTheta_merit)(void* data_opaque, double* z, double* w, double* F_merit, double* workV, double* JacThetaF_merit, SolverOptions* options); /**< function to get the descent direction, used for instance in the Newton-Josephy method */
  void* (*get_set_from_problem_data)(void* problem); /**< Function returning the set description from the  */
} functions_LSA;

// id of the stat structure 
#define NEWTON_STATS_ITERATION 1

/** \struct newton_stats Newton_methods.h */
typedef struct {
  int id; /**< id of this structure */
  double merit_value; /**< value of the merit function at the end of the iteration */
  double alpha; /**< value of the LS parameter */
  unsigned int status; /**< status of this newton iteration */
} newton_stats;

/** \struct newton_LSA_param Newton_methods.h*/
typedef struct {
  double p; /**<  p value for the acceptance test of the direction solution of the linear system */
  double sigma; /**< ratio for the decrease in norm of the C-function (\f$gamma'\f$ in VFBLSA)*/
  double rho; /**< coefficient for the direction check*/
  bool keep_H; /**< keep the matrix H untouched. Only used in the dense case, where a copy of the matrix is factorized */
  bool check_dir_quality; /**< Check the quality of the descent direction (Eqn 9.1.6 p. 805 in Facchinei & Pamg)*/
} newton_LSA_param;

/** \struct newton_LSA_data Newton_methods.h*/
typedef struct {
  NumericsMatrix* H; /**< matrix */
} newton_LSA_data;

// status of the newton step
#define NEWTON_STATS_NEWTON_STEP 1
#define NEWTON_STATS_DESC_DIR 2

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Newton algorithm for finding the zero of a function with a line search.
   * Mainly used for equation-based reformulation of CP or VI.
   * \param n size of the problem
   * \param z variable
   * \param w value of F(z)
   * \param info solver-specific values
   * \param data opaque problem definition
   * \param options options for this solver
   * \param functions struct of function pointers to compute F, H and the error
   */
  void newton_LSA(unsigned n, double *z, double *w, int *info, void* data, SolverOptions* options, functions_LSA* functions);

  /** Set some default values in the SolverOption when the solver is based on newton_LSA()
   * \param options the struct to modify
   */
  void newton_lsa_default_SolverOption(SolverOptions* options);

  /** Set the functions to compute F and F_merit and all the other pointers to NULL
   * \param functions structure to fill
   * \param compute_F function to compute F
   * \param merit_function function to compute F_merit
   */
  static inline void init_lsa_functions(functions_LSA* functions, compute_F_ptr compute_F, compute_F_merit_ptr merit_function)
  {
    functions->compute_F = compute_F;
    functions->compute_F_merit = merit_function;
    functions->compute_H = NULL;
    functions->compute_error = NULL;
    functions->compute_RHS_desc = NULL;
    functions->compute_H_desc = NULL;
    functions->compute_descent_direction = NULL;
    functions->compute_JacTheta_merit = NULL;
    functions->get_set_from_problem_data = NULL;
  }

  /** Set the parameters and data for newton_LSA
   * \param options the solver option
   * \param mat the */
 void set_lsa_params_data(SolverOptions* options, NumericsMatrix* mat);

 /** Check whether the solver uses the Newton_LSA framework or not
  * \param solverId the solver id
  * \return true if the solver is using newton_LSA, false otherwise
  */
 bool newton_LSA_check_solverId(int solverId);

 /** clear the solver-specific data
  * \param options the SolverOption structure
  */
 void newton_LSA_free_solverOptions(SolverOptions* options);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
