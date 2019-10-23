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
#ifndef SolverOptions_H
#define SolverOptions_H

/*!\file SolverOptions.h
  Structure used to send options (name, parameters and so on) to a specific solver-driver (mainly from Kernel to Numerics).
*/
#include <stdio.h> // for size_t
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep
#include "NumericsFwd.h"  // for SolverOptions

/** \struct Callback SolverOptions.h
Structure used to store user callbacks inside solvers
*/
typedef struct
{
  void *env; /**< general user environment */
  void (*collectStatsIteration)(void *env, int size, double*reaction,
                       double*velocity, double error, void* extra_data);/**< pointer on a function
* Its signature is: user env, problem size, reaction,
* velocity, error at end of solver iteration (when this makes sense) and an
* extra data structure */
} Callback;


/** \struct SolverOptions_ SolverOptions.h
    Structure used to send options (name, parameters and so on) to a specific solver (mainly from Kernel to Numerics).
*/
struct SolverOptions
{
  int solverId;                            /**< id number of the solver. */
  int isSet;                               /**< true(1) if the structure is ready to be used by a numerics driver. */
  int iSize;                               /**< iSize size of vector iparam */
  int * iparam;                            /**< list of solver parameters (integer type); Check solvers doc for details. */
  int dSize;                               /**< size of vector dparam */
  double * dparam;                         /**< list of solver parameters (double type); Check solvers doc for details. */
  int filterOn;                             /**< if true (1), check solution validity after the driver call. Default = 1. 
                                              For example if filterOn = 1 for a LCP, lcp_compute_error() 
                                              will be called at the end of the process). */
  size_t dWorkSize;                        /**< size of double type internal work array.*/
  double * dWork;                          /**< internal (double type) work array.*/
  size_t iWorkSize;                        /**< size of integer type internal work array.*/
  int * iWork;                          /**< internal (integer type) work array.*/
  size_t numberOfInternalSolvers;          /**< the number of internal or local 'sub-solvers' used by the solver.*/
  SolverOptions * internalSolvers;         /**< pointer to sub-solvers*/
  Callback * callback;                     /**< pointer to user-defined callback*/
  void * solverParameters;                 /**< additional parameters specific to the solver (GAMS and NewtonMethod only) */
  void * solverData;                       /**< additional data specific to the solver */
};

enum SICONOS_NUMERICS_PROBLEM_TYPE
{
  SICONOS_NUMERICS_PROBLEM_LCP = 0,
  SICONOS_NUMERICS_PROBLEM_MLCP = 1,
  SICONOS_NUMERICS_PROBLEM_EQUALITY = 2,
  SICONOS_NUMERICS_PROBLEM_FC2D = 3,
  SICONOS_NUMERICS_PROBLEM_FC3D = 4,
  SICONOS_NUMERICS_PROBLEM_NCP = 5,
  SICONOS_NUMERICS_PROBLEM_MCP = 6,
  SICONOS_NUMERICS_PROBLEM_VI = 7,
  SICONOS_NUMERICS_PROBLEM_AVI = 8,
  SICONOS_NUMERICS_PROBLEM_RELAY = 9,
};


/** Some value for iparam index */
enum SICONOS_IPARAM
{
  SICONOS_IPARAM_MAX_ITER = 0,
  SICONOS_IPARAM_ITER_DONE = 1
  ,
  SICONOS_IPARAM_PREALLOC =2
};
/** Some values for dparam index */
enum SICONOS_DPARAM
{
  SICONOS_DPARAM_TOL = 0,
  SICONOS_DPARAM_RESIDU = 1
};

/** for pivot based algorithm should be moved in a enum */
#define SICONOS_IPARAM_PIVOT_RULE 3
#define SICONOS_IPARAM_PATHSEARCH_STACKSIZE 5


extern const char* const SICONOS_NUMERICS_PROBLEM_LCP_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_MLCP_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_NCP_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_MCP_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_EQUALITY_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_FC2D_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_FC3D_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_VI_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_AVI_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_RELAY_STR;


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** screen display of solver parameters
      \param options the structure to be displayed
  */
  void solver_options_print(SolverOptions* options);

  /** Clear and free all pointer members of the structure.
   *   \param options the structure to be cleared.
   */
  void solver_options_delete(SolverOptions * options);

  /* Set all pointer fields to NULL, except iparam and dparam
   * \param options the struct to initialize
   */
  void solver_options_nullify(SolverOptions* options);

  /** fill a SolverOptions struct: set fields, allocate memory and set common
   * values
   * \param options struct to fill
   * \param solverId identity of the solver
   * \param iSize size of the iparam field (integer parameters)
   * \param dSize size of the dparam field (double parameters)
   * \param iter_max maximum number of iterations before the solver stops
   * if this does not make sense or is unwanted, give inf as value
   * \param tol tolerance for the solution.
   * if this does not make sense or is unwanted, give inf as value
   */
  void solver_options_fill(SolverOptions* options, int solverId, int iSize, int dSize, int iter_max, double tol);

  /** set parameters in SolverOption. This function should be used instead of
   * rewrittent each time a new function for setting the parameters
   * \param options the struct to set
   * \param solverId the id of the solver
   */
  void solver_options_set(SolverOptions* options, int solverId);

  /** return the id of a solver based on its name
   * \param pName the name of the solver
   * \return the id of the solver or 0 if it failed
   */
  int solver_options_name_to_id(const char * pName);

  /** return the name of a solver given its id
   * \param Id the id of the solver
   * \return the name of the solver
   */
  const char * solver_options_id_to_name(int Id);

  /** return the name of a problem type (LCP, NCP, VI, ...) based on its id
   * \param id the id of the problem
   * \return the name of the problem
   */
  const char * ns_problem_id_to_name(int id);

  /** free the solverData structure
   * \param options the structure to free
   */
  void solver_options_free_solver_specific_data(SolverOptions* options);

  /** copy SolverOptions
   * \param options_ori the structure to copy
   * \param options the output structure 
   */
  void solver_options_copy(SolverOptions* options_ori, SolverOptions* options);

  SolverOptions * solver_options_get_internal_solver(SolverOptions * options, size_t n);
  
  
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif



#endif
