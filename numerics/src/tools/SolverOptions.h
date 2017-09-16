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



#ifndef SolverOptions_H
#define SolverOptions_H

/*!\file SolverOptions.h
  Structure used to send options (name, parameters and so on) to a specific solver-driver (mainly from Kernel to Numerics).
  \author Franck Perignon
*/
#include "SiconosConfig.h"
#include "NumericsFwd.h"

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
  int solverId;                            /**< solverId Id of the solver (see ) */
  int isSet;                               /**< isSet int equal to false(0) if the parameters below have not been set (ie need to read default values) else true(1)*/
  int iSize;                               /**< iSize size of vector iparam */
  int * iparam;                            /**< iparam a list of int parameters (depends on each solver, see solver doc)*/
  int dSize;                               /**< dSize size of vector dparam */
  double * dparam;                         /**< dparam a list of double parameters (depends on each solver, see solver doc)*/
  int filterOn;                            /**< filterOn 1 to check solution validity after the driver call, else 0. Default = 1. (For example if
                                            * filterOn = 1 for a LCP, lcp_compute_error() will be called at the end of the process) */
  int dWorkSize;                           /**< dWorkSize size of vector iWork */
  double * dWork;                          /**< dWork is a pointer on a working memory zone (for doubles) reserved for the solver .*/
  int iWorkSize;                           /**< iWorkSize size of vector iWork */
  int * iWork;                             /**< iWork is a pointer on a working memory zone (for integers) reserved for the solver .*/
  int numberOfInternalSolvers;             /**< numberOfInternalSolvers the number of internal or local 'sub-solvers' used by the solver*/
  struct SolverOptions * internalSolvers; /**< internalSolvers pointer to sub-solvers*/
  Callback * callback;                     /**< callback a pointer to user Callback*/

  void * solverParameters;                 /**< additional parameters specific to the solver */

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
/* #define SICONOS_IPARAM_MAX_ITER 0 */
/* #define SICONOS_IPARAM_ITER_DONE 1 */
/* #define SICONOS_IPARAM_PREALLOC 2 */


enum SICONOS_IPARAM
{
  SICONOS_IPARAM_MAX_ITER = 0,
  SICONOS_IPARAM_ITER_DONE = 1,
  SICONOS_IPARAM_PREALLOC = 2
};

/** for pivot based algorithm */
#define SICONOS_IPARAM_PIVOT_RULE 3

/** pathsearch specific data */
#define SICONOS_IPARAM_PATHSEARCH_STACKSIZE 5

/** line search based algo use this */
#define SICONOS_IPARAM_LSA_NONMONOTONE_LS 3
#define SICONOS_IPARAM_LSA_NONMONOTONE_LS_M 4
#define SICONOS_IPARAM_LSA_FORCE_ARCSEARCH 5
#define SICONOS_IPARAM_LSA_SEARCH_CRITERION 6
#define SICONOS_IPARAM_GOLDSTEIN_ITERMAX 4

/** non-monotone specific part */
#define SICONOS_IPARAM_NMS_WATCHDOG_TYPE 7
#define SICONOS_IPARAM_NMS_PROJECTED_GRADIENT_TYPE 8
#define SICONOS_IPARAM_NMS_N_MAX 9

/** Some values for dparam index */
enum SICONOS_DPARAM
{
  SICONOS_DPARAM_TOL = 0,
  SICONOS_DPARAM_RESIDU = 1
};
/** line-search */
#define SICONOS_DPARAM_LSA_ALPHA_MIN 2
#define SICONOS_DPARAM_GOLDSTEIN_C 3
#define SICONOS_DPARAM_GOLDSTEIN_ALPHAMAX 4

/** non-monotone specific part */
#define SICONOS_DPARAM_NMS_DELTA 2
#define SICONOS_DPARAM_NMS_DELTA_VAR 3
#define SICONOS_DPARAM_NMS_SIGMA 4
#define SICONOS_DPARAM_NMS_ALPHA_MIN_WATCHDOG 5
#define SICONOS_DPARAM_NMS_ALPHA_MIN_PGRAD 6
#define SICONOS_DPARAM_NMS_MERIT_INCR 7


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


#include "SolverOptions_helpers.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** screen display of solver parameters
      \param options the structure to be displayed
  */
  void solver_options_print(SolverOptions* options);

  /** free some SolverOptions fields;
   *   \param options the structure to clean
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
  int solver_options_name_to_id(char * pName);

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

  SolverOptions * solver_options_get_internal_solver(SolverOptions * options, int n);
  
  
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif



#endif
