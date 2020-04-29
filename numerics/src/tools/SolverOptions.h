/* Siconos is a program dedicated to modeling, simulation and 
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include <stdbool.h> // for boolean type
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


// length of iparam/dparam arrays in solver options.
#define OPTIONS_PARAM_SIZE 20

/** \struct SolverOptions_ SolverOptions.h
    Structure used to send options (name, parameters and so on) to a specific solver (mainly from Kernel to Numerics).
    
    \rst
    Creation, update and destruction:

    :func:`solver_options_create`

    :func:`solver_options_update_internal`

    :func:`solver_options_delete`

    Details: :ref:`solver_options`
    \endrst

*/
struct SolverOptions
{
  int solverId;                            /**< id number of the solver. */
  bool isSet;                              /**< true(1) if the structure is ready to be used by a numerics driver. */
  int iSize;                               /**< iSize size of vector iparam */
  int * iparam;                            /**< list of solver parameters (integer type); Check solvers doc for details. */
  int dSize;                               /**< size of vector dparam */
  double * dparam;                         /**< list of solver parameters (double type); Check solvers doc for details. */
  bool filterOn;                           /**< if true (1), check solution validity after the driver call. Default = 1. 
                                              For example if filterOn = 1 for a LCP, lcp_compute_error() 
                                              will be called at the end of the process). */
  size_t dWorkSize;                        /**< size of double type internal work array.*/
  double * dWork;                          /**< internal (double type) work array.*/
  size_t iWorkSize;                        /**< size of integer type internal work array.*/
  int * iWork;                             /**< internal (integer type) work array.*/
  size_t numberOfInternalSolvers;          /**< the number of internal or local 'sub-solvers' used by the solver 
                                              (size of internalSolvers) .*/
  struct SolverOptions ** internalSolvers;  /**< list of internal solver options*/
  Callback * callback;                     /**< pointer to user-defined callback*/
  void * solverParameters;                 /**< additional parameters specific to the solver (GAMS and NewtonMethod only) */
  void * solverData;                       /**< additional data specific to the solver */
} ;


/** Some value for iparam index */
enum SICONOS_IPARAM
{
  SICONOS_IPARAM_MAX_ITER = 0,
  SICONOS_IPARAM_ITER_DONE = 1,
  SICONOS_IPARAM_PREALLOC = 2,
  SICONOS_IPARAM_NSGS_SHUFFLE = 5,
  SICONOS_IPARAM_ERROR_EVALUATION = 3, // CHECK IF THERE ARE NO CONFLICT WITH THIS !!
  SICONOS_IPARAM_PATHSEARCH_STACKSIZE = 19
};

/** allowed values for iparam[SICONOS_IPARAM_ERROR_EVALUATION */
enum SICONOS_IPARAM_ERROR_EVALUATION_ENUM
  {
   /** Complete error computation, including v computation*/
   SICONOS_ERROR_FULL_EVALUATION = 0,
   /** Light error computation with incremental values on r verification of absolute error at the end */
   SICONOS_ERROR_LIGHT_EVALUATION = 1,
   /**  only light error computation, do not update v unknown) */
   SICONOS_ERROR_LIGHT_EVALUATION_NO_UPDATE = 2
  };


/** Some values for dparam index */
enum SICONOS_DPARAM
  {
   SICONOS_DPARAM_TOL = 0,
   SICONOS_DPARAM_RESIDU = 1,
  };

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** screen display of solver parameters
      \param options the structure to be displayed
  */
  void solver_options_print(SolverOptions* options);

  /** Clear and free all pointer members of the structure, then
      release memory
      \param options the structure to be cleared.
  */
  void solver_options_delete(SolverOptions * options);

  /** Create and initialize a SolverOptions struct:
      allocate internal memories, set default values 
      depending on the id.
      \param id solver id number 
      \rst 
      It must belong to one of the available ids defined for each formulation,
      see :ref:`problems_and_solvers` for details.
      \endrst
      \return a pointer to options set, ready to use by a driver.
   */
  SolverOptions * solver_options_create(int solverId);

  /** Copy an existing set of options, to create a new one.

      Warning : callback, solverData and solverParameters of
      the new structure are pointer links to those of the original one!
      
      \param source an existing solver options structure
      \return a pointer to options set, ready to use by a driver.
   */
  SolverOptions * solver_options_copy(SolverOptions* source);

  /** Change one of the internal solver of a previously defined SolverOptions set.
      
      Allocate internal memories and set default values for the internal solver.
      
      Warning : the actual internal solver in position internal_solver_number and all its content will be destroyed
      and replaced by a new one.
      
      \param parent the top-level SolverOptions which contains the internal solver to be updated
      \param internal_solver_number number of the internal solver to be update (warning : this is the position
      in the list of internal solvers, not the id!)
      \param solver_id id number of the new internal solver to be created/updated
   */
  void solver_options_update_internal(SolverOptions* parent, size_t internal_solver_number, int solver_id);
  
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
  
  /** return the an internal solver options set
      \param options parent options
      \param number of the targeted solver
      \return a pointer to the internal solver options set 
   */
  SolverOptions * solver_options_get_internal_solver(SolverOptions * options, size_t n);
  
  
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif



#endif
