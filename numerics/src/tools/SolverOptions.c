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

#include "SolverOptions.h"
#include <assert.h>                         // for assert
#include <float.h>                          // for DBL_EPSILON
#include <math.h>                           // for INFINITY
#include <stdio.h>                          // for NULL, size_t, printf
#include <stdlib.h>                         // for free, calloc, malloc
#include <string.h>                         // for strcmp
#include "AVI_cst.h"                        // for SICONOS_AVI_CAOFERRIS_STR
#include "ConvexQP_Solvers.h"               // for convexQP_ADMM_set_default
#include "ConvexQP_cst.h"                   // for SICONOS_CONVEXQP_ADMM_STR
#include "Friction_cst.h"                   // for SICONOS_FRICTION_2D_CPG_STR
#include "GenericMechanical_Solvers.h"      // for gmp_set_default
#include "GenericMechanical_cst.h"          // for SICONOS_GENERIC_MECHANICA...
#include "LCP_Solvers.h"                    // for lcp_pivot_set_default
#include "MCP_cst.h"                        // for SICONOS_MCP_NEWTON_FB_FBL...
#include "MLCP_Solvers.h"                   // for mlcp_direct_set_default
#include "NCP_cst.h"                        // for SICONOS_NCP_NEWTON_FB_FBL...
#include "Newton_methods.h"                 // for newton_lsa_set_default
#include "NonSmoothNewton.h"                // for nonSmoothNewton_set_default
#include "PathSearch.h"                     // for pathsearch_set_default
#include "SOCLCP_Solvers.h"                 // for soclcp_nsgs_set_default
#include "SOCLCP_cst.h"                     // for SICONOS_SOCLCP_NSGS_STR
#include "SiconosNumerics_Solvers.h"        // for SICONOS_REGISTER_SOLVERS
#include "VI_cst.h"                         // for SICONOS_VI_BOX_AVI_LSA_STR
#include "VariationalInequality_Solvers.h"  // for variationalInequality_BOX...
#include "fc2d_Solvers.h"                   // for fc2d_nsgs_set_default
#include "fc3d_Solvers.h"                   // for fc3d_nsgs_set_default
#include "gfc3d_Solvers.h"                  // for gfc3d_aclmfp_set_default
#include "grfc3d_Solvers.h"      // for grfc3d_IPM_set_default
#include "lcp_cst.h"                        // for SICONOS_LCP_AVI_CAOFERRIS...
#include "mlcp_cst.h"                       // for SICONOS_MLCP_DIRECT_ENUM_STR
#include "numerics_verbose.h"               // for numerics_printf, numerics...
#include "relay_cst.h"                      // for SICONOS_RELAY_AVI_CAOFERR...
#include "rolling_fc_Solvers.h"           // for rfc3d_poc_set_default

/** Create a struct SolverOptions and initialize its content.

    Allocate memory for internal parameters array, ensure that pointers
    are properly set/nullified and tolerance, max iter and common default values.
    All other specific parameters must be set explicitely for each solver by
    its internal set default function.

    \param options struct to fill
    \param solver_id id number of the solver
    \param iter_max maximum number of iterations allowed for this solver
    \param tol tolerance for the solution.
*/
static SolverOptions* solver_options_initialize(int solver_id, int iter_max, double tol, size_t number_of_internal_solvers)
{

  // Warning : the C structure members after the first instanciation might be anything.
  // (C structure members can not have default values).
  // Maybe we should write a function to create so, like so * = so_new ?
  SolverOptions * options = (SolverOptions *)malloc(sizeof(SolverOptions));

  options->solverId = solver_id;
  options->iSize = 20;
  options->dSize = 20;
  options->iparam = (int *)calloc(20, sizeof(int));
  options->dparam = (double *)calloc(20, sizeof(double));

  // The content of iparam and dparam for indices in SICONOS_IPARAM and SICONOS_DPARAM enums
  // enum must be set in this function with a default value.

  options->iparam[SICONOS_IPARAM_MAX_ITER] = iter_max;
  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->iparam[SICONOS_IPARAM_PREALLOC] = 0;
  options->iparam[SICONOS_IPARAM_NSGS_SHUFFLE] = 0;  // not common to all formulations ?
  options->iparam[SICONOS_IPARAM_ERROR_EVALUATION] = 0; // not common to all formulations ?
  options->iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE] = 0;  // not common to all formulations ?

  options->dparam[SICONOS_DPARAM_TOL] = tol;
  options->dparam[SICONOS_DPARAM_RESIDU] = INFINITY;

  options->filterOn = true;

  // Everything regarding dwork, iwork, solverData and solverParameters is problem/formulation specific
  // and must be handled later, locally, by the driver when the problem is known.
  //
  options->dWorkSize = 0;
  options->dWork = NULL;
  options->dWorkSize = 0;
  options->iWork = NULL;
  options->callback = NULL;
  options->numberOfInternalSolvers = number_of_internal_solvers;
  options->internalSolvers = calloc(options->numberOfInternalSolvers, sizeof(SolverOptions*));
  options->solverData = NULL;
  options->solverParameters = NULL;

  options->isSet = true;
  return options;
}

/** */
static void recursive_solver_options_print(SolverOptions* options, int level)
{
  char* marge;
  marge = (char*) malloc((size_t)(level + 1) * sizeof(char));
  for(int i = 0; i < level; i++)
    marge[i] = ' ';
  marge[level] = '\0';

  numerics_printf("%s========== solver parameters: ", marge);
  if(!options->isSet)
    numerics_printf("%sThe solver parameters have not been set. \t options->isSet = %i ", marge, options->isSet);
  else
  {
    numerics_printf("%sThe solver parameters below have  been set \t options->isSet = %i", marge, options->isSet);
    numerics_printf("%sId of the solver\t\t\t\t options->solverId = %d ", marge, options->solverId);
    numerics_printf("%sName of the solver\t\t\t\t %s ", marge, solver_options_id_to_name(options->solverId));
    if(options->iparam != NULL)
    {
      numerics_printf("%ssize of the int parameters\t\t\t options->iSize = %i", marge, options->iSize);
      numerics_printf("%snon zero int parameters in options->iparam:", marge);
      for(int i = 0; i < options->iSize; ++i)
      {
        if(options->iparam[i]) numerics_printf("%s\t\t\t\t\t\t options->iparam[%i] = %d", marge, i, options->iparam[i]);
      }
    }
    if(options->dparam != NULL)
    {
      numerics_printf("%sdouble parameters \t\t\t\t options->dparam", marge);
      numerics_printf("%ssize of the double parameters\t\t\t options->dSize = %i", marge, options->dSize);
      numerics_printf("%snon zero double parameters in options->dparam:", marge);
      for(int i = 0; i < options->dSize; ++i)
      {
        if(options->dparam[i]>0.) numerics_printf("%s\t\t\t\t\t\t options->dparam[%i] = %.6le", marge, i, options->dparam[i]);
      }
    }
  }
  if(options->iWork == NULL)
  {
    numerics_printf("%sinteger work array have not been allocated. \t options->iWork = NULL ", marge);
  }
  else
  {
    numerics_printf("%sinteger work array have been allocated. \t options->iWork = %p ", marge, options->iWork);
    numerics_printf("%sinteger work array size \t\t\t options->iSize = %i ", marge, options->iSize);
  }
  if(options->dWork == NULL)
  {
    numerics_printf("%sdouble work array have not been allocated. \t options->dWork = NULL ", marge);
  }
  else
  {
    numerics_printf("%sdouble work array have been allocated. \t options->dWork = %p ", marge, options->dWork);
    numerics_printf("%sdouble work array size \t\t\t options->dSize = %i ", marge, options->dSize);
  }




  numerics_printf("%sSee %s documentation for parameters definition)", marge, solver_options_id_to_name(options->solverId));


  numerics_printf("%snumber of internal (or local) solvers \t\t options->numberOfInternalSolvers = %i", marge, options->numberOfInternalSolvers);
  for(size_t i = 0; i < options->numberOfInternalSolvers; i++)
  {
    recursive_solver_options_print(options->internalSolvers[i], level + 1);
  }
  free(marge);

}

void solver_options_print(SolverOptions* options)
{
  recursive_solver_options_print(options, 0);
}



void solver_options_delete(SolverOptions* op)
{
  if(op)
  {
    // Clear solverParameters and solverData, before anything.
    // Remark : these are specific data. And so, alloc/release
    // memory operations should be handled inside each
    // solver, at the end of the driver call.
    if(op->solverParameters)
      free(op->solverParameters);
    op->solverParameters = NULL;

    if(op->solverData)
      free(op->solverData);
    op->solverData = NULL;

    // Clear callback
    if(op->callback)
      free(op->callback);

    op->callback = NULL;

    // Clear internal solver(s)
    if(op->internalSolvers)
    {
      for(size_t i = 0; i < op->numberOfInternalSolvers; i++)
        solver_options_delete(op->internalSolvers[i]);
      op->numberOfInternalSolvers = 0;
      free(op->internalSolvers);
    }
    op->internalSolvers = NULL;

    // working arrays
    // Solver-specific data : see remark on top of this function,
    // regarding solverData.
    if(op->iWork)
      free(op->iWork);
    op->iWork = NULL;
    if(op->dWork)
      free(op->dWork);
    op->dWork = NULL;

    // int and double parameters arrays
    if(op->iparam)
      free(op->iparam);
    op->iparam = NULL;
    if(op->dparam)
      free(op->dparam);
    op->dparam = NULL;

    op->isSet = false;
  }
//  double free or corruption
//  free(op);
}

SolverOptions * solver_options_copy(SolverOptions* source)
{
  // Create a new solver options, with default setup
  SolverOptions * options = solver_options_create(source->solverId);

  for(size_t i=0; i < OPTIONS_PARAM_SIZE; ++i)
  {
    options->iparam[i] = source->iparam[i];
    options->dparam[i] = source->dparam[i];
  }

  if(source->dWork)
  {
    options->dWork = calloc(source->dWorkSize, sizeof(double));
    options->dWorkSize = source->dWorkSize;
    for(size_t i=0; i<options->dWorkSize; ++i)
      options->dWork[i] = source->dWork[i];
  }
  if(source->iWork)
  {
    options->iWork = calloc(source->iWorkSize, sizeof(int));
    options->iWorkSize = source->iWorkSize;
    for(size_t i=0; i<options->iWorkSize; ++i)
      options->iWork[i] = source->iWork[i];
  }
  assert(options->numberOfInternalSolvers == source->numberOfInternalSolvers);
  // this assert should be ensured by solver_options_create and initialize.

  for(size_t i=0; i<options->numberOfInternalSolvers; ++i)
    options->internalSolvers[i] = solver_options_copy(source->internalSolvers[i]);

  // Warning pointer links!
  if(source->callback)
    options->callback = source->callback; // Note FP: is it really safe to create pointer link here?

  if(source->solverData)
    options->solverData =source->solverData;

  if(source->solverParameters)
    options->solverParameters =source->solverParameters;

  return options;
}


SolverOptions * solver_options_create(int solverId)
{
  // This function must be the unique way for users to create a SolverOptions.
  // It ensures that the object is ready to use by a driver (pointers ready,
  // memory allocated, minimum default values set ...)
  // Any further modification of the default values must be done
  // by setting iparam or dparam explicitely of with solver_options_update_internal
  // function when it turns to internal solvers.

  SolverOptions * options = NULL;

  // In each case :
  // - first call a default setup, common to all solvers, with tolerance and max iter value,
  //   to ensure that pointers are ready, and minimum default values set;
  // - then call <formulation>_<solver_name>_set_default(options) to set value specific to
  //   each solver, if required (if not, the function does not exist).
  switch(solverId)
  {
  // --- VI Solvers ---
  // ref list : enum VI_SOLVER
  case SICONOS_VI_EG:
  case SICONOS_CONVEXQP_VI_EG:
  case SICONOS_SOCLCP_VI_EG:
  case SICONOS_FRICTION_3D_VI_EG:
  case SICONOS_GLOBAL_FRICTION_3D_VI_EG:
  {
    options = solver_options_initialize(solverId, 20000, 1e-3, 0);
    variationalInequality_ExtraGradient_set_default(options);
    break;
  }
  case SICONOS_VI_FPP:
  case SICONOS_CONVEXQP_VI_FPP:
  case SICONOS_SOCLCP_VI_FPP:
  case SICONOS_FRICTION_3D_VI_FPP:
  case SICONOS_FRICTION_3D_VI_FPP_Cylinder:
  case SICONOS_GLOBAL_FRICTION_3D_VI_FPP:
  {
    options = solver_options_initialize(solverId, 20000, 1e-3, 0);
    variationalInequality_FixedPointProjection_set_default(options);
    break;
  }
  case SICONOS_VI_HP:
  {
    options = solver_options_initialize(solverId, 20000, 1e-3, 0);
    variationalInequality_HyperplaneProjection_set_default(options);
    break;
  }
  case SICONOS_VI_BOX_QI:
  {
    options = solver_options_initialize(solverId, 1000, 1e-10, 0);
    variationalInequality_BOX_QI_set_default(options);
    break;
  }

  case SICONOS_VI_BOX_AVI_LSA:
  {
    options = solver_options_initialize(solverId, 100, 1e-12, 1);
    variationalInequality_BOX_AVI_set_default(options);
    break;
  }

  case SICONOS_VI_BOX_PATH:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    break;
  }

  // --- QP Solvers ---
  // ref list : enum CONVEXQP_SOLVER in ConvexQP_cst.h
  case SICONOS_LCP_CONVEXQP_PG:
  case SICONOS_CONVEXQP_PG:
  case SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER:
  {
    options = solver_options_initialize(solverId, 20000, 1e-6, 0);
    convexQP_ProjectedGradient_set_default(options);
    break;
  }

  case SICONOS_CONVEXQP_ADMM:
  {
    options = solver_options_initialize(solverId, 20000, 1e-6, 0);
    convexQP_ADMM_set_default(options);
    break;
  }
  // --- LCP and Relay Solvers ---
  // ref list : enum LCP_SOLVER in LCP_cst.h
  // ref list : enum RELAY_SOLVER in Relay_cst.h
  case SICONOS_LCP_LEMKE:
  case SICONOS_RELAY_LEMKE:
  case SICONOS_FRICTION_2D_LEMKE:
  {
    options = solver_options_initialize(solverId, 10000, 1e-6, 0);
    lcp_lexicolemke_set_default(options);
    break;
  }
  case SICONOS_LCP_NSGS_SBM:
  {
    options = solver_options_initialize(solverId, 1000, 1e-6, 1);
    lcp_nsgs_sbm_set_default(options);
    break;
  }
  case SICONOS_LCP_PGS:
  case SICONOS_RELAY_PGS:
  {
    options = solver_options_initialize(solverId, 10000, 1e-6, 0);
    break;
  }

  case SICONOS_LCP_CPG:
  {
    options = solver_options_initialize(solverId, 10000, 1e-6, 0);
    break;
  }
  case SICONOS_LCP_LATIN:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 0);
    lcp_latin_set_default(options);
    break;
  }

  case SICONOS_LCP_LATIN_W:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 0);
    lcp_latin_w_set_default(options);
    break;
  }

  case SICONOS_LCP_QP:
  case SICONOS_LCP_NSQP:
  {
    options = solver_options_initialize(solverId, 0, 1e-6, 0);
    break;
  }

  case SICONOS_LCP_NEWTONMIN:
  {
    options = solver_options_initialize(solverId, 1000, 1e-6, 0);
    break;
  }

  case SICONOS_LCP_NEWTON_FB_FBLSA:
  case SICONOS_LCP_NEWTON_MIN_FBLSA:
  {
    options = solver_options_initialize(solverId, 1000, 1e-10, 0);
    lcp_newton_FB_set_default(options);
    break;
  }

  case SICONOS_LCP_PSOR:
  {
    options = solver_options_initialize(solverId, 1000, 1e-6, 0);
    lcp_psor_set_default(options);
    break;
  }

  case SICONOS_LCP_RPGS:
  {
    options = solver_options_initialize(solverId, 10000, 1e-6, 0);
    lcp_rpgs_set_default(options);
    break;
  }

  case SICONOS_LCP_PATH:
  case SICONOS_RELAY_PATH:
  case SICONOS_MLCP_PATH:
  {
    options = solver_options_initialize(solverId, 1000, 1e-12, 0);
    break;
  }

  case SICONOS_LCP_ENUM:
  case SICONOS_RELAY_ENUM:
  case SICONOS_FRICTION_2D_ENUM:
  {
    options = solver_options_initialize(solverId, 0, 1e-6, 0);
    lcp_enum_set_default(options);
    break;
  }

  case SICONOS_LCP_PIVOT:
  {
    options = solver_options_initialize(solverId, 10000, 100 * DBL_EPSILON, 0);
    lcp_pivot_set_default(options);
    break;
  }

  case SICONOS_LCP_BARD:
  {
    options = solver_options_initialize(solverId, 10000, 100 * DBL_EPSILON, 0);
    lcp_pivot_set_default(options);
    // iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = SICONOS_LCP_PIVOT_BARD; set in lcp_driver
    break;
  }
  case SICONOS_LCP_MURTY:
  {
    options = solver_options_initialize(solverId, 10000, 100 * DBL_EPSILON, 0);
    lcp_pivot_set_default(options);
    // iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = SICONOS_LCP_PIVOT_LEAST_INDEX set in lcp_driver
    break;
  }

  case SICONOS_LCP_PATHSEARCH:
  {
    options = solver_options_initialize(solverId, 10000, 100 * DBL_EPSILON, 0);
    lcp_pathsearch_set_default(options);
    break;
  }

  case SICONOS_LCP_PIVOT_LUMOD:
  {
    options = solver_options_initialize(solverId, 10000, 100 * DBL_EPSILON, 0);
    lcp_pivot_lumod_set_default(options); // same as lcp_pivot
    break;
  }

  case SICONOS_LCP_GAMS:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    break;
  }

  // --- AVI Solvers ---
  // ref list : enum AVI_SOLVER in AVI_cst.h
  case SICONOS_AVI_CAOFERRIS:
  case SICONOS_AVI_PATHAVI:
  case SICONOS_RELAY_AVI_CAOFERRIS:
  case SICONOS_RELAY_AVI_CAOFERRIS_TEST:
  case SICONOS_LCP_AVI_CAOFERRIS:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    break;
  }

  // --- SOCLCP Solvers ---
  // ref list : enum SOCLCP_SOLVER in SOCLCP_cst.h
  case SICONOS_SOCLCP_NSGS:
  case SICONOS_FRICTION_3D_SOCLCP:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 1);
    soclcp_nsgs_set_default(options);
    break;
  }

  case SICONOS_SOCLCP_ProjectionOnConeWithLocalIteration:
  case SICONOS_SOCLCP_ProjectionOnCone:
  case SICONOS_SOCLCP_ProjectionOnConeWithRegularization:
  {
    options = solver_options_initialize(solverId, 1000, 1e-16, 0);
    soclcp_projection_set_default(options);
    break;
  }

  // --- MLCP Solvers ---
  // ref list : enum MLCP_SOLVER in mlcp_cst.h
  case SICONOS_MLCP_PGS:
  {
    options = solver_options_initialize(solverId, 1000, 1e-6, 0);
    mlcp_pgs_set_default(options);
    break;
  }
  case SICONOS_MLCP_PGS_SBM:
  {
    options = solver_options_initialize(solverId, 1000, 1e-6, 1);
    mlcp_pgs_sbm_set_default(options);
    break;
  }
  case SICONOS_MLCP_RPGS:
  {
    options = solver_options_initialize(solverId, 1000, 1e-6, 0);
    mlcp_rpgs_set_default(options);
    break;
  }
  case SICONOS_MLCP_PSOR:
  {
    options = solver_options_initialize(solverId, 1000, 1e-6, 0);
    mlcp_psor_set_default(options);
    break;
  }
  case SICONOS_MLCP_RPSOR:
  {
    options = solver_options_initialize(solverId, 1000, 1e-6, 0);
    mlcp_rpsor_set_default(options);
    break;
  }
  case SICONOS_MLCP_ENUM:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    mlcp_enum_set_default(options);
    break;
  }
  case SICONOS_MLCP_DIRECT_ENUM:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    mlcp_direct_enum_set_default(options);
    break;
  }
  case SICONOS_MLCP_SIMPLEX:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    break;
  }

  case SICONOS_MLCP_DIRECT_SIMPLEX:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    mlcp_direct_set_default(options);
    break;
  }
  case SICONOS_MLCP_PATH_ENUM:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    mlcp_enum_set_default(options);
    break;
  }
  case SICONOS_MLCP_DIRECT_PATH:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    mlcp_direct_set_default(options);
    break;
  }
  case SICONOS_MLCP_DIRECT_PATH_ENUM:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    mlcp_enum_set_default(options);
    mlcp_direct_set_default(options);
    break;
  }
  case SICONOS_MLCP_FB:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    break;
  }
  case SICONOS_MLCP_DIRECT_FB:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    mlcp_direct_set_default(options);
    break;
  }
  case SICONOS_MLCP_LCP_LEMKE:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    mlcp_direct_set_default(options);
    break;
  }

  // --- NCP Solvers ---
  // ref list : enum NCP_SOLVER in NCP_cst.h
  case SICONOS_NCP_NEWTON_FB_FBLSA:
  case SICONOS_NCP_NEWTON_MIN_FBLSA:
  {
    options = solver_options_initialize(solverId, 1000, 1e-12, 0);
    newton_lsa_set_default(options);
    break;
  };
  case SICONOS_NCP_PATHSEARCH:
  {
    options = solver_options_initialize(solverId, 100, 1e-12, 1);
    pathsearch_set_default(options);
    break;
  }
  case SICONOS_NCP_PATH:
  {
    options = solver_options_initialize(solverId, 10000, 1e-12, 0);
    break;
  }

  // --- MCP Solvers ---
  // ref list : enum MCP_SOLVER in MCP_cst.h
  case SICONOS_MCP_NEWTON_FB_FBLSA:
  case SICONOS_MCP_NEWTON_MIN_FBLSA:
  {
    options = solver_options_initialize(solverId, 1000, 1e-10, 0);
    newton_lsa_set_default(options);
    break;
  }
  case SICONOS_MCP_OLD_FB:
  {
    options = solver_options_initialize(solverId, 10, 1e-7, 0);
    nonSmoothNewton_set_default(options);
    break;
  }


  // --- Friction Solvers ---
  // ref list : enum FRICTION_SOLVER in Friction_cst.h
  case SICONOS_FRICTION_2D_NSGS:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 0);
    fc2d_nsgs_set_default(options);
    break;
  }
  case SICONOS_FRICTION_2D_CPG:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 0);
    break;
  }
  case SICONOS_FRICTION_3D_NSGS:
  case SICONOS_GLOBAL_FRICTION_3D_NSGS:
  case SICONOS_GLOBAL_FRICTION_3D_NSGS_WR:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 1);
    fc3d_nsgs_set_default(options);
    break;
  }

  case SICONOS_ROLLING_FRICTION_3D_NSGS:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 1);
    rfc3d_nsgs_set_default(options);
    break;
  }
  case SICONOS_ROLLING_FRICTION_3D_ADMM:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 0);
    rolling_fc3d_admm_set_default(options);
    break;
  }
  case SICONOS_ROLLING_FRICTION_2D_NSGS:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 1);
    rfc2d_nsgs_set_default(options);
    break;
  }

  case SICONOS_FRICTION_3D_NSGSV:
  case SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 1);
    fc3d_nsgs_velocity_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_PROX:
  case SICONOS_GLOBAL_FRICTION_3D_PROX_WR:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 1);
    fc3d_proximal_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_TFP:
  case SICONOS_GLOBAL_FRICTION_3D_TFP_WR:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 1);
    fc3d_tfp_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_NSN_AC:
  case SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR:
  {
    options = solver_options_initialize(solverId, 200, 1e-3, 0);
    fc3d_nsn_ac_set_default(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_NSN_AC:
  {
    options = solver_options_initialize(solverId, 200, 1e-10, 0);
    gfc3d_nsn_ac_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_NSN_AC_TEST:
  {
    options = solver_options_initialize(solverId, 1000, 1e-10, 0);
    newton_lsa_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_DSFP:
  case SICONOS_GLOBAL_FRICTION_3D_DSFP_WR:
  {
    options = solver_options_initialize(solverId, 20000, 1e-3, 0);
    fc3d_dsfp_set_default(options);
    break;
  }

  case SICONOS_FRICTION_3D_HP:
  {
    options = solver_options_initialize(solverId, 2000000, 1e-3, 0);
    fc3d_hp_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_FPP:
  {
    options = solver_options_initialize(solverId, 20000, 1e-3, 0);
    fc3d_fpp_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_EG:
  {
    options = solver_options_initialize(solverId, 20000, 1e-3, 0);
    fc3d_eg_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_NSN_FB:
  {
    options = solver_options_initialize(solverId, 200, 1e-3, 0);
    fc3d_nsn_fb_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_GAMS_PATH:
  case SICONOS_FRICTION_3D_GAMS_PATHVI:
  case SICONOS_FRICTION_3D_GAMS_LCP_PATH:
  case SICONOS_FRICTION_3D_GAMS_LCP_PATHVI:
  case SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH:
  case SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI:
  {
    options = solver_options_initialize(solverId, 10000, 1e-9, 0);
    break;
  }

  case SICONOS_FRICTION_3D_ACLMFP:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 1);
    fc3d_aclmfp_set_default(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_ACLMFP:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 1);
    gfc3d_aclmfp_set_default(options);
    break;
  }

  case SICONOS_FRICTION_3D_NSN_NM:
  {
    options = solver_options_initialize(solverId, 200, 1e-3, 0);
    fc3d_nsn_nm_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_PFP:
  {
    options = solver_options_initialize(solverId, 1000, 1e-4, 2);
    fc3d_pfp_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_ADMM:
  case SICONOS_GLOBAL_FRICTION_3D_ADMM_WR:
  {
    options = solver_options_initialize(solverId, 20000, 1e-6, 0);
    fc3d_admm_set_default(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_ADMM:
  {
    options = solver_options_initialize(solverId, 20000, 1e-6, 0);
    gfc3d_admm_set_default(options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_IPM:
  {
    options = solver_options_initialize(solverId, 20000, 1e-6, 0);
    gfc3d_ipm_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_NSN:
  {
    options = solver_options_initialize(solverId, 10, 1e-14, 0);
    fc3d_onecontact_nsn_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_GP:
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID:
  {
    options = solver_options_initialize(solverId, 10, 1e-14, 0);
    fc3d_onecontact_nsn_gp_set_default(options);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone:
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity:
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration:
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization:
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization:
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder:
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration:
  {
    options = solver_options_initialize(solverId, 1000, 1e-14, 0);
    fc3d_poc_set_default(options);
    break;
  }
  case SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone:
  {
    options = solver_options_initialize(solverId, 1000, 1e-12, 0);
    rfc3d_poc_set_default(options);
    break;
  }

  case SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration:
  {
    options = solver_options_initialize(solverId, 1000, 1e-12, 0);
    rfc3d_poc_withLocalIteration_set_default(options);
    break;
  }
  case SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone:
  {
    options = solver_options_initialize(solverId, 1000, 1e-12, 0);
    rfc2d_poc_set_default(options);
    break;
  }

  case SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration:
  {
    options = solver_options_initialize(solverId, 1000, 1e-12, 0);
    rfc2d_poc_withLocalIteration_set_default(options);
    break;
  }
  case SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR:
  {
    options = solver_options_initialize(solverId, 1000, 1e-12, 1);
    rfc3d_nsgs_set_default(options);
    break;
  }
  case SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM:
  {
    options = solver_options_initialize(solverId, 20000, 1e-6, 0);
    grfc3d_IPM_set_default(options);
    break;
  }

  case SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint:
  case SICONOS_FRICTION_3D_NCPGlockerFBNewton:
  case SICONOS_FRICTION_3D_NCPGlockerFBPATH:
  {
    options = solver_options_initialize(solverId, 1000, 1e-12, 0);
    break;
  }
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC:
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU:
  {
    options = solver_options_initialize(solverId, 1000, 1e-12, 0);
    break;
  }

  // --- GMP Solver ---
  case SICONOS_GENERIC_MECHANICAL_NSGS:
  {
    options = solver_options_initialize(solverId, 10000, 1e-4, 4);
    gmp_set_default(options);
    break;
  }
  default:
    numerics_error("solver_options_create", "Unknown solver : %s (%i) !", solver_options_id_to_name(solverId), solverId);

  }


  return options;
}

void solver_options_update_internal(SolverOptions* parent, size_t internal_solver_number, int solver_id)
{
  // Avoid access to a non existing internal solver.
  assert(parent->numberOfInternalSolvers > internal_solver_number);

  // Destroy current internal solver and create a new one with the new id.
  solver_options_delete(parent->internalSolvers[internal_solver_number]);
  parent->internalSolvers[internal_solver_number] = solver_options_create(solver_id);
}


const char * solver_options_id_to_name(int Id)
{
  switch(Id)
  {

#undef SICONOS_SOLVER_MACRO
#define SICONOS_SOLVER_MACRO(X) case X: return X ## _STR ;
    SICONOS_REGISTER_SOLVERS()
  default:
    return SICONOS_NONAME_STR;
  }
}

int solver_options_name_to_id(const char * pName)
{
#undef SICONOS_SOLVER_MACRO
#define SICONOS_SOLVER_MACRO(X) if (strcmp(X ## _STR, pName) == 0) return X;
  SICONOS_REGISTER_SOLVERS()
  return 0;
}

SolverOptions * solver_options_get_internal_solver(SolverOptions * options, size_t n)
{
  if(n+1 > options->numberOfInternalSolvers)
  {
    printf("solver_options_get_internal_solver : the index must be between 0 and  options->numberOfInternalSolvers -1 ");
    return NULL;
  }
  else
    return options->internalSolvers[n];
}

void solver_options_set_internal_solver(SolverOptions * options, size_t n, SolverOptions* NSO)
{
  if(n+1 > options->numberOfInternalSolvers)
  {
    printf("solver_options_set_internal_solver : the index must be between 0 and  options->numberOfInternalSolvers -1 ");
  }
  else
  {
    options->internalSolvers[n] = NSO;
  }
}

