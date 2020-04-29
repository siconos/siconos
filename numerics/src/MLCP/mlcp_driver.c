/* Siconos is a program dedicated to modeling, simulation and control
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

#include <stdio.h>                              // for NULL, fprintf, stderr
#include <stdlib.h>                             // for exit, EXIT_FAILURE
#include "MLCP_Solvers.h"                       // for mlcp_FB, mlcp_direct_FB
#include "MixedLinearComplementarityProblem.h"  // for mixedLinearComplement...
#include "NumericsFwd.h"                        // for SolverOptions, MixedL...
#include "NumericsMatrix.h"                     // for NumericsMatrix
#include "SolverOptions.h"                      // for SolverOptions
#include "mlcp_FB.h"                            // for mlcp_FB_getNbDWork
#include "mlcp_cst.h"                           // for SICONOS_MLCP_DIRECT_ENUM
#include "mlcp_direct.h"                     // for mlcp_direct_FB_getNbD...
#include "mlcp_direct_FB.h"                     // for mlcp_direct_FB_getNbD...
#include "mlcp_direct_enum.h"                   // for mlcp_direct_enum_getN...
#include "mlcp_direct_path.h"                   // for mlcp_direct_path_getN...
#include "mlcp_direct_path_enum.h"              // for mlcp_direct_path_enum
#include "mlcp_direct_simplex.h"                // for mlcp_direct_simplex_g...
#include "mlcp_enum.h"                          // for mlcp_enum_alloc_worki...
#include "mlcp_path_enum.h"                     // for mlcp_path_enum, mlcp_...
#include "mlcp_simplex.h"                       // for mlcp_simplex_init
#include "numerics_verbose.h"                   // for numerics_error, verbose

#ifndef MEXFLAG
#include "NonSmoothDrivers.h"    // for mlcp_driver
#endif

const char* const   SICONOS_NONAME_STR = "NONAME";
const char* const   SICONOS_MLCP_PGS_STR = "MLCP_PGS";
const char* const   SICONOS_MLCP_RPGS_STR = "MLCP_RPGS";
const char* const   SICONOS_MLCP_PSOR_STR = "MLCP_PSOR";
const char* const   SICONOS_MLCP_RPSOR_STR = "MLCP_RPSOR";
const char* const   SICONOS_MLCP_PATH_STR = "MLCP_PATH";
const char* const   SICONOS_MLCP_ENUM_STR = "MLCP_ENUM";
const char* const   SICONOS_MLCP_SIMPLEX_STR = "MLCP_SIMPLEX";
const char* const   SICONOS_MLCP_DIRECT_ENUM_STR = "MLCP_DIRECT_ENUM";
const char* const   SICONOS_MLCP_PATH_ENUM_STR = "MLCP_PATH_ENUM";
const char* const   SICONOS_MLCP_DIRECT_SIMPLEX_STR = "MLCP_DIRECT_SIMPLEX";
const char* const   SICONOS_MLCP_DIRECT_PATH_STR = "MLCP_DIRECT_PATH";
const char* const   SICONOS_MLCP_DIRECT_PATH_ENUM_STR = "MLCP_DIRECT_PATH_ENUM";
const char* const   SICONOS_MLCP_FB_STR = "MLCP_FB";
const char* const   SICONOS_MLCP_DIRECT_FB_STR = "MLCP_DIRECT_FB";
const char* const   SICONOS_MLCP_PGS_SBM_STR = "MLCP_PGS_SBM";


/** Compute the size of internal work arrays (for SolverOptions struct) and allocate them.

    static function : it does not need to be in the interface.

    \param[in] problem the MixedLinearComplementarityProblem structure which handles the problem (M,q)
    \param[in] options structure used to define the solver(s) and their parameters
    \param[out] iwize integer work array size
    \param[out] dwsize double work array size
*/
static void mlcp_driver_allocate_work_arrays(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  size_t iwsize, dwsize;
  // compute internal work arrays size, depending on the chosen solver.
  switch(options->solverId)
  {
  case SICONOS_MLCP_DIRECT_ENUM:
    iwsize = mlcp_direct_getNbIWork(problem, options) + mlcp_enum_getNbIWork(problem, options);
    dwsize = mlcp_direct_getNbDWork(problem, options) + mlcp_enum_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_PATH_ENUM:
    iwsize = mlcp_enum_getNbIWork(problem, options);
    dwsize = mlcp_enum_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_PATH_ENUM:
    iwsize =mlcp_direct_getNbIWork(problem, options) + mlcp_enum_getNbIWork(problem, options);
    dwsize =mlcp_direct_getNbDWork(problem, options) + mlcp_enum_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_ENUM:
    iwsize = mlcp_enum_getNbIWork(problem, options);
    dwsize = mlcp_enum_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_SIMPLEX:
    iwsize = mlcp_direct_getNbIWork(problem, options);
    dwsize = mlcp_direct_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_PATH:
    iwsize = mlcp_direct_getNbIWork(problem, options);
    dwsize = mlcp_direct_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_FB:
    iwsize = mlcp_FB_getNbIWork(problem, options);
    dwsize = mlcp_FB_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_FB:
    iwsize = mlcp_direct_getNbIWork(problem, options) + mlcp_FB_getNbIWork(problem, options);
    dwsize = mlcp_direct_getNbDWork(problem, options) + mlcp_FB_getNbDWork(problem, options);
    break;

  default :
    iwsize = 0;
    dwsize = 0;
  }
  // allocate solver options working arrays.
  options->iWorkSize = iwsize;
  options->dWorkSize = dwsize;
  if(options->iWorkSize)
    options->iWork = (int*)calloc(options->iWorkSize, sizeof(int));
  if(options->dWorkSize)
    options->dWork = (double*)calloc(options->dWorkSize, sizeof(double));
}

void mlcp_driver_init(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  // iWork, dWork arrays memory allocation
  mlcp_driver_allocate_work_arrays(problem, options);

  switch(options->solverId)
  {
  case SICONOS_MLCP_DIRECT_ENUM :
    mlcp_direct_enum_init(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_PATH_ENUM :
    mlcp_direct_path_enum_init(problem, options);
    break;
  case SICONOS_MLCP_PATH_ENUM :
    mlcp_path_enum_init(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_SIMPLEX :
    mlcp_direct_simplex_init(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_PATH :
    mlcp_direct_path_init(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_FB :
    mlcp_direct_FB_init(problem, options);
    break;
  case SICONOS_MLCP_SIMPLEX :
    mlcp_simplex_init(problem, options);
    break;
  case SICONOS_MLCP_FB :
    mlcp_FB_init(problem, options);
  default:
    ;/*Nothing to do*/
  }

}

void mlcp_driver_reset(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{


  switch(options->solverId)
  {
  case SICONOS_MLCP_DIRECT_ENUM :
    mlcp_direct_enum_reset();
    break;
  case SICONOS_MLCP_DIRECT_PATH_ENUM :
    mlcp_direct_path_enum_reset();
    break;
  case SICONOS_MLCP_PATH_ENUM :
    mlcp_path_enum_reset();
    break;
  case SICONOS_MLCP_DIRECT_SIMPLEX :
    mlcp_direct_simplex_reset();
    break;
  case SICONOS_MLCP_DIRECT_PATH :
    mlcp_direct_path_reset();
    break;
  case SICONOS_MLCP_DIRECT_FB :
    mlcp_direct_FB_reset();
    break;
  case SICONOS_MLCP_SIMPLEX :
    mlcp_simplex_reset();
    break;
  case SICONOS_MLCP_FB :
    mlcp_FB_reset();
    break;

  default:
    ;/*Nothing to do*/
  }

}

int mlcp_driver(MixedLinearComplementarityProblem* problem, double *z, double *w, SolverOptions* options)
{

  /* verbose=1; */

  if(options == NULL)
    numerics_error("mlcp_driver ", "null input for solver options.\n");

  /* Checks inputs */
  if(problem == NULL || z == NULL || w == NULL)
    numerics_error("mlcp_driver", "null input for MixedLinearComplementarityProblem and/or unknowns (z,w)");
  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;
//  if(verbose)
//    mixedLinearComplementarity_display(problem);
  if(verbose)
    solver_options_print(options);

  /* Switch to DenseMatrix or SparseBlockMatrix solver according to the type of storage for M */
  /* Storage type for the matrix M of the LCP */
  int storageType = problem->M->storageType;

  /* Sparse Block Storage */
  if(storageType == NM_SPARSE_BLOCK)
  {
    numerics_error("mlcp_driver", "not yet implemented for sparse block storage (NM_SPARSE_BLOCK)");
  }
  // else

  /*************************************************
   *  2 - Call specific solver (if no trivial sol.)
   *************************************************/

  /*  if(verbose==1){
    printf(" ========================== Call %s solver ==========================\n", name);
    printf("Initial z value:\n");
    for (i=0;i<problem->n+problem->m;i++)
      printf("z[%d]=%.32e\n",i,z[i]);

      }*/

  // Note FP : shouldn't we call mlcp_driver_init here in order to avoid explicit call by user ?
  if(!options->dWork)
    mlcp_driver_allocate_work_arrays(problem, options);

  switch(options->solverId)
  {
  case  SICONOS_MLCP_PGS:/****** PGS algorithm ******/
    mlcp_pgs(problem, z, w, &info, options);
    break;
  case  SICONOS_MLCP_PGS_SBM:/****** PGS algorithm ******/
    mlcp_pgs_SBM(problem, z, w, &info, options);
    break;
  case SICONOS_MLCP_RPGS:
    /****** RPGS algorithm ******/
    mlcp_rpgs(problem, z, w, &info, options);
    break;

  /****** PSOR algorithm ******/
  case SICONOS_MLCP_PSOR:
    mlcp_psor(problem, z, w, &info, options);
    break;

  /****** RPSOR algorithm ******/
  case SICONOS_MLCP_RPSOR:
    mlcp_rpsor(problem, z, w, &info, options);
    break;

  /****** PATH algorithm ******/
  case SICONOS_MLCP_PATH:
    mlcp_path(problem, z, w, &info, options);
    break;

  /****** ENUM algorithm ******/
  case  SICONOS_MLCP_ENUM:
  {
    numerics_printf(" ========================== Call ENUM solver for Mixed Linear Complementarity Problem (MLCP )===========\n");
    mlcp_enum(problem, z, w, &info, options);
    break;
  }
  /****** SIMPLEX algorithm ******/
  case SICONOS_MLCP_SIMPLEX:
    mlcp_simplex(problem, z, w, &info, options);
    break;

  /****** DIRECT ENUM algorithm ******/
  case SICONOS_MLCP_DIRECT_ENUM:
    mlcp_direct_enum(problem, z, w, &info, options);
    break;
  case SICONOS_MLCP_PATH_ENUM:
    mlcp_path_enum(problem, z, w, &info, options);
    break;

  /****** DIRECT SIMPLEX algorithm ******/
  case SICONOS_MLCP_DIRECT_SIMPLEX:
    mlcp_direct_simplex(problem, z, w, &info, options);
    break;

  /****** DIRECT PATH algorithm ******/
  case SICONOS_MLCP_DIRECT_PATH:
    mlcp_direct_path(problem, z, w, &info, options);
    break;
  case SICONOS_MLCP_DIRECT_PATH_ENUM:
    mlcp_direct_path_enum(problem, z, w, &info, options);
    break;

  /****** FB algorithm ******/
  case SICONOS_MLCP_FB :
    mlcp_FB(problem, z, w, &info, options);
    break;
  /****** DIRECT FB algorithm ******/
  case  SICONOS_MLCP_DIRECT_FB :
    mlcp_direct_FB(problem, z, w, &info, options);
    break;
  // need a svn add mlcp_GaussSeidel_SBM ...
  //  else if( strcmp( name , SICONOS_MLCP_MLCP_SBM ) == 0 )
  //    mlcp_GaussSeidel_SBM( problem, z , w , &info , options,1);

  /*error */
  default:
  {
    fprintf(stderr, "mlcp_driver error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);
  }
  }
  return info;
}



