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

#include "MLCP_Solvers.h"
#include "SiconosCompat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include "mlcp_cst.h"
#include "numerics_verbose.h"
#include "NumericsMatrix.h"

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


int mlcp_alloc_working_memory(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  switch (options->solverId)
  {
  case SICONOS_MLCP_ENUM :
    return mlcp_enum_alloc_working_memory(problem, options);
  default:

    return 0;/*Nothing to do*/
  }
}
void mlcp_free_working_memory(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  switch (options->solverId)
  {
  case SICONOS_MLCP_ENUM :
    mlcp_enum_free_working_memory(problem, options);
    break;
  default:
    ;/*Nothing to do*/
  }
}
void mlcp_driver_init(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  //const char* const  name = options->solverName;

  switch (options->solverId)
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


  switch (options->solverId)
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
int mlcp_driver_get_iwork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{

  switch (options->solverId)
  {
  case SICONOS_MLCP_DIRECT_ENUM:
    return  mlcp_direct_enum_getNbIWork(problem, options);
  case SICONOS_MLCP_PATH_ENUM:
    return  mlcp_path_enum_getNbIWork(problem, options);
  case SICONOS_MLCP_DIRECT_PATH_ENUM:
    return mlcp_direct_path_enum_getNbIWork(problem, options);
  case SICONOS_MLCP_ENUM:
    return  mlcp_enum_getNbIWork(problem, options);
  case SICONOS_MLCP_DIRECT_SIMPLEX:
    return  mlcp_direct_simplex_getNbIWork(problem, options);
  case SICONOS_MLCP_DIRECT_PATH:
    return  mlcp_direct_path_getNbIWork(problem, options);
  case SICONOS_MLCP_FB:
    return  mlcp_FB_getNbIWork(problem, options);
  case SICONOS_MLCP_DIRECT_FB:
    return  mlcp_direct_FB_getNbIWork(problem, options);
  default :
    return 0;
  }
}
int mlcp_driver_get_dwork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{


  switch (options->solverId)
  {
  case SICONOS_MLCP_DIRECT_ENUM:
    return  mlcp_direct_enum_getNbDWork(problem, options);
  case SICONOS_MLCP_PATH_ENUM:
    return  mlcp_path_enum_getNbDWork(problem, options);
  case SICONOS_MLCP_DIRECT_PATH_ENUM:
    return mlcp_direct_path_enum_getNbDWork(problem, options);
  case SICONOS_MLCP_ENUM:
    return  mlcp_enum_getNbDWork(problem, options);
  case SICONOS_MLCP_DIRECT_SIMPLEX:
    return  mlcp_direct_simplex_getNbDWork(problem, options);
  case SICONOS_MLCP_DIRECT_PATH:
    return  mlcp_direct_path_getNbDWork(problem, options);
  case SICONOS_MLCP_FB:
    return  mlcp_FB_getNbDWork(problem, options);
  case SICONOS_MLCP_DIRECT_FB:
    return  mlcp_direct_FB_getNbDWork(problem, options);
  default :
    return 0;
  }

  /* const char* const  name = options->solverName; */
  /* if (strcmp(name , "DIRECT_ENUM") == 0) */
  /*   return  mlcp_direct_enum_getNbDWork(problem, options); */
  /* else if (strcmp(name , "PATH_ENUM") == 0) */
  /*   return  mlcp_path_enum_getNbDWork(problem, options); */
  /* else if (strcmp(name , "DIRECT_PATH_ENUM") == 0) */
  /*   return  mlcp_direct_path_enum_getNbDWork(problem, options); */
  /* else if (strcmp(name , "ENUM") == 0) */
  /*   return  mlcp_enum_getNbDWork(problem, options); */
  /* else if (strcmp(name , "DIRECT_SIMPLEX") == 0) */
  /*   return  mlcp_direct_simplex_getNbDWork(problem, options); */
  /* else if (strcmp(name , "DIRECT_PATH") == 0) */
  /*   return  mlcp_direct_path_getNbDWork(problem, options); */
  /* else if (strcmp(name , "FB") == 0) */
  /*   return  mlcp_FB_getNbDWork(problem, options); */
  /* else if (strcmp(name , "DIRECT_FB") == 0) */
  /*   return  mlcp_direct_FB_getNbDWork(problem, options); */
  /* return 0; */
}

int mlcp_driver(MixedLinearComplementarityProblem* problem, double *z, double *w, SolverOptions* options)
{


  if (options == NULL)
    numerics_error("mlcp_driver ", "null input for solver options.\n");

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numerics_error("mlcp_driver", "null input for MixedLinearComplementarityProblem and/or unknowns (z,w)");
  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;
  if (verbose)
    mixedLinearComplementarity_display(problem);
  /* Switch to DenseMatrix or SparseBlockMatrix solver according to the type of storage for M */
  /* Storage type for the matrix M of the LCP */
  int storageType = problem->M->storageType;

  /* Sparse Block Storage */
  if (storageType == 1)
  {
    numerics_error("mlcp_driver", "not yet implemented for sparse storage.");
  }
  // else

  /*************************************************
   *  2 - Call specific solver (if no trivial sol.)
   *************************************************/

  /* Solver name */
  //  const char* const  name = options->solverName;

  /*  if(verbose==1){
    printf(" ========================== Call %s solver ==========================\n", name);
    printf("Initial z value:\n");
    for (i=0;i<problem->n+problem->m;i++)
      printf("z[%d]=%.32e\n",i,z[i]);

      }*/
  switch (options->solverId)
  {
  case  SICONOS_MLCP_PGS:/****** PGS algorithm ******/
    mlcp_pgs(problem, z , w , &info , options);
    break;
  case  SICONOS_MLCP_PGS_SBM:/****** PGS algorithm ******/
    mlcp_pgs_SBM(problem, z , w , &info , options);
    break;
  case SICONOS_MLCP_RPGS:
    /****** RPGS algorithm ******/
    mlcp_rpgs(problem, z , w , &info , options);
    break;

    /****** PSOR algorithm ******/
  case SICONOS_MLCP_PSOR:
    mlcp_psor(problem, z , w , &info , options);
    break;

    /****** RPSOR algorithm ******/
  case SICONOS_MLCP_RPSOR:
    mlcp_rpsor(problem, z , w , &info , options);
    break;

    /****** PATH algorithm ******/
  case SICONOS_MLCP_PATH:
    mlcp_path(problem, z , w , &info , options);
    break;

    /****** ENUM algorithm ******/
  case  SICONOS_MLCP_ENUM:
    mlcp_enum(problem, z , w , &info , options);
    break;

    /****** SIMPLEX algorithm ******/
  case SICONOS_MLCP_SIMPLEX:
    mlcp_simplex(problem, z , w , &info , options);
    break;

    /****** DIRECT ENUM algorithm ******/
  case SICONOS_MLCP_DIRECT_ENUM:
    mlcp_direct_enum(problem, z , w , &info , options);
    break;
  case SICONOS_MLCP_PATH_ENUM:
    mlcp_path_enum(problem, z , w , &info , options);
    break;

    /****** DIRECT SIMPLEX algorithm ******/
  case SICONOS_MLCP_DIRECT_SIMPLEX:
    mlcp_direct_simplex(problem, z , w , &info , options);
    break;

    /****** DIRECT PATH algorithm ******/
  case SICONOS_MLCP_DIRECT_PATH:
    mlcp_direct_path(problem, z , w , &info , options);
    break;
  case SICONOS_MLCP_DIRECT_PATH_ENUM:
    mlcp_direct_path_enum(problem, z , w , &info , options);
    break;

    /****** FB algorithm ******/
  case SICONOS_MLCP_FB :
    mlcp_FB(problem, z , w , &info , options);
    break;
    /****** DIRECT FB algorithm ******/
  case  SICONOS_MLCP_DIRECT_FB :
    mlcp_direct_FB(problem, z , w , &info , options);
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



