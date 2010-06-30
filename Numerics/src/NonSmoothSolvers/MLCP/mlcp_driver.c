/* Siconos-Numerics, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/

#include "MLCP_Solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include "mlcp_cst.h"

char  SICONOS_NONAME_STR[] = "NONAME";
char  SICONOS_MLCP_PGS_STR[] = "MLCP_PGS";
char  SICONOS_MLCP_RPGS_STR[] = "MLCP_RPGS";
char  SICONOS_MLCP_PSOR_STR[] = "MLCP_PSOR";
char  SICONOS_MLCP_RPSOR_STR[] = "MLCP_RPSOR";
char  SICONOS_MLCP_PATH_STR[] = "MLCP_PATH";
char  SICONOS_MLCP_ENUM_STR[] = "MLCP_ENUM";
char  SICONOS_MLCP_SIMPLEX_STR[] = "MLCP_SIMPLEX";
char  SICONOS_MLCP_DIRECT_ENUM_STR[] = "MLCP_DIRECT_ENUM";
char  SICONOS_MLCP_PATH_ENUM_STR[] = "MLCP_PATH_ENUM";
char  SICONOS_MLCP_DIRECT_SIMPLEX_STR[] = "MLCP_DIRECT_SIMPLEX";
char  SICONOS_MLCP_DIRECT_PATH_STR[] = "MLCP_DIRECT_PATH";
char  SICONOS_MLCP_DIRECT_PATH_ENUM_STR[] = "MLCP_DIRECT_PATH_ENUM";
char  SICONOS_MLCP_FB_STR[] = "MLCP_FB";
char  SICONOS_MLCP_DIRECT_FB_STR[] = "MLCP_DIRECT_FB";


void mlcp_driver_init(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  //char * name = options->solverName;

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
    break;
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
    break;
  case SICONOS_MLCP_PATH_ENUM:
    return  mlcp_path_enum_getNbIWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_PATH_ENUM:
    return mlcp_direct_path_enum_getNbIWork(problem, options);
    break;
  case SICONOS_MLCP_ENUM:
    return  mlcp_enum_getNbIWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_SIMPLEX:
    return  mlcp_direct_simplex_getNbIWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_PATH:
    return  mlcp_direct_path_getNbIWork(problem, options);
    break;
  case SICONOS_MLCP_FB:
    return  mlcp_FB_getNbIWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_FB:
    return  mlcp_direct_FB_getNbIWork(problem, options);
    break;
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
    break;
  case SICONOS_MLCP_PATH_ENUM:
    return  mlcp_path_enum_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_PATH_ENUM:
    return mlcp_direct_path_enum_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_ENUM:
    return  mlcp_enum_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_SIMPLEX:
    return  mlcp_direct_simplex_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_PATH:
    return  mlcp_direct_path_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_FB:
    return  mlcp_FB_getNbDWork(problem, options);
    break;
  case SICONOS_MLCP_DIRECT_FB:
    return  mlcp_direct_FB_getNbDWork(problem, options);
    break;
  default :
    return 0;
  }

  /* char * name = options->solverName; */
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

int mlcp_driver(MixedLinearComplementarityProblem* problem, double *z, double *w, SolverOptions* options, NumericsOptions* global_options)
{


  if (options == NULL || global_options == NULL)
    numericsError("mlcp_driver", "null input for solver and/or global options");

  /* Set global options */
  //  setNumericsOptions(global_options);

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("mlcp_driver", "null input for LinearComplementarityProblem and/or unknowns (z,w)");
  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;
  if (verbose)
    displayMLCP(problem);
  /* Switch to DenseMatrix or SparseBlockMatrix solver according to the type of storage for M */
  /* Storage type for the matrix M of the LCP */
  int storageType = problem->M->storageType;

  /* Sparse Block Storage */
  if (storageType == 1)
  {
    numericsError("mlcp_driver", "not yet implemented for sparse storage.");
  }
  // else

  /*************************************************
   *  2 - Call specific solver (if no trivial sol.)
   *************************************************/

  /* Solver name */
  //  char * name = options->solverName;

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



