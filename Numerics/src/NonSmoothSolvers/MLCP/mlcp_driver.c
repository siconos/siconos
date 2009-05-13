/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
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

void mlcp_driver_init(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  char * name = options->solverName;
  if (strcmp(name , "DIRECT_ENUM") == 0)
    mlcp_direct_enum_init(problem, options);
  else if (strcmp(name , "DIRECT_SIMPLEX") == 0)
    mlcp_direct_simplex_init(problem, options);
  else if (strcmp(name , "DIRECT_PATH") == 0)
    mlcp_direct_path_init(problem, options);
  else if (strcmp(name , "DIRECT_FB") == 0)
    mlcp_direct_FB_init(problem, options);
  else if (strcmp(name , "SIMPLEX") == 0)
    mlcp_simplex_init(problem, options);
  else if (strcmp(name , "FB") == 0)
    mlcp_FB_init(problem, options);
}
void mlcp_driver_reset(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  char * name = options->solverName;
  if (strcmp(name , "DIRECT_ENUM") == 0)
    mlcp_direct_enum_reset();
  else if (strcmp(name , "DIRECT_SIMPLEX") == 0)
    mlcp_direct_simplex_reset();
  else if (strcmp(name , "DIRECT_PATH") == 0)
    mlcp_direct_path_reset();
  else if (strcmp(name , "DIRECT_FB") == 0)
    mlcp_direct_FB_reset();
  else if (strcmp(name , "SIMPLEX") == 0)
    mlcp_simplex_reset();
}
int mlcp_driver_get_iwork(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  char * name = options->solverName;
  if (strcmp(name , "DIRECT_ENUM") == 0)
    return  mlcp_direct_enum_getNbIWork(problem, options);
  else if (strcmp(name , "ENUM") == 0)
    return  mlcp_enum_getNbIWork(problem, options);
  else if (strcmp(name , "DIRECT_SIMPLEX") == 0)
    return  mlcp_direct_simplex_getNbIWork(problem, options);
  else if (strcmp(name , "DIRECT_PATH") == 0)
    return  mlcp_direct_path_getNbIWork(problem, options);
  else if (strcmp(name , "FB") == 0)
    return  mlcp_FB_getNbIWork(problem, options);
  else if (strcmp(name , "DIRECT_FB") == 0)
    return  mlcp_direct_FB_getNbIWork(problem, options);
  return 0;
}
int mlcp_driver_get_dwork(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  char * name = options->solverName;
  if (strcmp(name , "DIRECT_ENUM") == 0)
    return  mlcp_direct_enum_getNbDWork(problem, options);
  else if (strcmp(name , "ENUM") == 0)
    return  mlcp_enum_getNbDWork(problem, options);
  else if (strcmp(name , "DIRECT_SIMPLEX") == 0)
    return  mlcp_direct_simplex_getNbDWork(problem, options);
  else if (strcmp(name , "DIRECT_PATH") == 0)
    return  mlcp_direct_path_getNbDWork(problem, options);
  else if (strcmp(name , "FB") == 0)
    return  mlcp_FB_getNbDWork(problem, options);
  else if (strcmp(name , "DIRECT_FB") == 0)
    return  mlcp_direct_FB_getNbDWork(problem, options);
  return 0;
}

int mlcp_driver(MixedLinearComplementarity_Problem* problem, double *z, double *w, Solver_Options* options, Numerics_Options* global_options)
{


  if (options == NULL || global_options == NULL)
    numericsError("mlcp_driver", "null input for solver and/or global options");

  /* Set global options */
  setNumericsOptions(global_options);

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("mlcp_driver", "null input for LinearComplementarity_Problem and/or unknowns (z,w)");
  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;
  //  displayMLCP(problem);
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
  char * name = options->solverName;

  /*  if(verbose==1){
    printf(" ========================== Call %s solver ==========================\n", name);
    printf("Initial z value:\n");
    for (i=0;i<problem->n+problem->m;i++)
      printf("z[%d]=%.32e\n",i,z[i]);

      }*/
  /****** PGS algorithm ******/
  if (strcmp(name , "PGS") == 0)
    mlcp_pgs(problem, z , w , &info , options);

  /****** RPGS algorithm ******/
  else if (strcmp(name , "RPGS") == 0)
    mlcp_rpgs(problem, z , w , &info , options);

  /****** PSOR algorithm ******/
  else if (strcmp(name , "PSOR") == 0)
    mlcp_psor(problem, z , w , &info , options);

  /****** RPSOR algorithm ******/
  else if (strcmp(name , "RPSOR") == 0)
    mlcp_rpsor(problem, z , w , &info , options);

  /****** PATH algorithm ******/
  else if (strcmp(name , "PATH") == 0)
    mlcp_path(problem, z , w , &info , options);

  /****** ENUM algorithm ******/
  else if (strcmp(name , "ENUM") == 0)
    mlcp_enum(problem, z , w , &info , options);

  /****** SIMPLEX algorithm ******/
  else if (strcmp(name , "SIMPLEX") == 0)
    mlcp_simplex(problem, z , w , &info , options);

  /****** DIRECT ENUM algorithm ******/
  else if (strcmp(name , "DIRECT_ENUM") == 0)
    mlcp_direct_enum(problem, z , w , &info , options);

  /****** DIRECT SIMPLEX algorithm ******/
  else if (strcmp(name , "DIRECT_SIMPLEX") == 0)
    mlcp_direct_simplex(problem, z , w , &info , options);

  /****** DIRECT PATH algorithm ******/
  else if (strcmp(name , "DIRECT_PATH") == 0)
    mlcp_direct_path(problem, z , w , &info , options);

  /****** FB algorithm ******/
  else if (strcmp(name , "FB") == 0)
    mlcp_FB(problem, z , w , &info , options);
  /****** DIRECT FB algorithm ******/
  else if (strcmp(name , "DIRECT_FB") == 0)
    mlcp_direct_FB(problem, z , w , &info , options);

  /*error */
  else
  {
    fprintf(stderr, "mlcp_driver error: unknown solver named: %s\n", name);
    exit(EXIT_FAILURE);
  }

  return info;
}



