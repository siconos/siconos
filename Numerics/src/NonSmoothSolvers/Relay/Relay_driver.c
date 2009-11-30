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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include <time.h>

int relay_driver(Relay_Problem* problem, double *z , double *w, Solver_Options* options, int numberOfSolvers, Numerics_Options* global_options)
{


  //Relay_display(problem);

  if (options == NULL || global_options == NULL)
    numericsError("Relay_driver", "null input for solver and/or global options");

  /* Set global options */
  setNumericsOptions(global_options);

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("Relay_driver", "null input for Relay_Problem and/or unknowns (z,w)");

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;

  /* Switch to DenseMatrix or SparseBlockMatrix solver according to the type of storage for M */
  /* Storage type for the matrix M of the LCP */
  int storageType = problem->M->storageType;

  /* Sparse Block Storage */
  if (storageType == 1)
  {
    numericsError("Relay_driver", "not yet implemented for sparse storage.");
  }
  // else

  /*************************************************
   *  2 - Call specific solver (if no trivial sol.)
   *************************************************/

  /* Solver name */
  char * name = options->solverName;

  if (verbose == 1)
    printf(" ========================== Call %s solver for Relayproblem ==========================\n", name);

  /****** NLGS algorithm ******/
  if (strcmp(name , "PGS") == 0)
    relay_pgs(problem, z , w , &info , options);

  else if (strcmp(name , "NLGS") == 0)
    fprintf(stderr, "Relay_driver error: NLGS solver obsolete use PGS:\n");
  else if ((strcmp(name , "Lemke") == 0) || (strcmp(name , "ENUM") == 0))
  {
    fprintf(stderr, "Relay_driver : Lemke and ENUM solver  not yet compeleted.  works only for ub=1, lb =-1:\n");
    // conversion into LCP
    LinearComplementarity_Problem* lcp_problem = (LinearComplementarity_Problem*)malloc(sizeof(LinearComplementarity_Problem));
    lcp_problem->size = 2 * problem->size ;
    lcp_problem->M = (NumericsMatrix *)malloc(sizeof(NumericsMatrix));
    lcp_problem->M->size0 = 2 * problem->size ;
    lcp_problem->M->size1 = 2 * problem->size ;

    lcp_problem->M->storageType = 0;
    lcp_problem->M->matrix1 = NULL;
    lcp_problem->M->matrix0 = (double*)malloc(lcp_problem->size * lcp_problem->size * sizeof(double));;
    lcp_problem->q = (double*)malloc(lcp_problem->size * sizeof(double));

    double *zlcp = (double*)malloc(lcp_problem->size * sizeof(double));
    double *wlcp = (double*)malloc(lcp_problem->size * sizeof(double));

    int i, j;

    for (i = 0; i < problem->size; i++)
    {
      for (j = 0; j < problem->size; j++)
      {
        lcp_problem->M->matrix0[i + j * lcp_problem->size] =  problem->M->matrix0[i + j * problem->size];
      }
    }
    for (i = 0; i < problem->size; i++)
    {
      for (j = problem->size; j < 2 * problem->size; j++)
      {
        lcp_problem->M->matrix0[i + j * lcp_problem->size] =  0.0;
      }
      lcp_problem->M->matrix0[i + (i + problem->size)*lcp_problem->size] =  1.0;
    }
    for (i = problem->size; i < 2 * problem->size; i++)
    {
      for (j = 0; j < 2 * problem->size; j++)
      {
        lcp_problem->M->matrix0[i + j * lcp_problem->size] =  0.0;
      }
      lcp_problem->M->matrix0[i + (i - problem->size)*lcp_problem->size] =  -1.0;
    }

    for (i = 0; i < problem->size; i++)
    {
      lcp_problem->q[i] = problem->q[i];
      lcp_problem->q[i + problem->size] = 2.0;
      for (j = 0; j < problem->size; j++)
      {
        lcp_problem->q[i] -= problem->M->matrix0[i + j * (problem->size)];
      }
    }
    int nbSolvers = 1;
    // Call the lcp_solver
    if ((strcmp(name , "ENUM") == 0))
    {
      lcp_enum_init(lcp_problem, options, 1);


    }
    info = lcp_driver(lcp_problem, zlcp , wlcp, options, nbSolvers, global_options);
    if ((strcmp(name , "ENUM") == 0))
    {
      lcp_enum_reset(lcp_problem, options, 1);

    }
    /*       printf("\n"); */

    // Conversion of result
    for (i = 0; i < problem->size; i++)
    {
      z[i] = 1.0 / 2.0 * (zlcp[i] - wlcp[i + problem->size]);
      //   printf("z[ %i]=%12.10e\n", i, z[i]);

      w[i] = wlcp[i] - zlcp[i + problem->size];
      //printf("w[ %i]=%12.10e\n", i, w[i]);
    }

    /*        for (i=0; i< lcp_problem->size; i++){  */
    /*     printf("zlcp[ %i]=%12.10e,\t wlcp[ %i]=%12.10e \n", i, zlcp[i],i, wlcp[i]); */
    /*        } */
    /*        printf("\n"); */

    /*        for (i=0; i< problem->size; i++){  */
    /*     printf("z[ %i]=%12.10e,\t w[ %i]=%12.10e\n", i, z[i],i, w[i]); */
    /*        } */


    /*        printf("\n"); */



    free(zlcp);
    free(wlcp);
    free(lcp_problem->q);
    freeNumericsMatrix(lcp_problem->M);
    free(lcp_problem->M);
    free(lcp_problem);


  }
  else if (strcmp(name , "Path") == 0)
  {
    relay_path(problem, z , w , &info , options);
  }
  /*error */
  else
  {
    fprintf(stderr, "Relay_driver error: unknown solver name: %s\n", name);
    exit(EXIT_FAILURE);
  }
  if (options[0].filterOn > 0)
    info = relay_compute_error(problem, z, w, options[0].dparam[0], &(options[0].dparam[1]));

  return info;
}

