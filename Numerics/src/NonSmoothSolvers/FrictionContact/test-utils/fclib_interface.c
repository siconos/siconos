/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "NumericsConfig.h"


#ifdef WITH_FCLIB
#include <fclib.h>

FrictionContactProblem* frictionContact_fclib_read(const char *path)
{

  struct fclib_local   *fclib_problem;

  fclib_problem = fclib_read_local(path);

  if (fclib_problem == NULL)
  {
    return NULL;
  }

  FrictionContactProblem* problem;

  problem = malloc(sizeof(FrictionContactProblem));

  problem->dimension = fclib_problem->spacedim;
  problem->mu = fclib_problem->mu;
  problem->q = fclib_problem->q;

  problem->numberOfContacts = fclib_problem->W->m / fclib_problem->spacedim; /* cf fclib spec */

  problem->M = (NumericsMatrix *)malloc(sizeof(NumericsMatrix));

  problem->M->storageType = 1; /* sparse */
  problem->M->size0 = fclib_problem->W->m;
  problem->M->size1 = fclib_problem->W->n;

  problem->M->matrix0 = NULL;
  problem->M->matrix1 = (SparseBlockStructuredMatrix*)malloc(sizeof(SparseBlockStructuredMatrix));

  problem->M->matrix1->block = NULL;
  problem->M->matrix1->index1_data = NULL;
  problem->M->matrix1->index2_data = NULL;

  sparseToSBM(problem->dimension, (SparseMatrix*)fclib_problem->W, problem->M->matrix1);

  return problem;

}

int frictionContact_fclib_write(FrictionContactProblem* problem, char * title, char * description, char * math_info,
                                const char *path)
{
  int info = 0;

  struct fclib_local   *fclib_problem;

  fclib_problem = malloc(sizeof(struct fclib_local));

  fclib_problem->spacedim = problem->dimension;
  fclib_problem->mu =  problem->mu;
  fclib_problem->q =  problem->q;



  fclib_problem->s =  NULL;

  fclib_problem->info = malloc(sizeof(struct fclib_info)) ;
  fclib_problem->info->title = title;
  fclib_problem->info->description = description;
  fclib_problem->info->math_info = math_info;

  fclib_problem->W = malloc(sizeof(struct fclib_matrix));
  fclib_problem->R = NULL;
  fclib_problem->V = NULL;

  fclib_problem->W->m = problem->M->size0;
  fclib_problem->W->n = problem->M->size1;
  fclib_problem->W->nz = -2;

  SparseMatrix * spmat ;
  if (problem ->M->storageType == 0) /* Dense Matrix */
  {
    fclib_problem->W->nzmax = problem->M->size0 * problem->M->size1;
    fclib_problem->W->p = (int*)malloc((fclib_problem->W->m + 1) * sizeof(int));
    fclib_problem->W->i = (int*)malloc((fclib_problem->W->nzmax) * sizeof(int));
    fclib_problem->W->x = (double*)malloc((fclib_problem->W->nzmax) * sizeof(double));
    for (int i = 0; i <  problem ->M->size0 ; i++)
    {
      fclib_problem->W->p[i] = i * problem ->M->size1;
      for (int j = 0; j <  problem ->M->size1 ; j++)
      {
        fclib_problem->W->x[i * problem ->M->size1 + j ] = problem ->M->matrix0[j * problem ->M->size0 + i  ];
      }
    }
    fclib_problem->W->p[fclib_problem->W->m + 1] = (fclib_problem->W->m + 1) * problem ->M->size1;

  }
  else if (problem ->M->storageType == 1) /* Sparse block storage */
  {
    spmat = malloc(sizeof(SparseMatrix));
    int res = SBMtoSparseInitMemory(problem ->M->matrix1, spmat);
    res = SBMtoSparse(problem->M->matrix1, spmat);
    fclib_problem->W->nzmax = spmat->nzmax;
    fclib_problem->W->x = spmat->x;
    fclib_problem->W->p = spmat->p;
    fclib_problem->W->i = spmat->i;
  }
  else
  {
    fprintf(stderr, "frictionContact_fclib_write, unknown storage type for A.\n");
    exit(EXIT_FAILURE); ;
  }

  info = fclib_write_local(fclib_problem, path);

  /*   fclib_delete_local (fclib_problem); */

  if (problem ->M->storageType == 0) /* Dense Matrix */
  {
    free(fclib_problem->W->p);
    free(fclib_problem->W->i);
    free(fclib_problem->W->x);
  }
  else if (problem ->M->storageType == 1)
  {
    freeSparse(spmat);
  }
  free(fclib_problem->W);
  free(fclib_problem->info);
  free(fclib_problem);

  return info;

}

#endif
