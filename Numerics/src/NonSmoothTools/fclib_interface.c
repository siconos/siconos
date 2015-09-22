/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#define DEBUG_MESSAGES

#include "debug.h"
#ifdef WITH_FCLIB
#include "fclib_interface.h"

FrictionContactProblem* from_fclib_local(const struct fclib_local* fclib_problem)
{
  FrictionContactProblem* problem;

  problem = malloc(sizeof(FrictionContactProblem));

  problem->dimension = fclib_problem->spacedim;
  problem->mu = fclib_problem->mu;
  problem->q = fclib_problem->q;

  problem->numberOfContacts = fclib_problem->W->m / fclib_problem->spacedim; /* cf fclib spec */

  problem->M = newNumericsMatrix();

  problem->M->storageType = 1; /* sparse */
  problem->M->size0 = fclib_problem->W->m;
  problem->M->size1 = fclib_problem->W->n;

  problem->M->matrix0 = NULL;
  problem->M->matrix1 = newSBM();

  problem->M->matrix1->block = NULL;
  problem->M->matrix1->index1_data = NULL;
  problem->M->matrix1->index2_data = NULL;

  sparseToSBM(problem->dimension, (CSparseMatrix*)fclib_problem->W, problem->M->matrix1);

  return problem;

}


FrictionContactProblem* frictionContact_fclib_read(const char *path)
{

  struct fclib_local   *fclib_problem;

  fclib_problem = fclib_read_local(path);

  if (!fclib_problem)
  {
    return NULL;
  }

  return from_fclib_local(fclib_problem);
}

int frictionContact_fclib_write(FrictionContactProblem* problem, char * title, char * description, char * mathInfo,
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
  fclib_problem->info->math_info = mathInfo;

  fclib_problem->W = malloc(sizeof(struct fclib_matrix));
  fclib_problem->R = NULL;
  fclib_problem->V = NULL;

  fclib_problem->W->m = problem->M->size0;
  fclib_problem->W->n = problem->M->size1;
  fclib_problem->W->nz = -2;

  CSparseMatrix * spmat = NULL;
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
    spmat = malloc(sizeof(CSparseMatrix));
    int MAYBE_UNUSED res = SBMtoSparseInitMemory(problem ->M->matrix1, spmat);
    res = SBMtoSparse(problem->M->matrix1, spmat);
    fclib_problem->W->nzmax = (int) spmat->nzmax;
    fclib_problem->W->x = spmat->x;
    fclib_problem->W->p = spmat->p;
    fclib_problem->W->i = spmat->i;
    fclib_problem->W->info = NULL;
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
    cs_spfree(spmat);
  }
  free(fclib_problem->W);
  free(fclib_problem->info);
  free(fclib_problem);

  return info;

}


GlobalFrictionContactProblem* from_fclib_global(const struct fclib_global* fclib_problem)
{
  GlobalFrictionContactProblem* problem;

  problem = malloc(sizeof(GlobalFrictionContactProblem));

  problem->dimension = fclib_problem->spacedim;
  problem->mu = fclib_problem->mu;
  problem->q = fclib_problem->f;
  problem->b = fclib_problem->w;

  problem->numberOfContacts = fclib_problem->H->n / fclib_problem->spacedim; /* cf fclib spec */

  problem->M = newNumericsMatrix();
  problem->M->storageType = 2; /* sparse */
  problem->M->size0 = fclib_problem->M->m;
  problem->M->size1 = fclib_problem->M->n;
  problem->M->matrix2 = newNumericsSparseMatrix();
  problem->M->matrix2->triplet=(CSparseMatrix*)fclib_problem->M;
  problem->M->matrix0 = NULL;
  problem->M->matrix1 = NULL;

  problem->H = newNumericsMatrix();
  problem->H->storageType = 2; /* sparse */
  problem->H->size0 = fclib_problem->H->m;
  problem->H->size1 = fclib_problem->H->n;
  problem->H->matrix2 = newNumericsSparseMatrix();
  problem->H->matrix2->triplet=(CSparseMatrix*)fclib_problem->H;
  problem->H->matrix0 = NULL;
  problem->H->matrix1 = NULL;

  return problem;

}


GlobalFrictionContactProblem* globalFrictionContact_fclib_read(const char *path)
{

  struct fclib_global   *fclib_problem;

  fclib_problem = fclib_read_global(path);

  if (!fclib_problem)
  {
    return NULL;
  }

  return from_fclib_global(fclib_problem);
}


int globalFrictionContact_fclib_write(
  GlobalFrictionContactProblem* problem,
  char * title,
  char * description,
  char * mathInfo,
  const char *path);
int globalFrictionContact_fclib_write(
  GlobalFrictionContactProblem* problem,
  char * title,
  char * description,
  char * mathInfo,
  const char *path)
{
  int rinfo = 0;


  globalFrictionContact_display(problem);
  FILE * file  =  fopen("toto.dat", "w");
  globalFrictionContact_printInFile(problem, file);
  DEBUG_PRINT("construcion of fclib_problem\n");
  struct fclib_global *fclib_problem;
  fclib_problem = malloc(sizeof(struct fclib_global));

  fclib_problem->info = malloc(sizeof(struct fclib_info)) ;

  fclib_problem->info->title = title;
  fclib_problem->info->description = description;
  fclib_problem->info->math_info = mathInfo;



  fclib_problem->spacedim = problem->dimension;
  fclib_problem->mu = problem->mu;
  fclib_problem->w =  problem->b;
  fclib_problem->f =  problem->q;

  /* only coordinates */
  assert(problem->M->matrix2);
  assert(problem->H->matrix2);

  fclib_problem->M = malloc(sizeof(struct fclib_matrix));
  fclib_problem->M->n = (int) problem->M->matrix2->triplet->n;
  fclib_problem->M->m = (int) problem->M->matrix2->triplet->m;
  fclib_problem->M->nzmax= (int) problem->M->matrix2->triplet->nzmax;
  fclib_problem->M->p= problem->M->matrix2->triplet->p;
  fclib_problem->M->i= problem->M->matrix2->triplet->i;
  fclib_problem->M->x= problem->M->matrix2->triplet->x;
  fclib_problem->M->nz= (int) problem->M->matrix2->triplet->nz;
  fclib_problem->M->info=NULL;
  fclib_problem->H = malloc(sizeof(struct fclib_matrix));
  fclib_problem->H->n = (int) problem->H->matrix2->triplet->n;
  fclib_problem->H->m = (int) problem->H->matrix2->triplet->m;
  fclib_problem->H->nzmax= (int) problem->H->matrix2->triplet->nzmax;
  fclib_problem->H->p= problem->H->matrix2->triplet->p;
  fclib_problem->H->i= problem->H->matrix2->triplet->i;
  fclib_problem->H->x= problem->H->matrix2->triplet->x;
  fclib_problem->H->nz= (int) problem->H->matrix2->triplet->nz;
  fclib_problem->H->info=NULL;


  fclib_problem->G = NULL;
  fclib_problem->b = NULL;



  DEBUG_PRINT("write in fclib of fclib_problem\n");

  rinfo = fclib_write_global(fclib_problem, path);
  DEBUG_PRINT("end of write in fclib of fclib_problem\n");
  /* free(fclib_problem->H); */
  /* free(fclib_problem->M); */
  /* free(fclib_problem->info); */
  /* free(fclib_problem); */

  return rinfo;

}

#endif
