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



void int_to_csi(int* o, csi* d, unsigned int n);
void int_to_csi(int* o, csi* d, unsigned int n)
{
  for(unsigned int i=0; i<n; ++i)
  {
    d[i] = (csi) o[i];
  }
}

void csi_to_int(csi* o, int* d, unsigned int n);
void csi_to_int(csi* o, int* d, unsigned int n)
{
  for(unsigned int i=0; i<n; ++i)
  {
    d[i] = (int) o[i];
  }
}



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

  CSparseMatrix W;

  W.nzmax = (csi) fclib_problem->W->nzmax;
  W.m = (csi) fclib_problem->W->m;
  W.n = (csi) fclib_problem->W->n;

  if (fclib_problem->W->nz == -1 || fclib_problem->W->nz == -2 ) /* We assume that the matrix is symmetric and square */
  {
    /* compressed colums */
    W.p = (csi*) malloc(sizeof(csi)*(W.n+1));
    int_to_csi(fclib_problem->W->p, W.p, (unsigned) (W.n+1));
  }
  /* else if (fclib_problem->W->nz == -2) */
  /* { */
  /*   /\* compressed rows *\/ */
  /*   fprintf(stderr, "from_fclib_local not implemented for csr matrices.\n"); */
  /*   exit(EXIT_FAILURE); ; */
  /* } */
  else
  {
    /* triplet */
    W.p = (csi*) malloc(sizeof(csi)*W.nzmax);
    int_to_csi(fclib_problem->W->p, W.p, (unsigned) W.nzmax);
  }

  W.i = (csi*) malloc(sizeof(csi)*W.nzmax);
  int_to_csi(fclib_problem->W->i, W.i, (unsigned) W.nzmax);

  W.x = fclib_problem->W->x;

  W.nz = fclib_problem->W->nz;

  sparseToSBM(problem->dimension, &W, problem->M->matrix1);

  free(W.p);
  free(W.i);

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
                                const char *path, int ndof)
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
    fclib_problem->W->m = (int) spmat->m;
    fclib_problem->W->n = (int) spmat->n;
    fclib_problem->W->x = spmat->x;



    if (spmat->nz == -1)
    {
      fclib_problem->W->p = (int*) malloc(sizeof(int)*(spmat->n+1));
      csi_to_int(spmat->p, fclib_problem->W->p, (unsigned) (spmat->n+1));
    }
    else if (spmat->nz == -2)
    {
      fclib_problem->W->p = (int*) malloc(sizeof(int)*(spmat->m+1));
      csi_to_int(spmat->p, fclib_problem->W->p, (unsigned) (spmat->m+1));
    }
    else
    {
      fclib_problem->W->p = (int*) malloc(sizeof(int)*spmat->nzmax);
      csi_to_int(spmat->p, fclib_problem->W->p, (unsigned) spmat->nzmax);
    }

    fclib_problem->W->i = (int*) malloc(sizeof(int)*spmat->nzmax);
    csi_to_int(spmat->i, fclib_problem->W->i, (unsigned) spmat->nzmax);

    fclib_problem->W->info = NULL;
  }
  else
  {
    fprintf(stderr, "frictionContact_fclib_write, unknown storage type for A.\n");
    exit(EXIT_FAILURE); ;
  }

  info = fclib_write_local(fclib_problem, path);

  info = fclib_create_int_attributes_in_info(path, "numberOfDegreeOfFreedom",
                                             ndof);



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
  free(fclib_problem->W->p);
  free(fclib_problem->W->i);
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
  problem->M->matrix0 = NULL;
  problem->M->matrix1 = NULL;

  CSparseMatrix * M = (CSparseMatrix*)malloc(sizeof(CSparseMatrix));
  M->nzmax = (csi) fclib_problem->M->nzmax;
  M->m = (csi) fclib_problem->M->m;
  M->n = (csi) fclib_problem->M->n;
  M->nz = (csi) fclib_problem->M->nz;
  M->x =  fclib_problem->M->x;
  
  if (fclib_problem->M->nz == -1)
  {
    /* compressed colums */
    problem->M->matrix2->csc= M;
    problem->M->matrix2->triplet=NULL;
    M->p = (csi*) malloc(sizeof(csi)*(M->n+1));
    int_to_csi(fclib_problem->M->p, M->p, (unsigned) (M->n+1));
  }
  else if (fclib_problem->M->nz == -2)
  {
    /* compressed rows */

    fprintf(stderr, "from_fclib_local not implemented for csr matrices.\n");
    exit(EXIT_FAILURE); ;
  }
  else
  {
    /* triplet */
    problem->M->matrix2->triplet=M;
    problem->M->matrix2->csc=NULL;
    M->p = (csi*) malloc(sizeof(csi)*M->nzmax);
    int_to_csi(fclib_problem->M->p, M->p, (unsigned) M->nzmax);
  }
  M->i = (csi*) malloc(sizeof(csi)*M->nzmax);
  int_to_csi(fclib_problem->M->i, M->i, (unsigned) M->nzmax);

  
  problem->H = newNumericsMatrix();
  problem->H->storageType = 2; /* sparse */
  problem->H->size0 = fclib_problem->H->m;
  problem->H->size1 = fclib_problem->H->n;
  problem->H->matrix2 = newNumericsSparseMatrix();
  problem->H->matrix0 = NULL;
  problem->H->matrix1 = NULL;
  
  CSparseMatrix * H = (CSparseMatrix*)malloc(sizeof(CSparseMatrix));;

  H->nzmax = (csi) fclib_problem->H->nzmax;
  H->m = (csi) fclib_problem->H->m;
  H->n = (csi) fclib_problem->H->n;
  H->nz = (csi) fclib_problem->H->nz;
  H->x =  fclib_problem->H->x;

  if (fclib_problem->H->nz == -1)
  {
    /* compressed colums */
    problem->H->matrix2->csc= H;
    problem->H->matrix2->triplet=NULL;
    H->p = (csi*) malloc(sizeof(csi)*(H->n+1));
    int_to_csi(fclib_problem->H->p, H->p, (unsigned) (H->n+1));
  }
  else if (fclib_problem->H->nz == -2)
  {
    /* compressed rows */
    fprintf(stderr, "from_fclib_local not implemented for csr matrices.\n");
    exit(EXIT_FAILURE); ;
  }
  else
  {
    /* triplet */
    problem->H->matrix2->triplet=H;
    problem->H->matrix2->csc=NULL;
    H->p = (csi*) malloc(sizeof(csi)*H->nzmax);
    int_to_csi(fclib_problem->H->p, H->p, (unsigned) H->nzmax);
  }

  H->i = (csi*) malloc(sizeof(csi)*H->nzmax);
  int_to_csi(fclib_problem->H->i, H->i, (unsigned) H->nzmax);



  
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

  fclib_problem->M->p= (int*) malloc(sizeof(int)*(problem->M->matrix2->triplet->nzmax));
  csi_to_int(problem->M->matrix2->triplet->p, fclib_problem->M->p,
             (unsigned) problem->M->matrix2->triplet->nzmax);
  fclib_problem->M->i= (int*) malloc(sizeof(int)*(problem->M->matrix2->triplet->nzmax));
  csi_to_int(problem->M->matrix2->triplet->i, fclib_problem->M->i,
             (unsigned) problem->M->matrix2->triplet->nzmax);
  fclib_problem->M->x= problem->M->matrix2->triplet->x;
  fclib_problem->M->nz= (int) problem->M->matrix2->triplet->nz;
  fclib_problem->M->info=NULL;
  fclib_problem->H = malloc(sizeof(struct fclib_matrix));
  fclib_problem->H->n = (int) problem->H->matrix2->triplet->n;
  fclib_problem->H->m = (int) problem->H->matrix2->triplet->m;
  fclib_problem->H->nzmax= (int) problem->H->matrix2->triplet->nzmax;
  fclib_problem->H->p= (int*) malloc(sizeof(int)*problem->H->matrix2->triplet->nzmax);
  csi_to_int(problem->H->matrix2->triplet->p, fclib_problem->H->p,
             (unsigned) problem->H->matrix2->triplet->nzmax);
  fclib_problem->H->i= (int*) malloc(sizeof(int)*problem->H->matrix2->triplet->nzmax);
  csi_to_int(problem->H->matrix2->triplet->i, fclib_problem->H->i,
             (unsigned) problem->H->matrix2->triplet->nzmax);
  fclib_problem->H->x= problem->H->matrix2->triplet->x;
  fclib_problem->H->nz= (int) problem->H->matrix2->triplet->nz;
  fclib_problem->H->info=NULL;


  fclib_problem->G = NULL;
  fclib_problem->b = NULL;



  DEBUG_PRINT("write in fclib of fclib_problem\n");

  rinfo = fclib_write_global(fclib_problem, path);
  DEBUG_PRINT("end of write in fclib of fclib_problem\n");

  free(fclib_problem->M->p);
  free(fclib_problem->M->i);
  free(fclib_problem->H->p);
  free(fclib_problem->H->i);
  free(fclib_problem->H);
  free(fclib_problem->M);
  free(fclib_problem->info);
  free(fclib_problem);

  return rinfo;

}

#endif
