/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <stdlib.h>
#include <assert.h>
#include "FrictionContactProblem.h"
#include "NumericsMatrix.h"
#include <stdio.h>
#include "numerics_verbose.h"
#include "SparseBlockMatrix.h"
#include <math.h>
#include <string.h>
#include "SiconosLapack.h"
//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

void frictionContact_display(FrictionContactProblem* problem)
{

  assert(problem);
  int n = problem->dimension * problem->numberOfContacts;
  printf("FrictionContact Display :\n-------------\n");
  printf("dimension :%d \n", problem->dimension);
  printf("numberOfContacts:%d \n", problem->numberOfContacts);

  if (problem->M)
  {
    printf("M matrix:\n");
    NM_display(problem->M);
  }
  else
    printf("No M matrix:\n");

  if (problem->q)
  {
    printf("q vector:\n");
    NM_vector_display(problem->q,n);
  }
  else
    printf("No q vector:\n");

  if (problem->mu)
  {
    printf("mu vector:\n");
    NM_vector_display(problem->mu,problem->numberOfContacts);
  }
  else
    printf("No mu vector:\n");

}





int frictionContact_printInFile(FrictionContactProblem*  problem, FILE* file)
{
  if (! problem)
  {
    fprintf(stderr, "Numerics, FrictionContactProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;

  int d  = problem->dimension;
  fprintf(file, "%d\n", d);
  int nc = problem->numberOfContacts;
  fprintf(file, "%d\n", nc);
  NM_write_in_file(problem->M, file);
  for (i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  for (i = 0; i < nc; i++)
  {
    fprintf(file, "%32.24e ", problem->mu[i]);
  }
  fprintf(file, "\n");
  return 0;
}

int frictionContact_printInFilename(FrictionContactProblem* problem, char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "w");

  if (!file)
  {
    return errno;
  }

  info = frictionContact_printInFile(problem, file);

  fclose(file);
  return info;
}

int frictionContact_newFromFile(FrictionContactProblem* problem, FILE* file)
{
  assert(file);
  DEBUG_PRINT("Start -- int frictionContact_newFromFile(FrictionContactProblem* problem, FILE* file)\n");
  int nc = 0, d = 0;
  int i;
  CHECK_IO(fscanf(file, "%d\n", &d));
  problem->dimension = d;
  DEBUG_PRINTF("problem->dimension = %i \n",problem->dimension );
  CHECK_IO(fscanf(file, "%d\n", &nc));
  problem->numberOfContacts = nc;
  problem->M =  NM_new_from_file(file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for (i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])));
  }

  problem->mu = (double *) malloc(nc * sizeof(double));
  for (i = 0; i < nc; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->mu[i])));
  }
  DEBUG_PRINT("End --  int frictionContact_newFromFile(FrictionContactProblem* problem, FILE* file)\n");

  return 0;
}

int frictionContact_newFromFilename(FrictionContactProblem* problem, char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "r");

  if (!file)
  {
    return errno;
  }

  info = frictionContact_newFromFile(problem, file);

  fclose(file);
  return info;
}

void frictionContactProblem_free(FrictionContactProblem* problem)
{
  assert(problem);
  if (problem->M)
  {
    NM_clear(problem->M);
    free(problem->M);
    problem->M = NULL;
  }

  if (problem->mu)
  {
  free(problem->mu);
  problem->mu = NULL;
  }

  if (problem->q)
  {
    free(problem->q);
    problem->q = NULL;
  }

  free(problem);

}

FrictionContactProblem* frictionContactProblem_new(void)
{
  FrictionContactProblem* fcp = (FrictionContactProblem*) malloc(sizeof(FrictionContactProblem));
  fcp->dimension = 0;
  fcp->numberOfContacts = 0;
  fcp->M = NULL;
  fcp->q = NULL;
  fcp->mu = NULL;

  return fcp;
}

FrictionContactProblem* frictionContactProblem_new_with_data(int dim, int nc, NumericsMatrix* M, double* q, double* mu)
{
  FrictionContactProblem* fcp = (FrictionContactProblem*) malloc(sizeof(FrictionContactProblem));

  fcp->dimension = dim;
  fcp->numberOfContacts = nc;
  fcp->M = M;
  fcp->q = q;
  fcp->mu = mu;

  return fcp;
}

//#define SN_SBM_TO_DENSE


void createSplittedFrictionContactProblem(FrictionContactProblem* problem,  SplittedFrictionContactProblem * splitted_problem)
{
  /* Number of contacts */
  int nc = problem->numberOfContacts;

  splitted_problem->fc3d = problem;

  /* Splitting vector q */

  splitted_problem->q_n = (double *) malloc(nc * sizeof(double));
  splitted_problem->q_t = (double *) malloc(2* nc * sizeof(double));

  for (int contact = 0 ; contact < nc; contact ++)
  {
    splitted_problem->q_n[contact] = problem->q[contact*3];
    splitted_problem->q_t[2*contact] = problem->q[contact*3+1];
    splitted_problem->q_t[2*contact+1] = problem->q[contact*3+2];
  }

  /* Splitting matrix M */

  int storageType = problem->M->storageType;
  NumericsMatrix * M =  problem->M;

  splitted_problem->M_nn =  NM_create(problem->M->storageType, nc, nc);
  splitted_problem->M_nt =  NM_create(problem->M->storageType, nc, 2*nc);
  splitted_problem->M_tn =  NM_create(problem->M->storageType, 2*nc, nc);
  splitted_problem->M_tt =  NM_create(problem->M->storageType, 2*nc, 2*nc);

  NumericsMatrix * M_nn =  splitted_problem->M_nn;
  NumericsMatrix * M_tn =  splitted_problem->M_tn;
  NumericsMatrix * M_nt =  splitted_problem->M_nt;
  NumericsMatrix * M_tt =  splitted_problem->M_tt;

  switch (storageType)
  {
  case NM_SPARSE_BLOCK:
  {
    DEBUG_PRINT("NM_SPARSE_BLOCK case\n");
    unsigned int row_components_nn[1] = {0};
    unsigned int row_components_size_nn =1;
    unsigned int col_components_nn[1] = {0};
    unsigned int col_components_size_nn =1;
    SBM_extract_component_3x3(M->matrix1, M_nn->matrix1,
                              row_components_nn, row_components_size_nn, col_components_nn, col_components_size_nn   );

    unsigned int row_components_nt[1] = {0};
    unsigned int row_components_size_nt =1;
    unsigned int col_components_nt[2] = {1,2};
    unsigned int col_components_size_nt =2;
    SBM_extract_component_3x3(M->matrix1, M_nt->matrix1,
                              row_components_nt, row_components_size_nt, col_components_nt, col_components_size_nt   );


    unsigned int row_components_tn[2] = {1,2};
    unsigned int row_components_size_tn =2;
    unsigned int col_components_tn[1] = {0};
    unsigned int col_components_size_tn =1;
    SBM_extract_component_3x3(M->matrix1, M_tn->matrix1,
                              row_components_tn, row_components_size_tn, col_components_tn, col_components_size_tn   );


    unsigned int row_components_tt[2] = {1,2};
    unsigned int row_components_size_tt =2;
    unsigned int col_components_tt[2] = {1,2};
    unsigned int col_components_size_tt =2;
    SBM_extract_component_3x3(M->matrix1, M_tt->matrix1,
                              row_components_tt, row_components_size_tt, col_components_tt, col_components_size_tt   );

#ifdef SN_SBM_TO_DENSE
    M_nn->matrix0 =(double *)malloc(M_nn->size0*M_nn->size1*sizeof(double));
    SBM_to_dense(M_nn->matrix1, M_nn->matrix0);
    M_nn->storageType=NM_DENSE;
    M_nt->matrix0 =(double *)malloc(M_nt->size0*M_nt->size1*sizeof(double));
    SBM_to_dense(M_nt->matrix1, M_nt->matrix0);
    M_nt->storageType=NM_DENSE;
    M_tn->matrix0 =(double *)malloc(M_tn->size0*M_tn->size1*sizeof(double));
    SBM_to_dense(M_tn->matrix1, M_tn->matrix0);
    M_tn->storageType=NM_DENSE;
    M_tt->matrix0 =(double *)malloc(M_tt->size0*M_tt->size1*sizeof(double));
    SBM_to_dense(M_tt->matrix1, M_tt->matrix0);
    M_tt->storageType=NM_DENSE;
#endif

    break;
  }
  default:
    numerics_error("createSplittedFrictionContactProblem", "storageType value %d not implemented yet !", storageType);
  }



}
void frictionContactProblem_compute_statistics(FrictionContactProblem* problem,
                                               double * reaction,
                                               double * velocity,

                                               double tol,
                                               int do_print)
{
  /* NumericsMatrix* M = problem->M; */
  /* double* q = problem->q; */
  double* mu = problem->mu;

  /* Number of contacts */
  int nc = problem->numberOfContacts;
  /* int n = problem->M->size0; */
  /* int m = 3 * nc; */


  int contact =0, pos =0;

  int take_off_count=0;
  int closed_count=0;

  int sticking_count=0;
  int sliding_count=0;
  int ambiguous_take_off_count=0;
  int ambiguous_closed_count=0;

  for (contact =0; contact < nc; contact ++)
  {
    pos=contact*3;
    if (velocity[pos] > tol)
    {
      take_off_count++;
      if (fabs(reaction[pos]) > tol)
        ambiguous_take_off_count++;
    }
    else if (reaction[pos] > tol)
    {
      closed_count++;
      if (fabs(velocity[pos]) > tol)
        ambiguous_closed_count++;

      double criteria = reaction[pos+1]*reaction[pos+1] + reaction[pos+2]*reaction[pos+2] - mu[contact]*mu[contact]*reaction[pos]*reaction[pos];
      if (criteria < 0 )
      {
        sticking_count++;
      }
      else
      {
        sliding_count++;
      }
    }
  }

  numerics_printf_verbose(do_print, "---- FC3D - STAT  - Number of contact = %i\t,  take off = %i\t, closed = %i\t, sliding = %i\t,sticking = %i\t, ambiguous take off = %i\t, ambiguous closed = %i",nc ,take_off_count,closed_count,sliding_count,sticking_count, ambiguous_take_off_count,ambiguous_closed_count);

}
FrictionContactProblem* frictionContact_copy(FrictionContactProblem* problem)
{
  assert(problem);

  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  FrictionContactProblem* new = (FrictionContactProblem*) malloc(sizeof(FrictionContactProblem));
  new->dimension = problem->dimension;
  new->numberOfContacts = problem->numberOfContacts;
  new->M = NM_new();
  NM_copy(problem->M, new->M);
  new->q = (double*)malloc(n*sizeof(double));
  memcpy(new->q, problem->q, n*sizeof(double));
  new->mu = (double*)malloc(nc*sizeof(double));
  memcpy(new->mu, problem->mu, nc*sizeof(double));
  return new;
}
void frictionContact_rescaling(
  FrictionContactProblem* problem,
  double alpha,
  double gamma)
{
  int n = problem->M->size0;
  /* scaling of M */
  NM_scal(alpha*gamma*gamma, problem->M);
  /* scaling of q */
  cblas_dscal(n,alpha*gamma,problem->q,1);

}
