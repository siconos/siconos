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
#include "RollingFrictionContactProblem.h"
#include "NumericsMatrix.h"
#include <stdio.h>
#include "numerics_verbose.h"
#include "SparseBlockMatrix.h"
#include <math.h>
//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

void rollingFrictionContact_display(RollingFrictionContactProblem* problem)
{

  assert(problem);
  int n = problem->dimension * problem->numberOfContacts;
  printf("RollingFrictionContact Display :\n-------------\n");
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
  
  if (problem->mu_r)
  {
    printf("mu_R vector:\n");
    NM_vector_display(problem->mu_r,problem->numberOfContacts);
  }
  else
    printf("No mu_R vector:\n");

}





int rollingFrictionContact_printInFile(RollingFrictionContactProblem*  problem, FILE* file)
{
  if (! problem)
  {
    fprintf(stderr, "Numerics, RollingFrictionContactProblem printInFile failed, NULL input.\n");
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
  for (i = 0; i < nc; i++)
  {
    fprintf(file, "%32.24e ", problem->mu_r[i]);
  }
  fprintf(file, "\n");
  return 0;
}

int rollingFrictionContact_printInFilename(RollingFrictionContactProblem* problem, char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "w");

  if (!file)
  {
    return errno;
  }

  info = rollingFrictionContact_printInFile(problem, file);

  fclose(file);
  return info;
}

int rollingFrictionContact_newFromFile(RollingFrictionContactProblem* problem, FILE* file)
{
  assert(file);
  DEBUG_PRINT("Start -- int rollingFrictionContact_newFromFile(RollingFrictionContactProblem* problem, FILE* file)\n");
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
  problem->mu_r = (double *) malloc(nc * sizeof(double));
  for (i = 0; i < nc; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->mu_r[i])));
  }
  DEBUG_PRINT("End --  int rollingFrictionContact_newFromFile(RollingFrictionContactProblem* problem, FILE* file)\n");

  return 0;
}

int rollingFrictionContact_newFromFilename(RollingFrictionContactProblem* problem, char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "r");

  if (!file)
  {
    return errno;
  }

  info = rollingFrictionContact_newFromFile(problem, file);

  fclose(file);
  return info;
}

void rollingFrictionContactProblem_free(RollingFrictionContactProblem* problem)
{
  assert(problem);
  if (problem->M)
  {
    NM_free(problem->M);
    free(problem->M);
    problem->M = NULL;
  }

  if (problem->mu)
  {
  free(problem->mu);
  problem->mu = NULL;
  }
  if (problem->mu_r)
  {
  free(problem->mu_r);
  problem->mu_r = NULL;
  }

  if (problem->q)
  {
    free(problem->q);
    problem->q = NULL;
  }

  free(problem);

}

RollingFrictionContactProblem* rollingFrictionContactProblem_new(void)
{
  RollingFrictionContactProblem* fcp = (RollingFrictionContactProblem*) malloc(sizeof(RollingFrictionContactProblem));
  fcp->dimension = 0;
  fcp->numberOfContacts = 0;
  fcp->M = NULL;
  fcp->q = NULL;
  fcp->mu = NULL;
  fcp->mu_r = NULL;

  return fcp;
}

RollingFrictionContactProblem* rollingFrictionContactProblem_new_with_data(int dim, int nc, NumericsMatrix* M, double* q, double* mu, double* mu_r)
{
  RollingFrictionContactProblem* fcp = (RollingFrictionContactProblem*) malloc(sizeof(RollingFrictionContactProblem));

  fcp->dimension = dim;
  fcp->numberOfContacts = nc;
  fcp->M = M;
  fcp->q = q;
  fcp->mu = mu;
  fcp->mu_r = mu_r;

  return fcp;
}

//#define SN_SBM_TO_DENSE

void rollingFrictionContactProblem_compute_statistics(RollingFrictionContactProblem* problem,
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


  int contact =0;

  int take_off_count=0;
  int closed_count=0;

  int sticking_count=0;
  int sliding_count=0;
  int ambiguous_take_off_count=0;
  int ambiguous_closed_count=0;

  for (contact =0; contact < nc; contact ++)
  {
    int pos=0;
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
