/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "GlobalRollingFrictionContactProblem.h"
#include "GlobalFrictionContactProblem.h"
#include "RollingFrictionContactProblem.h"
#include "FrictionContactProblem.h"
#include <assert.h>            // for assert
#include <math.h>              // for fabs
#include <stdio.h>             // for printf, fprintf, fscanf, NULL, fclose
#include <stdlib.h>            // for free, malloc, exit, EXIT_FAILURE
#include <sys/errno.h>         // for errno
#include "SiconosBlas.h"         // for cblas_dscal, cblas_dcopy
#include "NumericsMatrix.h"    // for NM_vector_display, NM_display, NM_clear
#include "numerics_verbose.h"  // for CHECK_IO, numerics_printf_verbose


//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "siconos_debug.h"             // for DEBUG_PRINT, DEBUG_PRINTF


GlobalRollingFrictionContactProblem* globalRollingFrictionContactProblem_new(void)
{
  GlobalRollingFrictionContactProblem* fcp = (GlobalRollingFrictionContactProblem*) malloc(sizeof(GlobalRollingFrictionContactProblem));
  fcp->dimension = 0;
  fcp->numberOfContacts = 0;
  fcp->M = NULL;
  fcp->H = NULL;
  fcp->q = NULL;
  fcp->b = NULL;
  fcp->mu = NULL;
  fcp->mu_r = NULL;

  return fcp;
}



void globalRollingFrictionContact_display(GlobalRollingFrictionContactProblem* problem)
{

  assert(problem);
  int n = problem->dimension * problem->numberOfContacts;
  printf("GlobalRollingFrictionContact Display :\n-------------\n");
  printf("dimension :%d \n", problem->dimension);
  printf("numberOfContacts:%d \n", problem->numberOfContacts);

  if(problem->M)
  {
    printf("M matrix:\n");
    NM_display(problem->M);
  }
  else
    printf("No M matrix:\n");

  if(problem->H)
  {
    printf("H matrix:\n");
    NM_display(problem->H);
  }
  else
    printf("No H matrix:\n");

  if(problem->q)
  {
    printf("q vector:\n");
    NM_vector_display(problem->q,n);
  }
  else
    printf("No q vector:\n");
  if(problem->b)
  {
    printf("b vector:\n");
    NM_vector_display(problem->b,n);
  }
  else
    printf("No q vector:\n");

  if(problem->mu)
  {
    printf("mu vector:\n");
    NM_vector_display(problem->mu,problem->numberOfContacts);
  }
  else
    printf("No mu vector:\n");

  if(problem->mu_r)
  {
    printf("mu_R vector:\n");
    NM_vector_display(problem->mu_r,problem->numberOfContacts);
  }
  else
    printf("No mu_R vector:\n");

}





int globalRollingFrictionContact_printInFile(GlobalRollingFrictionContactProblem*  problem, FILE* file)
{
  if(! problem)
  {
    fprintf(stderr, "Numerics, GlobalRollingFrictionContactProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;

  int d  = problem->dimension;
  fprintf(file, "%d\n", d);
  int nc = problem->numberOfContacts;
  fprintf(file, "%d\n", nc);
  NM_write_in_file(problem->M, file);
  NM_write_in_file(problem->H, file);
  for(i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  for(i = 0; i < problem->H->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->b[i]);
  }
  fprintf(file, "\n");
  for(i = 0; i < nc; i++)
  {
    fprintf(file, "%32.24e ", problem->mu[i]);
  }
  fprintf(file, "\n");
  for(i = 0; i < nc; i++)
  {
    fprintf(file, "%32.24e ", problem->mu_r[i]);
  }
  fprintf(file, "\n");
  return 0;
}

int globalRollingFrictionContact_printInFilename(GlobalRollingFrictionContactProblem* problem, char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "w");

  if(!file)
  {
    return errno;
  }

  info = globalRollingFrictionContact_printInFile(problem, file);

  fclose(file);
  return info;
}

GlobalRollingFrictionContactProblem* globalRollingFrictionContact_newFromFile(FILE* file)
{
  GlobalRollingFrictionContactProblem* problem = globalRollingFrictionContactProblem_new();
  assert(file);
  DEBUG_PRINT("Start -- int globalRollingFrictionContact_newFromFile(GlobalRollingFrictionContactProblem* problem, FILE* file)\n");
  int nc = 0, d = 0;
  int i;
  CHECK_IO(fscanf(file, "%d\n", &d));
  problem->dimension = d;
  DEBUG_PRINTF("problem->dimension = %i \n",problem->dimension);
  CHECK_IO(fscanf(file, "%d\n", &nc));
  problem->numberOfContacts = nc;
  problem->M =  NM_new_from_file(file);
  problem->H =  NM_new_from_file(file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for(i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])));
  }
  problem->b = (double *) malloc(problem->H->size1 * sizeof(double));
  for(int i = 0; i < problem->H->size1; ++i)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->b[i])));
  }
  problem->mu = (double *) malloc(nc * sizeof(double));
  for(i = 0; i < nc; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->mu[i])));
  }
  problem->mu_r = (double *) malloc(nc * sizeof(double));
  for(i = 0; i < nc; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->mu_r[i])));
  }
  DEBUG_PRINT("End --  int globalRollingFrictionContact_newFromFile(GlobalRollingFrictionContactProblem* problem, FILE* file)\n");

  return problem;
}

GlobalRollingFrictionContactProblem* globalRollingFrictionContact_new_from_filename(const char* filename)
{
  GlobalRollingFrictionContactProblem* problem = NULL;
  FILE * file = fopen(filename, "r");
  if(!file)
    numerics_error("GlobalRollingFrictionContactProblem", "Can not open file ", filename);

  problem = globalRollingFrictionContact_newFromFile(file);
  fclose(file);
  return problem;
}

void globalRollingFrictionContactProblem_free(GlobalRollingFrictionContactProblem* problem)
{
  assert(problem);
  if(problem->M)
  {
    NM_clear(problem->M);
    free(problem->M);
    problem->M = NULL;
  }
  if(problem->H)
  {
    NM_clear(problem->H);
    free(problem->H);
    problem->H = NULL;
  }

  if(problem->mu)
  {
    free(problem->mu);
    problem->mu = NULL;
  }
  if(problem->mu_r)
  {
    free(problem->mu_r);
    problem->mu_r = NULL;
  }

  if(problem->q)
  {
    free(problem->q);
    problem->q = NULL;
  }
  if(problem->b)
  {
    free(problem->b);
    problem->b = NULL;
  }

  free(problem);

}


GlobalRollingFrictionContactProblem* globalRollingFrictionContactProblem_new_with_data(int dim, int nc, NumericsMatrix* M, double* q, double* mu, double* mu_r)
{
  GlobalRollingFrictionContactProblem* fcp = (GlobalRollingFrictionContactProblem*) malloc(sizeof(GlobalRollingFrictionContactProblem));

  fcp->dimension = dim;
  fcp->numberOfContacts = nc;
  fcp->M = M;
  fcp->q = q;
  fcp->mu = mu;
  fcp->mu_r = mu_r;

  return fcp;
}

int globalRollingFrictionContact_computeGlobalVelocity(
  GlobalRollingFrictionContactProblem* problem,
  double * reaction,
  double * globalVelocity)
{
  int info = -1;

  int n = problem->M->size0;
  int m = problem->H->size1;

  /* globalVelocity <- problem->q */
  cblas_dcopy(n,  problem->q, 1, globalVelocity, 1);

  // We compute only if the problem has contacts
  if(m>0)
  {
    /* globalVelocity <-  H*reaction + globalVelocity*/
    NM_gemv(1.0, problem->H, reaction, 1.0, globalVelocity);
    DEBUG_EXPR(NM_vector_display(reaction, m));
  }

  /* Compute globalVelocity <- M^(-1) globalVelocity*/
  // info = NM_gesv_expert(problem->M, globalVelocity, NM_PRESERVE)
  info = NM_LU_solve(NM_preserve(problem->M), globalVelocity, 1);
  DEBUG_EXPR(NM_vector_display(globalVelocity, n));

  return info;
}



RollingFrictionContactProblem * globalRollingFrictionContact_reformulation_RollingFrictionContact(GlobalRollingFrictionContactProblem* problem)
{
  /* we use the code in gfc3d_reformulation_local_problem */
  GlobalFrictionContactProblem* gfc3d = globalFrictionContactProblem_new();
  gfc3d->numberOfContacts = problem->numberOfContacts;
  gfc3d->dimension =  problem->dimension;
  gfc3d->M = problem->M;
  gfc3d->H = problem->H;
  gfc3d->q = problem->q;
  gfc3d->b = problem->b;
  gfc3d->mu = problem->mu;

  FrictionContactProblem* fc3d = globalFrictionContact_reformulation_FrictionContact(gfc3d);

  RollingFrictionContactProblem * localproblem =  rollingFrictionContactProblem_new();

  localproblem->numberOfContacts = problem->numberOfContacts;
  localproblem->dimension =  problem->dimension;

  localproblem->mu =  fc3d->mu;
  localproblem->M = fc3d->M;
  localproblem->q = fc3d->q;
  localproblem->mu = fc3d->mu;

  localproblem->mu_r = (double*)malloc(problem->numberOfContacts * sizeof(double));
  cblas_dcopy(problem->numberOfContacts, problem->mu_r, 1, localproblem->mu_r, 1);

  free(gfc3d);
  free(fc3d);
  return localproblem;

}
