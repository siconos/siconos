/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include "GlobalFrictionContactProblem.h"
#include <assert.h>            // for assert
#include <stdlib.h>            // for free, malloc, exit, EXIT_FAILURE
#include <sys/errno.h>         // for errno
#include <string.h>            // for memcpy
#include "SiconosBlas.h"         // for cblas_dscal, cblas_dcopy
#include "NumericsMatrix.h"    // for NumericsMatrix, NM_display, NM_clear
#include "numerics_verbose.h"  // for CHECK_IO
#include "io_tools.h"
#include "debug.h"
#if defined(WITH_FCLIB)
#include "fclib_interface.h"
#endif

GlobalFrictionContactProblem* globalFrictionContactProblem_new(void)
{
  GlobalFrictionContactProblem* problem = malloc(sizeof(GlobalFrictionContactProblem));
  problem->M = NULL;
  problem->H = NULL;
  problem->q = NULL;
  problem->b = NULL;
  problem->mu = NULL;
  problem->env = NULL;
  problem->numberOfContacts = 0;
  problem->dimension = 0;
  return problem;
}

int globalFrictionContact_printInFile(GlobalFrictionContactProblem*  problem, FILE* file)
{
  if(! problem)
  {
    fprintf(stderr, "Numerics, GlobalFrictionContactProblem printInFile failed, NULL input.\n");
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
  return 0;
}

GlobalFrictionContactProblem* globalFrictionContact_newFromFile(FILE* file)
{
  GlobalFrictionContactProblem* problem = globalFrictionContactProblem_new();
  int nc = 0, d = 0;
  int info = 0;
  CHECK_IO(fscanf(file, "%d\n", &d), &info);
  problem->dimension = d;
  CHECK_IO(fscanf(file, "%d\n", &nc), &info);
  problem->numberOfContacts = nc;
  problem->M = NM_new_from_file(file);

  problem->H =  NM_new_from_file(file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for(int i = 0; i < problem->M->size1; ++i)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])), &info);
  }
  problem->b = (double *) malloc(problem->H->size1 * sizeof(double));
  for(int i = 0; i < problem->H->size1; ++i)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->b[i])), &info);
  }

  problem->mu = (double *) malloc(nc * sizeof(double));
  for(int i = 0; i < nc; ++i)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->mu[i])), &info);
  }

  if(info)
  {
    problem = NULL;
  }
  return problem;
}

GlobalFrictionContactProblem* globalFrictionContact_new_from_filename(const char* filename)
{
  GlobalFrictionContactProblem* problem = NULL;
  int is_hdf5 = check_hdf5_file(filename);
  if(is_hdf5)
  {
#if defined(WITH_FCLIB)
    problem = globalFrictionContact_fclib_read(filename);
#else
    numerics_error("GlobalFrictionContactProblem",
                   "Try to read an hdf5 file, while fclib interface is not active. Recompile Siconos with fclib.",
                   filename);
#endif
  }
  else
  {
    FILE * file = fopen(filename, "r");
    if(!file)
      numerics_error("GlobalFrictionContactProblem", "Can not open file ", filename);

    problem = globalFrictionContact_newFromFile(file);
    fclose(file);
  }
  return problem;

}

int globalFrictionContact_printInFileName(GlobalFrictionContactProblem* problem, const char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "w");

  if(!file)
  {
    return errno;
  }

  info = globalFrictionContact_printInFile(problem, file);

  fclose(file);
  return info;
}

void globalFrictionContact_free(GlobalFrictionContactProblem* problem)
{

  if(problem->M)
  {
    NM_clear(problem->M);
    free(problem->M);
  }
  problem->M = NULL;

  if(problem->H)
  {
    NM_clear(problem->H);
    free(problem->H);
  }
  problem->H = NULL;

  if(problem->mu)
  {
    free(problem->mu);
  }
  problem->mu = NULL;

  if(problem->q)
  {
    free(problem->q);
  }
  problem->q = NULL;

  if(problem->b)
  {
    free(problem->b);
  }
  problem->b = NULL;

  if(problem->env) assert(0 && "globalFrictionContact_free :: problem->env != NULL, don't know what to do");

  free(problem);
}

void globalFrictionContact_display(GlobalFrictionContactProblem* problem)
{

  assert(problem);
  int i, n = problem->dimension * problem->numberOfContacts;
  printf("GlobalFrictionContact Display :\n-------------\n");
  printf("dimension :%d \n", problem->dimension);
  printf("numberOfContacts:%d \n", problem->numberOfContacts);
  int m = problem->M->size0;
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
    for(i = 0; i < m; i++) printf("q[ %i ] = %12.8e\n", i, problem->q[i]);
  }
  else
    printf("No q vector:\n");

  if(problem->b)
  {
    printf("b vector:\n");
    for(i = 0; i < n; i++) printf("b[ %i ] = %12.8e\n", i, problem->b[i]);
  }
  else
    printf("No q vector:\n");

  if(problem->mu)
  {
    printf("mu vector:\n");
    for(i = 0; i < problem->numberOfContacts; i++) printf("mu[ %i ] = %12.8e\n", i, problem->mu[i]);
  }
  else
    printf("No mu vector:\n");

}

GlobalFrictionContactProblem* globalFrictionContact_copy(GlobalFrictionContactProblem* problem)
{
  assert(problem);

  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;

  GlobalFrictionContactProblem* new = (GlobalFrictionContactProblem*) malloc(sizeof(GlobalFrictionContactProblem));
  new->dimension = problem->dimension;
  new->numberOfContacts = problem->numberOfContacts;
  new->M = NM_new();
  NM_copy(problem->M, new->M);
  new->H = NM_new();
  NM_copy(problem->H, new->H);
  new->q = (double*)malloc(n*sizeof(double));
  memcpy(new->q, problem->q, n*sizeof(double));
  new->b = (double*)malloc(m*sizeof(double));
  memcpy(new->b, problem->b, m*sizeof(double));
  new->mu = (double*)malloc(nc*sizeof(double));
  memcpy(new->mu, problem->mu, nc*sizeof(double));
  new->env = NULL;
  return new;
}

void globalFrictionContact_rescaling(
  GlobalFrictionContactProblem* problem,
  double alpha,
  double beta,
  double gamma)
{
  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;

  /* scaling of M */
  NM_scal(alpha*gamma*gamma, problem->M);
  /* scaling of H */
  NM_scal(beta*gamma, problem->H);
  /* scaling of q */
  cblas_dscal(n,alpha*gamma,problem->q,1);
  /* scaling of b */
  cblas_dscal(m,beta,problem->b,1);

}
void globalFrictionContact_balancing_M_H(
  GlobalFrictionContactProblem* problem,
  BalancingMatrices * B_for_M,
  BalancingMatrices * B_for_H)
{


  assert(B_for_M);
  assert(B_for_H);

  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;
  /* scaling of M */
  /* NM_scal(alpha*gamma*gamma, problem->M); */
  NumericsMatrix* M_tmp = NM_create(NM_SPARSE, n, n);
  NM_triplet_alloc(M_tmp, n);
  NM_gemm(1.0,problem->M,B_for_M->D2,0.0, M_tmp);
  NM_gemm(1.0,B_for_M->D2,M_tmp, 0.0, problem->M);
  NM_gemm(1.0,B_for_M->D1,problem->M, 0.0, M_tmp);
  NM_copy(M_tmp, problem->M);

  /* scaling of H */
  /* NM_scal(beta*gamma, problem->H);*/


  /* Warning the matrix H must be scaled such
   * that the cone constraint is respected */
  NumericsMatrix* H_tmp = NM_create(NM_SPARSE, n, m);
  NM_triplet_alloc(H_tmp, n);
  NM_gemm(1.0, B_for_M->D2, problem->H, 0.0, H_tmp);
  NM_gemm(1.0, B_for_H->D1, H_tmp, 0.0, problem->H);

  /* scaling of q */
  /* cblas_dscal(n,alpha*gamma,problem->q,1); */
  double * q_tmp = (double *) malloc(n*sizeof(double));
  NM_gemv(1.0, B_for_M->D2, problem->q, 0.0, q_tmp);
  NM_gemv(1.0, B_for_M->D1, q_tmp, 0.0, problem->q);

  /* scaling of b */
  /* cblas_dscal(m,beta,problem->b,1); */
  double * b_tmp = (double *) malloc(m*sizeof(double));
  NM_gemv(1.0, B_for_H->D1, problem->b, 0.0, b_tmp);
  cblas_dcopy(m, b_tmp, 1, problem->b, 1);


  NM_clear(M_tmp);
  free(M_tmp);
  NM_clear(H_tmp);
  free(H_tmp);

  free(q_tmp);
  free(b_tmp);
}


void globalFrictionContact_balancing_M(
  GlobalFrictionContactProblem* problem,
  BalancingMatrices * B_for_M)
{
  assert(B_for_M);

  NM_compute_balancing_matrices(problem->M, 1e-03, 100, B_for_M);

  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;
  /* scaling of M */
  /* NM_scal(alpha*gamma*gamma, problem->M); */
  NumericsMatrix* M_tmp = NM_create(NM_SPARSE, n, n);
  NM_triplet_alloc(M_tmp, n);
  NM_gemm(1.0,problem->M,B_for_M->D2,0.0, M_tmp);
  NM_gemm(1.0,B_for_M->D1,M_tmp, 0.0, problem->M);

  /* scaling of H */
  /* NM_scal(beta*gamma, problem->H);*/
  /* Warning the matrix H must be scaled such
   * that the cone constraint is respected */
  NumericsMatrix* H_tmp = NM_create(NM_SPARSE, n, m);
  NM_triplet_alloc(H_tmp, n);
  NM_gemm(1.0, B_for_M->D2, problem->H, 0.0, H_tmp);
  NM_copy(H_tmp, problem->H);

  /* scaling of q */
  /* cblas_dscal(n,alpha*gamma,problem->q,1); */
  double * q_tmp = (double *) malloc(n*sizeof(double));
  NM_gemv(1.0, B_for_M->D2, problem->q, 0.0, q_tmp);
  cblas_dcopy(n, q_tmp, 1, problem->q, 1);

  /* scaling of b */
  /* cblas_dscal(m,beta,problem->b,1); */

  NM_clear(M_tmp);
  free(M_tmp);
  NM_clear(H_tmp);
  free(H_tmp);

  free(q_tmp);
  //free(b_tmp);
}

int globalFrictionContact_computeGlobalVelocity(
  GlobalFrictionContactProblem* problem,
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
