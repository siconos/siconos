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

#include <stdio.h>                        // for printf, NULL
#include <stdlib.h>                       // for free, malloc, calloc
#include "MCP_cst.h"                      // for SICONOS_MCP_NEWTON_FB_FBLSA
#include "MixedComplementarityProblem.h"  // for MixedComplementarityProblem
#include "NonSmoothDrivers.h"             // for mcp_driver
#include "NumericsFwd.h"                  // for MixedComplementarityProblem
#include "NumericsMatrix.h"               // for NM_create, NM_DENSE, Numeri...
#include "NumericsVerbose.h"              // for numerics_set_verbose
#include "SolverOptions.h"                // for solver_options_id_to_name
static double * M;
static double * q;

void testF(void* env, int size, double *z, double * F);
void testF(void* env, int size, double *z, double * F)
{
  /* printf("call to MCP function F(z) ...\n");   */
  /* for (int i =0 ; i <size; i++) */
  /* { */
  /*   printf("z[%i]= %lf\t",i,z[i]); */
  /* } */
  /* printf("\n"); */
  for(int i = 0 ; i < size; i++)
  {
    F[i] = q[i];
    for(int j = 0 ; j < size ;  j++)
    {
      F[i] += M[i + j * size] * z[j] ;
    }
  }
  /* for (int i =0 ; i <size; i++) */
  /* { */
  /*   printf("F[%i]= %lf\t",i,F[i]); */
  /* } */
  /* printf("\n"); */
  /* printf("End call to MCP function F(z) ...\n");   */
}

void testNablaF(void * env, int size, double *z, NumericsMatrix* nablaF);
void testNablaF(void * env, int size, double *z, NumericsMatrix* nablaF)
{
  /* printf("call to MCP function nablaF(z) ...\n"); */

  for(int i = 0 ; i < size; i++)
  {
    for(int j = 0 ; j < size ;  j++)
    {
      nablaF->matrix0[i + j * size] = M[i + j * size];
    }
  }
}

static MixedComplementarityProblem * create_mcp_1(void)
{
  /* Create a MixedComplementarityProblem */
  MixedComplementarityProblem* problem = (MixedComplementarityProblem *)malloc(sizeof(MixedComplementarityProblem));

  int n=10;

  problem->n1 = n-5;
  problem->n2 = 5;
  problem->compute_Fmcp = &testF ;
  problem->compute_nabla_Fmcp = &testNablaF ;
  problem->nabla_Fmcp =  NM_create(NM_DENSE, n, n);
  problem->env = NULL;

  M = (double *) calloc(n*n,sizeof(double));
  q = (double *) calloc(n,sizeof(double));

  for(int i =0; i< n; i++)
  {
    q[i] = - i + 7.;
    M[i+i*n] = 2.0 ;
    if(i < n-1)
      M[i+(i+1)*n] =1.0;
    if(i >0)
      M[i+(i-1)*n] =1.0;
  }
  //NM_dense_display(M,n,n,n);
  return problem;
}

static void free_mcp_1(MixedComplementarityProblem * problem)
{
  free(M);
  free(q);
  mixedComplementarityProblem_free(problem);
}

static int test_mcp_newton(int solverId)
{
  printf("test_mcp_newton() starts for solver %s.\n", solver_options_id_to_name(solverId));

  int info = 1 ;

  MixedComplementarityProblem* problem = create_mcp_1();

  /* Set solver options, FB solver */
  SolverOptions * options = solver_options_create(solverId);

  numerics_set_verbose(1);

  int size = problem->n1 + problem->n1 ;
  double * z = (double *)malloc(size * sizeof(double));
  double * w = (double *)malloc(size * sizeof(double));

  for(int i = 0 ; i < size; i++)
  {
    z[i] = 0.0;
    w[i] = 0.0;
  }

  info = mcp_driver(problem, z, w,  options);

  /// TODO : write a real test ... ////
  for(int i = 0 ; i < size; i++)
  {
    printf("z[%i]= %lf\t", i, z[i]);
  }
  printf("\n");
  for(int i = 0 ; i < size; i++)
  {
    printf("w[%i]= %lf\t", i, w[i]);
  }
  printf("\n");
  solver_options_delete(options);
  options = NULL;
  free(z);
  free(w);
  free_mcp_1(problem);
  printf("test_mcp_newton() starts for solver %s.\n\n", solver_options_id_to_name(solverId));
  return info;
}



int main(void)
{
  printf("Start tests for MCP solvers.\n\n");
  int info = 1 ;


  info = test_mcp_newton(SICONOS_MCP_NEWTON_MIN_FBLSA);
  info += test_mcp_newton(SICONOS_MCP_NEWTON_FB_FBLSA);
  printf("End of MCP solvers test. \n");
  return info;
}
