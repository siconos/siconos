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
#include <stdio.h>                         // for printf, FILE
#include <stdlib.h>                        // for calloc, free
#include <assert.h>                        // for assert
#include "GenericMechanicalProblem.h"      // for GenericMechanicalProblem
#include "GenericMechanical_Solvers.h"     // for gmp_compute_error, gmp_driver
#include "NumericsFwd.h"                   // for SolverOptions, GenericMech...
#include "SolverOptions.h"                 // for SolverOptions
#include "genericMechanical_test_utils.h"  // for genericMechanical_test_fun...
#include "test_utils.h"

int gmp_test_function(TestCase* current)
{

  int k, info = -1 ;
  GenericMechanicalProblem* problem = genericMechanicalProblem_new();
  info = genericMechanical_newFromFilename(problem, current->filename);
  double *reaction = (double*)calloc(problem->size, sizeof(double));
  double *velocity = (double*)calloc(problem->size, sizeof(double));
  info = gmp_driver(problem, reaction , velocity, current->options);
  double err = 0;
  gmp_compute_error(problem, reaction , velocity, current->options->dparam[SICONOS_DPARAM_TOL], current->options, &err);
  printf("\n");
  for (k = 0 ; k < problem->size; k++)
  {
    printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
  }
  printf("\n");
  
  if (!info)
  {
    if (err > current->options->dparam[SICONOS_DPARAM_TOL])
    {
      printf("test unsuccessful: err>tol\n");
      printf(" ---> info=%i err=%e and tol=%e\n", info, err, current->options->dparam[SICONOS_DPARAM_TOL]);

      return 1;
    }
    else
      printf("test successful: info=%i err=%e and tol=%e\n", info, err, current->options->dparam[SICONOS_DPARAM_TOL]);

  }
  else
  {
    printf("test unsuccessful.\n");
  }
  printf("GMP TEST: Nb GS it=%i\n",current->options->iparam[3]);
  free(reaction);
  free(velocity);
  //solver_options_delete(current->options); --> break tests. Should be investigated.
  genericMechanicalProblem_free(problem, NUMERICS_GMP_FREE_MATRIX);
  free(problem);

  return info;

}

void build_gmp_test(const char * filename,
                int solver_id, int* d_ind, double* dparam, int * i_ind, int* iparam,
                int internal_solver_id, int * i_d_ind, double * internal_dparam, int * i_i_ind, int * internal_iparam,
                TestCase* testname)
{
  // reference file name
  testname->filename = filename;
  // By default, test is expected to succeed.
  testname->will_fail = 0;

  // Set solver options to default.
  testname->options = (SolverOptions *)malloc(sizeof(SolverOptions));

  assert (internal_solver_id > 0); // internal solver is required

  gmp_setDefaultSolverOptions(testname->options, internal_solver_id);
  // Fill iparam and dparam in.
  if(iparam)
    for(int i=0; i<i_ind[0]; ++i)
      testname->options->iparam[i_ind[i+1]] = iparam[i];

  // dparam
  if(dparam)
    for(int i=0; i<d_ind[0]; ++i)
      testname->options->dparam[d_ind[i+1]] = dparam[i];

  // Internal solver setup
  if(internal_solver_id>0)
    {
      // internal iparam
      if(internal_iparam)
        for(int i=0; i<i_i_ind[0]; ++i)
          testname->options->internalSolvers[0].iparam[i_i_ind[i+1]] = internal_iparam[i];
      // internal dparam
      if(internal_dparam)
        for(int i=0; i<i_d_ind[0]; ++i)
          testname->options->internalSolvers[0].dparam[i_d_ind[i+1]] = internal_dparam[i];
    }
}


