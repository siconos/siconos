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
#include <stdio.h>                         // for printf
#include <stdlib.h>                        // for free, calloc
#include "GenericMechanicalProblem.h"      // for GenericMechanicalProblem
#include "GenericMechanical_Solvers.h"     // for gmp_compute_error, gmp_driver
#include "NumericsFwd.h"                   // for SolverOptions, GenericMech...
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_IPA...
#include "genericMechanical_test_utils.h"  // for gmp_test_function
#include "test_utils.h"                    // for TestCase
int gmp_test_function(TestCase* current)
{

  int k, info = -1 ;
  GenericMechanicalProblem* problem = genericMechanical_new_from_filename(current->filename);
  double *reaction = (double*)calloc(problem->size, sizeof(double));
  double *velocity = (double*)calloc(problem->size, sizeof(double));

  /* printf(" iparam: \n"); */
  /* for(size_t i=0;i<20;++i) */
  /*   printf(" %d \t ", current->options->internalSolvers[1]->iparam[i]); */
  /* printf(" dparam: \n"); */
  /* for(size_t i=0;i<20;++i) */
  /*   printf(" %g \t ", current->options->internalSolvers[1]->dparam[i]); */

  info = gmp_driver(problem, reaction, velocity, current->options);

  double err = 0;
  gmp_compute_error(problem, reaction, velocity, current->options->dparam[SICONOS_DPARAM_TOL], current->options, &err);
  printf("\n");
  for(k = 0 ; k < problem->size; k++)
  {
    printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k, reaction[k]);
  }
  printf("\n");


  if(!info)
  {
    if(err > current->options->dparam[SICONOS_DPARAM_TOL])
    {
      printf("test unsuccessful, residual = %g, info = %d, nb iter = %d\n", err, info, current->options->iparam[SICONOS_IPARAM_ITER_DONE]);
      info = 1;
    }
    else
      printf("test successful, residual = %g\t, number of iterations = %i \n", err, current->options->iparam[SICONOS_IPARAM_ITER_DONE]);

  }
  else
  {
    printf("test unsuccessful, residual = %g, info = %d, nb iter = %d\n", err, info, current->options->iparam[SICONOS_IPARAM_ITER_DONE]);
  }
  printf("GMP TEST: Nb GS it=%i\n",current->options->iparam[SICONOS_IPARAM_ITER_DONE]);
  free(reaction);
  free(velocity);
  genericMechanicalProblem_free(problem, NUMERICS_GMP_FREE_MATRIX);
  free(problem);
  return info;

}

