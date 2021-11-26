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
#include "test_utils.h"
#include <stdio.h>          // for printf, NULL
#include <stdlib.h>         // for free
#include "SolverOptions.h"  // for solver_options_id_to_name, SolverOptions
#include <time.h>
#include "numerics_verbose.h"
/* Auxiliary routine: printing a matrix */
void print_matrix(const char* desc, int m, int n, double* a, int lda)
{
  int i, j;
  printf("\n %s\n", desc);
  for(i = 0; i < m; i++)
  {
    for(j = 0; j < n; j++) printf(" %6.2f", a[i+j*lda]);
    printf("\n");
  }
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector(const char* desc, int n, int* a)
{
  int j;
  printf("\n %s\n", desc);
  for(j = 0; j < n; j++) printf(" %6i", a[j]);
  printf("\n");
}

void print_vector_norm(const char* desc, int m, int n, double* a, int lda)
{
  int i, j;
  printf("\n %s\n", desc);
  for(j = 0; j < n; j++)
  {
    double norm = 0.0;
    for(i = 0; i < m; i++) norm += a[i+j*lda] * a[i+j*lda];
    printf(" %6.2f", norm);
  }
  printf("\n");
}

void free_test_collection(TestCase* collection, int nb_tests)
{
  for(int i=0; i < nb_tests; ++i)
  {
    solver_options_delete(collection[i].options);
    free(collection[i].options);
    collection[i].options = NULL;
  }

  free(collection);
}

void print_test_info(int test_id, TestCase* current, const char* msg)
{
  printf("test %i (%s, on %s) %s.\n",
         test_id,
         solver_options_id_to_name(current->options->solverId),
         current->filename, msg);
  printf("\n");
}



void print_tests_collection_report(TestCase * collection, int n_failed, int * failed_tests,
                                   int n_succeeded, int * succeeded_tests)
{


  printf("\n Succeeded tests ids: [ ");
  for(int t =0; t < n_succeeded; t++)
    printf("%i, ", succeeded_tests[t]);
  printf("] :\n");

  printf("\n Failed tests ids: [ ");
  for(int t =0; t < n_failed; t++)
    printf("%i, ", failed_tests[t]);
  printf("] :\n");

  for(int t =0; t < n_failed; t++)
  {
    // was it expected to fail?
    if(collection[failed_tests[t]].will_fail == 1)
      print_test_info(failed_tests[t], &collection[failed_tests[t]], " was expected to fail");
    else if(collection[failed_tests[t]].will_fail == 0)  // Or not ...
      print_test_info(failed_tests[t], &collection[failed_tests[t]], " was expected to succeed");
    else if(collection[failed_tests[t]].will_fail == 2)  // Or is unstable.
      print_test_info(failed_tests[t], &collection[failed_tests[t]], " is unstable and has failed");
  }

  for(int t =0; t < n_succeeded; t++)
  {
    if(collection[succeeded_tests[t]].will_fail == 1)
      print_test_info(succeeded_tests[t], &collection[succeeded_tests[t]], " was expected to fail but has succeeded");
    else if(collection[failed_tests[t]].will_fail == 2)  // Or is unstable.
      print_test_info(failed_tests[t], &collection[failed_tests[t]], " is unstable and has suceeded");
  }
}

TestCase * build_test_collection_generic(int n_data, const char ** data_collection, int n_solvers, int * solvers_ids)
{
  int number_of_tests = n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc(number_of_tests * sizeof(TestCase));

  for(int s=0; s <n_solvers; ++s)
  {
    int current = s * n_data;
    for(int d=0; d <n_data; ++d)
    {
      collection[current].filename = data_collection[d];
      collection[current].will_fail = 0;
      collection[current].options = solver_options_create(solvers_ids[s]);
      current ++;
    }
  }
  return collection;
}

int run_test_collection(TestCase * collection, int number_of_tests, int (*test_function)(TestCase*))
{
  // arrays to save indices of failed and succeeded tests.
  int * failed_tests = (int*)calloc(number_of_tests, sizeof(int));
  int * succeeded_tests = (int*)calloc(number_of_tests, sizeof(int));
  int out = 0;
  int n_failed = 0;
  int n_succeeded = 0;
  long clk_tck = CLOCKS_PER_SEC;  

  // Loop through tests collection
  for(int test_num=0; test_num<number_of_tests; ++test_num)
  {
    // print solver details
    printf("\n################# start of test # %i #######################\n", test_num);
    printf("Solver : %s (id: %d) \n", solver_options_id_to_name(collection[test_num].options->solverId), collection[test_num].options->solverId);
    /* verbose=1; */
    solver_options_print(collection[test_num].options);
    /* verbose=0; */
    for(size_t i=0; i<collection[test_num].options->numberOfInternalSolvers; ++i)
    {
      int sid = collection[test_num].options->internalSolvers[i]->solverId;
      const char * internal_name = solver_options_id_to_name(sid);
      printf("Internal solver : %s (id: %d) \n", internal_name, sid);
    }
    printf("Data file : %s \n", collection[test_num].filename);
    clock_t t1 = clock();
    // Execute a single test
    int info = test_function(&collection[test_num]);
    clock_t t2 = clock();
    (void)printf("time (s) : %lf \n", (double)(t2-t1)/(double)clk_tck);
    
    // Update failed/succeeded lists
    if(info)  // info != 0 --> test is unsuccesful
    {
      failed_tests[n_failed++]  = test_num;
      // was it expected to succeed?
      if(collection[test_num].will_fail == 0)
        out = 1;
    }
    else
      succeeded_tests[n_succeeded++]  = test_num;
    printf("\n################# end of  test # %i #######################\n", test_num);
  }

  // report (on screen)
  print_tests_collection_report(collection, n_failed, failed_tests, n_succeeded, succeeded_tests);

  // tests status.
  free(failed_tests);
  free(succeeded_tests);
  return out;
}

