#ifndef TEST_UTILS_H
#define TEST_UTILS_H
#include "NumericsFwd.h"    // for SolverOptions
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

typedef struct
{
  /** name of the data file used for the test*/
  const char * filename;
  /** indicates if the test is supposed to fail (1) or not (0)*/
  int will_fail;
  /** solver options*/
  SolverOptions * options;
} TestCase;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Reads the list of reference data files
      returns an array of char
  */
  const char ** data_collection(void);

  void print_matrix( const char* desc, int m, int n, double* a, int lda );

  void print_int_vector( const char* desc, int n, int* a );

  void print_vector_norm( const char* desc, int m, int n, double* a, int lda );

  void print_test_info(int test_id, TestCase* current, const char* msg);

  void free_test_collection(TestCase* collection, int nb_tests);

  void print_tests_collection_report(TestCase * collection, int, int * failed_tests, int, int * succeeded_test);

  TestCase*  build_test_collection_generic(int, const char **, int, int*);

  int run_test_collection(TestCase * collection, int number_of_tests, int (*test_function)(TestCase*));

  /* void print_problem_data_in_Matlab_file(GlobalFrictionContactProblem * problem, FILE * file)   */


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
