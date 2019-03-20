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
#include "SiconosConfig.h"

#ifdef HAVE_MPI
#include "NumericsFwd.h"
#include <mpi.h>

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

MPI_Comm NM_MPI_comm(NumericsMatrix* A);
void NM_MPI_set_comm(NumericsMatrix* A, MPI_Comm comm);
#include <stdio.h>
#define CHECK_MPI(COMM, EXPR)                                           \
  do                                                                    \
  {                                                                     \
    int error_code = EXPR;                                              \
    MPI_Errhandler_set(COMM, MPI_ERRORS_RETURN);                        \
    if (error_code != MPI_SUCCESS) {                                    \
      char error_string[1024];                                          \
      int length_of_error_string, error_class;                          \
      MPI_Error_class(error_code, &error_class);                        \
      MPI_Error_string(error_class, error_string, &length_of_error_string); \
      fprintf(stderr, "%3d: %s\n", 0, error_string);                    \
      MPI_Error_string(error_code, error_string, &length_of_error_string); \
      fprintf(stderr, "%3d: %s\n", 0, error_string);                    \
      MPI_Abort(COMM, error_code);                                      \
    };                                                                  \
  } while(0)
#endif

  int NM_MPI_rank(NumericsMatrix* A);

  void NM_MPI_copy(const NumericsMatrix* A, NumericsMatrix*  B);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

