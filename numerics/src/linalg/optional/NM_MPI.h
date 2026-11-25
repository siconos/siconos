/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifdef SICONOS_HAS_MPI
#include <mpi.h>

#include "NumericsFwd.h"  // for NumericsMatrix

#if defined(__cplusplus)
extern "C" {
#endif

MPI_Comm NM_MPI_comm(NumericsMatrix* A);
void NM_MPI_set_comm(NumericsMatrix* A, MPI_Comm comm);
#include <stdio.h>
#define CHECK_MPI(COMM, EXPR)                                               \
  do {                                                                      \
    int error_code = EXPR;                                                  \
    MPI_Comm_set_errhandler(COMM, MPI_ERRORS_RETURN);                       \
    if (error_code != MPI_SUCCESS) {                                        \
      char error_string[1024];                                              \
      int length_of_error_string, error_class;                              \
      MPI_Error_class(error_code, &error_class);                            \
      MPI_Error_string(error_class, error_string, &length_of_error_string); \
      fprintf(stderr, "%3d: %s\n", 0, error_string);                        \
      MPI_Error_string(error_code, error_string, &length_of_error_string);  \
      fprintf(stderr, "%3d: %s\n", 0, error_string);                        \
      MPI_Abort(COMM, error_code);                                          \
    };                                                                      \
  } while (0)

int NM_MPI_rank(NumericsMatrix* A);

void NM_MPI_copy(const NumericsMatrix* A, NumericsMatrix* B);

#if defined(__cplusplus)
}
#endif
#endif
