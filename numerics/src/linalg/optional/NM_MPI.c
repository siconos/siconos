// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifdef SICONOS_HAS_MPI
#include "NM_MPI.h"

#include <assert.h>

#include "NumericsMatrix.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"
MPI_Comm NM_MPI_comm(NumericsMatrix* A) {
  if (!A) return MPI_COMM_NULL;
  if (NM_internalData(A)->mpi_comm == MPI_COMM_NULL) {
    if (verbose) {
      fprintf(stderr,
              "siconos/numerics: warning, MPI communicator has not been initialized,\n");
      fprintf(stderr, "siconos/numerics: MPI_COMM_WORLD will be used.\n");
    }
    NM_internalData(A)->mpi_comm = MPI_COMM_WORLD;
  }
  return NM_internalData(A)->mpi_comm = MPI_COMM_WORLD;
}

void NM_MPI_set_comm(NumericsMatrix* A, MPI_Comm comm) {
  assert(A);
  NM_internalData(A)->mpi_comm = comm;
}

int NM_MPI_rank(NumericsMatrix* A) {
  CHECK_NULL(A);
  int myid;
  CHECK_MPI(NM_MPI_comm(A), MPI_Comm_rank(NM_MPI_comm(A), &myid));
  return myid;
}

void NM_MPI_copy(const NumericsMatrix* A, NumericsMatrix* B) {
  assert(A);
  assert(B);
  if (A->internalData && A->internalData->mpi_comm) {
    NM_MPI_set_comm(B, A->internalData->mpi_comm);
  }
}
#endif /* SICONOS_HAS_MPI */
