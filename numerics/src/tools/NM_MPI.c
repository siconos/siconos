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
#include "NM_MPI.h"
#include <assert.h>
#ifdef HAVE_MPI
#include "NumericsMatrix.h"
#include "numerics_verbose.h"
MPI_Comm NM_MPI_comm(NumericsMatrix* A)
{
  assert(A);

  if (NM_internalData(A)->mpi_comm == MPI_COMM_NULL)
  {
    if (verbose)
    {
      fprintf(stderr, "siconos/numerics: warning, MPI communicator has not been initialized,\n");
      fprintf(stderr, "siconos/numerics: MPI_COMM_WORLD will be used.\n");
    }
    NM_internalData(A)->mpi_comm = MPI_COMM_WORLD;
  }
  return NM_internalData(A)->mpi_comm = MPI_COMM_WORLD;
}

void NM_MPI_set_comm(NumericsMatrix* A, MPI_Comm comm)
{
  assert(A);
  NM_internalData(A)->mpi_comm = comm;
}


#endif /* WITH_MPI */

int NM_MPI_rank(NumericsMatrix* A)
{
  assert(A);
  int myid;
#ifdef HAVE_MPI
  CHECK_MPI(NM_MPI_comm(A), MPI_Comm_rank(NM_MPI_comm(A), &myid));
#else
  myid = 0;
#endif
  return myid;
}

void NM_MPI_copy(const NumericsMatrix* A, NumericsMatrix* B)
{
  assert(A);
  assert(B);
#ifdef HAVE_MPI
  if (A->internalData && A->internalData->mpi_comm)
  {
    NM_MPI_set_comm(B, A->internalData->mpi_comm);
  }
#endif
}
