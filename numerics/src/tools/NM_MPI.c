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
NM_MPI_comm_t NM_MPI_comm(NumericsMatrix* A)
{
  assert(A);
  NumericsMatrixInternalData* p = NM_internalData(A);
  if (p->mpi_comm == MPI_COMM_NULL)
  {
    int myid;
    int argc = 0;
    /* C99 requires that argv[argc] == NULL. With openmpi 1.8, we get a
     * segfault if this is not true */
    char *argv0 = 0;
    char **argv = &argv0;
    CHECK_MPI(MPI_Init(&argc, &argv));
    CHECK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &myid));

    p->mpi_comm = MPI_COMM_WORLD;
  }
  return p->mpi_comm;

}

#endif /* WITH_MPI */

int NM_MPI_rank(NumericsMatrix* A)
{
  assert(A);
  NumericsMatrixInternalData* p = NM_internalData(A);
  int myid;
#ifdef HAVE_MPI
  CHECK_MPI(MPI_Comm_rank(p->mpi_comm, &myid));
#else
  myid = 0;
#endif
  return myid;
}
