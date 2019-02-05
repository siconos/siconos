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

/* thread_local madness for the MPI communicator */
#include "tlsdef.h"


/* MPI_INIT should be called only once. Therefore we have to remember if this
 * was already done or not. Using TLS seems to be the best option here --xhub  */
tlsvar MPI_Comm NM_mpi_com = MPI_COMM_NULL;


#if defined(_WIN32) && defined(_MSC_VER)
BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpReserved)
{
  if (NM_mpi_com != MPI_COMM_NULL)
  {
    MPI_Finalize();
    NM_mpi_com = MPI_COMM_NULL;
  }
  return true;
}
#elif defined(__GNUC__) & !defined(__APPLE__)
/* MB: called at program exit */
static DESTRUCTOR_ATTR void cleanup_MPI(void)
{
  if (NM_mpi_com != MPI_COMM_NULL)
  {
    MPI_Finalize();
    NM_mpi_com = MPI_COMM_NULL;
  }
}
#endif /* MB: so it is ok for gnu or win/visual,
       /* what about the rest of the world ? */



#ifndef NULL
#define NULL 0
#endif

MPI_Comm NM_MPI_com(MPI_Comm m)
{
  assert(m);
  if (NM_mpi_com == MPI_COMM_NULL)
  {
    if (m != MPI_COMM_NULL)
    {
      NM_mpi_com = m;
    }
    else
    {
      int myid;
      int argc = 0;
      /* C99 requires that argv[argc] == NULL. With openmpi 1.8, we get a
       * segfault if this is not true */
      char *argv0 = NULL;
      char **argv = &argv0;
      CHECK_MPI(MPI_Init(&argc, &argv));
      CHECK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &myid));

      NM_mpi_com = MPI_COMM_WORLD;
    }
  }

  return NM_mpi_com;

}

#endif /* WITH_MPI */

int NM_MPI_rank(NM_MPI_comm com)
{
  int myid;
#ifdef HAVE_MPI
  CHECK_MPI(MPI_Comm_rank(com, &myid));
#else
  myid = 0;
#endif
  return myid;
}
