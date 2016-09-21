/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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


#include "NumericsMatrix_private.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"

#include "debug.h"
#include "numerics_verbose.h"

#ifdef WITH_MUMPS

#ifdef HAVE_MPI

/* thread_local madness for the MPI communicator */

#ifndef __cplusplus

  #if __STDC_VERSION__ >= 201112L
    #include <threads.h>
    #define tlsvar thread_local
  #else

    #if defined(__GNUC__)
      #define tlsvar __thread 
    #else
      #error "Don't know how to create a thread-local variable"
    #endif
  #endif

#else

  #if SICONOS_CXXVERSION >= 201103L
    #define tlsvar thread_local
  #else
    #if defined(__GNUC__)
      #define tlsvar __thread
    #elif defined(_MSC_VER)
      #define tlsvar __declspec(thread)
    #else
      #error "Don't know how to create a thread-local variable"
    #endif
  #endif

#endif



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
#elif defined(__GNU__) & !defined(__APPLE__)

  void cleanup_MPI(void) __attribute__ ((destructor))
  {
    if (NM_mpi_com != MPI_COMM_NULL)
    {
      MPI_Finalize();
      NM_mpi_com = MPI_COMM_NULL;
    }
  }

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

int* NM_MUMPS_irn(NumericsMatrix* A)
{

  if (NM_sparse(A)->triplet)
  {
    CSparseMatrix* triplet = NM_sparse(A)->triplet;
    csi nz = triplet->nz;

    int* iWork = NM_iWork(A, (int) (2*nz) + 1);

    for (int k=0; k<nz; ++k)
    {
      iWork [k + nz] = (int) (triplet->p [k]) + 1;
      iWork [k]      = (int) (triplet->i [k]) + 1;
    }

    iWork [2*nz] = (int) nz;
  }
  else
  {
    CSparseMatrix* csc = NM_sparse(A)->csc;
    csi nzmax = csc->nzmax ;

    int* iWork = NM_iWork(A, (int) (2*nzmax) + 1);

    csi n = csc->n ;
    csi nz = 0;
    csi* csci = csc->i ;
    csi* cscp = csc->p ;

    for (csi j=0; j<n; ++j)
    {
      for (csi p = cscp [j]; p < cscp [j+1]; ++p)
      {
        assert (csc->x [p] != 0.);
        nz++;
        iWork [p + nzmax] = (int) j;
        iWork [p]         = (int) csci [p];
      }
    }

    iWork [2*nzmax] = (int) nz;
  }

  return NM_iWork(A, 0);
}


int* NM_MUMPS_jcn(NumericsMatrix* A)
{
  if (NM_sparse(A)->triplet)
  {
    return NM_iWork(A, 0) + NM_sparse(A)->triplet->nz;
  }
  else
  {
    csi nzmax = NM_csc(A)->nzmax;
    return NM_iWork(A, 0) + nzmax;
  }
}


DMUMPS_STRUC_C* NM_MUMPS_id(NumericsMatrix* A)
{
  NumericsSparseLinearSolverParams* params = NM_linearSolverParams(A);
  DMUMPS_STRUC_C* mumps_id;

  if (!params->solver_data)
  {
    /* valgrind reports some conditional move on initialized data in MUMPS
     * --xhub */
    params->solver_data = calloc(1, sizeof(DMUMPS_STRUC_C));

    mumps_id = (DMUMPS_STRUC_C*) params->solver_data;

    // Initialize a MUMPS instance. Use MPI_COMM_WORLD.
    mumps_id->job = JOB_INIT;
    mumps_id->par = 1;
    mumps_id->sym = 0;

#ifdef HAVE_MPI
    if (NM_MPI_com(MPI_COMM_NULL) == MPI_COMM_WORLD)
    {
      mumps_id->comm_fortran = (MUMPS_INT) USE_COMM_WORLD;
    }
    else
    {
      mumps_id->comm_fortran = (MUMPS_INT) MPI_Comm_c2f(NM_MPI_com(MPI_COMM_NULL));
    }
#endif /* WITH_MPI */

    dmumps_c(mumps_id);

    if (verbose == 1)
    {
      mumps_id->ICNTL(1) = -1; // Error messages, standard output stream.
      mumps_id->ICNTL(2) = -1; // Diagnostics,    standard output stream.
      mumps_id->ICNTL(3) = -1; // Global infos,   standard output stream.

    }
    else if (verbose == 2)
    {
      mumps_id->ICNTL(1) = -1; // Error messages, standard output stream.
      mumps_id->ICNTL(2) = -1; // Diagnostics,    standard output stream.
      mumps_id->ICNTL(3) = 6; // Global infos,   standard output stream.

//      mumps_id->ICNTL(4) = 4; // Errors, warnings and information on
                              // input, output parameters printed.

      mumps_id->ICNTL(11) = 2; // Error analysis
    }
    else if (verbose >= 3)
    {
      mumps_id->ICNTL(1) = 6; // Error messages, standard output stream.
      mumps_id->ICNTL(2) = 6; // Diagnostics,    standard output stream.
      mumps_id->ICNTL(3) = 6; // Global infos,   standard output stream.

//      mumps_id->ICNTL(4) = 4; // Errors, warnings and information on
                              // input, output parameters printed.

//      mumps_id->ICNTL(10) = 1; // One step of iterative refinment
      mumps_id->ICNTL(11) = 1; // Error analysis
    }
    else
    {
      mumps_id->ICNTL(1) = -1;
      mumps_id->ICNTL(2) = -1;
      mumps_id->ICNTL(3) = -1;
    }

    mumps_id->ICNTL(24) = 1; // Null pivot row detection see also CNTL(3) & CNTL(5)
//      mumps_id->ICNTL(10) = -2; // One step of iterative refinment
    // ok for a cube on a plane & four contact points
    // computeAlartCurnierSTD != generated in this case...

    //mumps_id->CNTL(3) = ...;
    //mumps_id->CNTL(5) = ...;

    mumps_id = (DMUMPS_STRUC_C*) params->solver_data;
    mumps_id->n = (int) NM_triplet(A)->n;
    mumps_id->irn = NM_MUMPS_irn(A);
    mumps_id->jcn = NM_MUMPS_jcn(A);

    int nz;
    if (NM_sparse(A)->triplet)
    {
      nz = (int) NM_sparse(A)->triplet->nz;
      mumps_id->nz = nz;
      mumps_id->a = NM_sparse(A)->triplet->x;
    }
    else
    {
      nz = NM_linearSolverParams(A)->iWork[2 * NM_csc(A)->nzmax];
      mumps_id->nz = nz;
      mumps_id->a = NM_sparse(A)->csc->x;
    }
  }
  else
  {
    mumps_id = (DMUMPS_STRUC_C*) params->solver_data;
    DEBUG_EXPR_WE(double data_ptr = NULL;
        if (NM_sparse(A)->triplet) { data_ptr = NM_sparse(A)->triplet->x; }
        else { data_ptr = NM_sparse(A)->csc->x; }
        if (data_ptr != mumps_id->a) { fprintf(stderr, "the data array and mumps_id->a don't match: %p != %p\n", data_ptr, mumps_id-a); } )

  }

  return mumps_id;
}


void NM_MUMPS_free(void* p)
{
  NumericsSparseLinearSolverParams* params = (NumericsSparseLinearSolverParams*) p;
  DMUMPS_STRUC_C* mumps_id = (DMUMPS_STRUC_C*) params->solver_data;

  /* clean the mumps instance */
  mumps_id->job = -2;
  dmumps_c(mumps_id);

  /* Here we free mumps_id ...  */
  free(params->solver_data);
  params->solver_data = NULL;

}

void NM_MUMPS_extra_display(DMUMPS_STRUC_C* mumps_id)
{
  if (mumps_id->ICNTL(11) == 2 || mumps_id->ICNTL(11) == 1)
  {
    printf("MUMPS : inf norm of A is %g\n", mumps_id->RINFOG(4));
    printf("MUMPS : inf norm of x is %g\n", mumps_id->RINFOG(5));
    printf("MUMPS : component wise scaled residual %g\n", mumps_id->RINFOG(6));
    printf("MUMPS : backward error estimate omega1 %g\n", mumps_id->RINFOG(7));
    printf("MUMPS : backward error estimate omega2 %g\n", mumps_id->RINFOG(8));
    printf("MUMPS : \n");
  }

  if (mumps_id->ICNTL(11) == 1)
  {
    printf("MUMPS : estimate for error in solution %g\n", mumps_id->RINFOG(10));
    printf("MUMPS : condition number 1 %g\n", mumps_id->RINFOG(10));
    printf("MUMPS : condition number 2 %g\n", mumps_id->RINFOG(11));
  }
}
#endif

