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
#include "NM_MUMPS.h"
#ifdef WITH_MUMPS
#include "CSparseMatrix_internal.h"
#include "NumericsMatrix_internal.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "numerics_verbose.h"
#include "NM_MPI.h"

/* #define DEBUG_MESSAGES */
#include "debug.h"

void NM_MUMPS_free(void* p)
{
  NSM_linear_solver_params* params = (NSM_linear_solver_params*) p;
  free(params->linear_solver_data);
  params->linear_solver_data = NULL;
}

void NM_MUMPS_set_irn_jcn(NumericsMatrix* A)
{
  /* MUMPS works on triplet format. */

  CSparseMatrix* triplet = NM_triplet(A);
  CS_INT nz = triplet->nz;
  assert(nz > 0);

  /* TODO: do not allocate when sizeof(MUMPS_INT) == sizeof(CS_INT),
   * just do triplet->p[k]++*/
  MUMPS_INT* iWork = (MUMPS_INT*)NM_iWork(A, (size_t)(2*nz) + 1, sizeof(MUMPS_INT));

  for (size_t k=0 ; k < (size_t)nz; ++k)
  {
    iWork [k + nz] = (MUMPS_INT) (triplet->p [k]) + 1;
    iWork [k]      = (MUMPS_INT) (triplet->i [k]) + 1;
  }

  iWork [2*nz] = (MUMPS_INT) nz;

  NM_MUMPS_id(A)->irn = (MUMPS_INT*)NM_iWork(A, 0, 0);
  NM_MUMPS_id(A)->jcn = &((MUMPS_INT*)NM_iWork(A, 0, 0))[numericsSparseMatrix(A)->triplet->nz];
}



DMUMPS_STRUC_C* NM_MUMPS_id(NumericsMatrix* A)
{
  NSM_linear_solver_params* params = NSM_linearSolverParams(A);
  DMUMPS_STRUC_C* mumps_id;

  if (!params->linear_solver_data)
  {
    /* valgrind reports some conditional move on initialized data in MUMPS
     * --xhub */
    params->linear_solver_data = calloc(1, sizeof(DMUMPS_STRUC_C));
    mumps_id = (DMUMPS_STRUC_C*) params->linear_solver_data;
    mumps_id->job = 0;
  }
  mumps_id = (DMUMPS_STRUC_C*) params->linear_solver_data;
  return mumps_id;
}

void NM_MUMPS_set_id(NumericsMatrix* A, DMUMPS_STRUC_C* id)
{
  NSM_linearSolverParams(A)->linear_solver_data = id;
}

void NM_MUMPS(NumericsMatrix* A, int job)
{
#ifdef HAVE_MPI
  if (NM_MPI_rank(A)==0)
  {
    NM_MUMPS_id(A)->job = job;
    /* we send the job number for listening processes */
    DEBUG_PRINTF("MPI: %d sending job %d\n", NM_MPI_rank(A),job);
    CHECK_MPI(NM_MPI_comm(A), MPI_Bcast(&job, 1, MPI_INT, 0, NM_MPI_comm(A)));
    if (job)
    {
      dmumps_c(NM_MUMPS_id(A));
    }
  }
  else
  {
    int ijob=-1;
    /* Loop until the end of everything. Note: job = -2 is not the end
       of everything as some job = -1 may be done later. Here we have
       a convention: job = 0 (which is not part of MUMPS) means that
       we will never call again MUMPS and that we are going to stop
       all the processes. */
    while(ijob)
    {
      /* I am a listening process, I will stop when ijob = 0 */
      CHECK_MPI(NM_MPI_comm(A), MPI_Bcast(&ijob, 1, MPI_INT, 0, NM_MPI_comm(A)));
      DEBUG_PRINTF("MPI: %d receiving job %d\n", NM_MPI_rank(A),ijob);
      if (ijob)
      {
        NM_MUMPS_id(A)->job = ijob;
        dmumps_c(NM_MUMPS_id(A));
      }
    }
  }
#else
  NM_MUMPS_id(A)->job = job;
  dmumps_c(NM_MUMPS_id(A));
#endif
}

void NM_MUMPS_set_control_params(NumericsMatrix* A)
{
  DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);
  mumps_id->par = 1; /* host (rank=0) is also involved in computations */
  mumps_id->sym = 0; /* unsymmetric */

#ifdef HAVE_MPI
  if (NM_MPI_comm(A) == MPI_COMM_WORLD)
  {
    mumps_id->comm_fortran = (MUMPS_INT) USE_COMM_WORLD;
  }
  else
  {
    mumps_id->comm_fortran = (MUMPS_INT) MPI_Comm_c2f(NM_MPI_comm(A));
  }
#endif /* WITH_MPI */


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


}

void NM_MUMPS_set_problem(NumericsMatrix* A, double *b)
{
  /* numerics matrices are not distributed */
  if (NM_MPI_rank(A) == 0)
  {
    DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);
    mumps_id->n = (MUMPS_INT) NM_triplet(A)->n;
    NM_MUMPS_set_irn_jcn(A);

    MUMPS_INT nz;
    nz = (MUMPS_INT) NM_triplet(A)->nz;
    mumps_id->nz = nz;
    mumps_id->a = NM_triplet(A)->x;
    mumps_id->rhs = b;
  }
}

void NM_MUMPS_extra_display(NumericsMatrix* A)
{
  DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

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

