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

#ifdef WITH_MUMPS
#include "CSparseMatrix_internal.h"
#include "NumericsMatrix_internal.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"

#include "debug.h"
#include "numerics_verbose.h"


#include "NM_MPI.h"

MUMPS_INT* NM_MUMPS_irn(NumericsMatrix* A)
{

  if (numericsSparseMatrix(A)->triplet)
  {
    CSparseMatrix* triplet = numericsSparseMatrix(A)->triplet;
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
  }
  else
  {
    fprintf(stderr, "NM_MUMPS_irn :: xhub doubt this code is correct");
    exit(EXIT_FAILURE);

#if 0
    CSparseMatrix* csc = numericsSparseMatrix(A)->csc;
    CS_INT nzmax = csc->nzmax ;

    MUMPS_INT* iWork = NM_iWork(A, (int) (2*nzmax) + 1);

    CS_INT n = csc->n ;
    CS_INT nz = 0;
    CS_INT* csci = csc->i ;
    CS_INT* cscp = csc->p ;

    for (CS_INT j=0; j<n; ++j)
    {
      for (CS_INT p = cscp [j]; p < cscp [j+1]; ++p)
      {
        assert (csc->x [p] != 0.);
        nz++;
        iWork [p + nzmax] = (MUMPS_INT) j;
        iWork [p]         = (MUMPS_INT) csci [p];
      }
    }

    iWork [2*nzmax] = (MUMPS_INT) nz;
#endif
  }

  return (MUMPS_INT*)NM_iWork(A, 0, 0);
}


MUMPS_INT* NM_MUMPS_jcn(NumericsMatrix* A)
{
  if (numericsSparseMatrix(A)->triplet)
  {
    return &((MUMPS_INT*)NM_iWork(A, 0, 0))[numericsSparseMatrix(A)->triplet->nz];
  }
  else
  {
    fprintf(stderr, "NM_MUMPS_irn :: xhub doubt this code is correct");
    exit(EXIT_FAILURE);
#if 0
    CS_INT nzmax = NM_csc(A)->nzmax;
    return NM_iWork(A, 0, 0) + nzmax;
#endif
  }
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

    // Initialize a MUMPS instance. Use MPI_COMM_WORLD.
    mumps_id->job = JOB_INIT;
    mumps_id->par = 1;
    mumps_id->sym = 0;

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

    /* numerics matrices are not distributed */
    if (NM_MPI_rank(A) == 0)
    {
      mumps_id->n = (MUMPS_INT) NM_triplet(A)->n;
      mumps_id->irn = NM_MUMPS_irn(A);
      mumps_id->jcn = NM_MUMPS_jcn(A);

      MUMPS_INT nz;
      nz = (MUMPS_INT) NM_triplet(A)->nz;
      mumps_id->nz = nz;
      mumps_id->a = NM_triplet(A)->x;
      }
    }
/* MB: mumps does not accept csc format for input
/*      else
      {
        fprintf(stderr, "NM_MUMPS_irn :: xhub doubt this code is correct");
        exit(EXIT_FAILURE);
#if 0
        nz = NSM_linearSolverParams(A)->iWork[2 * NM_csc(A)->nzmax];
        mumps_id->nz = nz;
        mumps_id->a = numericsSparseMatrix(A)->csc->x;
#endif
      }
*/

  else
  {
    mumps_id = (DMUMPS_STRUC_C*) params->linear_solver_data;
    DEBUG_EXPR_WE(double data_ptr = NULL;
                  if (numericsSparseMatrix(A)->triplet) { data_ptr = numericsSparseMatrix(A)->triplet->x; }
                  else { data_ptr = numericsSparseMatrix(A)->csc->x; }
                  if (data_ptr != mumps_id->a) { fprintf(stderr, "the data array and mumps_id->a don't match: %p != %p\n", data_ptr, mumps_id-a); } )

      }

  return mumps_id;
}


void NM_MUMPS_free(void* p)
{
  NSM_linear_solver_params* params = (NSM_linear_solver_params*) p;
  DMUMPS_STRUC_C* mumps_id = (DMUMPS_STRUC_C*) params->linear_solver_data;

  /* clean the mumps instance */
  mumps_id->job = -2;
  dmumps_c(mumps_id);

  /* Here we free mumps_id ...  */
  free(params->linear_solver_data);
  params->linear_solver_data = NULL;

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

