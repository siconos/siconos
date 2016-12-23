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

#ifndef NumericsMatrix_private_H
#define NumericsMatrix_private_H

/*!\file NumericsMatrix_private.h
 * \brief non-public functions and data structures
 */

#include "SiconosConfig.h"
#include "NumericsMatrix.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Free the internalData of a NumericsMatrix
   * \param m the matrix */
  void NM_internalData_free(NumericsMatrix* m);


#ifdef WITH_MUMPS

#ifdef HAVE_MPI
#include <mpi.h>
#endif /* HAVE_MPI */

#include <dmumps_c.h>

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1]
#define CNTL(I) cntl[(I)-1]
#define RINFOG(I) rinfog[(I)-1]

#ifdef HAVE_MPI
  /** Get the MPI communicator. Call MPI_Init if needed.
   * \param[in] m an MPI communicator
   * \return the MPI communicator.
   */
  MPI_Comm NM_MPI_com(MPI_Comm m);
#endif /*  HAVE_MPI */

  int* NM_MUMPS_irn(NumericsMatrix* A);
  int* NM_MUMPS_jcn(NumericsMatrix* A);

  /** Get (and create if necessary) the working data for MUMPS
   * \param A the matrix to be factorized
   */
  DMUMPS_STRUC_C* NM_MUMPS_id(NumericsMatrix* A);

  /** Free the working data for MUMPS
   * \param p a NumericsSparseLinearSolverParams object holding the data
   */
  void NM_MUMPS_free(void* p);

  /** Display extra information about the solve
   * \param mumps_id the working space of MUMPS
   */
  void NM_MUMPS_extra_display(DMUMPS_STRUC_C* mumps_id);

#endif

#ifdef WITH_UMFPACK
#include <umfpack.h>


#ifdef SICONOS_INT64
#define UMFPACKPREFIX(X) umfpack_dl ## X
#else
#define UMFPACKPREFIX(X) umfpack_di ## X
#endif

#define UMFPACK_FN(X) UMFPACKPREFIX(_ ## X)

/** \struct NM_UMFPACK_WS NumericsMatrix_private.h
 * Structure for holding the data UMFPACK needs
 */
typedef struct {
  void* symbolic; /**< for the symbolic analysis */
  void* numeric;  /**< for the numerical factorization */
  double control[UMFPACK_CONTROL]; /**< control parameters */
  double info[UMFPACK_INFO]; /**< informations from UMFPACK */
  csi* wi; /**< integer workspace, size n */
  double* wd; /**< double workspace, size: with iterative refinement: 5n, without: n */
  double* x; /**< solution of the problem, size n */
} NM_UMFPACK_WS;

  /** Get (and create if necessary) the working data for UMPFACK
   * \param A the matrix to be factorized
   */
  NM_UMFPACK_WS* NM_UMFPACK_ws(NumericsMatrix* A);

  /** Free the working data for UMFPACK
   * \param p a NumericsSparseLinearSolverParams object holding the data
   */
  void NM_UMFPACK_free(void* p);

  /** Display extra information about the solve
   * \param umfpack_ws the working space of UMFPACK
   */
  void NM_UMFPACK_extra_display(NM_UMFPACK_WS* umfpack_ws);

  /** Factorize a matrix
   * \param A the matrix to factorize
   * \return the workspace containing the factorized form and other infos
   */
  NM_UMFPACK_WS* NM_UMFPACK_factorize(NumericsMatrix* A);

#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
