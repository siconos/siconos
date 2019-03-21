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

#ifndef NM_MUMPS_h
#define NM_MUMPS_h

#include "SiconosConfig.h"
#include "NumericsFwd.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

#ifdef WITH_MUMPS
#include <dmumps_c.h>

#ifndef MUMPS_INT
#define MUMPS_INT int
#endif

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1]
#define CNTL(I) cntl[(I)-1]
#define RINFOG(I) rinfog[(I)-1]

  void NM_MUMPS_set_irn_jcn(NumericsMatrix* A);

  /** Get (and create if necessary) the config data for MUMPS
   * \param A the matrix to be factorized
   */
  DMUMPS_STRUC_C* NM_MUMPS_id(NumericsMatrix* A);

  /** Set the config data for MUMPS
   * \param A the matrix
   * \param id the config data
   */
  void NM_MUMPS_set_id(NumericsMatrix* A, DMUMPS_STRUC_C* id);

  /** Free the config data for MUMPS
   * \param param p a pointer on the linear solver parameters
   */
  void NM_MUMPS_free(void* p);

  /** Set control parameters (mpi communicator, verbosity)
   * \param A the matrix holding the MUMPS config
   */
  void NM_MUMPS_set_control_params(NumericsMatrix* A);

  /** Set linear problem
   * \param A the matrix holding the MUMPS config
   * \param b a pointer on double values
   */
  void NM_MUMPS_set_problem(NumericsMatrix* A, double *b);

  /** Set parameters for the solve
   * \param A the matrix holding the MUMPS config
  void NM_MUMPS_set_params(NumericsMatrix* A);

  /** Display extra information about the solve
   * \param A the matrix holding the MUMPS config
   */
  void NM_MUMPS_extra_display(NumericsMatrix* A);

  /** call MUMPS. If several MPI process are running, only rank=0
   * returns after sending the job to all the process, the others will
   * wait for other job, until job=0.
   * \param A the matrix holding the MUMPS config
   * \param job the int that specify the MUMPS job.  job=0 makes all process return.
   */
  void NM_MUMPS(NumericsMatrix* A, int job);

#endif /* WITH_MUMPS */
  /** copy MUMPS id if compiled WITH_MUMPS, otherwise do nothing.
   * \param A, the source NumericsMatrix
   * \param B, the destination NumericsMatrix
   */
  void NM_MUMPS_copy(const NumericsMatrix* A, NumericsMatrix* B);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
