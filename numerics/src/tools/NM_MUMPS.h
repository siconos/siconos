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

  /** Get (and create if necessary) the config data for MUMPS.
   * \param A, the matrix to be factorized.
   */
  DMUMPS_STRUC_C* NM_MUMPS_id(NumericsMatrix* A);

  /** Set the config data for MUMPS.
   * \param A the matrix,
   * \param id the config data,
   * \return the DMUMPS_STRUC_C config.
   */
  void NM_MUMPS_set_id(NumericsMatrix* A, DMUMPS_STRUC_C* id);

  /** Set control parameters (mpi communicator, verbosity).
   * \param A, the matrix holding the MUMPS config.
   */
  void NM_MUMPS_set_control_params(NumericsMatrix* A);

  /** Set icntl control.
   * \param A, the matrix holding the MUMPS config,
   * \param index, the fortran index in the icntl array,
   * \param val, the new integer value.
   */
  void NM_MUMPS_set_icntl(NumericsMatrix* A, unsigned int index, int val);

  /** Return icntl value.
   * \param A, the matrix holding the MUMPS config,
   * \param index, the FORTRAN index in the icntl array,
   * \return the ICNTL(index) value.
   */
  int NM_MUMPS_icntl(NumericsMatrix* A, unsigned int index);

  /** Set cntl control.
   * \param A, the matrix holding the MUMPS config,
   * \param index, the FORTRAN index in the cntl array,
   * \param val the new double value.
   */
  void NM_MUMPS_set_cntl(NumericsMatrix* A, unsigned int index, double val);

   /** Return cntl value.
   * \param A, the matrix holding the MUMPS config,
   * \param index, the FORTRAN index in the cntl array.
   * \return the CNTL(index) value.
   */
  double NM_MUMPS_cntl(NumericsMatrix* A, unsigned int index);

  /** Set linear problem
   * \param A, the matrix holding the MUMPS config,
   * \param b, a pointer on double values.
   */
  void NM_MUMPS_set_problem(NumericsMatrix* A, double *b);

  /** Set MUMPS verbosity.
   * \param A, the matrix holding the MUMPS config,
   * \param verbosity, an integer: 0 for silent or 1 for verbose.
   */
  void NM_MUMPS_set_verbosity(NumericsMatrix* A, unsigned int verbosity);

  /** Set some default parameters for the solve.
   * \param A the matrix holding the MUMPS config
   */
  void NM_MUMPS_set_default_params(NumericsMatrix* A);

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

  /** Free the config data for MUMPS
   * \param param p a pointer on the linear solver parameters
   */
  void NM_MUMPS_free(void* p);


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
