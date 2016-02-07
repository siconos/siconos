/* Siconos-Numerics, Copyright INRIA 2005-2015
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#ifndef NumericsMatrix_private_H
#define NumericsMatrix_private_H

#include "SiconosConfig.h"
#include "NumericsMatrix.h"

#ifdef WITH_MUMPS
#include <mpi.h>
#include <dmumps_c.h>

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1]
#define CNTL(I) cntl[(I)-1]
#define RINFOG(I) rinfog[(I)-1]

  /** Get the MPI communicator. Call MPI_Init if needed.
   * \param[in] m an MPI communicator
   * \return the MPI communicator.
   */
  MPI_Comm NM_MPI_com(MPI_Comm m);

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

#endif
