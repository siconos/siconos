/* Siconos-Numerics, Copyright INRIA 2005-2011.
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

#ifndef NCP_H
#define NCP_H

/*! \page NCProblem Nonlinear Complementarity Problems (NCP)

  \section ncpIntro The problem
  Find \f$z \in \mathcal{R}^n\f$ such that:\n\n
  \f$
  0 \le F(z) \perp z \ge 0 \\
  \f$

  \section ncpSolvers Available solvers/formulations:
   - Fischer-Burmeister (see NCP_FischerBurmeister.h)\n\n
   - NCP_Path() Interface to Path (Ferris)

*/

/*!\file NCP_Solvers.h
  \brief Functions related to NCP formulation and solvers.
  \author Franck Perignon, last modification: 21/05/2008
*/

#include "SparseBlockMatrix.h"
#include "NCP_Path.h"
#include "NCP_FixedP.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * This function checks the validity of the vector z as a solution \n
   * of the NCP : \n
   * \f$
   *    0 \le z \perp Mz + q \ge 0
   * \f$
   * \author Houari Khenous
   \warning temporary function - To be reviewed
  */
  void NCP_compute_error(int n, double *vec , double *q , double *z , int verbose, double *w, double *err);

  /** This function adapts the NCP_compute_error routine for M saved as a SparseBlockStructuredMatrix.
   * TEMPORARY FUNCTION, used to compare pfc3D and FrictionContact3D functions.
   * \author Franck Perignon
   */
  void NCP_block_compute_error(int n, SparseBlockStructuredMatrix *M , double *q , double *z , int verbose, double *w, double *err);

#ifdef __cplusplus
}
#endif

#endif
