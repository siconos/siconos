/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

#ifndef FC3D2NCP_GLOCKER_H
#define FC3D2NCP_GLOCKER_H

/*!\file FrictionContact3D2NCP_Glocker.h
  interface to functions used to write the Friction-Contact 3D problem
  as a NCP, using Glocker formulation.

  The input problem looks like:
  \f[
  velocity = M.reaction +q
  \f]
  velocity, reaction (the unknowns) and q are vectors of size n. M is a nXn matrix.

  All details concerning this formulation are given in:
  Acary, V. and B. Brogliato (2008). Numerical Methods for Nonsmooth Dynamical Systems: Applications
  in Mechanics and Electronics. Vol. 35 of LNACM. Springer Verlag.
  Chapter 13, part 4.3 p 450.

*/
#include "SparseBlockMatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

  /** Initialize some parameters required for the formulation
     \param size of the problem (number of contacts X 3), n
     \param the global matrix M
     \param the global vector q
     \param the global vector mu of the friction coefficients (size = n/3)
  */
  void NCPGlocker_initialize(int, const double*const, const double*const, const double*const);

  /**  */
  void NCPGlocker_initialize_SBS(int n0, const SparseBlockStructuredMatrix*const M0, const double*const q0, const double*const mu0);

  /** Pick the required sub-blocks in q, M ... according to the considered contact and write the
     operators required for the Glocker formulation
     \param the number of the considered contact (its position in global M)
     \param the global reaction vector (size n)
  */
  void NCPGlocker_update(int, double*);

  /** Retrieve global reaction values after solving, from computed "reactionGlocker".
     \param the number of the considered contact
     \param the global reaction (in-out parameter)
  */
  void NCPGlocker_post(int, double *);

  /** To compute F
     \param the resulting FOut, in-out parameter (warning: must be null on input)
     \param up2Date, bool to avoid recomputation of some parameters: true if F or jacobianF has been computed and if
     the considered local problem has not changed, else false.
  */
  void computeFGlocker(double **, int);

  /** To compute jacobianF
     \param the resulting jacobianFOut, in-out parameter (warning: must be null on input)
     \param up2Date, bool to avoid recomputation of some parameters: true if F or jacobianF has been computed and if
     the considered local problem has not changed, else false.
  */
  void computeJacobianFGlocker(double **, int);

  /** free memory for friction contact to NCP-Glocker */
  void NCPGlocker_free();

#ifdef __cplusplus
}
#endif

#endif
