/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#ifndef FC3D2NCP_GLOCKER_H
#define FC3D2NCP_GLOCKER_H

/*!\file FrictionContact3D2NCP_Glocker.h
  \brief interface to functions used to write the Friction-Contact 3D problem
  as a NCP, using Glocker formulation.

  The idea is to write the Friction Contact problem as a NCP:
  The input problem looks like:
  \f[
  velocity = M.reaction +q
  \f]
  velocity, reaction (the unknowns) and q are vectors of size n. M is a nXn matrix.

  and is formulate as:
  Find \f$reaction_G \in \mathcal{R}^5\f$ such that:\n\n
  \f$
  0 \le F_G(reaction_G) \perp reaction_G \ge 0 \\
  \f$
  with
  \f$ F_G(reaction_G) = M_G.reaction_G + g(reaction_G) + q_G \f$
  index \f$ G \f$ stands for "Glocker" in all related operators.

  The relations between \$_G\$ operators and input Friction-Contact problem, plus all other details, are given in:
  Acary, V. and B. Brogliato (2008). Numerical Methods for Nonsmooth Dynamical Systems: Applications
  in Mechanics and Electronics. Vol. 35 of LNACM. Springer Verlag.
  Chapter 13, part 4.3 p 420.



*/
#include "NumericsMatrix.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Initialize some parameters required for the formulation as
     the size of the problem (number of contacts X 3), n
     the the global matrix M
     the the global vector q
     the global vector mu of the friction coefficients (size = n/3)
     \param problem the global problem
     \param localproblem the local problem
  */
  void NCPGlocker_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem);

  /** Pick the required sub-blocks in q, M ... according to the considered contact and write the
     operators required for the Glocker formulation
     \param problem the global problem
     \param localproblem the local problem
     \param pos the number of the considered contact (its position in global M)
     \param options
  */
  void NCPGlocker_update(int, FrictionContactProblem* problem, FrictionContactProblem* localproblem,  double* pos, SolverOptions* options);

  /** Retrieve global reaction values after solving, from computed "reactionGlocker".
     \param contactnumber the number of the considered contact
     \param[in,out] reaction he global reaction (in-out parameter)
  */
  void NCPGlocker_post(int contactnumber, double * reaction);

  /** To compute F
     \param[in,out] FOut the resulting FOut (warning: must be null on input)
     \param up2Date  boolean variable to avoid recomputation of some parameters: true if F or jacobianF has been computed and if
     the considered local problem has not changed, else false.
  */
  void computeFGlocker(double ** FOut, int up2Date);

  /** To compute jacobianF
     \param[in,out] jacobianFOut the resulting (warning: must be null on input)
     \param up2Date Boolean variable to avoid recomputation of some parameters: true if F or jacobianF has been computed and if
     the considered local problem has not changed, else false.
  */
  void computeJacobianFGlocker(double ** jacobianFOut, int up2Date);

  /** compute NCP error for Fischer-Burmeister formulation **/
  double Compute_NCP_error1(int contact, double error);
  double Compute_NCP_error2(int contact, double error);

  /** compute Fixed Point Solution for the NCP formulation  **/
  void compute_Z_GlockerFixedP(int contact, double *reactionstep);

  /** free memory for friction contact to NCP-Glocker */
  void NCPGlocker_free();

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
