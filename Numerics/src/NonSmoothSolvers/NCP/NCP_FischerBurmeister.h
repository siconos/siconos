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

#ifndef NCPFB_H
#define NCPFB_H

/*!\file NCP_FischerBurmeister.h

  \brief Fischer Burmeister function for NCP

   It provides routines to compute Fischer-Burmeister function and its jacobian,
   written as a NCP-function:

    \f[
      \phi(z,F(z)) = \sqrt( z^2 + F(z)^2) - z - F(z)
    \f]

   For details see the paper of Kanzow and Kleinmichel, "A New Class of Semismooth Newton-type Methods for Nonlinear
   Complementarity Problems", Computational Optimization and Applications 11, 227-251 (1998).

   The notations below are more or less those of this paper.

   Functions:

   phi_FB(int size, double* z, double* F, double* phiVector)

   jacobianPhi_FB(int size, double* z, double* F, double* jacobianF, double* phiVector)

   \author Houari Khenous, Franck Perignon last modification (13/12/2007)

*/

#ifdef __cplusplus
extern "C" {
#endif

  /** NCP Fischer Burmeister function, \f$ \phi(z,F(z)) \f$
      \param size of vector z
      \param vector z
      \param vector F(z)
      \param vector \f$ \phi(z,F(z)) \f$, in-out arg.
  */
  void phi_FB(int, double*, double*, double*);

  /** Jacobian of NCP Fischer Burmeister function, \f$ \nabla_z \phi(z,F(z)) \f$
      \param size of vector z
      \param vector z
      \param vector F(z)
      \param \f$ \nabla_z F(z) \f$
      \param \f$ \nabla_z \phi(z,F(z)) \f$, in-out arg.
  */
  void jacobianPhi_FB(int, double*, double*, double*, double*);


#ifdef __cplusplus
}
#endif

#endif
