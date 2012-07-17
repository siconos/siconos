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

#ifndef FB_H
#define FB_H

/*!\file FischerBurmeister.h

  \brief Fischer Burmeister functions.

  A set of routines to compute the Fischer-Burmeister function and its jacobian.

  The Fischer-Burmeister function is defined as :
  \f[
  \phi(z,w) = \sqrt( z^2 + w^2) - z - w
  \f]

  This function is used to solve MLCP, MCP and NCP. The inequalities are rewritten using Fischer function with \f$ w = F(z) \f$ and solved with a semi-smooth Newton algorithm.

  For "mixed" problems (i.e. including equality constraints), the Fischer function is defined as :

  \f[ \phi_{mixed}(z,F(z)) =
  \left\lbrace \begin{array}{c}
  F_e(z) \\
  \sqrt( z^2 + F_i(z)^2) - z - F_i(z) \end{array}\right. \f]

  where index "e" stands for equalities part in F and "i" for inequalities.

  For details see the paper of Kanzow and Kleinmichel, "A New Class of Semismooth Newton-type Methods for Nonlinear
  Complementarity Problems", Computational Optimization and Applications 11, 227-251 (1998).

  The notations below are more or less those of this paper.

  \author Houari Khenous, Franck Pérignon last modification (17/07/2012)

*/

#ifdef __cplusplus
extern "C"
{
#endif

  /** Fischer Burmeister function, \f$ \phi(z,F(z)) \f$
      \param[in] : size of vector z
      \param[in] : vector z
      \param[in] : vector F(z)
      \param[in,out] : phi vector \f$ \phi(z,F(z)) \f$
  */
  void phi_FB(int size, double* z, double* F, double* phi);

  /** Jacobian of the Fischer Burmeister function, \f$ \nabla_z \phi(z,F(z)) \f$
      \param[in] : size of vector z
      \param[in] : vector z
      \param[in] : vector F(z)
      \param[in] : \f$ \nabla_z F(z) \f$
      \param[in,out] : \f$ \nabla_z \phi(z,F(z)) \f$.
  */
  void jacobianPhi_FB(int size, double* z, double* F, double* jacobianF, double* jacobianPhi);

  /** Mixed Fischer Burmeister function,
      \f[ \phi(z,F(z)) = \left\lbrace \begin{array}{c} F(z) \\ \sqrt( z^2 + F(z)^2) - z - F(z) \end{array}\right. \f], the upper for equalities and the rest for inequalities.
      \param[in] : number of equality constraints.
      \param[in] : number of complementarity constraints.
      \param[in] : vector z (size = sizeEq + sizeIneq)
      \param[in] : vector F(z)
      \param[in,out] : \f$ \phi(z,F(z)) \f$.
  */
  void phi_Mixed_FB(int sizeEq, int sizeIneq, double* z, double* F, double* phi);

  /** Jacobian of the mixed Fischer Burmeister function, \f$ \nabla_z \phi(z,F(z)) \f$
      \param[in] : number of equality constraints.
      \param[in] : number of complementarity constraints.
      \param[in] : vector z
      \param[in] : vector F(z)
      \param[in] : \f$ \nabla_z F(z) \f$
      \param[in,out] : \f$ \nabla_z \phi(z,F(z)) \f$ .
  */
  void jacobianPhi_Mixed_FB(int sizeEq, int sizeIneq, double* z, double* F, double* jacobianF, double* jacobianPhi);


#ifdef __cplusplus
}
#endif

#endif
