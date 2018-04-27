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

#include "SiconosConfig.h"
#include "NumericsMatrix.h"

#ifdef __cplusplus
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Fischer Burmeister function, \f$ \phi(z,F(z)) \f$
      \param[in] size of vector z
      \param[in] z vector \f$ z \f$
      \param[in] F vector \f$ F(z) \f$
      \param[in,out]  phi vector \f$ \phi(z,F(z)) \f$
  */
  void phi_FB(int size, double* restrict z, double* restrict F, double* restrict phi);

  /** Jacobian of the Fischer Burmeister function, \f$ \nabla_z \phi(z,F(z)) \f$
      \param[in] size of vector \f$ z \f$
      \param[in] z vector \f$ z \f$
      \param[in] F vector \f$ F(z) \f$
      \param[in] jacobianF \f$ \nabla_z F(z) \f$
      \param[in,out]  jacobianPhi \f$ \nabla_z \phi(z,F(z)) \f$.
      \warning this function looks broken !
  */
  void jacobianPhi_FB(int size, double* z, double* F, double* jacobianF, double* jacobianPhi);

  /** Mixed Fischer Burmeister function,
      \f[ \phi(z,F(z)) = \left\lbrace \begin{array}{c} F(z) \\ \sqrt( z^2 + F(z)^2) - z - F(z) \end{array}\right. \f], the upper for equalities and the rest for inequalities.
      \param[in] sizeEq number of equality constraints.
      \param[in] sizeIneq number of complementarity constraints.
      \param[in] z vector z (size = sizeEq + sizeIneq)
      \param[in] F vector F(z)
      \param[in,out] phi \f$ \phi(z,F(z)) \f$.
  */
  void phi_Mixed_FB(int sizeEq, int sizeIneq, double* restrict z, double* restrict F, double* restrict phi);

  /** Jacobian of the mixed Fischer Burmeister function, \f$ \nabla_z \phi(z,F(z)) \f$
      \param[in] sizeEq number of equality constraints.
      \param[in] sizeIneq number of complementarity constraints.
      \param[in] z vector \f$z\f$
      \param[in] F vector \f$F(z)\f$
      \param[in] jacobianF \f$ \nabla_z F(z) \f$
      \param[in,out] jacobianPhi \f$ \nabla_z \phi(z,F(z)) \f$ .
      \warning this function looks broken !
  */
  void jacobianPhi_Mixed_FB(int sizeEq, int sizeIneq, double* z, double* F, double* jacobianF, double* jacobianPhi);

  /** Computes an element of \f$Jac \mathbf{F}_{\mathrm{FB}}\f$ (possibly mixed) Fischer-Burmeister function, see Facchinei--Pang (2003) p. 808
      \param[in] n1 number of equality constraints.
      \param[in] n2 number of complementarity constraints.
      \param[in] z vector \f$z\f$
      \param[in] F vector \f$F(z)\f$
      \param[in] workV1 work vector (value gets overwritten)
      \param[in] workV2 work vector (value gets overwritten)
      \param[in] nabla_F \f$ \nabla_z F(z) \f$
      \param[in,out] H element of Jac_F_merit
  */
void Jac_F_FB(int n1, int n2, double* restrict z, double* restrict F, double* restrict workV1, double* restrict workV2, NumericsMatrix* restrict nabla_F, NumericsMatrix* restrict H);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
