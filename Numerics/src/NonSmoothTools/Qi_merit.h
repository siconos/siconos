/* Siconos-Numerics, Copyright INRIA 2005-2014
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
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#ifndef QI_MERIT_H
#define QI_MERIT_H

/*!\file Qi_merit.h
  \brief functions related to the Qi C-functions used as a merit function for box VI problems

  Reference Facchinei--Pang pp.869 - 877.

  \author Olivier Huber, 21/05/2014

*/

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Evaluates the C function for a box-constrained VI
   * \param n size of the problem
   * \param[in] x box-constrained variable of the VI
   * \param[out] F value of the function
   * \param[out] Fbox value of the function
   * \param[in] lb lower bounds, that is lb <= x
   * \param[in] ub upper bounds, that is ub >= x
   * */
  void phi_Qi(int n, double* x, double* F, double* Fbox, double* lb, double* ub);

  /** Evaluates the Jacobian of the C function for a box-constrained VI
   * \param n size of the problem
   * \param[in] x box-constrained variable of the VI
   * \param[out] Fbox value of the function
   * \param workV1 work vector
   * \param workV2 work vector
   * \param[in] nabla_Fbox gradient of the C-function
   * \param[in] lb lower bounds, that is lb <= x
   * \param[in] ub upper bounds, that is ub >= x
   * \param[out] H an element of the Jacobian
   * */
  void Jac_F_Qi(int n, double* x, double* Fbox, double* workV1, double* workV2, double* nabla_Fbox, double* lb, double* ub, double* H);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
