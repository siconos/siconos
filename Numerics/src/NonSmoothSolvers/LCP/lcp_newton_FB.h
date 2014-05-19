/* Siconos-Numerics, Copyright INRIA 2005-2014.
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

#ifndef LCP_NEWTON_FB
#define LCP_NEWTON_FB

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Compute F(z) = Mz + q
   * \param data_opaque a LinearComplementarityProblem but casted
   * \param[in] z non-basic variable
   * \param[out] w basic variable (result)
   */
  void FB_compute_F_lcp(void* data_opaque, double* z, double* w);

  /** Compute an element of JacF_FB, see Facchinei--Pang p. 808
   * \param data_opaque a LinearComplementarityProblem but casted
   * \param[in] z non-basic variable
   * \param[in] w basic variable
   * \param workV1 work vector which contains "z"
   * \param workV2 work vector
   * \param[out] H an element of JacF_FB
   */
  void FB_compute_H_lcp(void* data_opaque, double* z, double* w, double* workV1, double* workV2, double* H);

  /** Compute the error for termination, here lcp_compute_error
   * \param data_opaque a LinearComplementarityProblem but casted
   * \param[in] z non-basic variable
   * \param[in] w basic variable
   * \param notused not used here
   * \param[in] tol the tolerance
   * \param[out] err the error on the LCP (not FB)
   */
  void FB_compute_error_lcp(void* data_opaque, double* z, double* w, double* notused, double tol, double* err);

  /** Compute F_FB : \f${F_FB}_i = \sqrt(z_i^2 + F_i^2) - (z_i + F_i)\f$
   * \param data_opaque a LinearComplementarityProblem but casted
   * \param[in] z non-basic variable
   * \param[in] w basic variable
   * \param[out] F_FB value of the function
   */
  void lcp_FB(void* data_opaque, double* z, double* F, double* F_FB);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
