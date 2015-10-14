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
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#ifndef MIN_MERIT
#define MIN_MERIT
/*!\file min_merit.h

  \brief Functions for the min based merit function

  A set of routines used in the min reformulation of a CP

  The min function is:
  \f{equation*}
  \mathbf{F}_{\mathrm{min}}(z) = \min( z, F(z))
  \f}
*/

#include "SiconosConfig.h"
#include "NumericsMatrix.h"

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Compute \f$\mathbf{F}_{\mathrm{min}}(z)\f$, possibly in the mixed case
   * \param n1 number of equality constraints
   * \param n2 size of the complementary variables
   * \param[in] z input vector
   * \param[in] F value of F
   * \param[out] Fmin returned vector
   */
  void F_min(int n1, int n2, double* restrict z, double* restrict F, double* restrict Fmin);

  /** Compute an element of Jac F_min
   * \param n1 number of equality constraints
   * \param n2 size of the complementarity variables
   * \param[in] z input vector
   * \param[in] F value of F
   * \param[in] nabla_F value of nabla_F
   * \param[out] H returned vector
   */
  void Jac_F_min(int n1, int n2, double* restrict z, double* restrict F, NumericsMatrix* nabla_F, NumericsMatrix* H);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
