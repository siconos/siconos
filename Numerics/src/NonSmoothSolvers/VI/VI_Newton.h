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

#ifndef VI_NEWTON_H
#define VI_NEWTON_H

/*!\file VI_Newton.h
 * \brief Functions for solving VI using Newton method
 */

#include "SiconosConfig.h"
#include "NumericsMatrix.h"

#if defined(__cplusplus)
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  void VI_compute_F(void* data_opaque, double* x, double* F);
  void VI_compute_error_box(void* data_opaque, double* x, double* F, double* Jac_F_merit, double tol, double* err);
  void VI_compute_F_box_Qi(void* data_opaque, double* x, double* F, double* Fbox);
  void VI_compute_H_box_Qi(void* data_opaque, double* x, double* F, double* workV1, double* workV2, NumericsMatrix* H);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
