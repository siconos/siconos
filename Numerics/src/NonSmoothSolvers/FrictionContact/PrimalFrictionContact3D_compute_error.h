/* Siconos-Numerics, Copyright INRIA 2005-2010.
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

#ifndef PrimalFrictionContact3D_compute_error_H
#define PrimalFrictionContact3D_compute_error_H
#include "PrimalFrictionContactProblem.h"
/*!\file FrictionContact3D_compute_error.h
  \brief functions related to error computation for friction-contact problems

  \author Vincent Acary, 26/05/2008

*/

#ifdef __cplusplus
extern "C"
{
#endif

  /** Error computation for friction-contact 3D problem
      \param problem the structure which defines the friction-contact problem
      \param reaction
      \param velocity
      \param globalVelocity
      \param tolerance value for error computation
      \param[in,out] error value
   */
  void PrimalFrictionContact3D_compute_error(PrimalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, double tolerance, double * error);

#ifdef __cplusplus
}
#endif

#endif
