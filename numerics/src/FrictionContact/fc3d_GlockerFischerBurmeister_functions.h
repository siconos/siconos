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
#ifndef FRICTIONCONTACT3D_GlockerFischerBurmeister_functions_H
#define FRICTIONCONTACT3D_GlockerFischerBurmeister_functions_H

/*!\file fc3d_GlockerFischerBurmeister_function.h

  \brief Typedef and functions declarations related GlockerFischerBurmeister
 */
#include "FrictionContactProblem.h"
#include "SparseBlockMatrix.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
typedef void (*UpdateSolverPtr)(int, double*);


  void F_GlockerFischerBurmeister(int sizeF, double* reaction, double* FVector, int up2Date);


  void jacobianF_GlockerFischerBurmeister(int sizeF, double* reaction, double* jacobianFMatrix, int up2Date);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
