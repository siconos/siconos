/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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

/*! \file PluginTypes.hpp
  \brief list of typedef for pointers to functions used in plugin mechanism.
*/

#ifndef PLUGINTYPES_HPP
#define PLUGINTYPES_HPP

/** Pointer to function used for plug-in for matrix-type operators that depends only on time */
typedef void (*MatrixFunctionOfTime)(double, unsigned int, unsigned int, double*, unsigned int, double*);

/** Pointer to function used for plug-in for vector-type operators that depends only on time */
typedef void (*VectorFunctionOfTime)(double, unsigned int, double*, unsigned int, double*);

#include "SimpleMatrix.h"
#include "PluggedObject.hpp"

/** Matrix plugged to a MatrixFunctionOfTime */
typedef PluggedObject<MatrixFunctionOfTime, SimpleMatrix> Plugged_Matrix_FTime;

/** Vector plugged to a VectorFunctionOfTime */
typedef PluggedObject<VectorFunctionOfTime, SimpleVector> Plugged_Vector_FTime;

/** Smart pointer to a Plugged_Matrix_FTime */
TYPEDEF_SPTR(Plugged_Matrix_FTime);

/** Smart pointer to a Plugged_Vector_FTime */
TYPEDEF_SPTR(Plugged_Vector_FTime);

/** */
typedef void (*FPtr1)(double, unsigned int, const double*, double*, unsigned int, double*);

/** */
typedef void (*FPtr2)(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*);

/** */
typedef void (*FPtr3)(unsigned int, const double*, unsigned int, double*, unsigned int, double*);

/** */
typedef void (*FPtr4)(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*);

/** */
typedef void (*FPtr5)(unsigned int, const double*, const double*, double*, unsigned int, double*);

/** */
typedef void (*FPtr6)(double, unsigned int, const double*, const double*, double*, unsigned int, double*);

/** */
typedef void (*FPtr7)(unsigned int, const double*, double*, unsigned int, double*);

#endif
