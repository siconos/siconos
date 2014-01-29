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
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef FRICTIONCONTACTPROBLEM_AS_VI_H
#define FRICTIONCONTACTPROBLEM_AS_VI_H

/*! \page fcProblem Friction-contact problems (2D or 3D) as VI
 *
 *
 */


/*!\file FrictionContactProblem_as_VI.h
  \brief Definition of a structure to handle with friction-contact (2D or 3D) problems.
*/

#include "NumericsMatrix.h"
#include "FrictionContactProblem.h"
#include "VariationalInequality.h"

/** \struct FrictionContactProblem_as_VI FrictionContactProblem_as_VI.h
 *
 */
typedef struct
{
  /* the VI associated with the FC3D probelem */
  VariationalInequality * vi;
  /* the FC3D associated with the VI  */
  FrictionContactProblem * fc3d;
} FrictionContactProblem_as_VI;



#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  void Function_VI_FC3D(void * self, double *x, double *F);

  void Projection_VI_FC3D(void *viIn, double *x, double *PX);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
