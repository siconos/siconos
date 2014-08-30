/* Siconos-Numerics, Copyright INRIA 2005-2013.
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
#ifndef BOX_H
#define BOX_H

/*!\file box.h
 * \brief Some helpers for dealing with polyhedra and polytopes
 *
 * \author Olivier Huber
*/

typedef struct {
  int id;
} generic_set;

/** \struct box box.h
 * Definition of a polytope in terms of (H,K) representation
 */
typedef struct box_constraints
{
  int id; /**< id of the structure, usually solver specific */
  double* lb; /**< lower bounds */
  double* ub; /**< upper bounds */
} box_constraints;

#include "NumericsConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** free a box struct
  * \param b the box struct to free
  */
  void free_box(box_constraints* b);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
