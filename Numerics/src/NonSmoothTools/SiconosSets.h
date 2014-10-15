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

#ifndef SICONOS_SETS_H
#define SICONOS_SETS_H

/*!\file SiconosSets.h
 * \brief Sets structures used in Siconos:
 * - box constraints \f$ K = \{ x \in \mathbb{R}^n | lb_i \leq x_i \leq ub_i\quad i = 1 .. n\}\f$
 * - polytopes and polyhedra \f$ \{ x \in \mathbb{R}^n | Hx\leq K\}\f$
 * - sets defined by a set of inequalities \f$\{g_i(x)\leq 0\}\f$ ( work in progress)
 *
 * \author Olivier Huber
*/

/** \struct generic_set SiconosSets.h
 * Generic set (can be seen as a kind of ``base class''). Mainly used to infer
 * the type of set (box, polytope, ...) to properly operate on it
 */
typedef struct {
  int id; /**< type of the set */
} generic_set;

/** The positive orthant does not need to contain much info */
typedef generic_set positive_orthant;

/** \struct box_constraints SiconosSets.h
 * Definition of a rectangular set, also known as box 
 */
typedef struct
{
  int id; /**< id of the structure, usually solver specific */
  double* lb; /**< lower bounds */
  double* ub; /**< upper bounds */
} box_constraints;

/** \struct Polyhedron SiconosSets.h
 * Definition of a polytope in terms of (H,K) representation
 */
typedef struct
{
  int id; /**< id of the structure, usually solver specific */
  unsigned size_ineq; /**< number of inequalities */
  unsigned size_eq; /**< number of equalities */
  double* H; /**< H matrix in an (H,K) representation of a polytope H x <= K */
  double* K; /**< K vector in an (H,K) representation of a polytope H x <= K */
  double* Heq; /**< Heq matrix for the equality constraints Heq x = Keq */
  double* Keq; /**< Keq vector for the equality constraints Heq x = Keq */
} polyhedron;

enum SICONOS_SET_ID { SICONOS_SET_POSITIVE_ORTHANT, SICONOS_SET_BOX, SICONOS_SET_POLYHEDRON };

#include "NumericsConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** project the point x on a set
   * \param n the size of x
   * \param x the point to project
   * \param set the set on which we project x
   */
  void project_on_set(int n, double* x, void* set);

  /** free a set
   * \param set struct to be freed
   */
  void free_siconos_set(void* set);

  /** free a box struct
  * \param b the box struct to free
  */
  void free_box(box_constraints* b);

  /** free a Polyhedron struct
  * \param poly the Polyhedron struct to free
  */
  void free_polyhedron(polyhedron* poly);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif


