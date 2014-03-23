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
#ifndef POLYHEDRON_H
#define POLYHEDRON_H

/*!\file polyhedron.h
 * \brief Some helpers for dealing with polyhedra and polytopes
 *
 * \author Olivier Huber
*/

/** \struct Polyhedron polyhedron.h
 * Definition of a polytope in terms of (H,K) representation
 */
typedef struct Polyhedron
{
  unsigned int size;
  double* H; /**< H matrix in an (H,K) representation of a polytope H x <= K */
  double* K; /**< K vector in an (H,K) representation of a polytope H x <= K */
  double* Heq; /**< H matrix in an (H,K) representation of a polytope Heq x = Keq */
  double* Keq; /**< K vector in an (H,K) representation of a polytope Heq x = Keq */
} Polyhedron;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** free a Polyhedron struct
  * \param poly the Polyhedron struct to free
  */
  void freePolyhedron(Polyhedron* poly);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
