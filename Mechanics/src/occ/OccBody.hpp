/* Siconos-Kernel  Copyright INRIA 2005-2014.
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
 * Contact: Vincent ACARY siconos-team@lists.gforge.inria.fr
 */
/** \file OccBody.hpp
    \brief A Siconos Newton Euler dynamical system with
    associated contact shapes
 */

#ifndef OccBody_hpp
#define OccBody_hpp

#include "MechanicsFwd.hpp"
#include <NewtonEulerDS.hpp>

class OccBody : public NewtonEulerDS
{
public:
  //  Note: with c++11 =>
  //  using NewtonEulerDS::NewtonEulerDS;


  /* Default constructor.
   */
  OccBody() : NewtonEulerDS() {};


  //! Constructor from a minimum set of data.
  //  \param position : initial coordinates of this DynamicalSystem.
  //  \param velocity: initial velocity of this DynamicalSystem.
  //  \param mass : the mass.
  //  \param inertia : the inertia matrix.
  //
  OccBody(SP::SiconosVector position,
          SP::SiconosVector velocity,
          double mass ,
          SP::SiconosMatrix inertia);

  /** Association of a contact shape.
   * \param shape : the contact shape.
   */
  void addContactShape(SP::OccContactShape shape);

  /** Update positions and orientations of contact shapes.
   */
  void updateContactShapes();

  /** Get an associated contact shapes by its rank of association.
      \param id : the number of the shape.
   */
  OccContactShape& contactShape(unsigned int id) const;

  ACCEPT_STD_VISITORS();

protected:
  SP::ContactShapes _contactShapes;

  ACCEPT_SERIALIZATION(OccBody);
};

#endif
