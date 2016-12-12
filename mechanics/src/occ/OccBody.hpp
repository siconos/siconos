/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
/** \file OccBody.hpp
    \brief A Siconos Newton Euler dynamical system with
    associated contact shapes
 */

#ifndef OccBody_hpp
#define OccBody_hpp

#include "MechanicsFwd.hpp"

#include "SiconosFwd.hpp"
#include "NewtonEulerDS.hpp"

class OccBody : public NewtonEulerDS
{
public:
  //  Note: with c++11 =>
  //  using NewtonEulerDS::NewtonEulerDS;


  /* Default constructor.
   */
  OccBody() : NewtonEulerDS() {};


  //! Constructor from a minimum set of data.
  //  \param position initial coordinates of this DynamicalSystem.
  //  \param velocity initial velocity of this DynamicalSystem.
  //  \param mass the mass.
  //  \param inertia the inertia matrix.
  //
  OccBody(SP::SiconosVector position,
          SP::SiconosVector velocity,
          double mass ,
          SP::SiconosMatrix inertia);

  /** Association of a contact shape.
   * \param shape the contact shape.
   * \param position relative position (x, y, z).
   * \param orientation relative orientation quaternion w, x, y, z
   * \param group contact group default 0
   */
  void addContactShape(SP::OccContactShape shape,
                       SP::SiconosVector position = SP::SiconosVector(),
                       SP::SiconosVector orientation = SP::SiconosVector(),
                       unsigned int group=0);


  /** Association of a shape without contact.
   * \param shape the shape
   * \param position relative position (x, y, z).
   * \param orientation relative orientation quaternion w, x, y, z
   */
  void addShape(SP::TopoDS_Shape shape,
                SP::SiconosVector position = SP::SiconosVector(),
                SP::SiconosVector orientation = SP::SiconosVector());

  /** Update positions and orientations of contact shapes.
   */
  void updateContactShapes();

  /** Update positions and orientations of shapes.
   */
  void updateShapes();

  /** Get an associated contact shape by its rank of association.
   *  \param id the number of the shape.
   */
  const OccContactShape& contactShape(unsigned int id) const;

  /** Get an associated shape by its rank of association.
   *  \param id the number of the shape.
   */
  const TopoDS_Shape& shape(unsigned int id) const;

  ACCEPT_BASE_STD_VISITORS(NewtonEulerDS);

protected:
  SP::ContactShapes _contactShapes;
  SP::TopoDS_Shapes _shapes;

  ACCEPT_SERIALIZATION(OccBody);
};

#endif
