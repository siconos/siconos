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
/** \file OccContactShape.hpp
    \brief OpenCascade contact shape
 */
#ifndef OccContactShape_hpp
#define OccContactShape_hpp

#include "MechanicsFwd.hpp"

#include <SiconosFwd.hpp>
#include <SiconosVisitor.hpp>
#include <string>

DEFINE_SPTR(TopoDS_Shape);
DEFINE_SPTR(TopoDS_Edge);
DEFINE_SPTR(TopoDS_Face);


struct OccContactShape
{

  /** Default constructor.
   */
  OccContactShape();

  /** Constructor from OpenCascade object.
      \param shape
   */
  OccContactShape(TopoDS_Shape& shape)
    : _shape(createSPtrTopoDS_Shape(shape))
  {};

  /** Constructor from const OpenCascade object : remove constness.
      \param shape
   */
  OccContactShape(const TopoDS_Shape& shape)
    : _shape(createSPtrTopoDS_Shape(const_cast<TopoDS_Shape&>(shape)))
  {};

  /** Constructor from OccContactShape
      \param shape : a SP::OccContactShape
   */
  OccContactShape(const OccContactShape& shape)
    : _shape(shape.shape()) {};

  /** Destructor.
   */
  virtual ~OccContactShape() {};

  /** Return shared pointer on data
   */
  SP::TopoDS_Shape shape() const { return _shape;};

  /** Return OpenCascade data.
   */
  TopoDS_Shape& data() const { return *_shape;};

  /** Known contacts.
   */
  enum ContactTypeValue { Face, Edge, Unknown };

  /** Get the contact type.
   */
  ContactTypeValue contactType() const;

  /** Get a face from its index.
      \param index : the index of the face.
  */
  SPC::TopoDS_Face face(unsigned int index) const;

  /** Get an edge from its index.
      \param index : the index of the face.
  */
  SPC::TopoDS_Edge edge(unsigned int index) const;

  /** Export the contact shape into a string.
   */
  std::string exportBRepToString() const;

  /** Import the contact shape from a string.
   *  \param brepstr : the string containing the whole brep.
   */
  void importBRepFromString(const std::string& brepstr);

  /** Compute and store UV bounds of the shape.
   */
  virtual void computeUVBounds();


  /** Set id.
   *  \param id the new id
  */
  void setId(unsigned int id) { this->_id = id; }

  /** Get id.
   * \return an unsigned int
   */
  unsigned int id() { return this->_id; }

  /** Set shape position and orientation.
      \param q : NewtonEulerDS state
  */
  virtual void move(const SiconosVector& q);


  /** Distance to a general contact shape.
      \param sh2 : the contact shape.
      \param normalFromFace1 : normal on first contact shape, default on second.
      \return the distance, contact points and normal in ContactShapeDistance
   */
  virtual SP::ContactShapeDistance distance(
    const OccContactShape& sh2, bool normalFromFace1=false) const;

  /** Distance to a contact face.
      \param sh2 : the contact face.
      \param normalFromFace1 : normal on first contact shape, default on second.
      \return the distance, contact points and normal in ContactShapeDistance
   */
  virtual SP::ContactShapeDistance distance(
    const OccContactFace& sh2, bool normalFromFace1=false) const;

  /** Distance to a contact edge.
      \param sh2 : the contact edge.
      \param normalFromFace1 : normal on first contact shape, default on second.
      \return the distance, contact points and normal in ContactShapeDistance
   */
  virtual SP::ContactShapeDistance distance(
    const OccContactEdge& sh2, bool normalFromFace1=false) const;

  /** Computed UV bounds.
   * @{
   */
  double bsup1[2];
  double binf1[2];
  double bsup2[2];
  double binf2[2];
  /**
   * @}
   */

  /** Contact group. */
  unsigned int contactGroup;

  SP::TopoDS_Shape _shape;

  unsigned int _id;

  /** visitors hook
   */
  VIRTUAL_ACCEPT_VISITORS(OccContactShape);

};

#endif
