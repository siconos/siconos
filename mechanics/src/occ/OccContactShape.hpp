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
      \param shape : an OccContactShape
   */
  OccContactShape(const OccContactShape& shape)
    : _shape(shape.shape()) {};

  /** Destructor.
   */
  virtual ~OccContactShape() {};

  /** Return shared pointer on OpenCascade data
   */
  SP::TopoDS_Shape shape() const { return this->_shape;};

  /** Return reference on OpenCascade data.
   */
  TopoDS_Shape& data() const { return *this->_shape;};

  /** Set OpenCascade data from a shared pointer.
   */
  void setShape(SP::TopoDS_Shape shape)
  {
    this->_shape = shape;
  }

  /** Set OpenCascade data.
   */
  void setData(TopoDS_Shape& data)
  {
    this->_shape = createSPtrTopoDS_Shape(data);
  }

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
