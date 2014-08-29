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

#include <string>
#include <SiconosFwd.hpp>
#include <MechanicsFwd.hpp>
DEFINE_SPTR(TopoDS_Shape);

struct OccContactShape
{

  /** Default constructor.
   */
  OccContactShape() {};

  /** Constructor from OpenCascade object.
      \param shape : FACE or EDGE shape
   */
  OccContactShape(TopoDS_Shape& shape)
    : _shape(createSPtrTopoDS_Shape(shape))
  {};

  /** Constructor from const OpenCascade object : remove constness.
      \param shape : FACE or EDGE shape
   */
  OccContactShape(const TopoDS_Shape& shape)
    : _shape(createSPtrTopoDS_Shape(const_cast<TopoDS_Shape&>(shape)))
  {};

  /** Return OpenCascade data with default initialization if needed.
   */
  TopoDS_Shape& data();

  /** Return OpenCascade data.
   */
  TopoDS_Shape& data() const {return *_shape;};

  /** Known contacts.
   */
  enum ContactTypeValue { Face, Edge, Unknown };

  /** Get the contact type.
   */
  ContactTypeValue contactType() const;

  /** Export the contact shape into a string.
   */
  std::string exportBRepAsString() const;

  /** import the contact shape from a string.
   *  \param brepstr : the string containing the whole brep.
   */
  void importBRepFromString(const std::string& brepstr);

  /** Compute and store UV bounds of the shape.
   */
  void computeUVBounds();

  /** Set shape position and orientation
      \param q : NewtonEulerDS state
  */
  void move(const SiconosVector& q);


  /** Distance to another contact shape.
      \param sh2 : the other contact shape.
   */
  void distance(
    const OccContactShape& sh2,
    double& X1, double& Y1, double& Z1,
    double& X2, double& Y2, double& Z2,
    double& nX, double& nY, double& nZ,
    bool normalFromFace1,double& MinDist) const;

  /** computed UV bounds
   * @{
   */
  double bsup1[2];
  double binf1[2];
  double bsup2[2];
  double binf2[2];
  /**
   * @}
   */

  /** contact group */
  unsigned int contactGroup;

  SP::TopoDS_Shape _shape;

};

#endif
