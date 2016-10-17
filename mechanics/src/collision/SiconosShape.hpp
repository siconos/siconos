/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

/*! \file SiconosShape.hpp
  \brief Definition of an abstract rigid shape
*/


#ifndef SiconosShape_h
#define SiconosShape_h

#include "MechanicsFwd.hpp"
#include <SiconosVisitor.hpp>
#include <SiconosSerialization.hpp>
#include <SiconosVector.hpp>
#include <SiconosMatrix.hpp>

class SiconosShape
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosShape);

  double _inside_margin;
  double _outside_margin;
  unsigned int _group;   // the collision group to identify non-smooth law
  unsigned int _version; // version number tracks changes to shape properties

  SiconosShape()
    : _inside_margin(0.1)
    , _outside_margin(0.1)
    , _group(0)
    , _version(0)
    {}

public:

  virtual ~SiconosShape() {}

  void setInsideMargin (double margin)
  {
    _inside_margin = margin;
    _version ++;
  }

  void setOutsideMargin(double margin)
  {
    _outside_margin = margin;
    _version ++;
  }

  double insideMargin() { return _inside_margin; }

  double outsideMargin() { return _outside_margin; }

  void setGroup(int group) { _group = group; }

  unsigned int group() { return _group; }

  unsigned int version() const { return _version; }

  /** visitors hook
   */
  VIRTUAL_ACCEPT_VISITORS();
};

class SiconosPlane : public SiconosShape,
                     public std11::enable_shared_from_this<SiconosPlane>
{
public:
  SiconosPlane() : SiconosShape() {}

  virtual ~SiconosPlane() {}

  /** visitors hook
   */
  ACCEPT_VISITORS();
};

class SiconosSphere : public SiconosShape,
                      public std11::enable_shared_from_this<SiconosSphere>
{
protected:
  float _radius;

public:
  SiconosSphere(float radius)
    : SiconosShape(), _radius(radius) {}

  virtual ~SiconosSphere() {}

  float radius() const { return _radius; }
  void setRadius(float r) { _radius = r; _version ++; }

  /** visitors hook
   */
  ACCEPT_VISITORS();
};

class SiconosBox : public SiconosShape,
                   public std11::enable_shared_from_this<SiconosBox>
{
protected:
  SP::SiconosVector _dimensions;

public:
  SiconosBox(float width, float height, float depth)
    : SiconosShape(), _dimensions(new SiconosVector(3))
  {
    (*_dimensions)(0) = width;
    (*_dimensions)(1) = height;
    (*_dimensions)(2) = depth;
  }

  SiconosBox(SP::SiconosVector dimensions)
    : SiconosShape(), _dimensions(dimensions) {}

  virtual ~SiconosBox() {}

  SP::SiconosVector dimensions() const { return _dimensions; }

  void setDimensions(SP::SiconosVector dim)
  {
    (*_dimensions)(0) = (*dim)(0);
    (*_dimensions)(1) = (*dim)(1);
    (*_dimensions)(2) = (*dim)(2);
    _version ++;
  }

  /** visitors hook
   */
  ACCEPT_VISITORS();
};

class SiconosConvexHull : public SiconosShape,
                          public std11::enable_shared_from_this<SiconosConvexHull>
{
protected:
  SP::SiconosMatrix _vertices;

public:
  SiconosConvexHull(SP::SiconosMatrix vertices)
    : SiconosShape(), _vertices(vertices)
  {
    if (_vertices && _vertices->size(1) != 3)
      throw SiconosException("Convex hull vertices matrix must have 3 columns.");
  }

  virtual ~SiconosConvexHull() {}

  SP::SiconosMatrix vertices() const { return _vertices; }

  void setVertices(SP::SiconosMatrix vertices)
  {
    _vertices = vertices;
    _version ++;
  }

  /** visitors hook
   */
  ACCEPT_VISITORS();
};

#endif /* SiconosShape_h */
