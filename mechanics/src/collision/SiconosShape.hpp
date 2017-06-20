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
  unsigned int _version; // version number tracks changes to shape properties

  SiconosShape()
    : _inside_margin(0.1)
    , _outside_margin(0.0)
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

  unsigned int version() const { return _version; }

  /** visitors hook
   */
  VIRTUAL_ACCEPT_VISITORS();
};

class SiconosPlane : public SiconosShape,
                     public std11::enable_shared_from_this<SiconosPlane>
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosPlane);

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
private:
  SiconosSphere() : SiconosShape() {};

protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosSphere);
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
private:
  SiconosBox() : SiconosShape() {};

protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosBox);
  SP::SiconosVector _dimensions;

public:
  SiconosBox(double width, double height, double depth)
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

  void setDimensions(double width, double height, double depth)
  {
    (*_dimensions)(0) = width;
    (*_dimensions)(1) = height;
    (*_dimensions)(2) = depth;
    _version ++;
  }

  void setDimensions(SP::SiconosVector dim)
  {
    _dimensions = dim;
    _version ++;
  }

  void setDimensions(const SiconosVector& dim)
  {
    (*_dimensions)(0) = dim(0);
    (*_dimensions)(1) = dim(1);
    (*_dimensions)(2) = dim(2);
    _version ++;
  }

  /** visitors hook
   */
  ACCEPT_VISITORS();
};

class SiconosCylinder : public SiconosShape,
                        public std11::enable_shared_from_this<SiconosCylinder>
{
private:
  SiconosCylinder() : SiconosShape() {};

protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosCylinder);
  double _radius;
  double _length;

public:
  SiconosCylinder(float radius, float length)
    : SiconosShape(), _radius(radius), _length(length)
  {
  }

  virtual ~SiconosCylinder() {}

  void setRadius(double radius)
  {
    _radius = radius;
    _version ++;
  }

  double radius() { return _radius; }

  void setLength(double length)
  {
    _length = length;
    _version ++;
  }

  double length() { return _length; }

  /** visitors hook
   */
  ACCEPT_VISITORS();
};

class SiconosConvexHull : public SiconosShape,
                          public std11::enable_shared_from_this<SiconosConvexHull>
{
private:
  SiconosConvexHull() : SiconosShape() {};

protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosConvexHull);
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

typedef std::vector<unsigned int> VUInt;
TYPEDEF_SPTR(VUInt)

class SiconosMesh : public SiconosShape,
                    public std11::enable_shared_from_this<SiconosMesh>
{
private:
  SiconosMesh() : SiconosShape() {};

protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosMesh);
  SP::VUInt _indexes;
  SP::SiconosMatrix _vertices;

public:
  SiconosMesh(SP::VUInt indexes,
              SP::SiconosMatrix vertices)
    : SiconosShape(), _indexes(indexes), _vertices(vertices)
  {
    if (!_indexes || (_indexes->size() % 3) != 0)
      throw SiconosException("Mesh indexes size must be divisible by 3.");
    if (!_vertices || _vertices->size(1) != 3)
      throw SiconosException("Mesh vertices matrix must have 3 columns.");
  }

  SP::VUInt indexes() { return _indexes; }
  SP::SiconosMatrix vertices() { return _vertices; }

  virtual ~SiconosMesh() {}

  /** visitors hook
   */
  ACCEPT_VISITORS();
};

class SiconosHeightMap : public SiconosShape,
                         public std11::enable_shared_from_this<SiconosHeightMap>
{
private:
  SiconosHeightMap() : SiconosShape() {};

protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosHeightMap);
  SP::SiconosMatrix _height_data;
  double _length_x;
  double _length_y;

public:
  SiconosHeightMap(SP::SiconosMatrix height_data,
                   double length_x, double length_y)
    : SiconosShape(), _height_data(height_data),
      _length_x(length_x), _length_y(length_y)
  {
  }

  SP::SiconosMatrix height_data() { return _height_data; }
  double length_x() { return _length_x; }
  double length_y() { return _length_y; }

  virtual ~SiconosHeightMap() {}

  /** visitors hook
   */
  ACCEPT_VISITORS();
};

#endif /* SiconosShape_h */
