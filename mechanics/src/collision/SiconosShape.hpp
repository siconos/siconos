/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

  /** Set the inside margin of the shape.  This is a distance that the
   *  contour should be shrunk to improve contact detection robustness.
   *  It will have an effect on the roundness of corners. */
  void setInsideMargin (double margin)
  {
    _inside_margin = margin;
    _version ++;
  }

  /** Set the outside margin of the shape.  This is the distance from
   *  the contact shell to an external shell used to detect contacts
   *  in advance.  The implementation will detect contact points on
   *  the external shell and project them back to the contact shell.
   *  Note: Currently not working in Bullet implementation!  Better to
   *  leave at zero. */
  void setOutsideMargin(double margin)
  {
    _outside_margin = margin;
    _version ++;
  }

  double insideMargin() { return _inside_margin; }

  double outsideMargin() { return _outside_margin; }

  unsigned int version() const { return _version; }

  VIRTUAL_ACCEPT_VISITORS();
};

class SiconosPlane : public SiconosShape,
                     public std::enable_shared_from_this<SiconosPlane>
{
protected:
  ACCEPT_SERIALIZATION(SiconosPlane);

public:
  SiconosPlane() : SiconosShape() {}

  virtual ~SiconosPlane() {}

  ACCEPT_VISITORS();
};

class SiconosSphere : public SiconosShape,
                      public std::enable_shared_from_this<SiconosSphere>
{
private:
  SiconosSphere() : SiconosShape() {};

protected:
  ACCEPT_SERIALIZATION(SiconosSphere);
  float _radius;

public:
  SiconosSphere(float radius)
    : SiconosShape(), _radius(radius) {}

  virtual ~SiconosSphere() {}

  float radius() const { return _radius; }
  void setRadius(float r) { _radius = r; _version ++; }

  ACCEPT_VISITORS();
};

class SiconosBox : public SiconosShape,
                   public std::enable_shared_from_this<SiconosBox>
{
private:
  SiconosBox() : SiconosShape() {};

protected:
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

  ACCEPT_VISITORS();
};

class SiconosCylinder : public SiconosShape,
                        public std::enable_shared_from_this<SiconosCylinder>
{
private:
  SiconosCylinder() : SiconosShape() {};

protected:

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

  ACCEPT_VISITORS();
};

class SiconosCone : public SiconosShape,
                        public std::enable_shared_from_this<SiconosCone>
{
private:
  SiconosCone() : SiconosShape() {};

protected:
  ACCEPT_SERIALIZATION(SiconosCone);
  double _radius;
  double _length;

public:
  SiconosCone(float radius, float length)
    : SiconosShape(), _radius(radius), _length(length)
  {
  }

  virtual ~SiconosCone() {}

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

  ACCEPT_VISITORS();
};

class SiconosCapsule : public SiconosShape,
                        public std::enable_shared_from_this<SiconosCapsule>
{
private:
  SiconosCapsule() : SiconosShape() {};

protected:

  ACCEPT_SERIALIZATION(SiconosCapsule);
  double _radius;
  double _length;

public:
  SiconosCapsule(float radius, float length)
    : SiconosShape(), _radius(radius), _length(length)
  {
  }

  virtual ~SiconosCapsule() {}

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

  ACCEPT_VISITORS();
};



class SiconosConvexHull : public SiconosShape,
                          public std::enable_shared_from_this<SiconosConvexHull>
{
private:
  SiconosConvexHull() : SiconosShape() {};

protected:

  ACCEPT_SERIALIZATION(SiconosConvexHull);
  SP::SiconosMatrix _vertices;

public:
  SiconosConvexHull(SP::SiconosMatrix vertices)
    : SiconosShape(), _vertices(vertices)
  {
    if (_vertices && _vertices->size(1) != 3)
      THROW_EXCEPTION("Convex hull vertices matrix must have 3 columns.");
  }

  virtual ~SiconosConvexHull() {}

  SP::SiconosMatrix vertices() const { return _vertices; }

  void setVertices(SP::SiconosMatrix vertices)
  {
    _vertices = vertices;
    _version ++;
  }

  ACCEPT_VISITORS();
};

typedef std::vector<unsigned int> VUInt;
TYPEDEF_SPTR(VUInt)

class SiconosMesh : public SiconosShape,
                    public std::enable_shared_from_this<SiconosMesh>
{
private:
  SiconosMesh() : SiconosShape() {};

protected:

  ACCEPT_SERIALIZATION(SiconosMesh);
  SP::VUInt _indexes;
  SP::SiconosMatrix _vertices;

public:
  SiconosMesh(SP::VUInt indexes,
              SP::SiconosMatrix vertices)
    : SiconosShape(), _indexes(indexes), _vertices(vertices)
  {
    if (!_indexes || (_indexes->size() % 3) != 0)
      THROW_EXCEPTION("Mesh indexes size must be divisible by 3.");
    if (!_vertices || _vertices->size(0) != 3)
      THROW_EXCEPTION("Mesh vertices matrix must have 3 columns.");
  }

  SP::VUInt indexes() { return _indexes; }
  SP::SiconosMatrix vertices() { return _vertices; }

  virtual ~SiconosMesh() {}

  ACCEPT_VISITORS();
};

class SiconosHeightMap : public SiconosShape,
                         public std::enable_shared_from_this<SiconosHeightMap>
{
private:
  SiconosHeightMap() : SiconosShape() {};

protected:

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

  ACCEPT_VISITORS();
};


class SiconosDisk : public SiconosShape,
                    public std::enable_shared_from_this<SiconosDisk>
{
private:
  SiconosDisk() : SiconosShape() {};

protected:
  ACCEPT_SERIALIZATION(SiconosDisk);
  float _radius;

public:
  SiconosDisk(float radius)
    : SiconosShape(), _radius(radius) {}

  virtual ~SiconosDisk() {}

  float radius() const { return _radius; }
  void setRadius(float r) { _radius = r; _version ++; }

  ACCEPT_VISITORS();
};

class SiconosBox2d : public SiconosShape,
                   public std::enable_shared_from_this<SiconosBox2d>
{
private:
  SiconosBox2d() : SiconosShape() {};

protected:

  ACCEPT_SERIALIZATION(SiconosBox2d);
  SP::SiconosVector _dimensions;

public:
  SiconosBox2d(double width, double height)
    : SiconosShape(), _dimensions(new SiconosVector(2))
  {
    (*_dimensions)(0) = width;
    (*_dimensions)(1) = height;
  }

  SiconosBox2d(SP::SiconosVector dimensions)
    : SiconosShape(), _dimensions(dimensions) {}

  virtual ~SiconosBox2d() {}

  SP::SiconosVector dimensions() const { return _dimensions; }

  void setDimensions(double width, double height)
  {
    (*_dimensions)(0) = width;
    (*_dimensions)(1) = height;
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
    _version ++;
  }

  ACCEPT_VISITORS();
};

class SiconosConvexHull2d : public SiconosShape,
                          public std::enable_shared_from_this<SiconosConvexHull2d>
{
private:
  SiconosConvexHull2d() : SiconosShape() {};

protected:

  ACCEPT_SERIALIZATION(SiconosConvexHull2d);
  SP::SiconosMatrix _vertices;

  /* boolean to use the normal to the selected edge of the convexhull
     to avoid contact with vertex */
  bool _avoidInternalEdgeContact;


public:
  /* index of the first point of the selected edge to compute the normal edge (default=0) */
  int _normal_edge_pointA;
  /* index of the first point of the selected edge to compute the normal edge (default=1) */
  int _normal_edge_pointB;


  SiconosConvexHull2d(SP::SiconosMatrix vertices)
    : SiconosShape(), _vertices(vertices), _avoidInternalEdgeContact(false), _normal_edge_pointA(0), _normal_edge_pointB(1)
  {
    if (_vertices && _vertices->size(1) != 2)
      THROW_EXCEPTION("Convex hull vertices matrix must have 2 columns in 2d.");
  }

  virtual ~SiconosConvexHull2d() {}

  SP::SiconosMatrix vertices() const { return _vertices; }

  void setVertices(SP::SiconosMatrix vertices)
  {
    _vertices = vertices;
    _version ++;
  }
  bool avoidInternalEdgeContact() const {return _avoidInternalEdgeContact;}

  void setAvoidInternalEdgeContact(bool value)
  {
    _avoidInternalEdgeContact = value;
  }

  ACCEPT_VISITORS();
};

#endif /* SiconosShape_h */
