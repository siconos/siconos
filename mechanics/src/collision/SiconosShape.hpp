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
  \brief Definition of rigid shapes (derivated from top level abstract class SiconosShape)
  derivated
*/

#ifndef SiconosShape_h
#define SiconosShape_h

#include <SiconosSerialization.hpp>
#include <memory>
#include <vector>

namespace siconos::algebra {
class SiconosVector;
class SiconosMatrix;
}  // namespace siconos::algebra

namespace siconos::collision {

namespace internal {
class ShapeVisitor;
}

class SiconosShape {
 protected:
  ACCEPT_SERIALIZATION(SiconosShape);

  double _inside_margin{0.1};
  double _outside_margin{0.};
  unsigned int _version{0};
  // version number tracks changes to shape properties

  // Rule of five
  SiconosShape() = default;
  SiconosShape(const SiconosShape&) = delete;
  SiconosShape(SiconosShape&&) = delete;
  SiconosShape& operator=(const SiconosShape&) = delete;
  SiconosShape& operator=(SiconosShape&&) = delete;

 public:
  virtual ~SiconosShape() noexcept = default;

  /** Set the inside margin of the shape.  This is a distance that the
   *  contour should be shrunk to improve contact detection robustness.
   *  It will have an effect on the roundness of corners. */
  void setInsideMargin(double margin) {
    _inside_margin = margin;
    _version++;
  }

  /** Set the outside margin of the shape.  This is the distance from
   *  the contact shell to an external shell used to detect contacts
   *  in advance.  The implementation will detect contact points on
   *  the external shell and project them back to the contact shell.
   *  Note: Currently not working in Bullet implementation!  Better to
   *  leave at zero. */
  void setOutsideMargin(double margin) {
    _outside_margin = margin;
    _version++;
  }

  double insideMargin() { return _inside_margin; }

  double outsideMargin() { return _outside_margin; }

  unsigned int version() const { return _version; }

  virtual void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) = 0;
};

class SiconosPlane : public SiconosShape, public std::enable_shared_from_this<SiconosPlane> {
 protected:
  ACCEPT_SERIALIZATION(SiconosPlane);

 public:
  SiconosPlane() = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;

  virtual ~SiconosPlane() noexcept = default;
};

class SiconosSphere : public SiconosShape, public std::enable_shared_from_this<SiconosSphere> {
 private:
  SiconosSphere() = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosSphere);
  float _radius{0.};

 public:
  SiconosSphere(float radius) : SiconosShape(), _radius(radius) {}

  virtual ~SiconosSphere() noexcept = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;

  float radius() const { return _radius; }
  void setRadius(float r) {
    _radius = r;
    _version++;
  }
};

class SiconosBox : public SiconosShape, public std::enable_shared_from_this<SiconosBox> {
 private:
  SiconosBox() = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosBox);
  std::shared_ptr<siconos::algebra::SiconosVector> _dimensions{nullptr};

 public:
  SiconosBox(double width, double height, double depth);

  SiconosBox(std::shared_ptr<siconos::algebra::SiconosVector> dimensions)
      : SiconosShape(), _dimensions(dimensions) {}

  virtual ~SiconosBox() noexcept = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;

  std::shared_ptr<siconos::algebra::SiconosVector> dimensions() const { return _dimensions; }

  void setDimensions(double width, double height, double depth);

  void setDimensions(std::shared_ptr<siconos::algebra::SiconosVector> dim) {
    _dimensions = dim;
    _version++;
  }

  void setDimensions(const siconos::algebra::SiconosVector& dim);
};

class SiconosCylinder : public SiconosShape,
                        public std::enable_shared_from_this<SiconosCylinder> {
 private:
  SiconosCylinder() = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosCylinder);
  double _radius{0.};
  double _length{0.};

 public:
  SiconosCylinder(float radius, float length)
      : SiconosShape(), _radius(radius), _length(length) {}

  virtual ~SiconosCylinder() noexcept = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;

  void setRadius(double radius) {
    _radius = radius;
    _version++;
  }

  double radius() { return _radius; }

  void setLength(double length) {
    _length = length;
    _version++;
  }

  double length() { return _length; }
};

class SiconosCone : public SiconosShape, public std::enable_shared_from_this<SiconosCone> {
 private:
  SiconosCone() = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosCone);
  double _radius{0.};
  double _length{0.};

 public:
  SiconosCone(float radius, float length) : SiconosShape(), _radius(radius), _length(length) {}

  virtual ~SiconosCone() noexcept = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;

  void setRadius(double radius) {
    _radius = radius;
    _version++;
  }

  double radius() { return _radius; }

  void setLength(double length) {
    _length = length;
    _version++;
  }

  double length() { return _length; }
};

class SiconosCapsule : public SiconosShape,
                       public std::enable_shared_from_this<SiconosCapsule> {
 private:
  SiconosCapsule() = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosCapsule);
  double _radius{0.};
  double _length{0.};

 public:
  SiconosCapsule(float radius, float length)
      : SiconosShape(), _radius(radius), _length(length) {}

  virtual ~SiconosCapsule() noexcept = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;

  void setRadius(double radius) {
    _radius = radius;
    _version++;
  }

  double radius() { return _radius; }

  void setLength(double length) {
    _length = length;
    _version++;
  }

  double length() { return _length; }
};

class SiconosConvexHull : public SiconosShape,
                          public std::enable_shared_from_this<SiconosConvexHull> {
 private:
  SiconosConvexHull() = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosConvexHull);
  std::shared_ptr<siconos::algebra::SiconosMatrix> _vertices{nullptr};

 public:
  SiconosConvexHull(std::shared_ptr<siconos::algebra::SiconosMatrix> vertices);

  virtual ~SiconosConvexHull() noexcept = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;

  std::shared_ptr<siconos::algebra::SiconosMatrix> vertices() const { return _vertices; }

  void setVertices(std::shared_ptr<siconos::algebra::SiconosMatrix> vertices) {
    _vertices = vertices;
    _version++;
  }
};

class SiconosMesh : public SiconosShape, public std::enable_shared_from_this<SiconosMesh> {
 private:
  SiconosMesh() = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosMesh);
  std::shared_ptr<std::vector<unsigned int>> _indexes{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _vertices{nullptr};

 public:
  SiconosMesh(std::shared_ptr<std::vector<unsigned int>> indexes,
              std::shared_ptr<siconos::algebra::SiconosMatrix> vertices);

  std::shared_ptr<std::vector<unsigned int>> indexes() { return _indexes; }
  std::shared_ptr<siconos::algebra::SiconosMatrix> vertices() { return _vertices; }

  virtual ~SiconosMesh() noexcept = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;
};

class SiconosHeightMap : public SiconosShape,
                         public std::enable_shared_from_this<SiconosHeightMap> {
 private:
  SiconosHeightMap() = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosHeightMap);
  std::shared_ptr<siconos::algebra::SiconosMatrix> _height_data{nullptr};
  double _length_x{0.};
  double _length_y{0.};

 public:
  SiconosHeightMap(std::shared_ptr<siconos::algebra::SiconosMatrix> height_data,
                   double length_x, double length_y)
      : SiconosShape(), _height_data(height_data), _length_x(length_x), _length_y(length_y) {}

  std::shared_ptr<siconos::algebra::SiconosMatrix> height_data() { return _height_data; }
  double length_x() { return _length_x; }
  double length_y() { return _length_y; }
  //
  virtual ~SiconosHeightMap() noexcept = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;
};

class SiconosDisk : public SiconosShape, public std::enable_shared_from_this<SiconosDisk> {
 private:
  SiconosDisk() = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosDisk);
  float _radius{0.};

 public:
  SiconosDisk(float radius) : SiconosShape(), _radius(radius) {}

  virtual ~SiconosDisk() noexcept = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;

  float radius() const { return _radius; }
  void setRadius(float r) {
    _radius = r;
    _version++;
  }
};

class SiconosBox2d : public SiconosShape, public std::enable_shared_from_this<SiconosBox2d> {
 private:
  SiconosBox2d() = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosBox2d);
  std::shared_ptr<siconos::algebra::SiconosVector> _dimensions{nullptr};

 public:
  SiconosBox2d(double width, double height);

  SiconosBox2d(std::shared_ptr<siconos::algebra::SiconosVector> dimensions)
      : SiconosShape(), _dimensions(dimensions) {}

  virtual ~SiconosBox2d() noexcept = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;

  std::shared_ptr<siconos::algebra::SiconosVector> dimensions() const { return _dimensions; }

  void setDimensions(double width, double height);

  void setDimensions(std::shared_ptr<siconos::algebra::SiconosVector> dim) {
    _dimensions = dim;
    _version++;
  }

  void setDimensions(const siconos::algebra::SiconosVector& dim);
};

class SiconosConvexHull2d : public SiconosShape,
                            public std::enable_shared_from_this<SiconosConvexHull2d> {
 private:
  SiconosConvexHull2d() = delete;
  ;

 protected:
  ACCEPT_SERIALIZATION(SiconosConvexHull2d);
  std::shared_ptr<siconos::algebra::SiconosMatrix> _vertices{nullptr};

  /* boolean to use the normal to the selected edge of the convexhull
     to avoid contact with vertex */
  bool _avoidInternalEdgeContact{false};

 public:
  /* index of the first point of the selected edge to compute the normal edge (default=0) */
  int _normal_edge_pointA{0};
  /* index of the first point of the selected edge to compute the normal edge (default=1) */
  int _normal_edge_pointB{1};

  SiconosConvexHull2d(std::shared_ptr<siconos::algebra::SiconosMatrix> vertices);

  virtual ~SiconosConvexHull2d() noexcept = default;
  void accept(std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) override;

  std::shared_ptr<siconos::algebra::SiconosMatrix> vertices() const { return _vertices; }

  void setVertices(std::shared_ptr<siconos::algebra::SiconosMatrix> vertices) {
    _vertices = vertices;
    _version++;
  }
  bool avoidInternalEdgeContact() const { return _avoidInternalEdgeContact; }

  void setAvoidInternalEdgeContact(bool value) { _avoidInternalEdgeContact = value; }
};
}  // namespace siconos::collision
#endif /* SiconosShape_h */
