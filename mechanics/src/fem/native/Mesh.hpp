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

/*! \file Mesh.hpp

 */
#ifndef MESH_H
#define MESH_H

// #include <iostream>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "FETypes.hpp"

namespace siconos::mechanics::fem {

class MeshElement;

/** a Mesh vertex */
class MeshVertex {
 private:
  /** Vertex Number */
  size_t num_ = 0;

  /* Vextex Coordinate */
  double _x = 0.;
  double _y = 0.;
  double _z = 0.;

  /** elements to which the node belongs **/
  std::vector<std::shared_ptr<MeshElement>> elements_ = {};

  // Rule of five
  MeshVertex() = delete;
  MeshVertex(MeshVertex&) = delete;
  MeshVertex& operator=(const MeshVertex&) = delete;
  MeshVertex(MeshVertex&&) = delete;
  MeshVertex& operator=(MeshVertex&&) = delete;

 public:
  /* Constructor from data */
  MeshVertex(size_t num, double x, double y, double z) : num_{num}, _x{x}, _y{y}, _z{z} {};

  virtual ~MeshVertex() noexcept = default;

  double x() const { return _x; };
  double y() const { return _y; };
  double z() const { return _z; };

  size_t num() const { return num_; }

  void link_element(std::shared_ptr<MeshElement> elem) { elements_.push_back(elem); }
  const std::span<const std::shared_ptr<MeshElement>> elements() const { return elements_; }

  virtual void display() const;
};

class MeshBeamVertex : public MeshVertex {
 private:
  /* Vextex Rotation Coordinate */
  double _qx = 0.;
  double _qy = 0.;
  double _qz = 0.;

 public:
  /* Constructor from data */
  MeshBeamVertex(size_t num, double x, double y, double z, double qx, double qy, double qz)
      : MeshVertex(num, x, y, z) {
    _qx = qx;
    _qy = qy;
    _qz = qz;
  };

  virtual ~MeshBeamVertex() noexcept = default;

  auto qx() { return _qx; };
  auto qy() { return _qy; };
  auto qz() { return _qz; };
  void display() const override;
};

/** a mesh element */
class MeshElement {
 protected:
  /** Element number */
  size_t num_{0};

  /** type (following gmsh convention), default = 2, 3-node triangle */
  FiniteElementType type_{FiniteElementType::T3};

  /** vertices **/
  std::vector<std::shared_ptr<MeshVertex>> vertices_ = {};

  /** tags **/
  std::vector<int> tags_ = {};

  /** Rule of five */
  MeshElement() = delete;
  MeshElement(MeshElement&) = delete;
  MeshElement& operator=(const MeshElement&) = delete;
  MeshElement(MeshElement&&) = delete;
  MeshElement& operator=(MeshElement&&) = delete;

 public:
  MeshElement(size_t num) {};
  MeshElement(size_t num, FiniteElementType type,
              std::vector<std::shared_ptr<MeshVertex>> vertices, std::vector<int> tags = {})
      : num_{num}, type_{type}, vertices_{std::move(vertices)}, tags_{std::move(tags)} {};

  //   MeshElement(size_t num, FiniteElementType type,
  //               std::vector<std::shared_ptr<MeshVertex>> vertices)
  //       : MeshElement(num, type, vertices, {0}) {};

  ~MeshElement() noexcept = default;

  FiniteElementType type() const { return type_; }
  size_t num() const { return num_; }
  const std::span<const std::shared_ptr<MeshVertex>> vertices() const { return vertices_; }
  auto tags(size_t n) const { return tags_[n]; }
  void display() const;
};

/** a Mesh container */
class Mesh {
 protected:
  /** Space dimension */
  int dimension_ = 2;

  /** vertices */
  std::vector<std::shared_ptr<MeshVertex>> vertices_ = {};

  /** elements */
  std::vector<std::shared_ptr<MeshElement>> elements_ = {};

  /** Physical entities
      Connection between the tag and the name of the corresponding physical entity
      (like Dirichlet BC, Applied force ...)
   */
  std::vector<std::tuple<int, std::string>> physical_entities_ = {};

  /** Rule of five */
  Mesh() = delete;
  Mesh(Mesh&) = delete;
  Mesh& operator=(const Mesh&) = delete;
  Mesh(Mesh&&) = delete;
  Mesh& operator=(Mesh&&) = delete;

 public:
  /** constructor
      \param dim dimension
      \param vertices a vector of vertices
      \param elements a vector of elements
      \param physical_entities connection between tags and entities
  */
  Mesh(int dim, std::vector<std::shared_ptr<MeshVertex>> vertices,
       std::vector<std::shared_ptr<MeshElement>> elements,
       std::vector<std::tuple<int, std::string>> physical_entities);

  /** constructor
      \param dim dimension
      \param vertices a vector of vertices
      \param elements a vector of elements
   */
  Mesh(int dim, std::vector<std::shared_ptr<MeshVertex>> vertices,
       std::vector<std::shared_ptr<MeshElement>> elements)
      : Mesh(dim, vertices, elements, {}) {};
  ;

  /** destructor */
  ~Mesh() noexcept = default;

  int dim() const { return dimension_; };

  const std::span<const std::shared_ptr<MeshVertex>> vertices() const { return vertices_; }

  const std::span<const std::shared_ptr<MeshElement>> elements() const { return elements_; }

  std::vector<std::tuple<int, std::string>> physical_entities() { return physical_entities_; }

  /** print the data of the Mesh

      \param brief true to toggle verbose mode off, default = true
   */

  void display(bool brief = true) const;
};

}  // namespace siconos::mechanics::fem
#endif  // MESH_H
