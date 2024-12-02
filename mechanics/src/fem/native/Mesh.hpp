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
#include <string>
#include <vector>

#include "FETypes.hpp"

namespace siconos::mechanics::fem {

class MElement;

/** a Mesh vertex */
class MVertex {
 private:
  /** Vertex Number */
  size_t _num = 0;

  /* Vextex Coordinate */
  double _x = 0.;
  double _y = 0.;
  double _z = 0.;

  /** elements to which the node belongs **/
  std::vector<std::shared_ptr<MElement>> _elements = {};

  // Rule of five
  MVertex() = delete;
  MVertex(MVertex &) = delete;
  MVertex &operator=(const MVertex &) = delete;
  MVertex(MVertex &&) = delete;
  MVertex &operator=(MVertex &&) = delete;

 public:
  /* Constructor from data */
  MVertex(size_t num, double x, double y, double z) : _num{num}, _x{x}, _y{y}, _z{z} {};

  ~MVertex() noexcept = default;

  auto x() { return _x; };
  auto y() { return _y; };
  auto z() { return _z; };

  auto num() { return _num; }

  auto &elements() { return _elements; };

  void display();
};

/** a mesh element */
class MElement {
 protected:
  /** Element number */
  size_t _num{0};

  /** type (following gmsh convention), default = 2, 3-node triangle */
  FiniteElementType _type{FiniteElementType::T3};

  /** vertices **/
  std::vector<std::shared_ptr<MVertex>> _vertices = {};

  /** tags **/
  std::vector<int> _tags = {};

  /** Rule of five */
  MElement() = delete;
  MElement(MElement &) = delete;
  MElement &operator=(const MElement &) = delete;
  MElement(MElement &&) = delete;
  MElement &operator=(MElement &&) = delete;

 public:
  MElement(size_t num, FiniteElementType type, std::vector<std::shared_ptr<MVertex>> &vertices,
           std::vector<int> tags)
      : _num{num}, _type{type}, _vertices{vertices}, _tags{tags} {};

  MElement(size_t num, FiniteElementType type, std::vector<std::shared_ptr<MVertex>> vertices)
      : MElement(num, type, vertices, {0}) {};

  ~MElement() noexcept = default;

  auto type() { return _type; }
  auto num() { return _num; }
  auto &vertices() { return _vertices; }
  auto tags(int n) { return _tags[n]; }
  void display();
};

/** a Mesh container */
class Mesh {
 protected:
  /** Space dimension */
  int _dim = 2;

  /** vertices */
  std::vector<std::shared_ptr<MVertex>> _vertices = {};

  /** elements */
  std::vector<std::shared_ptr<MElement>> _elements = {};

  /** Physical entities
      Connection between the tag and the name of the corresponding physical entity
      (like Dirichlet BC, Applied force ...)
   */
  std::vector<std::tuple<int, std::string>> _physical_entities = {};

  /** Rule of five */
  Mesh() = delete;
  Mesh(Mesh &) = delete;
  Mesh &operator=(const Mesh &) = delete;
  Mesh(Mesh &&) = delete;
  Mesh &operator=(Mesh &&) = delete;

 public:
  /** constructor
      \param dim dimension
      \param vertices a vector of vertices
      \param elements a vector of elements
      \param physical_entities connection between tags and entities
  */
  Mesh(int dim, std::vector<std::shared_ptr<MVertex>> vertices,
       std::vector<std::shared_ptr<MElement>> elements,
       std::vector<std::tuple<int, std::string>> physical_entities);

  /** constructor
      \param dim dimension
      \param vertices a vector of vertices
      \param elements a vector of elements
   */
  Mesh(int dim, std::vector<std::shared_ptr<MVertex>> vertices,
       std::vector<std::shared_ptr<MElement>> elements)
      : Mesh(dim, vertices, elements, {}) {};
  ;

  /** destructor */
  ~Mesh() noexcept = default;

  auto dim() { return _dim; };

  auto &vertices() { return _vertices; }

  auto &elements() { return _elements; }

  auto physical_entities() { return _physical_entities; }

  /** print the data of the Mesh

      \param brief true to toggle verbose mode off, default = true
   */

  void display(bool brief = true) const;
};

}  // namespace siconos::mechanics::fem
#endif  // MESH_H
