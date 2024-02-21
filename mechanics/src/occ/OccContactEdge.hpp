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
#ifndef OccContactEdge_hpp
#define OccContactEdge_hpp

#include "OccContactShape.hpp"

// OpenCASCADE classes
class TopoDS_Edge;

namespace siconos::mechanics::occ {

// template <typename T>
// class OccContactShape1

// {
//  public:
//   std::shared_ptr<TopoDS_Shape> _shape{nullptr};
//   std::shared_ptr<const T> _edge{nullptr};
//   unsigned int _index{0};
//   unsigned int contactGroup{0};
//   unsigned int _id{0};
//   std::array<double, 8> boundaries{0., 0., 0., 0., 0., 0., 0., 0.};

//   OccContactShape1(const OccContactShape& shape, unsigned int index);

//   ~OccContactShape1() noexcept = default;

//   const std::shared_ptr<const T> contact() const;

//   void computeUVBounds();

//   std::shared_ptr<T> edge_or_face(int index) const;
// };

struct OccContactEdge : public OccContactShape {
  OccContactEdge(const OccContactShape& shape, unsigned int index);

  virtual ~OccContactEdge() noexcept = default;

  virtual const std::shared_ptr<const TopoDS_Edge> contact() const;

  virtual void computeUVBounds();

  std::shared_ptr<const TopoDS_Edge> _edge{nullptr};
  unsigned int _index{0};
};
}  // namespace siconos::mechanics::occ
#endif
