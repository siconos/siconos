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

struct OccContactEdge : public OccContactShape {
  OccContactEdge(const OccContactShape& shape, int index);

  ~OccContactEdge() noexcept = default;

  virtual std::shared_ptr<TopoDS_Edge> contact() const;

  virtual void computeUVBounds() override;

  std::shared_ptr<const TopoDS_Edge> _edge{nullptr};
  int _index{0};
};
}  // namespace siconos::mechanics::occ
#endif
