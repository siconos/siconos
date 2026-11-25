// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
