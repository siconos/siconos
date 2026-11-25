// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef OccContactFace_hpp
#define OccContactFace_hpp

#include "OccContactShape.hpp"

// OpenCASCADE classes
class TopoDS_Face;

namespace siconos::mechanics::occ {
struct OccContactFace : public OccContactShape {
  OccContactFace(const OccContactShape& shape, int index);
  ~OccContactFace() noexcept = default;

  virtual std::shared_ptr<TopoDS_Face> contact() const;

  virtual void computeUVBounds() override;

  int _index{0};
  std::shared_ptr<const TopoDS_Face> _face{nullptr};
};
}  // namespace siconos::mechanics::occ
#endif
