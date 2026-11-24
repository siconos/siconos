/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/** \file ContactPoint.hpp
    \brief Implementation details for OccR class.
 */
#ifndef ContactPoint_hpp
#define ContactPoint_hpp

#include "OccContactShape.hpp"

namespace siconos::mechanics::occ {

struct OccContactShape;

struct ContactPoint {
  ContactPoint(const OccContactShapeV& shape) : contactShape{shape} {};

  // Rule of 5
  ContactPoint() = delete;
  ContactPoint(const ContactPoint&) = delete;
  ContactPoint(ContactPoint&&) = delete;
  ContactPoint& operator=(const ContactPoint&) = delete;
  ContactPoint& operator=(ContactPoint&&) = delete;
  ~ContactPoint() noexcept = default;

  // const OccContactShape& contactShape() const { return _shape; };

  // by value: the variant holds a shared pointer to the shape (a reference could be bound
  // to a temporary, e.g. created by the python bindings)
  const OccContactShapeV contactShape;

  double u{0.};
  double v{0.};
};
}  // namespace siconos::mechanics::occ
#endif
