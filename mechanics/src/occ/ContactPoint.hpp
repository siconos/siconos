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

  const OccContactShapeV& contactShape;

  double u{0.};
  double v{0.};
};
}  // namespace siconos::mechanics::occ
#endif
