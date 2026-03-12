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

/*! \file ShapeVisitor.hpp
  \brief Definition of classes used to visit different Siconos Shape types.
  Internal use only (SiconosBulletCollisionmanager_impl : must not be propagated to the user
  interface.
*/

#ifndef SHAPE_VISITORS_HPP
#define SHAPE_VISITORS_HPP

#include "SiconosShape.hpp"

namespace siconos::collision::internal {

class ShapeVisitor {
 protected:
  ShapeVisitor() = default;
  ShapeVisitor(const ShapeVisitor &) = delete;
  ShapeVisitor(ShapeVisitor &&) = delete;
  ShapeVisitor &operator=(ShapeVisitor &&) = delete;
  ShapeVisitor &operator=(const ShapeVisitor &) = delete;

  virtual ~ShapeVisitor() noexcept = default;

 public:
  virtual void visit(std::shared_ptr<siconos::collision::SiconosPlane> shape) {};

  virtual void visit(std::shared_ptr<siconos::collision::SiconosSphere> shape) {};

  virtual void visit(std::shared_ptr<siconos::collision::SiconosBox> shape) {};

  virtual void visit(std::shared_ptr<siconos::collision::SiconosCylinder> shape) {};

  virtual void visit(std::shared_ptr<siconos::collision::SiconosCone> shape) {};

  virtual void visit(std::shared_ptr<siconos::collision::SiconosCapsule> shape) {};

  virtual void visit(std::shared_ptr<siconos::collision::SiconosConvexHull> shape) {};

  virtual void visit(std::shared_ptr<siconos::collision::SiconosMesh> shape) {};

  virtual void visit(std::shared_ptr<siconos::collision::SiconosHeightMap> shape) {};

  virtual void visit(std::shared_ptr<siconos::collision::SiconosDisk> shape) {};

  virtual void visit(std::shared_ptr<siconos::collision::SiconosBox2d> shape) {};

  virtual void visit(std::shared_ptr<siconos::collision::SiconosConvexHull2d> shape) {};
};

}  // namespace siconos::collision::internal
#endif



