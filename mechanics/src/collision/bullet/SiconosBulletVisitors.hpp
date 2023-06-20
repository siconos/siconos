/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#ifndef SICONOSBULLET_VISITORS_HPP
#define SICONOSBULLET_VISITORS_HPP

#include "ShapeVisitor.hpp"
#include "SiconosBulletCollisionManager_impl.hpp"  // For BodyBoxRecord and so on
#include "SiconosShape.hpp"

namespace siconos::modeling {

class SecondOrderDS;
}

namespace siconos::collision {

class SiconosContactor;
class StaticBody;

}  // namespace siconos::collision

namespace siconos::collision::bullet::internal {

class SiconosBulletCollisionManager_impl;

class CreateCollisionObjectShapeVisitor : public siconos::collision::internal::ShapeVisitor {
 public:
  std::shared_ptr<SiconosBulletCollisionManager_impl> impl{nullptr};
  const std::shared_ptr<siconos::modeling::SecondOrderDS> ds{nullptr};
  const std::shared_ptr<siconos::algebra::SiconosVector> base{nullptr};
  std::shared_ptr<siconos::collision::SiconosContactor> contactor{nullptr};
  std::shared_ptr<siconos::collision::StaticBody> staticBody{nullptr};

  ~CreateCollisionObjectShapeVisitor() noexcept = default;

  CreateCollisionObjectShapeVisitor(
      std::shared_ptr<SiconosBulletCollisionManager_impl> _impl,
      const std::shared_ptr<siconos::modeling::SecondOrderDS> _ds,
      const std::shared_ptr<siconos::algebra::SiconosVector> _base,
      const std::shared_ptr<StaticBody> _staticBody)
      : impl(_impl), ds(_ds), base(_base), staticBody(_staticBody) {}

  virtual void visit(std::shared_ptr<siconos::collision::SiconosPlane> shape) override;

  virtual void visit(std::shared_ptr<siconos::collision::SiconosSphere> shape) override;

  virtual void visit(std::shared_ptr<siconos::collision::SiconosBox> shape) override;

  virtual void visit(std::shared_ptr<siconos::collision::SiconosCylinder> shape) override;

  virtual void visit(std::shared_ptr<siconos::collision::SiconosCone> shape) override;

  virtual void visit(std::shared_ptr<siconos::collision::SiconosCapsule> shape) override;

  virtual void visit(std::shared_ptr<siconos::collision::SiconosConvexHull> shape) override;

  virtual void visit(std::shared_ptr<siconos::collision::SiconosMesh> shape) override;

  virtual void visit(std::shared_ptr<siconos::collision::SiconosHeightMap> shape) override;

  virtual void visit(std::shared_ptr<siconos::collision::SiconosDisk> shape) override;

  virtual void visit(std::shared_ptr<siconos::collision::SiconosBox2d> shape) override;

  virtual void visit(std::shared_ptr<siconos::collision::SiconosConvexHull2d> shape) override;
};

}  // namespace siconos::collision::bullet::internal
#endif
