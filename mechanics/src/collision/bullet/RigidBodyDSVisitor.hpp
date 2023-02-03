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

#ifndef RIGIDBODYDS_VISITORS_HPP
#define RIGIDBODYDS_VISITORS_HPP

#include "SiconosVisitor.hpp"

namespace siconos::collision {
class RigidBody2dDS;
class RigidBodyDS;

}  // namespace siconos::collision

namespace siconos::collision::bullet::internal {
class SiconosBulletCollisionManager_impl;

class RigidBodyDSVisitor : public siconos::internal::SiconosVisitor {
 protected:
  std::shared_ptr<SiconosBulletCollisionManager_impl> impl{nullptr};

  using SiconosVisitor::visit;

  RigidBodyDSVisitor() = delete;
  RigidBodyDSVisitor(const RigidBodyDSVisitor &) = delete;
  RigidBodyDSVisitor(RigidBodyDSVisitor &&) = delete;
  RigidBodyDSVisitor &operator=(RigidBodyDSVisitor &&) = delete;
  RigidBodyDSVisitor &operator=(const RigidBodyDSVisitor &) = delete;

 public:
  RigidBodyDSVisitor(std::shared_ptr<SiconosBulletCollisionManager_impl> implv)
      : SiconosVisitor{}, impl{implv} {};

  virtual ~RigidBodyDSVisitor() noexcept = default;

  virtual void visit(std::shared_ptr<siconos::collision::RigidBodyDS> bds);

  virtual void visit(std::shared_ptr<siconos::collision::RigidBody2dDS> bds);
};

}  // namespace siconos::collision::bullet::internal
#endif
