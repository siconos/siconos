/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef Bullet2d3DR_hpp
#define Bullet2d3DR_hpp

#include "BulletDeclarations.h"
#include "Contact2d3DR.hpp"

namespace siconos::collision {
class RigidBody2dDS;
}

namespace siconos::collision::bullet {
class Bullet2d3DR : public siconos::collision::Contact2d3DR {
 private:
  ACCEPT_SERIALIZATION(Bullet2d3DR);

 public:
  virtual ~Bullet2d3DR() noexcept = default;

  /* For users that may require extra information about contacts. */
  std::shared_ptr<btCollisionObject> btObject[2] = {nullptr, nullptr};
  std::shared_ptr<btCollisionShape> btShape[2] = {nullptr, nullptr};

  virtual void updateContactPointsFromManifoldPoint(
      const btPersistentManifold& manifold, const btManifoldPoint& point, bool flip,
      double scaling, std::shared_ptr<siconos::collision::RigidBody2dDS> ds1,
      std::shared_ptr<siconos::collision::RigidBody2dDS> ds2);
};
}  // namespace siconos::collision::bullet
#endif
