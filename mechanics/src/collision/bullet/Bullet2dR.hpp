/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef Bullet2dR_hpp
#define Bullet2dR_hpp

#include "BulletDeclarations.h"
#include "Contact2dR.hpp"

namespace siconos::collision {
class RigidBody2dDS;
}
namespace siconos::collision::bullet {

class Bullet2dR : public siconos::collision::Contact2dR {
 private:
  ACCEPT_SERIALIZATION(Bullet2dR);

 public:
  virtual ~Bullet2dR() noexcept = default;

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
