// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef BulletVelocityAngularVelocityR_hpp
#define BulletVelocityAngularVelocityR_hpp

#include "BulletDeclarations.h"
#include "ContactVelocityAngularVelocityR.hpp"

namespace siconos::modeling {
class NewtonEulerDS;
}

namespace siconos::collision::bullet {
class BulletVelocityAngularVelocityR
    : public siconos::collision::ContactVelocityAngularVelocityR {
 private:
  ACCEPT_SERIALIZATION(BulletVelocityAngularVelocityR);

 public:
  virtual ~BulletVelocityAngularVelocityR() noexcept = default;

  /* For users that may require extra information about contacts. */
  std::shared_ptr<btCollisionObject> btObject[2] = {nullptr, nullptr};
  std::shared_ptr<btCollisionShape> btShape[2] = {nullptr, nullptr};

  virtual void updateContactPointsFromManifoldPoint(
      const btPersistentManifold& manifold, const btManifoldPoint& point, bool flip,
      double scaling, std::shared_ptr<siconos::modeling::NewtonEulerDS> ds1,
      std::shared_ptr<siconos::modeling::NewtonEulerDS> ds2);
};
}  // namespace siconos::collision::bullet
#endif
