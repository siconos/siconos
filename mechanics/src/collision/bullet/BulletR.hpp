/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BulletR_hpp
#define BulletR_hpp

#include "BulletDeclarations.h"
#include "ContactR.hpp"

namespace siconos::modeling {
class NewtonEulerDS;
}

namespace siconos::collision::bullet {

class BulletR : public siconos::collision::ContactR {
 private:
  ACCEPT_SERIALIZATION(BulletR);

 public:
  virtual ~BulletR() noexcept = default;

  /* For users that may require extra information about contacts. */
  std::shared_ptr<btCollisionObject> btObject[2] = {nullptr, nullptr};
  std::shared_ptr<btCollisionShape> btShape[2] = {nullptr, nullptr};

  virtual void updateContactPointsFromManifoldPoint(
      const btPersistentManifold& manifold, const btManifoldPoint& point, bool flip,
      double scaling, std::shared_ptr<siconos::modeling::NewtonEulerDS> ds1,
      std::shared_ptr<siconos::modeling::NewtonEulerDS> ds2);

  void updateRelativeContactPointsFromManifoldPoint(
      const btPersistentManifold& manifold, const btManifoldPoint& point, bool flip,
      double scaling, std::shared_ptr<siconos::modeling::NewtonEulerDS> ds1,
      std::shared_ptr<siconos::modeling::NewtonEulerDS> ds2);
  void display() const override;

  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::collision::bullet

#endif
