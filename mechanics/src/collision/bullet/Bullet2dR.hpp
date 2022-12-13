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
