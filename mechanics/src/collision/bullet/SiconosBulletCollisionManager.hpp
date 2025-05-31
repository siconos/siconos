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

/*! \file SiconosBulletCollisionManager.hpp
  \brief Definition of a Bullet-based interaction handler for contact
  detection.
*/

#ifndef SiconosBulletCollisionManager_h
#define SiconosBulletCollisionManager_h

#include "BulletDeclarations.h"
#include "SiconosBulletOptions.hpp"  // header only, for options and stats
#include "SiconosCollisionManager.hpp"

namespace siconos::collision {
class RigidBodyDS;
class RigidBody2dDS;
class SiconosShape;
}  // namespace siconos::collision

namespace siconos::collision::bullet {

class BulletR;
class Bullet2d3DR;
class Bullet5DR;
class Bullet2dR;

namespace internal {  // An "impl" class to hide implementation. See
                      // SiconosBulletCollisionManager.cpp
class SiconosBulletCollisionManager_impl;
}

class SiconosBulletCollisionManager : public siconos::collision::SiconosCollisionManager {
 protected:
  ACCEPT_SERIALIZATION(SiconosBulletCollisionManager);

  using BodiesVariant = std::variant<std::shared_ptr<siconos::collision::RigidBodyDS>,
                                     std::shared_ptr<siconos::collision::RigidBody2dDS>>;

  // This value is compared to the initial distance computed
  // at the creation of the interaction
  // if distance < - WARNING_TOLERANCE_AT_CREATION_INTERACTION
  // a warning is raised.
  static constexpr double WARNING_TOLERANCE_AT_CREATION_INTERACTION = 1e-5;

  std::shared_ptr<internal::SiconosBulletCollisionManager_impl> _impl{nullptr};

  static siconos::simulation::Simulation* gSimulation;

  bool _with_equality_constraints{false};
  std::shared_ptr<SiconosBulletOptions> _options{nullptr};
  SiconosBulletStatistics _stats{};

  void initialize_impl();

  // callback for contact point removal, and a global for context
  static bool bulletContactClear(void* userPersistentData);

  // callback to modify the contact point when it has just been added in the
  // manifold.
  static bool bulletContactAddedCallback(btManifoldPoint& cp,
                                         const btCollisionObjectWrapper* colObj0Wrap,
                                         int partId0, int index0,
                                         const btCollisionObjectWrapper* colObj1Wrap,
                                         int partId1, int index1);

  /** Provided so that creation of collision points can be overridden.
   *  See modify_normals.py in examples/Mechanics/Hacks */
  virtual std::shared_ptr<BulletR> makeBulletR(
      std::shared_ptr<siconos::collision::RigidBodyDS> ds1,
      std::shared_ptr<siconos::collision::SiconosShape> shape1,
      std::shared_ptr<siconos::collision::RigidBodyDS> ds2,
      std::shared_ptr<siconos::collision::SiconosShape> shape2, const btManifoldPoint&);

  /** Provided so that creation of collision points can be overridden.
   *  See modify_normals.py in examples/Mechanics/Hacks */
  virtual std::shared_ptr<Bullet5DR> makeBullet5DR(
      std::shared_ptr<siconos::collision::RigidBodyDS> ds1,
      std::shared_ptr<siconos::collision::SiconosShape> shape1,
      std::shared_ptr<siconos::collision::RigidBodyDS> ds2,
      std::shared_ptr<siconos::collision::SiconosShape> shape2, const btManifoldPoint&);

  /** Provided so that creation of collision points can be overridden.
   *  See modify_normals.py in examples/Mechanics/Hacks */
  virtual std::shared_ptr<Bullet2dR> makeBullet2dR(
      std::shared_ptr<siconos::collision::RigidBody2dDS> ds1,
      std::shared_ptr<siconos::collision::SiconosShape> shape1,
      std::shared_ptr<siconos::collision::RigidBody2dDS> ds2,
      std::shared_ptr<siconos::collision::SiconosShape> shape2, const btManifoldPoint&);

  /** Provided so that creation of collision points can be overridden.
   *  See modify_normals.py in examples/Mechanics/Hacks */
  virtual std::shared_ptr<Bullet2d3DR> makeBullet2d3DR(
      std::shared_ptr<siconos::collision::RigidBody2dDS> ds1,
      std::shared_ptr<siconos::collision::SiconosShape> shape1,
      std::shared_ptr<siconos::collision::RigidBody2dDS> ds2,
      std::shared_ptr<siconos::collision::SiconosShape> shape2, const btManifoldPoint&);

 public:
  SiconosBulletCollisionManager();
  SiconosBulletCollisionManager(std::shared_ptr<SiconosBulletOptions> options);
  virtual ~SiconosBulletCollisionManager() noexcept;

  /** Add a static body in the collision detector.
   */
  std::shared_ptr<siconos::collision::StaticBody> addStaticBody(
      std::shared_ptr<siconos::collision::SiconosContactorSet> cs,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& position =
          std::nullopt,
      int number = 0) override;

  /** Remove a body from the collision detector.
   */
  void removeStaticBody(const std::shared_ptr<siconos::collision::StaticBody>& body) override;

  /** Remove a body from the collision detector. This must be done
   *  after removing a body from the NonSmoothDynamicalSystem
   *  otherwise contact will occur with a non-graph body which results
   *  in failure. */
  void removeBody(const std::shared_ptr<siconos::modeling::SecondOrderDS>& body) override;

  void updateInteractions(
      std::shared_ptr<siconos::simulation::Simulation> simulation) override;

  std::vector<std::shared_ptr<siconos::collision::SiconosCollisionQueryResult>>
  lineIntersectionQuery(const siconos::algebra::SiconosVector& start,
                        const siconos::algebra::SiconosVector& end, bool closestOnly = false,
                        bool sorted = true) override;

  void clearOverlappingPairCache();

  auto options() const { return _options; }
  const SiconosBulletStatistics& statistics() const { return _stats; }
  void resetStatistics() { _stats = SiconosBulletStatistics(); }

  /**
     Set the usage of equality constraints. When the number
     of objects is huge as in granular material, the usage
     of equality constraint breaks scalability.
     This have to be fixed.

     \param choice a boolean, default is True.
  */
  void useEqualityConstraints(bool choice = true) { _with_equality_constraints = choice; };
};
}  // namespace siconos::collision::bullet
#endif /* SiconosBulletCollisionManager.hpp */
