/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include <MechanicsFwd.hpp>
#include <BulletSiconosFwd.hpp>

#include <SiconosCollisionManager.hpp>
#include <SiconosShape.hpp>
#include <SiconosContactor.hpp>

#include <map>

DEFINE_SPTR(SiconosBulletCollisionManager_impl);


enum SiconosBulletDimension
{
  SICONOS_BULLET_3D=0,
  SICONOS_BULLET_2D=1
};

struct SiconosBulletOptions
{
  SiconosBulletOptions();

  int dimension;
  double contactBreakingThreshold;
  double contactProcessingThreshold;
  double worldScale;
  bool useAxisSweep3;
  bool clearOverlappingPairCache;
  unsigned int perturbationIterations;
  unsigned int minimumPointsPerturbationThreshold;
  bool enableSatConvex;
  bool enablePolyhedralContactClipping;
};

struct SiconosBulletStatistics
{
  SiconosBulletStatistics()
    : new_interactions_created(0)
    , existing_interactions_processed(0)
    , interaction_warnings(0)
    {}
  int new_interactions_created;
  int existing_interactions_processed;
  int interaction_warnings;
};

class SiconosBulletCollisionManager : public SiconosCollisionManager
{
protected:
  SP::SiconosBulletCollisionManager_impl impl;

  void initialize_impl();

  // callback for contact point removal, and a global for context
  static bool bulletContactClear(void* userPersistentData);
  static Simulation *gSimulation;

public:
  SiconosBulletCollisionManager();
  SiconosBulletCollisionManager(const SiconosBulletOptions &options);
  virtual ~SiconosBulletCollisionManager();

protected:
  bool _with_equality_constraints;
  SiconosBulletOptions _options;
  SiconosBulletStatistics _stats;

  /** Provided so that creation of collision points can be overridden.
   * See modify_normals.py in examples/Mechanics/Hacks */
  virtual SP::BulletR makeBulletR(SP::RigidBodyDS ds1, SP::SiconosShape shape1,
                                  SP::RigidBodyDS ds2, SP::SiconosShape shape2,
                                  const btManifoldPoint &);
  
  /** Provided so that creation of collision points can be overridden.
   * See modify_normals.py in examples/Mechanics/Hacks */
  virtual SP::Bullet5DR makeBullet5DR(SP::RigidBodyDS ds1, SP::SiconosShape shape1,
                                      SP::RigidBodyDS ds2, SP::SiconosShape shape2,
                                      const btManifoldPoint &);

  /** Provided so that creation of collision points can be overridden.
   * See modify_normals.py in examples/Mechanics/Hacks */
  virtual SP::Bullet2dR makeBullet2dR(SP::RigidBody2dDS ds1, SP::SiconosShape shape1,
                                      SP::RigidBody2dDS ds2, SP::SiconosShape shape2,
                                      const btManifoldPoint &);

  /** Provided so that creation of collision points can be overridden.
   * See modify_normals.py in examples/Mechanics/Hacks */
  virtual SP::Bullet2d3DR makeBullet2d3DR(SP::RigidBody2dDS ds1, SP::SiconosShape shape1,
                                          SP::RigidBody2dDS ds2, SP::SiconosShape shape2,
                                          const btManifoldPoint &);

public:
  StaticContactorSetID insertStaticContactorSet(
    SP::SiconosContactorSet cs, SP::SiconosVector position = SP::SiconosVector());

  bool removeStaticContactorSet(StaticContactorSetID id);

  void removeBody(const SP::RigidBodyDS& body);
  void removeBody(const SP::RigidBody2dDS& body);

  void updateInteractions(SP::Simulation simulation);

  std::vector<SP::SiconosCollisionQueryResult>
  lineIntersectionQuery(const SiconosVector& start, const SiconosVector& end,
                        bool closestOnly=false, bool sorted=true);

  void clearOverlappingPairCache();

  const SiconosBulletOptions &options() const { return _options; }
  const SiconosBulletStatistics &statistics() const { return _stats; }
  void resetStatistics() { _stats = SiconosBulletStatistics(); }

  /** Set the usage of equality constraints. When the number
      of objects is huge as in granular material, the usage
      of equality constraint breaks scalability.
      This have to be fixed.
   * \param choice a boolean, default is True.
   */
  void useEqualityConstraints(bool choice=true)
  { _with_equality_constraints = choice; };
};

#endif /* SiconosBulletCollisionManager.hpp */
