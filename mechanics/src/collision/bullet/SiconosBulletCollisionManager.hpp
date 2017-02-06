/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#include <SiconosCollisionManager.hpp>
#include <SiconosShape.hpp>
#include <SiconosContactor.hpp>

#include <map>

DEFINE_SPTR(SiconosBulletCollisionManager_impl);

struct SiconosBulletOptions
{
  SiconosBulletOptions();

  double breakingThreshold;
  double worldScale;
  bool useAxisSweep3;
  unsigned int perturbationIterations;
  unsigned int minimumPointsPerturbationThreshold;
  double contactProcessingThreshold;
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
  SiconosBulletOptions _options;
  SiconosBulletStatistics _stats;

public:
  StaticContactorSetID insertStaticContactorSet(
    SP::SiconosContactorSet cs, SP::SiconosVector position = SP::SiconosVector());

  bool removeStaticContactorSet(StaticContactorSetID id);

  void updateInteractions(SP::Simulation simulation);

  const SiconosBulletOptions &options() const { return _options; }
  const SiconosBulletStatistics &statistics() const { return _stats; }
  void resetStatistics() { _stats = SiconosBulletStatistics(); }
};

#endif /* SiconosBulletCollisionManager.hpp */
