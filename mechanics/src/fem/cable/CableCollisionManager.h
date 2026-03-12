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
#pragma once
#include "InteractionManager.hpp"

/*! \file CableCollisionmanager.h
  CableCollisionManager class
*/

namespace siconos::modeling {
class Interaction;
}

namespace siconos::fem::cable {

class CableDS;
class Support;

/**
 *  A class which handles:
 *   - a dynamical system (the cable)
 *   - a vector of supports
 *   and manage their interactions.
 */
class CableCollisionManager : public siconos::simulation::InteractionManager {
 public:
  CableCollisionManager(const std::shared_ptr<CableDS> a_cableDS,
                        const std::vector<std::shared_ptr<Support>> &a_supports,
                        double a_tolContact = 1e-3)
      : cable_ds_{a_cableDS}, supports_{a_supports}, tolAtContact_{a_tolContact} {};

  virtual ~CableCollisionManager() noexcept = default;

  virtual void updateInteractions(
      std::shared_ptr<siconos::simulation::Simulation> simulation) override;

 protected:
  std::shared_ptr<CableDS> cable_ds_{nullptr};
  std::vector<std::shared_ptr<Support>> supports_ = {};
  double tolAtContact_{1e-3};

  using t_contacts = std::map<unsigned int, std::shared_ptr<siconos::modeling::Interaction>>;
  t_contacts contacts_map_;
};
}  // namespace siconos::fem::cable
