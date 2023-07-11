/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2023 INRIA.
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

  Manager for interactions in the cable setup.
 */

namespace siconos::modeling {
class Interaction;
}

namespace siconos::fem::cable {

class CableDS;
class Support;

class CableCollisionManager : public siconos::simulation::InteractionManager {
 public:
  CableCollisionManager(const std::shared_ptr<CableDS> a_model,
                        const std::vector<std::shared_ptr<Support>> &a_supports,
                        double a_tolContact = 1e-3)
      : m_model{a_model}, m_supports{a_supports}, m_tolContact{a_tolContact} {};

  virtual ~CableCollisionManager() noexcept = default;

  virtual void updateInteractions(
      std::shared_ptr<siconos::simulation::Simulation> simulation) override;

 protected:
  std::shared_ptr<CableDS> m_model{nullptr};
  std::vector<std::shared_ptr<Support>> m_supports = {};
  double m_tolContact{1e-3};

  using t_contacts = std::map<unsigned int, std::shared_ptr<siconos::modeling::Interaction>>;
  t_contacts m_contacts;
};
}  // namespace siconos::fem::cable
