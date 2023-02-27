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
#pragma once
#include <CableDS.hpp>
#include <InteractionManager.hpp>
#include <Support.h>

using namespace std;
using namespace siconos::mechanics::fem;

class CableCollisionManager : public InteractionManager {

public:
  CableCollisionManager(const std::shared_ptr<CableDS> a_model,
                        const vector<std::shared_ptr<Support>> &a_supports,
                        double a_tolContact = 1e-3);
  virtual ~CableCollisionManager();
  virtual void updateInteractions(std::shared_ptr<Simulation> simulation);

protected:
  std::shared_ptr<CableDS> m_model;
  vector<std::shared_ptr<Support>> m_supports;
  double m_tolContact;

  typedef std::map<unsigned int, std::shared_ptr<Interaction>> t_contacts;
  t_contacts m_contacts;
};
