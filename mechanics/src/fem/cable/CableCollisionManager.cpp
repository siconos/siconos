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
#include "CableCollisionManager.h"

#include "Cable2d3DR.hpp"
#include "CableDS.hpp"
#include "Interaction.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "Support.h"

/** Called by Simulation after updating positions prior to starting
 * the Newton loop. */
void siconos::fem::cable::CableCollisionManager::updateInteractions(
    std::shared_ptr<siconos::simulation::Simulation> simulation) {
  siconos::algebra::SiconosVector &q = *(m_model->q());
  size_t nb = q.size();

  unsigned int node_idx = 0;
  for (size_t i = 0; i < nb; i += 3, node_idx++) {  // boucle par 3 pour récupérer x,y,z

    auto contactItr = m_contacts.find(node_idx);
    std::shared_ptr<siconos::modeling::Interaction> contact = nullptr;
    if (contactItr != m_contacts.end()) {
      contact = contactItr->second;
    }
    auto pc1 = std::make_shared<siconos::algebra::SiconosVector>(3);
    pc1->setValue(0, q.getValue(i));
    pc1->setValue(1, q.getValue(i + 1));
    pc1->setValue(2, q.getValue(i + 2));

    for (auto &s : m_supports) {
      if (s->isContact(pc1, m_tolContact)) {
        // test si le point x,y,z est en contact avec le support
        // récupérer pc2 (projection du point sur l'obstacle), normal, tangent
        auto pc2 = s->pc2();
        auto normal = s->normal();
        auto tangent = s->tangent();

        if (contact) {
          auto relation = std::static_pointer_cast<Cable2d3DR>(contact->relation());
          relation->updateContactPoint(pc1, pc2, normal, tangent);
        } else {
          // create relation
          auto relation = std::make_shared<Cable2d3DR>(node_idx, pc2, normal, tangent);

          // create interaction
          auto inter = std::make_shared<siconos::modeling::Interaction>(s->nslaw(), relation);
          m_contacts[node_idx] = inter;

          // link the interaction and the dynamical system
          simulation->link(inter, m_model);
        }
      } else {
        // pas en contact
        if (contact) {
          // si existe un contact -> remove
          m_contacts.erase(node_idx);
        }
      }
    }
  }
  std::cout << "Fin UpInterac \n";
}
