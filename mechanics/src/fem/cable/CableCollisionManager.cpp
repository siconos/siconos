/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include "CableCollisionManager.h"

#include <memory>

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
  auto& q = cable_ds_->q_read();
  auto nb = q.size();

  unsigned int node_idx = 0;
  for (siconos::algebra::Index i = 0; i < nb; i += 3, node_idx++) {
    auto contactItr = contacts_map_.find(node_idx);
    std::shared_ptr<siconos::modeling::Interaction> contact = nullptr;
    if (contactItr != contacts_map_.end()) {
      contact = contactItr->second;
    }
    siconos::algebra::SiconosVector3 pc1;
    pc1(0) = q(i);
    pc1(1) = q(i + 1);
    pc1(2) = q(i + 2);

    for (auto& s : supports_) {
      if (s->isContact(pc1, tolAtContact_)) {
        // If the current point is in contact with the support then
        // we get it's projection (pc2) on the obstacle, the normal and the tangent.
        auto pc2 = s->pc2();
        auto normal = s->normal();
        auto tangent = s->tangent();

        if (contact) {  // The interaction already exists
          auto relation = std::dynamic_pointer_cast<Cable2d3DR>(contact->relation());
          assert(relation);
          relation->updateContactPoints(pc1, *pc2, *normal, *tangent);
        } else {
          // create relation
          auto relation = std::make_shared<Cable2d3DR>(node_idx, *pc2, *normal, *tangent);

          // create interaction
          auto inter = std::make_shared<siconos::modeling::Interaction>(s->nslaw(), relation);
          contacts_map_[node_idx] = inter;

          // link the interaction and the dynamical system
          simulation->link(inter, cable_ds_);
        }
      } else {
        // no contact
        if (contact) {
          // update the contacts map: no interaction for the current node
          contacts_map_.erase(node_idx);
        }
      }
    }
  }
}
