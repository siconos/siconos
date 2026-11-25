// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
                        const std::vector<std::shared_ptr<Support>>& a_supports,
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
