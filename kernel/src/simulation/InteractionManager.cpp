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

#include "InteractionManager.hpp"

#include <algorithm>

#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NonSmoothLaw.hpp"
#include "Relation.hpp"
#include "Simulation.hpp"
#include "SimulationGraphs.hpp"
#include "siconos_debug.h"
#include "SiconosMemory.hpp"

void siconos::simulation::InteractionManager::insertNonSmoothLaw(
    std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw, long unsigned int group1,
    long unsigned int group2) {
  // ublas::matrix size type is not the same on 32 bits and 64 bits
  auto maxgroup = std::max((siconos::modeling::NSLawMatrix::size_type)group1,
                           (siconos::modeling::NSLawMatrix::size_type)group2);
  _nslaws.resize(std::max(_nslaws.size1(), maxgroup + 1));
  _nslaws(group1, group2) = nslaw;
}

std::shared_ptr<siconos::modeling::NonSmoothLaw>
siconos::simulation::InteractionManager::nonSmoothLaw(long unsigned int group1,
                                                      long unsigned int group2) {
  if (group1 < _nslaws.size1() && group2 < _nslaws.size2())
    return _nslaws(group1, group2);
  else
    return std::shared_ptr<siconos::modeling::NonSmoothLaw>();  // ??
}

bool siconos::simulation::interactions_manager::is_interaction_present(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds1,
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds2,
    std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> DSG0) {
  for (auto [oei, oeiend] = DSG0->out_edges(DSG0->descriptor(ds1)); oei != oeiend; ++oei) {
    if (DSG0->bundle(DSG0->target(*oei)) == ds2) return true;
  }
  return false;
}

void siconos::simulation::interactions_manager::build_and_link_interaction(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds1,
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds2,
    std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> DSG0,
    std::shared_ptr<siconos::modeling::Relation> rel,
    std::shared_ptr<siconos::simulation::Simulation> sim,
    std::shared_ptr<siconos::simulation::InteractionManager> parent) {
  if (!is_interaction_present(ds1, ds2, DSG0)) {
    auto nslaw = parent->nonSmoothLaw(DSG0->groupId[DSG0->descriptor(ds1)],
                                      DSG0->groupId[DSG0->descriptor(ds2)]);
    assert(nslaw);
    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, rel);
    sim->nonSmoothDynamicalSystem()->link(inter, ds1, ds2);
  }
}

void siconos::simulation::interactions_manager::remove_interaction_if_exists(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds1,
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds2,
    std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> DSG0,
    std::shared_ptr<siconos::simulation::Simulation> sim) {
  for (auto [oei, oeiend] = DSG0->out_edges(DSG0->descriptor(ds1)); oei != oeiend; ++oei) {
    if (DSG0->bundle(DSG0->target(*oei)) == ds2) {
      DEBUG_PRINTF("remove interaction : %d\n", DSG0->bundle(*oei)->number());
      sim->unlink(DSG0->bundle(*oei));
      return;
    }
  }
}