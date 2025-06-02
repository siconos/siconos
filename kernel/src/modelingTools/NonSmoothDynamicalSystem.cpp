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
#include "NonSmoothDynamicalSystem.hpp"

#include "DynamicalSystem.hpp"
#include "Interaction.hpp"
#include "Relation.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "SimulationGraphs.hpp"
#include "Tools.hpp"  // for enum_to_string
#include "Topology.hpp"

// #include <SiconosConfig.h>
// #include <functional>
// using namespace std::placeholders;

// #include <limits>

// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "siconos_debug.h"

//  constructor
siconos::modeling::NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(double t0, double T)
    : _t0(t0), _T(T) {
  // === Builds an empty topology ===
  _topology = std::make_shared<siconos::simulation::Topology>();
  // we push a first element in the list to avoid acces to null when
  // we call --_changeLog.end();
  _changeLog.push_back(Change(ChangeType::clearTopology));
  DEBUG_EXPR(std::prev(_changeLog.end())->display());

  // see Simulation::initialize() for an explanation of why we
  // implement this changelog
};

siconos::modeling::NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem() noexcept { clear(); }

// changelog
void siconos::modeling::NonSmoothDynamicalSystem::Change::display() const {
  std::cout << "Changes display   " << this << std::endl;
  auto changeval = siconos::tools::enum_to_string(typeOfChange);
  if (typeOfChange == ChangeType::addDynamicalSystem) {
    std::cout << "typeOfChange : " << changeval << " : addDynamicalSystem\n";
  } else if (typeOfChange == ChangeType::rmDynamicalSystem) {
    std::cout << "typeOfChange : " << changeval << " : rmDynamicalSystem\n";
  } else if (typeOfChange == ChangeType::addInteraction) {
    std::cout << "typeOfChange : " << changeval << " : addInteraction\n";
  } else if (typeOfChange == ChangeType::rmInteraction) {
    std::cout << "typeOfChange : " << changeval << " : rmInteraction\n";
  } else if (typeOfChange == ChangeType::clearTopology) {
    std::cout << "typeOfChange : " << changeval << " : clearTopology\n";
  }
}

void siconos::modeling::NonSmoothDynamicalSystem::clearChangeLogTo(
    const std::list<Change>::const_iterator& it) {
  /* Given an interator into the changelog list, clear everything that
   * comes before it. User must be careful calling this if he has two
   * simulations, but in the one-simulation case (currently 100% of
   * cases), calling this will prevent changelog from building up
   * forever. Important especially for simulations using an
   * InteractionManager, e.g. mechanics_run.py. */
  _changeLog.erase(_changeLog.begin(), it);
}

// === DynamicalSystems management ===

void siconos::modeling::NonSmoothDynamicalSystem::display() const {
  std::cout << " ===== Non Smooth Dynamical System display =====\n ";
  std::cout << "---> isBVP = " << _BVP << std::endl;
  _topology->indexSet0()->display();
  std::cout << "---> last change :\n";
  std::prev(_changeLog.end())->display();
  std::cout << "===================================================\n";
}

void siconos::modeling::NonSmoothDynamicalSystem::insertDynamicalSystem(
    std::shared_ptr<DynamicalSystem> ds) {
  // some checks here ...
  if (!ds) {
    THROW_EXCEPTION(
        "siconos::modeling::NonSmoothDynamicalSystem::insertDynamicalSystem :: DS is nul");
  }

  // Do not insert the same ds several times : results in errors in initialisation process.
  if (!_topology->hasDynamicalSystem(ds)) {
    _topology->insertDynamicalSystem(ds);
    _changeLog.push_back(Change(ChangeType::addDynamicalSystem, ds));
    _mIsLinear = ((ds)->isLinear() && _mIsLinear);
  }
}

void siconos::modeling::NonSmoothDynamicalSystem::removeDynamicalSystem(
    std::shared_ptr<DynamicalSystem> ds) {
  _topology->removeDynamicalSystem(ds);
  _changeLog.push_back(Change(ChangeType::rmDynamicalSystem, ds));
}

void siconos::modeling::NonSmoothDynamicalSystem::removeInteraction(
    std::shared_ptr<Interaction> inter) {
  _topology->removeInteraction(inter);
  _changeLog.push_back(Change(ChangeType::rmInteraction, inter));
}

void siconos::modeling::NonSmoothDynamicalSystem::link(std::shared_ptr<Interaction> inter,
                                                       std::shared_ptr<DynamicalSystem> ds1,
                                                       std::shared_ptr<DynamicalSystem> ds2) {
  _mIsLinear = (inter->relation()->isLinear() && _mIsLinear);
  _topology->link(inter, ds1, ds2);
  _changeLog.push_back(Change(ChangeType::addInteraction, inter));
};

void siconos::modeling::NonSmoothDynamicalSystem::clear() {
  _topology->clear();
  _changeLog.push_back(Change(ChangeType::clearTopology));
}

void siconos::modeling::NonSmoothDynamicalSystem::setSymmetric(bool val) {
  _topology->setSymmetric(val);
}

void siconos::modeling::NonSmoothDynamicalSystem::resetNonSmoothPart(unsigned int level) {
  for (auto& vi : *dynamicalSystems()) {
    dynamicalSystems()->bundle(vi)->resetNonSmoothPart(level);
  }
}

void siconos::modeling::NonSmoothDynamicalSystem::swapInMemory() {
  // could be better to call bind method
  for (auto& vi : *dynamicalSystems()) {
    dynamicalSystems()->bundle(vi)->swapInMemory();
  }
}
void siconos::modeling::NonSmoothDynamicalSystem::pushInteractionsInMemory() {
  // Save Interactions state into Memory.

  if (_topology->indexSet0()->size() > 0) {
    // Temp FP : saveInOldVar was called for each osns and each osns call
    // swapInOldVar for all interactions in the nsds.
    // ==> let's do it only once, by the simu.

    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    auto indexSet0 = _topology->indexSet0();
    for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
      indexSet0->bundle(*ui)->swapInMemory();
    }
  }
}
void siconos::modeling::NonSmoothDynamicalSystem::updateDSPlugins(double time) {
  // could be better to call bind method
  for (auto& vi : *dynamicalSystems()) {
    dynamicalSystems()->bundle(vi)->updatePlugins(time);
  }
}
void siconos::modeling::NonSmoothDynamicalSystem::updateInput(double time,
                                                              unsigned int level) {
  DEBUG_BEGIN("Nonsmoothdynamicalsystem::updateInput(double time, unsigned int level)\n");
  DEBUG_PRINTF("with level = %i\n", level);

  // To compute input(level) (ie with lambda[level]) for all Interactions.
  //  assert(level>=0);

  // Set dynamical systems non-smooth part to zero.
  resetNonSmoothPart(level);

  // We compute input using lambda(level).
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  std::shared_ptr<siconos::modeling::Interaction> inter;
  auto indexSet0 = _topology->indexSet0();
  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForInput() <= level);
    assert(inter->upperLevelForInput() >= level);
    inter->computeInput(time, level);
  }

  DEBUG_END("Nonsmoothdynamicalsystem::updateInput(double time, unsigned int level)\n");
}

void siconos::modeling::NonSmoothDynamicalSystem::updateOutput(double time,
                                                               unsigned int level) {
  // To compute output(level) (ie with y[level]) for all Interactions.
  //  assert(level>=0);

  DEBUG_BEGIN(
      "siconos::modeling::NonSmoothDynamicalSystem::updateOutput(unsigned int level)\n");
  DEBUG_PRINTF("with level = %i\n", level);
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  std::shared_ptr<siconos::modeling::Interaction> inter;
  auto indexSet0 = _topology->indexSet0();
  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForOutput() <= level);
    assert(inter->upperLevelForOutput() >= level);
    inter->computeOutput(time, level);
  }
  DEBUG_END("siconos::modeling::NonSmoothDynamicalSystem::updateOutput(unsigned int level)\n");
}

void siconos::modeling::NonSmoothDynamicalSystem::updateOutput(double time,
                                                               unsigned int level_min,
                                                               unsigned int level_max) {
  // To compute output(level) (ie with y[level]) for all Interactions in I0
  // and for a range of levels in a single pass through I0.
  //  assert(level>=0);

  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  std::shared_ptr<siconos::modeling::Interaction> inter;
  auto indexSet0 = _topology->indexSet0();
  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForOutput() <= level_max);
    assert(inter->upperLevelForOutput() >= level_min);
    for (unsigned int level = level_min; level <= level_max; ++level)
      inter->computeOutput(time, level);
  }
}

void siconos::modeling::NonSmoothDynamicalSystem::computeInteractionJacobians(double time) {
  DEBUG_BEGIN(
      "siconos::modeling::NonSmoothDynamicalSystem::computeInteractionJacobians(double "
      "time)\n");
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  std::shared_ptr<siconos::modeling::Interaction> inter;
  auto indexSet0 = _topology->indexSet0();
  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    inter = indexSet0->bundle(*ui);
    inter->relation()->computeJach(time, *inter);
    inter->relation()->computeJacg(time, *inter);
  }
  DEBUG_END(
      "siconos::modeling::NonSmoothDynamicalSystem::computeInteractionJacobians(double "
      "time)\n");
}

void siconos::modeling::NonSmoothDynamicalSystem::computeInteractionJacobians(
    double time, siconos::graphs::InteractionsGraph& indexSet) {
  DEBUG_BEGIN(
      "siconos::modeling::NonSmoothDynamicalSystem::computeInteractionJacobians(double "
      "time)\n");
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  std::shared_ptr<siconos::modeling::Interaction> inter;
  for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui) {
    inter = indexSet.bundle(*ui);
    inter->relation()->computeJach(time, *inter);
    inter->relation()->computeJacg(time, *inter);
  }
  DEBUG_END(
      "siconos::modeling::NonSmoothDynamicalSystem::computeInteractionJacobians(double "
      "time)\n");
}

void siconos::modeling::NonSmoothDynamicalSystem::visitDynamicalSystems(
    siconos::modeling::dynamical_systems::Visitor& visitor) {
  auto& dsg = *dynamicalSystems();
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsiend;
  std::tie(dsi, dsiend) = dsg.vertices();
  for (; dsi != dsiend; ++dsi) {
    dsg.bundle(*dsi)->accept(visitor);
  }
}

size_t siconos::modeling::NonSmoothDynamicalSystem::getNumberOfDS() const {
  return _topology->dSG(0)->size();
}

const std::shared_ptr<siconos::graphs::DynamicalSystemsGraph>
siconos::modeling::NonSmoothDynamicalSystem::dynamicalSystems() const {
  return _topology->dSG(0);
}

std::shared_ptr<siconos::modeling::DynamicalSystem>
siconos::modeling::NonSmoothDynamicalSystem::dynamicalSystem(unsigned int nb) const {
  return _topology->getDynamicalSystem(nb);
}

void siconos::modeling::NonSmoothDynamicalSystem::displayDynamicalSystems() const {
  _topology->displayDynamicalSystems();
}

size_t siconos::modeling::NonSmoothDynamicalSystem::getNumberOfInteractions() const {
  return _topology->indexSet0()->size();
};

const std::shared_ptr<siconos::graphs::InteractionsGraph>
siconos::modeling::NonSmoothDynamicalSystem::interactions() const {
  return _topology->indexSet0();
};

std::shared_ptr<siconos::modeling::Interaction>
siconos::modeling::NonSmoothDynamicalSystem::interaction(unsigned int nb) const {
  return _topology->getInteraction(nb);
}

std::shared_ptr<siconos::modeling::Interaction>
siconos::modeling::NonSmoothDynamicalSystem::interaction(std::string name) const {
  return _topology->getInteraction(name);
}

void siconos::modeling::NonSmoothDynamicalSystem::setName(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds, const std::string& name) {
  _topology->setName(ds, name);
};

std::string siconos::modeling::NonSmoothDynamicalSystem::name(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  return _topology->name(ds);
}

void siconos::modeling::NonSmoothDynamicalSystem::setName(
    std::shared_ptr<siconos::modeling::Interaction> interaction, const std::string& name) {
  _topology->setName(interaction, name);
};

std::string siconos::modeling::NonSmoothDynamicalSystem::name(
    std::shared_ptr<siconos::modeling::Interaction> inter) {
  return _topology->name(inter);
}

void siconos::modeling::NonSmoothDynamicalSystem::setControlProperty(
    std::shared_ptr<siconos::modeling::Interaction> inter, const bool isControlInteraction) {
  _topology->setControlProperty(inter, isControlInteraction);
}

std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>>
siconos::modeling::NonSmoothDynamicalSystem::dynamicalSystemsVector() const {
  std::vector<std::shared_ptr<DynamicalSystem>> dynamicalSystemsVector;
  auto& dsg = *dynamicalSystems();
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsiend;
  std::tie(dsi, dsiend) = dsg.vertices();
  for (; dsi != dsiend; ++dsi) {
    dynamicalSystemsVector.push_back(dsg.bundle(*dsi));
  }

  return dynamicalSystemsVector;
}
std::vector<std::shared_ptr<siconos::modeling::Interaction>>
siconos::modeling::NonSmoothDynamicalSystem::InteractionsVector() const {
  std::vector<std::shared_ptr<siconos::modeling::Interaction>> interactionsVector;
  auto indexSet0 = _topology->indexSet0();
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    interactionsVector.push_back(indexSet0->bundle(*ui));
  }

  return interactionsVector;
}
