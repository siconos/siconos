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
#include "Simulation.hpp"

#include "DynamicalSystem.hpp"
#include "EventsManager.hpp"
#include "Interaction.hpp"
#include "InteractionManager.hpp"
#include "OneStepIntegrator.hpp"
#include "OneStepNSProblem.hpp"
#include "SiconosConst.hpp"  // siconos::internal::LEVELMAX
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "Tools.hpp"  // for enum_to_string
#include "Topology.hpp"
// for Debug
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::simulation::Simulation::Simulation(
    std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
    std::shared_ptr<TimeDiscretisation> td)
    : _nsds{nsds} {
  if (!td) THROW_EXCEPTION("Simulation constructor - timeDiscretisation == nullptr.");
  assert(_nsds);
  _allOSI =
      std::make_shared<std::set<std::shared_ptr<siconos::integrators::OneStepIntegrator>>>();
  _allNSProblems = std::make_shared<
      std::vector<std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>>>();
  _eventsManager = std::make_shared<EventsManager>(td);
  _eventsManager->updateT(_nsds->finalT());
  _nsdsChangeLogPosition = _nsds->changeLog().begin();
}

// --- Destructor ---
siconos::simulation::Simulation::~Simulation() noexcept {
  clear();
  // -> see shared ressources for this
  if (statOut.is_open()) statOut.close();
}

double siconos::simulation::Simulation::getTk() const { return _eventsManager->getTk(); }

double siconos::simulation::Simulation::getTkp1() const { return _eventsManager->getTkp1(); }

double siconos::simulation::Simulation::getTkp2() const { return _eventsManager->getTkp2(); }

double siconos::simulation::Simulation::currentTimeStep() const {
  return _eventsManager->currentTimeStep();
}

double siconos::simulation::Simulation::startingTime() const {
  return _eventsManager->startingTime();
}

double siconos::simulation::Simulation::nextTime() const { return _eventsManager->nextTime(); }

bool siconos::simulation::Simulation::hasNextEvent() const {
  return _eventsManager->hasNextEvent();
}

// clear all maps to break shared_ptr cycle
void siconos::simulation::Simulation::clear() {
  if (_allOSI) {
    _allOSI->clear();
  }
  if (_allNSProblems) {
    _allNSProblems->clear();
  }
}

// Getters/setters

void siconos::simulation::Simulation::insertIntegrator(
    std::shared_ptr<siconos::integrators::OneStepIntegrator> osi) {
  _allOSI->insert(osi);
}

void siconos::simulation::Simulation::associate(
    std::shared_ptr<siconos::integrators::OneStepIntegrator> osi,
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  _allOSI->insert(osi);

  _OSIDSmap[osi].push_back(ds);
}

std::shared_ptr<siconos::graphs::InteractionsGraph> siconos::simulation::Simulation::indexSet(
    unsigned int i) {
  return _nsds->topology()->indexSet(i);
}

std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>
siconos::simulation::Simulation::oneStepNSProblem(int Id) {
  if (!(*_allNSProblems)[Id])
    THROW_EXCEPTION(
        "Simulation - oneStepNSProblem(Id) - The One Step NS Problem is not in "
        "the "
        "simulation.");

  return (*_allNSProblems)[Id];
}

void siconos::simulation::Simulation::updateIndexSets() {
  DEBUG_BEGIN("siconos::simulation::Simulation::updateIndexSets()\n");
  // update I0 indices
  unsigned int nindexsets = _nsds->topology()->indexSetsSize();

  DEBUG_PRINTF("  nindexsets = %d\n", nindexsets);
  if (nindexsets > 1) {
    for (unsigned int i = 1; i < nindexsets; ++i) {
      updateIndexSet(i);
      _nsds->topology()->indexSet(i)->update_vertices_indices();
      _nsds->topology()->indexSet(i)->update_edges_indices();
    }
  }
  DEBUG_END("siconos::simulation::Simulation::updateIndexSets()\n");
}

void siconos::simulation::Simulation::updateDSPlugins(double time) {
  _nsds->updateDSPlugins(time);
}

void siconos::simulation::Simulation::insertNonSmoothProblem(
    std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osns, int Id) {
  if (_allNSProblems->size() > (unsigned int)Id) {
    if ((*_allNSProblems)[Id])
      THROW_EXCEPTION(
          "Simulation - insertNonSmoothProblem(osns), trying to insert a OSNSP "
          "already "
          "existing. ");
    (*_allNSProblems)[Id] = osns;
  } else {
    _allNSProblems->resize(Id + 1);
    (*_allNSProblems)[Id] = osns;
  }
}

void siconos::simulation::Simulation::initializeOSIAssociations() {
  // 1-  OneStepIntegrators initialization ===
  // we set the simulation pointer and the graph of DS in osi for all
  // integrators declared in the set.
  for (auto osi : *_allOSI) {
    if (!osi->isInitialized()) {
      DEBUG_PRINT("- 1 - set simulation pointer  and the graph of ds in osi\n");
      osi->setSimulationPtr(shared_from_this());
      // a subgraph has to be implemented.
      osi->setDynamicalSystemsGraph(_nsds->topology()->dSG(0));
    }
  }

  // 2 - we set the osi of DS that has been defined through associate(ds,osi)
  std::list<std::shared_ptr<siconos::modeling::DynamicalSystem>>::iterator itlist;
  for (auto& osi_it : _OSIDSmap) {
    DEBUG_PRINT(
        "- 2 - we set the osi of DS that has been defined through "
        "associate(ds,osi)\n");
    for (auto& ds : osi_it.second) {
      _nsds->topology()->setOSI(ds, osi_it.first);
    }
    // osi_it.second.clear();  // clear the tmp list
  }
  _OSIDSmap.clear();  // clear the tmp map.
}

void siconos::simulation::Simulation::applyNSDSChangelogForDS() {
  // 4- we initialize new  ds and interaction
  /* Changes to the NSDS are tracked by a changelog, making it fast
   * for the Simulation to scan any changes it has not yet seen and
   * initialize the associated ata structures.  It is just an
   * optimisation over scanning the whole NSDS for new elements at
   * each step. */
  auto DSG = _nsds->topology()->dSG(0);

  auto _nsdsChangeLogPosition_save = _nsdsChangeLogPosition;

  auto& itc = _nsdsChangeLogPosition;

  bool interactionInitialized = false;
  itc++;

  while (itc != _nsds->changeLog().end()) {
    DEBUG_PRINT("- 3 - we initialize new  ds and interaction \n");
    DEBUG_PRINT("The nsds has changed\n");
    const auto& change = *itc;
    itc++;

    DEBUG_EXPR(change.display());
    if (change.typeOfChange ==
        siconos::modeling::NonSmoothDynamicalSystem::ChangeType::addDynamicalSystem) {
      auto ds = change.ds;
      DEBUG_PRINTF("ds number : %i\n", ds->number());
      if (!DSG->properties(DSG->descriptor(ds)).osi) {
        if (_allOSI->size() == 0)
          THROW_EXCEPTION("Simulation::initialize - there is no osi in this Simulation !!");
        DEBUG_PRINTF("_allOSI->size() = %lu\n", _allOSI->size());
        auto osi_default = *_allOSI->begin();
        _nsds->topology()->setOSI(ds, osi_default);
        if (_allOSI->size() > 1) {
          std::cout << "Warning. The simulation has multiple OneStepIntegrators "
                       "(OSI) but the DS number "
                    << ds->number()
                    << " is not assigned to an "
                       "OSI. We assign the following OSI to this DS."
                    << std::endl;
        }
      }
      auto& osi = *DSG->properties(DSG->descriptor(ds)).osi;
      osi.initializeWorkVectorsForDS(getTk(), ds);
    }
    // else if(change.typeOfChange == NonSmoothDynamicalSystem::addInteraction)
    // {
    //   SP::Interaction inter = change.i;
    //   initializeInteraction(getTk(), inter);
    //   interactionInitialized = true;
    // }
    // else if(change.typeOfChange ==
    // NonSmoothDynamicalSystem::rmDynamicalSystem)
    // {
    //   // also need to force an update in this case since indexSet1 may
    //   // still have Interactions that refer to DSs that are not in graph
    //   interactionInitialized = true;
    // }
  }
  _nsdsChangeLogPosition = _nsdsChangeLogPosition_save;

  // _nsdsChangeLogPosition = _nsds->changeLogPosition();

  // // (re)initialize OneStepNSProblem(s) if necessary
  // if(interactionInitialized || !_isInitialized)
  // {
  //   DEBUG_PRINT("(re)Initialize OneStepNSProblem(s)\n");
  //   // Initialize OneStepNSProblem(s). Depends on the type of simulation.
  //   // Warning FP : must be done in any case, even if the interactions set
  //   // is empty.
  //   initializeOneStepNSProblem();

  //   // Since initializeOneStepNSProblem calls updateIndexSets() which resets
  //   the
  //   // topology->hasChanged() flag, it must be specified explicitly.
  //   // Otherwise OneStepNSProblem may fail to update its matrices.
  //   _nsds->topology()->setHasChanged(true);
  // }
}

void siconos::simulation::Simulation::initializeNSDSChangelog() {
  // 4- we initialize new  ds and interaction
  /* Changes to the NSDS are tracked by a changelog, making it fast
   * for the Simulation to scan any changes it has not yet seen and
   * initialize the associated ata structures.  It is just an
   * optimisation over scanning the whole NSDS for new elements at
   * each step. */
  auto DSG = _nsds->topology()->dSG(0);
  auto& itc = _nsdsChangeLogPosition;

  bool interactionInitialized = false;
  itc++;

  while (itc != _nsds->changeLog().end()) {
    DEBUG_PRINT("- 3 - we initialize new  ds and interaction \n");
    DEBUG_PRINT("The nsds has changed\n");
    const auto& change = *itc;
    itc++;

    DEBUG_EXPR(change.display());
    // if(change.typeOfChange == NonSmoothDynamicalSystem::addDynamicalSystem)
    // {
    //   SP::DynamicalSystem ds = change.ds;
    //   if(!DSG->properties(DSG->descriptor(ds)).osi)
    //   {
    //     if(_allOSI->size() == 0)
    //       THROW_EXCEPTION
    //       ("Simulation::initialize - there is no osi in this Simulation !!");
    //     DEBUG_PRINTF("_allOSI->size() = %lu\n", _allOSI->size());
    //     SP::OneStepIntegrator osi_default = *_allOSI->begin();
    //     _nsds->topology()->setOSI(ds, osi_default);
    //     if(_allOSI->size() > 1)
    //     {
    //       std::cout << "Warning. The simulation has multiple
    //       OneStepIntegrators "
    //                 "(OSI) but the DS number " << ds->number() << " is not
    //                 assigned to an " "OSI. We assign the following OSI to
    //                 this DS." << std::endl;
    //     }
    //   }
    //   OneStepIntegrator& osi = *DSG->properties(DSG->descriptor(ds)).osi;
    //   osi.initializeWorkVectorsForDS(getTk(),ds);
    // }
    if (change.typeOfChange ==
        siconos::modeling::NonSmoothDynamicalSystem::ChangeType::addInteraction) {
      auto inter = change.i;
      initializeInteraction(getTk(), inter);
      interactionInitialized = true;
    } else if (change.typeOfChange ==
               siconos::modeling::NonSmoothDynamicalSystem::ChangeType::rmDynamicalSystem) {
      // ---- A ds has been removed from the NDS ? ----
      // also need to force an update in this case since indexSet1 may
      // still have Interactions that refer to DSs that are not in graph
      interactionInitialized = true;
    }
  }
  _nsdsChangeLogPosition = std::prev(_nsds->changeLog().end());

  // (re)initialize OneStepNSProblem(s) if necessary
  if (interactionInitialized || !_isInitialized) {
    DEBUG_PRINT("(re)Initialize OneStepNSProblem(s)\n");
    // Initialize OneStepNSProblem(s). Depends on the type of simulation.
    // Warning FP : must be done in any case, even if the interactions set
    // is empty.
    updateIndexSets();
    initializeOneStepNSProblem();

    // Since updateIndexSets() call resets the
    // topology->hasChanged() flag, it must be specified explicitly.
    // Otherwise OneStepNSProblem may fail to update its matrices.
    _nsds->topology()->setHasChanged(true);
  }
}

void siconos::simulation::Simulation::initializeIndexSets() {
  // 3 - we finalize the initialization of osi

  // symmetry in indexSets Do we need it ?
  _nsds->topology()->setProperties();

  // === OneStepIntegrators initialization ===
  for (auto osi : *_allOSI) {
    if (!osi->isInitialized()) {
      DEBUG_PRINT("- 4 - we finalize the initialization of osi\n");
      DEBUG_PRINT("osi->initialize\n")
      osi->initialize();
      _numberOfIndexSets = std::max<int>(osi->numberOfIndexSets(), _numberOfIndexSets);
    }
  }

  auto topo = _nsds->topology();
  auto indxSize = topo->indexSetsSize();
  assert(_numberOfIndexSets > 0);
  if ((indxSize == siconos::internal::LEVELMAX) || (indxSize < _numberOfIndexSets)) {
    DEBUG_PRINT("Topology : a different number of indexSets has been found \n");
    DEBUG_PRINT("Topology :  we resize the number of index sets \n");
    topo->indexSetsResize(_numberOfIndexSets);
    // Init if the size has changed
    for (auto i = indxSize; i < topo->indexSetsSize(); i++)  // ++i ???
      topo->resetIndexSetPtr(i);                             // Creates the interaction graph
  }
}

void siconos::simulation::Simulation::firstInitialize() {
  if (!_isInitialized) {
    DEBUG_PRINT(" - 6 - First initialization of the simulation\n");
    _T = _nsds->finalT();

    // === Events manager initialization ===
    _eventsManager->initialize(_T);
    _tinit = _eventsManager->startingTime();

    // Process events at time _tinit. Useful to save values in memories
    // for example.  Warning: can not be called during
    // eventsManager->initialize, because it needs the initialization of
    // OSI, OSNS ...
    // _eventsManager->preUpdate(*this);

    _tend = _eventsManager->nextTime();

    // End of initialize:

    //  - all OSI and OSNS (ie DS and Interactions) states are computed
    //  - for time _tinit and saved into memories.
    //  - Sensors or related objects are updated for t=_tinit.
    //  - current time of the model is equal to t1, time of the first
    //  - event after _tinit.
    //  - currentEvent of the simu. corresponds to _tinit and nextEvent
    //  - to _tend.

    _isInitialized = true;
  }
}

void siconos::simulation::Simulation::initialize() {
  DEBUG_BEGIN("siconos::simulation::Simulation::initialize()");
  DEBUG_EXPR_WE(std::cout << "Simulation name :" << name() << std::endl;);

  // 1 - Process any pending OSI->DS associations
  initializeOSIAssociations();

  // 2 - Initialize index sets for OSIs
  initializeIndexSets();

  // 3 - initialize new ds
  applyNSDSChangelogForDS();

  // 4 - update the world from DS
  // for external contact detection library for instance
  updateWorldFromDS();

  // 5 - call the InteractionManager to add/remove interactions
  // Warning: this routine sometimes may update the external objects
  // from the state of ds without using updateWorldFromDS()
  updateInteractions();

  // 6 - initialize new interactions
  initializeNSDSChangelog();

  // 7 - First initialization of the simulation
  firstInitialize();

  DEBUG_END("Simulation::initialize()\n");
}

void siconos::simulation::Simulation::initializeInteraction(
    double time, std::shared_ptr<siconos::modeling::Interaction> inter) {
  DEBUG_BEGIN(
      "siconos::simulation::Simulation::initializeInteraction(double time, "
      "std::shared_ptr<siconos::modeling::Interaction> inter)\n");
  // Get the interaction properties from the topology for initialization.
  auto indexSet0 = _nsds->topology()->indexSet0();
  auto ui = indexSet0->descriptor(inter);

  // This calls computeOutput() and initializes qMemory and q_k.
  auto& DSG = *_nsds->topology()->dSG(0);

  // auto osi = indexSet0->properties(ui).osi;
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds1;
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds2;
  // --- Get the dynamical system(s) (edge(s)) connected to the current
  // interaction (vertex)
  // ---
  if (indexSet0->properties(ui).source != indexSet0->properties(ui).target) {
    DEBUG_PRINT("a two DS Interaction\n");
    ds1 = indexSet0->properties(ui).source;
    ds2 = indexSet0->properties(ui).target;
  } else {
    DEBUG_PRINT("a single DS Interaction\n");
    ds1 = indexSet0->properties(ui).source;
    ds2 = ds1;
    // \warning this looks like some debug code, but it gets executed even with
    // NDEBUG. may be compiler does something smarter, but still it should be
    // rewritten. --xhub
    siconos::graphs::InteractionsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = indexSet0->out_edges(ui); oei != oeiend; ++oei) {
      // note : at most 4 edges
      ds2 = indexSet0->bundle(*oei);
      if (ds2 != ds1) {
        assert(false);
        break;
      }
    }
  }
  assert(ds1);
  assert(ds2);

  auto& osi1 = *DSG.properties(DSG.descriptor(ds1)).osi;
  auto& osi2 = *DSG.properties(DSG.descriptor(ds2)).osi;
  auto& i_prop = indexSet0->properties(ui);
  i_prop.osi1 = DSG.properties(DSG.descriptor(ds1)).osi;
  i_prop.osi2 = DSG.properties(DSG.descriptor(ds2)).osi;

  if (&osi1 == &osi2) {
    osi1.initializeWorkVectorsForInteraction(*inter, i_prop, DSG);
    osi1.updateAndSwapAllOutput(*inter, time);
  } else {
    osi1.initializeWorkVectorsForInteraction(*inter, i_prop, DSG);
    osi1.updateAndSwapAllOutput(*inter, time);
    osi2.initializeWorkVectorsForInteraction(*inter, i_prop, DSG);
    osi2.updateAndSwapAllOutput(*inter, time);
  }
  DEBUG_END(
      "siconos::simulation::Simulation::initializeInteraction(double time, "
      "std::shared_ptr<siconos::modeling::Interaction> inter)\n");

  auto osi1Type = osi1.getType();
  auto osi2Type = osi2.getType();

  // Check consistency of the OneStepIntegrator
  // We assume that the osi of ds1 (osi1) is integrating the interaction (to be
  // reworked for more general case)
  if (osi1Type != osi2Type) {
    THROW_EXCEPTION(
        "Simulation::initializeInteraction: integration of Interaction not yet "
        "implemented "
        "for OSI1 and OSI2 of type " +
        siconos::tools::enum_to_string(osi1Type) + siconos::tools::enum_to_string(osi2Type));
  }
}

int siconos::simulation::Simulation::computeOneStepNSProblem(int Id) {
  DEBUG_BEGIN("siconos::simulation::Simulation::computeOneStepNSProblem(int Id)\n");
  DEBUG_PRINTF("with Id = %i\n", Id);

  if (!(*_allNSProblems)[Id])
    THROW_EXCEPTION(
        "Simulation - computeOneStepNSProblem, OneStepNSProblem == nullptr, "
        "Id: " +
        std::to_string(Id));

  // Before compute, inform all OSNSs if topology has changed
  if (_nsds->topology()->hasChanged()) {
    for (auto osns : *_allNSProblems) {
      osns->setHasBeenUpdated(false);
    }
  }

  int info = (*_allNSProblems)[Id]->compute(nextTime());

  DEBUG_END("siconos::simulation::Simulation::computeOneStepNSProblem(int Id)\n");
  return info;
}

std::shared_ptr<siconos::algebra::SiconosVector> siconos::simulation::Simulation::y(
    unsigned int level, unsigned int coor) {
  // return output(level) (ie with y[level]) for all Interactions.
  // assert(level>=0);

  DEBUG_BEGIN(
      "siconos::simulation::Simulation::output(unsigned int level, unsigned "
      "int coor)\n");
  DEBUG_PRINTF("with level = %i and coor = %i \n", level, coor);

  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  std::shared_ptr<siconos::modeling::Interaction> inter;
  auto indexSet0 = _nsds->topology()->indexSet0();

  auto y = std::make_shared<siconos::algebra::SiconosVector>(
      _nsds->topology()->indexSet0()->size());
  int i = 0;
  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForOutput() <= level);
    assert(inter->upperLevelForOutput() >= level);
    y->setValue(i, inter->y(level)->getValue(coor));
    i++;
  }
  DEBUG_END(
      "siconos::simulation::Simulation::output(unsigned int level, unsigned "
      "int coor)\n");
  return y;
}

std::shared_ptr<siconos::algebra::SiconosVector> siconos::simulation::Simulation::lambda(
    unsigned int level, unsigned int coor) {
  // return input(level) (ie with lambda[level]) for all Interactions.
  // assert(level>=0);

  DEBUG_BEGIN(
      "siconos::simulation::Simulation::input(unsigned int level, unsigned int "
      "coor)\n");
  DEBUG_PRINTF("with level = %i and coor = %i \n", level, coor);

  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  std::shared_ptr<siconos::modeling::Interaction> inter;
  auto indexSet0 = _nsds->topology()->indexSet0();

  auto lambda = std::make_shared<siconos::algebra::SiconosVector>(
      _nsds->topology()->indexSet0()->size());
  int i = 0;
  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForOutput() <= level);
    assert(inter->upperLevelForOutput() >= level);
    lambda->setValue(i, inter->lambda(level)->getValue(coor));
    i++;
  }
  DEBUG_END(
      "siconos::simulation::Simulation::input(unsigned int level, unsigned int "
      "coor)\n");
  return lambda;
}

void siconos::simulation::Simulation::run() {
  unsigned int count = 0;  // events counter.

  std::cout << " ==== Simulation is running - This may take a while ... ====\n";
  while (hasNextEvent()) {
    advanceToEvent();
    processEvents();
    count++;
  }
  std::cout << "===== End of simulation. " << count << " events have been processed. ==== \n";
}

void siconos::simulation::Simulation::processEvents() {
  DEBUG_BEGIN("void siconos::simulation::Simulation::processEvents()\n");
  _eventsManager->processEvents(*this);

  // if(_eventsManager->hasNextEvent())
  // {
  //   // For TimeStepping Scheme, need to update IndexSets, but not for
  //   EventDriven scheme if(Type::value(*this) != Type::EventDriven)
  //   {
  //     updateIndexSets();
  //   }
  // }
  DEBUG_END("void siconos::simulation::Simulation::processEvents()\n");
}

void siconos::simulation::Simulation::clearNSDSChangeLog() {
  _nsds->clearChangeLogTo(_nsdsChangeLogPosition);
}

void siconos::simulation::Simulation::updateT(double T) {
  _T = T;
  _eventsManager->updateT(T);
}

void siconos::simulation::Simulation::link(
    std::shared_ptr<siconos::modeling::Interaction> inter,
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds1,
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds2) {
  DEBUG_PRINTF("link interaction : %d\n", inter->number());

  nonSmoothDynamicalSystem()->link(inter, ds1, ds2);
}

void siconos::simulation::Simulation::unlink(
    std::shared_ptr<siconos::modeling::Interaction> inter) {
  nonSmoothDynamicalSystem()->removeInteraction(inter);
}

void siconos::simulation::Simulation::updateInteractions() {
  // Update interactions if a manager has been provided.  Changes will be
  // detected by siconos::simulation::Simulation::initialize() changelog code.
  if (_interman) _interman->updateInteractions(shared_from_this());
}

void siconos::simulation::Simulation::computeResidu() {
  DEBUG_BEGIN("siconos::simulation::Simulation::computeResidu()\n");
  for (auto osi : *_allOSI) {
    osi->computeResidu();
  }

  DEBUG_END("siconos::simulation::Simulation::computeResidu()\n");
}

void siconos::simulation::Simulation::updateAllInput() {
  DEBUG_BEGIN("Simulation::updateAllInput()\n");
  // nonSmoothDynamicalSystem()->resetNonSmoothPart();

  // 1 - compute input (lambda -> r)
  if (!_allNSProblems->empty()) {
    for (auto osi : *_allOSI) osi->resetAllNonSmoothParts();

    for (auto osi : *_allOSI) osi->updateInput(nextTime());
    //_nsds->updateInput(nextTime(),levelInput);
  }
  DEBUG_END("siconos::simulation::Simulation::updateAllInput()\n");
}

void siconos::simulation::Simulation::updateInput(unsigned int level) {
  DEBUG_BEGIN("Simulation::updateInput()\n");

  // nonSmoothDynamicalSystem()->resetNonSmoothPart(level);

  // 1 - compute input (lambda -> r)
  if (!_allNSProblems->empty()) {
    for (auto osi : *_allOSI) osi->resetNonSmoothPart(level);

    for (auto osi : *_allOSI) osi->updateInput(nextTime(), level);
  }
  DEBUG_END("Simulation::updateInput()\n");
}

void siconos::simulation::Simulation::updateState(unsigned int) {
  DEBUG_BEGIN("siconos::simulation::Simulation::updateState()\n");
  // 2 - compute state for each dynamical system
  for (auto osi : *_allOSI) {
    osi->updateState();
  }

  /*Because the dof of DS have been updated,
    the world (CAO for example) must be updated.*/
  updateWorldFromDS();

  DEBUG_END("siconos::simulation::Simulation::updateState()\n");
}

void siconos::simulation::Simulation::updateOutput(unsigned int) {
  DEBUG_BEGIN("siconos::simulation::Simulation::updateOutput()\n");

  // 3 - compute output ( x ... -> y)
  if (!_allNSProblems->empty()) {
    for (auto osi : *_allOSI) {
      osi->updateOutput(nextTime());
    }
  }
  DEBUG_END("siconos::simulation::Simulation::updateOutput()\n");
}
