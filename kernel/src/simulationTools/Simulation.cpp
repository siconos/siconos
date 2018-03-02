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
#include "Simulation.hpp"
#include "EventDriven.hpp"
#include "DynamicalSystem.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Topology.hpp"
#include "Interaction.hpp"
#include "Relation.hpp"
#include "EventsManager.hpp"
#include "LagrangianDS.hpp"

// One Step Integrators
#include "OneStepIntegrator.hpp"

// One Step Non Smooth Problems
#include "LCP.hpp"
#include "QP.hpp"
#include "Relay.hpp"
#include "NonSmoothLaw.hpp"
#include "TypeName.hpp"
// for Debug
//#define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"
#include <fstream>

// --- Constructor with a TimeDiscretisation (and thus a NonSmoothDynamicalSystem) and an
// --- id ---
Simulation::Simulation(SP::NonSmoothDynamicalSystem nsds, SP::TimeDiscretisation td):
  _name("unnamed"), _tinit(0.0), _tend(0.0), _tout(0.0),
  _nsds(nsds),
  _numberOfIndexSets(0),
  _tolerance(DEFAULT_TOLERANCE), _printStat(false),
  _staticLevels(false),_isInitialized(false)
{
  if (!td)
    RuntimeException::selfThrow("Simulation constructor - timeDiscretisation == NULL.");
  _useRelativeConvergenceCriterion = false;
  _relativeConvergenceCriterionHeld = false;
  _relativeConvergenceTol = 10e-3;

  // === indexSets will be updated during initialize() call ===

  _allOSI.reset(new OSISet());
  _allNSProblems.reset(new OneStepNSProblems());
  _eventsManager.reset(new EventsManager(td)); //
  _eventsManager->updateT(_nsds->finalT());

  _nsdsChangeLogPosition = nsds->changeLog().begin();
}



// --- Constructor with a TimeDiscretisation (and thus a Model) and an
// --- id ---
Simulation::Simulation(SP::TimeDiscretisation td):
  _name("unnamed"), _tinit(0.0), _tend(0.0), _tout(0.0),
  _numberOfIndexSets(0),
  _tolerance(DEFAULT_TOLERANCE), _printStat(false),
  _staticLevels(false),_isInitialized(false)
{
  if (!td)
    RuntimeException::selfThrow("Simulation constructor - timeDiscretisation == NULL.");
  _useRelativeConvergenceCriterion = false;
  _relativeConvergenceCriterionHeld = false;
  _relativeConvergenceTol = 10e-3;

  // === indexSets will be updated during initialize() call ===

  _allOSI.reset(new OSISet());
  _allNSProblems.reset(new OneStepNSProblems());
  _eventsManager.reset(new EventsManager(td)); //
  _eventsManager->updateT(_nsds->finalT());
}

// --- Destructor ---
Simulation::~Simulation()
{
  clear();
  // -> see shared ressources for this
  if (statOut.is_open()) statOut.close();
}

double Simulation::getTk() const
{
  return _eventsManager->getTk();
}

double Simulation::getTkp1() const
{
  return _eventsManager->getTkp1();
}

double Simulation::getTkp2() const
{
  return _eventsManager->getTkp2();
}

double Simulation::currentTimeStep() const
{
  return _eventsManager->currentTimeStep();
}

double Simulation::startingTime() const
{
  return _eventsManager->startingTime();
}

double Simulation::nextTime() const
{
  return _eventsManager->nextTime();
}

bool Simulation::hasNextEvent() const
{
  return _eventsManager->hasNextEvent();
}


// clear all maps to break shared_ptr cycle
void Simulation::clear()
{
  if (_allOSI)
  {
    _allOSI->clear();
  }
  if (_allNSProblems)
  {
    _allNSProblems->clear();
  }
}

// Getters/setters

void Simulation::insertIntegrator(SP::OneStepIntegrator osi)
{
  _allOSI->insert(osi);
}

void Simulation::associate(SP::OneStepIntegrator osi, SP::DynamicalSystem ds)
{
  _allOSI->insert(osi);

  _OSIDSmap[osi].push_back(ds);

}

SP::InteractionsGraph Simulation::indexSet(unsigned int i)
{
  return _nsds->topology()->indexSet(i) ;
}

SP::OneStepNSProblem Simulation::oneStepNSProblem(int Id)
{
  if (!(*_allNSProblems)[Id])
    RuntimeException::selfThrow("Simulation - oneStepNSProblem(Id) - The One Step NS Problem is not in the simulation.");

  return (*_allNSProblems)[Id];
}

void Simulation::updateIndexSets()
{

  DEBUG_BEGIN("Simulation::updateIndexSets()\n");
  // update I0 indices
  unsigned int nindexsets = _nsds->topology()->indexSetsSize();

  DEBUG_PRINTF("  nindexsets = %d\n", nindexsets   );
  if (nindexsets > 1)
  {
    for (unsigned int i = 1; i < nindexsets ; ++i)
    {
      updateIndexSet(i);
      _nsds->topology()->indexSet(i)->update_vertices_indices();
      _nsds->topology()->indexSet(i)->update_edges_indices();
    }
  }
  DEBUG_END("Simulation::updateIndexSets()\n");

}

void Simulation::insertNonSmoothProblem(SP::OneStepNSProblem osns, int Id)
{
  if (_allNSProblems->size() > (unsigned int)Id)
  {
    if ((*_allNSProblems)[Id])
      RuntimeException::selfThrow("Simulation - insertNonSmoothProblem(osns), trying to insert a OSNSP already existing. ");
    (*_allNSProblems)[Id] = osns;
  }
  else
  {
    _allNSProblems->resize(Id+1);
    (*_allNSProblems)[Id] = osns;
  }
}

void Simulation::initializeOSIAssociations()
{
  // 1-  OneStepIntegrators initialization ===
  // we set the simulation pointer and the graph of DS in osi
  for (OSIIterator itosi = _allOSI->begin();
       itosi != _allOSI->end(); ++itosi)
  {
    if (!(*itosi)->isInitialized()){
      DEBUG_PRINT("- 1 - set simulation pointer  and the graph of ds in osi\n");
      (*itosi)->setSimulationPtr(shared_from_this());
      // a subgraph has to be implemented.
      (*itosi)->setDynamicalSystemsGraph(_nsds->topology()->dSG(0));
    }
  }

  // 2 - we set the osi of DS that has been defined through associate(ds,osi)
  std::map< SP::OneStepIntegrator, std::list<SP::DynamicalSystem> >::iterator  it;
  std::list<SP::DynamicalSystem> ::iterator  itlist;
  for ( it = _OSIDSmap.begin();  it !=_OSIDSmap.end(); ++it)
  {
    DEBUG_PRINT("- 2 - we set the osi of DS that has been defined through associate(ds,osi)\n");
    for ( itlist = it->second.begin();  itlist !=it->second.end(); ++itlist)
    {
      SP::DynamicalSystem ds =  *itlist;
      SP::OneStepIntegrator osi =it->first;

      _nsds->topology()->setOSI( ds , osi);
    }
    it->second.clear();
  }
}

bool Simulation::initializeNSDSChangelog()
{
  // 3- we initialize new  ds and interaction
  /* Changes to the NSDS are tracked by a changelog, making it fast
   * for the Simulation to scan any changes it has not yet seen and
   * initialize the associated ata structures.  It is just an
   * optimisation over scanning the whole NSDS for new elements at
   * each step. */
  SP::DynamicalSystemsGraph DSG = _nsds->topology()->dSG(0);
  NonSmoothDynamicalSystem::ChangeLogIter itc = _nsdsChangeLogPosition;

  bool interactionInitialized = false;
  itc++;
  while(itc != _nsds->changeLog().end())
  {
    DEBUG_PRINT("- 3 - we initialize new  ds and interaction \n");
    DEBUG_PRINT("The nsds has changed\n")
    const NonSmoothDynamicalSystem::Changes& changes = *itc;
    itc++;

    DEBUG_EXPR(changes.display());
    if (changes.typeOfChange == NonSmoothDynamicalSystem::addDynamicalSystem)
    {
      SP::DynamicalSystem ds = changes.ds;
      if (!DSG->properties(DSG->descriptor(ds)).osi)
      {
        if (_allOSI->size() == 0)
        RuntimeException::selfThrow
            ("Simulation::initialize - there is no osi in this Simulation !!");
        DEBUG_PRINTF("_allOSI->size() = %lu\n", _allOSI->size());
        SP::OneStepIntegrator osi_default = *_allOSI->begin();
        _nsds->topology()->setOSI(ds, osi_default);
        if (_allOSI->size() > 1)
        {
          std::cout << "Warning. The simulation has multiple OneStepIntegrators "
            "(OSI) but the DS number " << ds->number() << " is not assigned to an "
            "OSI. We assign the following OSI to this DS." << std::endl;
        }
      }
      OneStepIntegrator& osi = *DSG->properties(DSG->descriptor(ds)).osi;
      osi.initializeWorkVectorsForDS(getTk(),ds);
    }
    else if (changes.typeOfChange == NonSmoothDynamicalSystem::addInteraction)
    {
      SP::Interaction inter = changes.i;
      initializeInteraction(getTk(), inter);
      interactionInitialized = true;
    }
  }
  _nsdsChangeLogPosition = _nsds->changeLogPosition();

  return interactionInitialized;
}

void Simulation::initializeIndexSets()
{
  // 4 - we finalize the initialization of osi

  // symmetry in indexSets Do we need it ?
  _nsds->topology()->setProperties();

  // === OneStepIntegrators initialization ===
  for (OSIIterator itosi = _allOSI->begin();
       itosi != _allOSI->end(); ++itosi)
  {
    if (!(*itosi)->isInitialized()){
      DEBUG_PRINT("- 4 - we finalize the initialization of osi\n");
      DEBUG_PRINT("osi->initialize\n")
      (*itosi)->initialize();
      _numberOfIndexSets = std::max<int>((*itosi)->numberOfIndexSets(), _numberOfIndexSets);
    }
  }

  SP::Topology topo = _nsds->topology();
  unsigned int indxSize = topo->indexSetsSize();
  assert (_numberOfIndexSets >0);
  if ((indxSize == LEVELMAX) || (indxSize < _numberOfIndexSets ))
  {
    DEBUG_PRINT("Topology : a different number of indexSets has been found \n");
    DEBUG_PRINT("Topology :  we resize the number of index sets \n");
    topo->indexSetsResize(_numberOfIndexSets);
    // Init if the size has changed
    for (unsigned int i = indxSize; i < topo->indexSetsSize(); i++) // ++i ???
      topo->resetIndexSetPtr(i);
  }
}

void Simulation::initialize()
{
  DEBUG_BEGIN("Simulation::initialize()");
  DEBUG_EXPR_WE(std::cout << "Simulation name :"<< name() << std::endl;);

  // 1 - Process any pending OSI->DS associations
  initializeOSIAssociations();

  // 2 - allow the InteractionManager to add/remove any interactions it wants
  updateWorldFromDS();
  updateInteractions();

  // 3 - initialize new ds and interactions
  bool interactionInitialized =
    initializeNSDSChangelog();

  // 4 - Initialize index sets for OSIs
  initializeIndexSets();

  // 5 - Initialize OneStepNSProblem(s)
  if (interactionInitialized || !_isInitialized)
  {
    DEBUG_PRINT( "- 5 -Initialize OneStepNSProblem(s)\n");
    // Initialize OneStepNSProblem(s). Depends on the type of simulation.
    // Warning FP : must be done in any case, even if the interactions set
    // is empty.
    initOSNS();

    // Since initOSNS calls updateIndexSets() which resets the
    // topology->hasChanged() flag, it must be specified explicitly.
    // Otherwise OneStepNSProblem may fail to update its matrices.
    _nsds->topology()->setHasChanged(true);
  }


  // 6 - First initialization of the simulation
  if(!_isInitialized)
  {
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

    _tend =  _eventsManager->nextTime();

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


  DEBUG_END("Simulation::initialize()\n");
}

void Simulation::initializeInteraction(double time, SP::Interaction inter)
{
  DEBUG_BEGIN("Simulation::initializeInteraction(double time, SP::Interaction inter)\n");
  // Get the interaction properties from the topology for initialization.
  SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet0();
  InteractionsGraph::VDescriptor ui = indexSet0->descriptor(inter);

  // This calls computeOutput() and initializes qMemory and q_k.
  DynamicalSystemsGraph &DSG = *_nsds->topology()->dSG(0);

  //SP::OneStepIntegrator osi = indexSet0->properties(ui).osi;
  SP::DynamicalSystem ds1;
  SP::DynamicalSystem ds2;
  // --- Get the dynamical system(s) (edge(s)) connected to the current interaction (vertex) ---
  if (indexSet0->properties(ui).source != indexSet0->properties(ui).target)
  {
    DEBUG_PRINT("a two DS Interaction\n");
    ds1 = indexSet0->properties(ui).source;
    ds2 = indexSet0->properties(ui).target;
  }
  else
  {
    DEBUG_PRINT("a single DS Interaction\n");
    ds1 = indexSet0->properties(ui).source;
    ds2 = ds1;
    // \warning this looks like some debug code, but it gets executed even with NDEBUG.
    // may be compiler does something smarter, but still it should be rewritten. --xhub
    InteractionsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = indexSet0->out_edges(ui);
         oei != oeiend; ++oei)
    {
      // note : at most 4 edges
      ds2 = indexSet0->bundle(*oei);
      if (ds2 != ds1)
      {
        assert(false);
        break;
      }
    }
  }
  assert(ds1);
  assert(ds2);

  OneStepIntegrator& osi1 = *DSG.properties(DSG.descriptor(ds1)).osi;
  OneStepIntegrator& osi2 = *DSG.properties(DSG.descriptor(ds2)).osi;

  InteractionProperties& i_prop = indexSet0->properties(ui);
  if (&osi1 == &osi2 )
    {
      osi1.initializeWorkVectorsForInteraction(*inter, i_prop,  DSG);
      osi1.update_interaction_output(*inter, time, i_prop);
    }
  else
    {
      osi1.initializeWorkVectorsForInteraction(*inter, i_prop,  DSG);
      osi1.update_interaction_output(*inter, time, i_prop);
      osi2.initializeWorkVectorsForInteraction(*inter, i_prop,  DSG);
      osi2.update_interaction_output(*inter, time, i_prop);
    }
  DEBUG_END("Simulation::initializeInteraction(double time, SP::Interaction inter)\n");
}



int Simulation::computeOneStepNSProblem(int Id)
{
  DEBUG_BEGIN("Simulation::computeOneStepNSProblem(int Id)\n");
  DEBUG_PRINTF("with Id = %i\n", Id);

  if (!(*_allNSProblems)[Id])
    RuntimeException::selfThrow("Simulation - computeOneStepNSProblem, OneStepNSProblem == NULL, Id: " + Id);

  // Before compute, inform all OSNSs if topology has changed
  if (_nsds->topology()->hasChanged())
  {
    for (OSNSIterator itOsns = _allNSProblems->begin();
         itOsns != _allNSProblems->end(); ++itOsns)
    {
      (*itOsns)->setHasBeenUpdated(false);
    }
  }

  int info = (*_allNSProblems)[Id]->compute(nextTime());
  
  DEBUG_END("Simulation::computeOneStepNSProblem(int Id)\n");
  return info;
}


SP::SiconosVector Simulation::y(unsigned int level, unsigned int coor)
{
  // return output(level) (ie with y[level]) for all Interactions.
  // assert(level>=0);

  DEBUG_BEGIN("Simulation::output(unsigned int level, unsigned int coor)\n");
  DEBUG_PRINTF("with level = %i and coor = %i \n", level,coor);

  InteractionsGraph::VIterator ui, uiend;
  SP::Interaction inter;
  SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet0();

  SP::SiconosVector y (new SiconosVector (_nsds->topology()->indexSet0()->size() ));
  int i=0;
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForOutput() <= level);
    assert(inter->upperLevelForOutput() >= level);
    y->setValue(i,inter->y(level)->getValue(coor));
    i++;
  }
  DEBUG_END("Simulation::output(unsigned int level, unsigned int coor)\n");
  return y;
}

SP::SiconosVector Simulation::lambda(unsigned int level, unsigned int coor)
{
  // return input(level) (ie with lambda[level]) for all Interactions.
  // assert(level>=0);

  DEBUG_BEGIN("Simulation::input(unsigned int level, unsigned int coor)\n");
  DEBUG_PRINTF("with level = %i and coor = %i \n", level,coor);

  InteractionsGraph::VIterator ui, uiend;
  SP::Interaction inter;
  SP::InteractionsGraph indexSet0 = _nsds->topology()->indexSet0();

  SP::SiconosVector lambda (new SiconosVector (_nsds->topology()->indexSet0()->size() ));
  int i=0;
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet0->bundle(*ui);
    assert(inter->lowerLevelForOutput() <= level);
    assert(inter->upperLevelForOutput() >= level);
    lambda->setValue(i,inter->lambda(level)->getValue(coor));
    i++;
  }
  DEBUG_END("Simulation::input(unsigned int level, unsigned int coor)\n");
  return lambda;
}


void Simulation::run()
{
  unsigned int count = 0; // events counter.

  std::cout << " ==== Start of " << Type::name(*this) << " simulation - This may take a while ... ====" <<std::endl;
  while (hasNextEvent())
  {
    advanceToEvent();
    processEvents();
    count++;
  }
  std::cout << "===== End of " << Type::name(*this) << " simulation. " << count << " events have been processed. ==== " <<std::endl;
}

void Simulation::processEvents()
{
  DEBUG_BEGIN("void Simulation::processEvents()\n");
  _eventsManager->processEvents(*this);

  if (_eventsManager->hasNextEvent())
  {
    // For TimeStepping Scheme, need to update IndexSets, but not for EventDriven scheme
    if (Type::value(*this) != Type::EventDriven)
    {
      updateIndexSets();
    }
  }
  DEBUG_END("void Simulation::processEvents()\n");
}

void Simulation::clearNSDSChangeLog()
{
  _nsds->clearChangeLogTo(_nsdsChangeLogPosition);
}

void Simulation::updateT(double T)
{
  _T = T;
  _eventsManager->updateT(T);
}


void Simulation::link(SP::Interaction inter,
                      SP::DynamicalSystem ds1,
                      SP::DynamicalSystem ds2)
{
  DEBUG_PRINTF("link interaction : %d\n", inter->number());

  nonSmoothDynamicalSystem()->link(inter, ds1, ds2);
}

void Simulation::unlink(SP::Interaction inter)
{
  nonSmoothDynamicalSystem()->removeInteraction(inter);
}

void Simulation::updateInteractions()
{
  // Update interactions if a manager was provided.  Changes will be
  // detected by Simulation::initialize() changelog code.
  if (_interman)
    _interman->updateInteractions(shared_from_this());
}

void Simulation::updateInput(unsigned int)
{
  DEBUG_BEGIN("Simulation::updateInput()\n");
  OSIIterator itOSI;
  // 1 - compute input (lambda -> r)
  if (!_allNSProblems->empty())
  {
    for (itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
      (*itOSI)->updateInput(nextTime());
    //_nsds->updateInput(nextTime(),levelInput);
  }
  DEBUG_END("Simulation::updateInput()\n");
}

void Simulation::updateState(unsigned int)
{
  DEBUG_BEGIN("Simulation::updateState()\n");
  OSIIterator itOSI;
  // 2 - compute state for each dynamical system
  for (itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
    (*itOSI)->updateState();
  /*Because the dof of DS have been updated,
    the world (CAO for example) must be updated.*/
  updateWorldFromDS();

  DEBUG_END("Simulation::updateState()\n");
}

void Simulation::updateOutput(unsigned int)
{
  DEBUG_BEGIN("Simulation::updateOutput()\n");

  // 3 - compute output ( x ... -> y)
  if (!_allNSProblems->empty())
  {
    OSIIterator itOSI;
    for (itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
      (*itOSI)->updateOutput(nextTime());
  }
  DEBUG_END("Simulation::updateOutput()\n");
}

// void Simulation::prepareIntegratorForDS(SP::OneStepIntegrator osi,
//                                         SP::DynamicalSystem ds,
//                                         SP::Model m, double time)
// {
//   assert(m && m->nonSmoothDynamicalSystem() && "Simulation::prepareIntegratorForDS requires a Model with an NSDS.");

//   /*
//    * Steps to be accomplished when adding a DS to a Model and
//    * Simulation:
//    *
//    * 1. Add the DS to model->_nsds (Model::insertDynamicalSystem(ds))
//    *    (assumed done before this function is called, everything else
//    *    done in this function)
//    *
//    * 2. Add the OSI to simulation->_allOSI (Simulation::insertIntegrator)
//    *
//    * 3. Assign the OSI to the DS via the pointer in
//    *   _nsds->_topology->_DSG properties for the DS (setOSI).  Since
//    *   _nsds is not necessarily available yet, so take it from Model.
//    *
//    * 4. If Simulation already initialized, then DS work vectors in
//    *    _dynamicalSystemsGraph properties for the DS must be
//    *    initialized (OSI::initializeWorkVectorsForDS), otherwise it will
//    *    be called later during Simulation::initialize().
//   */

//   // Keep OSI in the set, no effect if already present.
//   insertIntegrator(osi);

//   // Associate the OSI to the DS in the topology.
//   m->nonSmoothDynamicalSystem()->topology()->setOSI(ds, osi);

//   // Prepare work vectors, etc.
//   // If OSI has no DSG yet, assume DS will be initialized later.
//   // (Typically, during Simulation::initialize())
//   if (osi->dynamicalSystemsGraph())
//     osi->initializeWorkVectorsForDS(*m, time, ds);
// }
