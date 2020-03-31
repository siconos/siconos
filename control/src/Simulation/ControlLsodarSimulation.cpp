/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include "EventDriven.hpp"
#include "LsodarOSI.hpp"
#include "EventsManager.hpp"
#include "Event.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "SimulationGraphs.hpp"

#include "ControlLinearAdditionalTermsED.hpp"
#include "ControlManager.hpp"
#include "ControlLsodarSimulation.hpp"
#include "ControlSimulation_impl.hpp"
#include "Observer.hpp"
#include "Actuator.hpp"

#include <boost/timer/timer.hpp>


ControlLsodarSimulation::ControlLsodarSimulation(double t0, double T, double h):
  ControlSimulation(t0, T, h)
{
  _processIntegrator.reset(new LsodarOSI());
  _processSimulation.reset(new EventDriven(_nsds, _processTD, 0));
  _processSimulation->setName("plant simulation");
  _processSimulation->insertIntegrator(_processIntegrator);
  std::static_pointer_cast<LsodarOSI>(_processIntegrator)->setExtraAdditionalTerms(
    std::shared_ptr<ControlLinearAdditionalTermsED>(new ControlLinearAdditionalTermsED()));

  _DSG0 = _nsds->topology()->dSG(0);
  _IG0 = _nsds->topology()->indexSet0();

  // Control part
  _CM.reset(new ControlManager(_processSimulation));
}

void ControlLsodarSimulation::run()
{
  EventsManager& eventsManager = *_processSimulation->eventsManager();
  unsigned k = 0;
  boost::timer::cpu_timer time;
  EventDriven& sim = static_cast<EventDriven&>(*_processSimulation);

  while(sim.hasNextEvent())
  {
    if(eventsManager.needsIntegration())
    {
      sim.advanceToEvent();
    }
    sim.processEvents();
    Event& currentEvent = *eventsManager.currentEvent();
    if(currentEvent.getType() == ACTUATOR_EVENT)
    {
      // this is necessary since we changed the control input, hence the RHS
      sim.setIstate(1);
    }
    if(sim.hasNextEvent() && eventsManager.nextEvent()->getType() == TD_EVENT)  // We store only on TD_EVENT, this should be settable
    {
      (*_dataM)(k, 0) = sim.startingTime();
      storeData(k);
      ++k;
    }
  }

  /* saves last status */
  (*_dataM)(k, 0) = sim.startingTime();
  storeData(k);
  ++k;

  // Warning FP : with the new interface boost::timer, the
  // result of elpased is in ns while it was in seconds
  // in the old interface.
  // elapsed is a tuple with wall, user and system times.
  _elapsedTime = time.elapsed().user * 1e-9;
  _dataM->resize(k, _nDim + 1);
}
