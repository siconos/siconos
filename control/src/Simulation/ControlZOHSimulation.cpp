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

#include "TimeStepping.hpp"
#include "ZeroOrderHoldOSI.hpp"
#include "Model.hpp"
#include "EventsManager.hpp"
#include "Event.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "SimulationGraphs.hpp"

#include "ControlManager.hpp"
#include "ControlZOHAdditionalTerms.hpp"
#include "ControlZOHSimulation.hpp"
#include "ControlSimulation_impl.hpp"

#include <boost/progress.hpp>
#include <boost/timer.hpp>

ControlZOHSimulation::ControlZOHSimulation(double t0, double T, double h):
  ControlSimulation(t0, T, h)
{
  _processIntegrator.reset(new ZeroOrderHoldOSI());
  std11::static_pointer_cast<ZeroOrderHoldOSI>(_processIntegrator)->setExtraAdditionalTerms(
      std11::shared_ptr<ControlZOHAdditionalTerms>(new ControlZOHAdditionalTerms()));
  _processSimulation.reset(new TimeStepping(_processTD, 0));
  _processSimulation->setName("plant simulation");
  _processSimulation->insertIntegrator(_processIntegrator);

  _processSimulation->setNonSmoothDynamicalSystemPtr(
    _model->nonSmoothDynamicalSystem());

  _DSG0 = _model->nonSmoothDynamicalSystem()->topology()->dSG(0);
  _IG0 = _model->nonSmoothDynamicalSystem()->topology()->indexSet0();

  // Control part
  _CM.reset(new ControlManager(_processSimulation));
}

void ControlZOHSimulation::run()
{
  EventsManager& eventsManager = *_processSimulation->eventsManager();
  unsigned k = 0;
  boost::progress_display show_progress(_N);
  boost::timer time;
  time.restart();

  TimeStepping& sim = static_cast<TimeStepping&>(*_processSimulation);

  while (sim.hasNextEvent())
  {
    Event& nextEvent = *eventsManager.nextEvent();
    if (nextEvent.getType() == TD_EVENT)
    {
      sim.computeOneStep();
    }

    sim.nextStep();

    if (sim.hasNextEvent() && eventsManager.nextEvent()->getType() == TD_EVENT)  // We store only on TD_EVENT
    {
      (*_dataM)(k, 0) = sim.startingTime();
      storeData(k);
      ++k;
      if (!_silent)
      {
        ++show_progress;
      }
    }
  }

  /* saves last status */
  (*_dataM)(k, 0) = sim.startingTime();
  storeData(k);
  ++k;

  _elapsedTime = time.elapsed();
  _dataM->resize(k, _nDim + 1);
}
