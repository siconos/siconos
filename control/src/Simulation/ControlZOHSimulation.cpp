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

#include "ControlZOHSimulation.hpp"

#include <chrono>

#include "ControlManager.hpp"
#include "ControlZOHAdditionalTerms.hpp"
#include "Event.hpp"
#include "EventsManager.hpp"
#include "SiconosMatrix.hpp"
#include "TimeStepping.hpp"
#include "Topology.hpp"  //#define DEBUG_BEGIN_END_ONLY
#include "ZeroOrderHoldOSI.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"
#include "SiconosException.hpp"

siconos::control::ControlZOHSimulation::ControlZOHSimulation(double t0, double T, double h)
    : ControlSimulation(t0, T, h) {
  _processIntegrator = std::make_shared<siconos::integrators::ZeroOrderHoldOSI>();
  std::static_pointer_cast<siconos::integrators::ZeroOrderHoldOSI>(_processIntegrator)
      ->setExtraAdditionalTerms(std::make_shared<ControlZOHAdditionalTerms>());
  _processSimulation =
      std::make_shared<siconos::simulation::TimeStepping>(_nsds, _processTD, 0);
  _processSimulation->setName("plant simulation");
  _processSimulation->insertIntegrator(_processIntegrator);

  _DSG0 = _nsds->topology()->dSG(0);
  _IG0 = _nsds->topology()->indexSet0();

  // Control part
  _CM = std::make_shared<ControlManager>(_processSimulation);
}

void siconos::control::ControlZOHSimulation::run() {
  DEBUG_BEGIN("void siconos::control::ControlZOHSimulation::run()\n");
  auto& eventsManager = *_processSimulation->eventsManager();
  unsigned k = 0;
  auto start = std::chrono::system_clock::now();

  
  auto& sim = static_cast<siconos::simulation::TimeStepping&>(*_processSimulation);
  try {
    while (sim.hasNextEvent()) {
      auto& nextEvent = *eventsManager.nextEvent();
      if (nextEvent.getType() == siconos::simulation::EventType::TD) {
        sim.computeOneStep();
      }

      sim.nextStep();
      if (sim.hasNextEvent() &&
          eventsManager.nextEvent()->getType() ==
              siconos::simulation::EventType::TD)  // We store only on TD_EVENT
      {
        (*_dataM)(k, 0) = sim.startingTime();
        storeData(k);
        ++k;
      }
    }
  } catch (...) {
    siconos::exception::process();
  }
  /* saves last status */
  (*_dataM)(k, 0) = sim.startingTime();
  storeData(k);
  ++k;

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double, std::milli> fp_s = end - start;
  _elapsedTime = fp_s.count();

  _dataM->resize(k, _nDim + 1);
  DEBUG_END("void siconos::control::ControlZOHSimulation::run()\n");
}
