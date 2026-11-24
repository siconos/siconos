/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "ControlZOHSimulation.hpp"

#include <chrono>

#include "ControlManager.hpp"
#include "ControlZOHAdditionalTerms.hpp"
#include "Event.hpp"
#include "EventsManager.hpp"
#include "TimeStepping.hpp"
#include "ZeroOrderHoldOSI.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "SiconosException.hpp"
#include "siconos_debug.h"

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
  int k = 0;
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

  DEBUG_END("void siconos::control::ControlZOHSimulation::run()\n");
}
