/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "ControlLsodarSimulation.hpp"

#include <chrono>

#include "ControlLinearAdditionalTermsED.hpp"
#include "ControlManager.hpp"
#include "Event.hpp"
#include "EventDriven.hpp"
#include "EventsManager.hpp"
#include "LsodarOSI.hpp"

siconos::control::ControlLsodarSimulation::ControlLsodarSimulation(double t0, double T,
                                                                   double h)
    : ControlSimulation(t0, T, h) {
  _processIntegrator = std::make_shared<siconos::integrators::LsodarOSI>();
  _processSimulation =
      std::make_shared<siconos::simulation::EventDriven>(_nsds, _processTD, 0);
  _processSimulation->setName("plant simulation");
  _processSimulation->insertIntegrator(_processIntegrator);
  std::static_pointer_cast<siconos::integrators::LsodarOSI>(_processIntegrator)
      ->setExtraAdditionalTerms(std::make_shared<ControlLinearAdditionalTermsED>());

  _DSG0 = _nsds->topology()->dSG(0);
  _IG0 = _nsds->topology()->indexSet0();

  // Control part
  _CM = std::make_shared<ControlManager>(_processSimulation);
}

void siconos::control::ControlLsodarSimulation::run() {
  auto& eventsManager = *_processSimulation->eventsManager();
  int k = 0;
  auto start = std::chrono::system_clock::now();
  auto& sim = static_cast<siconos::simulation::EventDriven&>(*_processSimulation);

  while (sim.hasNextEvent()) {
    if (eventsManager.needsIntegration()) {
      sim.advanceToEvent();
    }
    sim.processEvents();
    auto& currentEvent = *eventsManager.currentEvent();
    if (currentEvent.getType() == siconos::simulation::EventType::Actuator) {
      // this is necessary since we changed the control input, hence the RHS
      sim.setIstate(1);
    }
    if (sim.hasNextEvent() &&
        eventsManager.nextEvent()->getType() ==
            siconos::simulation::EventType::TD)  // We store only on TD_EVENT, this should be
                                                 // settable
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

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double, std::milli> fp_s = end - start;
  _elapsedTime = fp_s.count();
  /// _dataM->resize(k, _nDim + 1);
}
