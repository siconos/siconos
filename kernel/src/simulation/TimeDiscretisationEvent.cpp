// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "TimeDiscretisationEvent.hpp"

#include "Simulation.hpp"
#include "TimeDiscretisation.hpp"

void siconos::simulation::TimeDiscretisationEvent::process(Simulation& simulation) {
  // Update y[i] values in Interactions with new DS states.
  // simulation->updateOutput(0, 1);
  // Save state(s) in Memories (DS and Interactions, through OSI and OSNS).
  if (noSaveInMemory_) return;

  simulation.nonSmoothDynamicalSystem()->swapInMemory();  // To save pre-impact values
  simulation.nonSmoothDynamicalSystem()
      ->pushInteractionsInMemory();  // To save pre-impact values
}

void siconos::simulation::TimeDiscretisationEvent::update(unsigned int k) {
  assert(k > _k &&
         "TimeDiscretisationEvent::update - next step has to be greater than the current one");
  if (_td)  // if no TimeDiscretisation, then do nothing
  {
    if (_td->hGmp())
      incrementTime(k - _k);
    else
      setTime(_td->getTk(k));

    _k = k;
  }
}
