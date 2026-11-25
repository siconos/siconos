// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file ControlSensor.cpp
 * A generic control sensor
 */

#include "ControlSensor.hpp"

void siconos::control::ControlSensor::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) {
  Sensor::initialize(nsds);
  //  if (_delay > 0)
  //  {
  //    if (_timeDiscretisation->getTDCase() != 2)
  //    {
  //       THROW_EXCEPTION("ControlSensor::initialize the timediscretization should be of type
  //       2");
  //    }
  //    else
  //    {
  //       double h = _timeDiscretisation->currentTimeStep();
  //       double shift = _delay;
  //       if (_delay >= h)
  //       {
  //         size_t size = ceil(_delay/h);
  //         shift = fmod(_delay, h);
  //         _bufferY.resize(size);
  //       }
  //      _timeDiscretisation->setT0(_timeDiscretisation->currentTime() + h - shift);
  //    }
  //  }
  //  else if (_delay < 0)
  //    THROW_EXCEPTION("ControlSensor::initialize the delay value should be >= 0");
}

auto siconos::control::ControlSensor::getYDim() const { return _storedY->size(); }
