/* Siconos-Kernel, Copyright INRIA 2005-2013.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file ControlSensor.cpp
 * A generic control sensor
*/

#include "ControlSensor.hpp"
#include <cmath>
#include "TimeDiscretisation.hpp"
#include "SiconosVector.hpp"

void ControlSensor::initialize(const Model& m)
{
  Sensor::initialize(m);
  if (_delay > 0)
  {
    if (_timeDiscretisation->getTDCase() != 2)
    {
       RuntimeException::selfThrow("ControlSensor::initialize the timediscretization should be of type 2");
    }
    else
    {
       double h = _timeDiscretisation->currentTimeStep();
       double shift = _delay;
       if (_delay >= h)
       {
         unsigned int size = ceil(_delay/h);
         shift = fmod(_delay, h);
         _bufferY.resize(size);
       }
      _timeDiscretisation->setT0(_timeDiscretisation->currentTime() + h - shift);
    }
  }
  else if (_delay < 0)
    RuntimeException::selfThrow("ControlSensor::initialize the delay value should be >= 0");
}

unsigned int ControlSensor::getYDim() const
{
  return _storedY->size();
}
