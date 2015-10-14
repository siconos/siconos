/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

#include "PID.hpp"

#include "FirstOrderLinearTIDS.hpp"
#include "SiconosVector.hpp"
#include "ActuatorFactory.hpp"
#include "TimeDiscretisation.hpp"
#include "ControlSensor.hpp"
#include "SimpleMatrix.hpp"
#include "Model.hpp"
#include "Simulation.hpp"
#include "EventsManager.hpp"

#include <iostream>

PID::PID(SP::ControlSensor sensor, SP::SimpleMatrix B): Actuator(PID_, sensor), _ref(0), _curDeltaT(0)
{
  _B = B;
}

PID::~PID()
{
}

void PID::initialize(const Model& m)
{
  _u.reset(new SiconosVector(1, 0));
  Actuator::initialize(m);

  _curDeltaT = m.simulation()->eventsManager()->timeDiscretisation()->currentTimeStep(0);

  // initialize _err
  _err.reset(new boost::circular_buffer<double> (3));
  for (unsigned int i = 0; i < 3; ++i)
    (*_err).push_front(0.0);
}

void PID::actuate()
{
  /** \todo We have to distinguish two cases : linear or nonlinear
   *  support the nonlinear case
   */

  // Compute the new error

  (*_err).push_front(_ref - _sensor->y()(0));

  // compute the new control and update it
  (*_u)(0) += ((*_K)(0) + (*_K)(2) / _curDeltaT + (*_K)(1) * _curDeltaT) * (*_err)[0] +
              (-(*_K)(0) - 2 * (*_K)(2) / _curDeltaT) * (*_err)[1] + (*_K)(2) / _curDeltaT * (*_err)[2];
}

void PID::setK(SP::SiconosVector K)
{
  // check dimensions ...
  if (K->size() != 3)
  {
    RuntimeException::selfThrow("PID::setK - the size of K should be 3");
  }
  else
  {
    _K = K;
  }
}

void PID::setTimeDiscretisation(const TimeDiscretisation& td)
{
}

void PID::display() const
{
  Actuator::display();
  std::cout << "current error vector: ";
  std::cout << (*_err)[0] << " "  << (*_err)[1] << " " << (*_err)[2] << std::endl;
}
AUTO_REGISTER_ACTUATOR(PID_, PID)
