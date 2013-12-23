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

#include "FirstOrderLinearTIDS.hpp"
#include "TimeStepping.hpp"
#include "Relay.hpp"
#include "EventsManager.hpp"

#include "LinearSMC.hpp"

#include "SiconosVector.hpp"
#include "ControlSensor.hpp"
#include "ZeroOrderHoldOSI.hpp"
#include "TimeDiscretisation.hpp"
#include "ActuatorFactory.hpp"

LinearSMC::LinearSMC(SP::ControlSensor sensor, unsigned int type):
  CommonSMC(type, sensor)
{
}

LinearSMC::LinearSMC(SP::ControlSensor sensor, SP::SiconosMatrix B, SP::SiconosMatrix D, unsigned int type):
  CommonSMC(type, sensor, B, D)
{
}

LinearSMC::~LinearSMC()
{
}

void LinearSMC::actuate()
{

  if (!_noUeq)
    computeUeq();

    *(_DS_SMC->x()) = _sensor->y();
    prod(*_B, *_ueq, *(_DS_SMC->b()));
  _simulationSMC->computeOneStep();
//  if (_indx > 0)
  {
    _simulationSMC->nextStep();
  }


  // discontinous part
  *_us = *_lambda;
  *_u = *_us;
  *_u += *_ueq;
  _indx++;

}

void LinearSMC::setD(const SiconosMatrix& D)
{
  _D.reset(new SimpleMatrix(D));
}

AUTO_REGISTER_ACTUATOR(LINEAR_SMC, LinearSMC)
