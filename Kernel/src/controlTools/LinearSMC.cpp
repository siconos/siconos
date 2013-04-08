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
#include "ZeroOrderHold.hpp"
#include "TimeDiscretisation.hpp"
#include "ActuatorFactory.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "RelayNSL.hpp"

using namespace ActuatorFactory;

LinearSMC::LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, int name):
  CommonSMC(name, t, ds)
{
}

LinearSMC::LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, SP::SiconosMatrix B, SP::SiconosMatrix D, int name):
  CommonSMC(name, t, ds, B, D)
{
}

LinearSMC::LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList, int name):
  CommonSMC(name, t, ds, sensorList)
{
}

LinearSMC::~LinearSMC()
{
}

void LinearSMC::actuate()
{

  if (!_noUeq)
    computeUeq();

    *(_DS_SMC->x()) = *(_sensor->y()); // XXX this is sooo wrong
    prod(*_B, *_ueq, *(_DS_SMC->b()));
  _simulationSMC->computeOneStep();
//  if (_indx > 0)
  {
    _simulationSMC->nextStep();
  }


  // discontinous part
  *_us = *_lambda;
  prod(1.0, *_B, *_us, *_sampledControl);
  prod(1.0, *_B, *_ueq, *_sampledControl, false);
  _indx++;

}

void LinearSMC::setD(const SiconosMatrix& D)
{
  _D.reset(new SimpleMatrix(D));
}

AUTO_REGISTER_ACTUATOR(LINEAR_SMC, LinearSMC)
