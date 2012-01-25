/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

#include "LinearChatteringSMC.hpp"
using namespace std;
using namespace ActuatorFactory;

LinearChatteringSMC::LinearChatteringSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds): CommonSMC(103, t, ds)
{
}

LinearChatteringSMC::LinearChatteringSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList): CommonSMC(103, t, ds, sensorList)
{
}

LinearChatteringSMC::~LinearChatteringSMC()
{
}

void LinearChatteringSMC::initialize(SP::Model m)
{
  CommonSMC::initialize(m);

  // We can only work with FirstOrderNonLinearDS, FirstOrderLinearDS and FirstOrderLinearTIDS
  // We can use the Visitor mighty power to check if we have the right type
  Type::Siconos dsType;
  dsType = Type::value(*_DS);
  if (dsType != Type::FirstOrderLinearDS && dsType != Type::FirstOrderLinearTIDS)
  {
    RuntimeException::selfThrow("LinearChatteringSMC::initialize - the control of nonlinear System is not yet implemented");
  }

  // Get the dimension of the output
  // XXX What if there is more than one sensor ...

  _sensor = dynamic_pointer_cast<ControlSensor>(*(_allSensors->begin()));
  if (_sensor == NULL)
  {
    RuntimeException::selfThrow("LinearChatteringSMC::initialize - the given sensor is not a ControlSensor");
  }
  else
  {
    _u.reset(new SimpleVector(_nDim, 0));

    // XXX really stupid stuff
    _DS->setzPtr(_u);
  }
  _indx = 0;
  _initDone = true;
  _s.reset(new SimpleVector(_sDim));
}

void LinearChatteringSMC::actuate()
{
  prod(*_Csurface, *(_sensor->y()), *_s);
  for (unsigned int i = 0; i < _sDim; i++)
  {
    if ((*_s)(i) > 0)
      (*_u)(_nDim - 1) = -2;
    else if ((*_s)(i) < 0)
      (*_u)(_nDim - 1) = 2;
    else
      (*_u)(_nDim - 1) = 0;
  }
  _indx++;
}

AUTO_REGISTER_ACTUATOR(103, LinearChatteringSMC);
