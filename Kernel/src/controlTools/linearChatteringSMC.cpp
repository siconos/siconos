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

#include "linearChatteringSMC.hpp"
using namespace std;
using namespace ActuatorFactory;

linearChatteringSMC::linearChatteringSMC(int name, SP::TimeDiscretisation t, SP::Model m): commonSMC(name, t, m)
{
}

linearChatteringSMC::linearChatteringSMC(int name, SP::TimeDiscretisation t, SP::Model m, const Sensors& sensorList): commonSMC(name, t, m)
{
}

linearChatteringSMC::~linearChatteringSMC()
{
}

void linearChatteringSMC::initialize()
{
  commonSMC::initialize();

  // We can only work with FirstOrderNonLinearDS, FirstOrderLinearDS and FirstOrderLinearTIDS
  // XXX see how to check this
  _DS = dynamic_pointer_cast<FirstOrderLinearDS>(_model->nonSmoothDynamicalSystem()->dynamicalSystemNumber(0));
  if (_DS == NULL)
  {
    RuntimeException::selfThrow("The control of nonlinear System is not yet implemented");
  }

  // Get the dimension of the output
  // XXX What if there is more than one sensor ...

  _sensor = dynamic_pointer_cast<controlSensor>(*(_allSensors->begin()));
  if (_sensor == NULL)
  {
    RuntimeException::selfThrow("linearChatteringSMC::initialize - the given sensor is not a controlSensor");
  }
  else
  {
    _nDim = _sensor->getYDim();

    _u.reset(new SimpleVector(_nDim, 0));

    // XXX really stupid stuff
    _DS->setzPtr(_u);
  }
  _indx = 0;
  _initDone = true;
}

void linearChatteringSMC::actuate()
{
  _s = inner_prod(*_Csurface, *(_sensor->y()));
  if (_s > 0)
    (*_u)(_nDim - 1) = -2;
  else if (_s < 0)
    (*_u)(_nDim - 1) = 2;
  else
    (*_u)(_nDim - 1) = 0;

  _indx++;
}

void linearChatteringSMC::setCsurface(const SimpleVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != _nDim)
  {
    RuntimeException::selfThrow("linearChatteringSMC::setCsurface - inconstency between the size of the dimension of the state space and Csurface");
  }
  else
  {
    if (_Csurface)
    {
      *_Csurface = newValue;
    }
    else
    {
      _Csurface.reset(new SimpleVector(newValue));
    }
  }
}

void linearChatteringSMC::setCsurfacePtr(SP::SimpleVector newPtr)
{
  // check dimensions ...
  if (newPtr->size() != _nDim)
  {
    RuntimeException::selfThrow("linearChatteringSMC::setCsurfacePtr - inconstency between the size of the dimension of the state space and Csurface");
  }
  else
  {
    _Csurface = newPtr;
  }
}

AUTO_REGISTER_ACTUATOR(103, linearChatteringSMC);
