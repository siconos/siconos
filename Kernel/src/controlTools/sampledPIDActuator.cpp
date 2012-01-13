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

#include "sampledPIDActuator.hpp"
using namespace std;
using namespace ActuatorFactory;

sampledPIDActuator::sampledPIDActuator(int name, SP::TimeDiscretisation t, SP::Model m): Actuator(name, t, m)
{
}

sampledPIDActuator::sampledPIDActuator(int name, SP::TimeDiscretisation t, SP::Model m, const Sensors& sensorList): Actuator(name, t, m)
{
}

sampledPIDActuator::~sampledPIDActuator()
{
}

void sampledPIDActuator::initialize()
{
  Actuator::initialize();

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
    RuntimeException::selfThrow("sampledPIDActuator::initialize - the given sensor is not a controlSensor");
  }
  else
  {
    _nDim = _DS->getN();

    // initialize _err
    _err.reset(new boost::circular_buffer<double> (3));
    for (unsigned int i = 0; i < 3; i++)
      (*_err).push_front(0.0);

    _u.reset(new SimpleVector(_nDim, 0));
    // to have immediatly the right SP
    _DS->setb(_u);
  }
}

void sampledPIDActuator::actuate()
{
  // We have to distinguish two cases : linear or nonlinear
  // TODO support the nonlinear case

  // Get DeltaT
  _curDeltaT = _timeDiscretisation->currentTimeStep();

  // Compute the new error

  (*_err).push_front(_ref - (*_sensor->y())(0));

  // compute the new control and update it
  (*_u)(1) += ((*_K)(0) + (*_K)(2) / _curDeltaT + (*_K)(1) * _curDeltaT) * (*_err)[0] +
              (-(*_K)(0) - 2 * (*_K)(2) / _curDeltaT) * (*_err)[1] + (*_K)(2) / _curDeltaT * (*_err)[2];
}

void sampledPIDActuator::setK(const SimpleVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != 3)
  {
    RuntimeException::selfThrow("sampledPIDActuator::setK - the size of K is not 3");
  }
  else
  {
    if (_K)
    {
      *_K = newValue;
    }
    else
    {
      _K.reset(new SimpleVector(newValue));
    }
  }
}

void sampledPIDActuator::setKPtr(SP::SimpleVector newPtr)
{
  // check dimensions ...
  if (newPtr->size() != 3)
  {
    RuntimeException::selfThrow("sampledPIDActuator::setKPtr - the size of K is not 3");
  }
  else
  {
    _K = newPtr;
  }
}

sampledPIDActuator* sampledPIDActuator::convert(Actuator* s)
{
  return dynamic_cast<sampledPIDActuator*>(s);
}

AUTO_REGISTER_ACTUATOR(100, sampledPIDActuator);
