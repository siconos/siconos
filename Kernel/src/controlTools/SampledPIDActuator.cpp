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

#include "SampledPIDActuator.hpp"
using namespace std;
using namespace ActuatorFactory;

SampledPIDActuator::SampledPIDActuator(SP::TimeDiscretisation t, SP::DynamicalSystem ds): Actuator(SAMPLED_PID_ACTUATOR, t, ds)
{
}

SampledPIDActuator::SampledPIDActuator(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList): Actuator(SAMPLED_PID_ACTUATOR, t, ds, sensorList)
{
}

SampledPIDActuator::~SampledPIDActuator()
{
}

void SampledPIDActuator::initialize(SP::Model m)
{
  Actuator::initialize(m);

  // We can only work with FirstOrderNonLinearDS, FirstOrderLinearDS and FirstOrderLinearTIDS
  // We can use the Visitor mighty power to check if we have the right type
  Type::Siconos dsType;
  dsType = Type::value(*_DS);
  if (dsType != Type::FirstOrderLinearDS && dsType != Type::FirstOrderLinearTIDS)
  {
    RuntimeException::selfThrow("SampledPIDActuator is not yet implemented for system of type" + dsType);
  }
  else
  {
    SP::FirstOrderLinearDS FOLDS = cpp11ns::static_pointer_cast<FirstOrderLinearDS>(_DS);


    // Get the dimension of the output
    // XXX What if there is more than one sensor ...

    _sensor = cpp11ns::dynamic_pointer_cast<ControlSensor>(*(_allSensors->begin()));
    if (_sensor == NULL)
    {
      RuntimeException::selfThrow("SampledPIDActuator::initialize - the given sensor is not a ControlSensor");
    }
    else
    {
      // initialize _err
      _err.reset(new boost::circular_buffer<double> (3));
      for (unsigned int i = 0; i < 3; i++)
        (*_err).push_front(0.0);

      _u.reset(new SiconosVector(_nDim, 0));
      // to have immediatly the right SP
      FOLDS->setb(_u);
    }
  }
}

void SampledPIDActuator::actuate()
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

void SampledPIDActuator::setK(const SiconosVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != 3)
  {
    RuntimeException::selfThrow("SampledPIDActuator::setK - the size of K should be 3");
  }
  else
  {
    if (_K)
    {
      *_K = newValue;
    }
    else
    {
      _K.reset(new SiconosVector(newValue));
    }
  }
}

void SampledPIDActuator::setKPtr(SP::SiconosVector newPtr)
{
  // check dimensions ...
  if (newPtr->size() != 3)
  {
    RuntimeException::selfThrow("SampledPIDActuator::setKPtr - the size of K should be 3");
  }
  else
  {
    _K = newPtr;
  }
}
AUTO_REGISTER_ACTUATOR(SAMPLED_PID_ACTUATOR, SampledPIDActuator)
