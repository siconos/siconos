/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "LinearSensor.hpp"

#include "SensorFactory.hpp"
#include "SiconosVector.hpp"
#include "DynamicalSystem.hpp"
#include "TimeDiscretisation.hpp"


LinearSensor::LinearSensor(SP::DynamicalSystem ds): ControlSensor(LINEAR_SENSOR, ds), _k(0), _nSteps(0)
{}

LinearSensor::LinearSensor(SP::DynamicalSystem ds, SP::SimpleMatrix matC, SP::SimpleMatrix matD):
  ControlSensor(LINEAR_SENSOR, ds), _k(0), _matC(matC), _matD(matD), _nSteps(0)
{}

LinearSensor::~LinearSensor()
{
}

void LinearSensor::initialize(const NonSmoothDynamicalSystem& nsds)
{
  // Call initialize of base class
  ControlSensor::initialize(nsds);

  // consistency checks
  if (!_matC)
  {
    RuntimeException::selfThrow("LinearSensor::initialize - no C matrix was given");
  }

  unsigned int colC = _matC->size(1);
  unsigned int rowC = _matC->size(0);
  // What happen here if we have more than one DS ?
  // This may be unlikely to happen.
  //  _DS = _model->nonSmoothDynamicalSystem()->dynamicalSystemNumber(0);
  if (colC != _DS->n())
  {
    RuntimeException::selfThrow(" LinearSensor::initialize - The number of column of the C matrix must be equal to the length of x");
  }
  if (_matD)
  {
    unsigned int rowD = _matD->size(0);
    if (rowC != rowD)
    {
      RuntimeException::selfThrow("C and D must have the same number of rows");
    }
  }

  // --- Get the values ---
  // -> saved in a matrix data
  // -> event
  _storedY.reset(new SiconosVector(rowC));
  //  (_data[_eSensor])["StoredY"] = storedY;
  // set the dimension of the output
  *_storedY = prod(*_matC, *_DSx);
}

void LinearSensor::capture()
{
  *_storedY = prod(*_matC, *_DSx);
  // untested
  if (_matD)
//    *_storedY += prod(*_matD, *_DS->z());
  //  _dataPlot->setSubRow(_k, 1, _storedY);
  _k++;

  if (_delay > 0)
  {
    _bufferY.push_back(_storedY);
  }
}
void  LinearSensor::setC(const SimpleMatrix& C)
{
  *_matC = C;
}

void  LinearSensor::setD(const SimpleMatrix& D)
{
  *_matD = D;
}

AUTO_REGISTER_SENSOR(LINEAR_SENSOR, LinearSensor)
