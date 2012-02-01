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
#include "LinearSensor.hpp"
#include "SensorFactory.hpp"
#include "ioMatrix.hpp"
#include "DynamicalSystem.hpp"
#include "Model.hpp"
#include "TimeDiscretisation.hpp"
#include "NonSmoothDynamicalSystem.hpp"

using namespace std;
using namespace SensorFactory;

LinearSensor::LinearSensor(SP::TimeDiscretisation t, SP::DynamicalSystem ds): ControlSensor(100, t, ds)
{}

LinearSensor::LinearSensor(SP::TimeDiscretisation t, SP::DynamicalSystem ds, SP::SimpleMatrix matC, SP::SimpleMatrix matD): ControlSensor(100, t, ds), _matC(matC), _matD(matD)
{}

LinearSensor::~LinearSensor()
{
  cout << "LinearSensor destructor called" << endl;
}

void LinearSensor::initialize(SP::Model m)
{
  // Call initialize of base class
  ControlSensor::initialize(m);

  // consistency checks
  unsigned int colC = _matC->size(1);
  unsigned int rowC = _matC->size(0);
  // What happen here if we have more than one DS ?
  // This may be unlikely to happen.
  //  _DS = _model->nonSmoothDynamicalSystem()->dynamicalSystemNumber(0);
  if (colC != _nDim)
  {
    char err[200];
    snprintf(err, 200, "The number of column of the C matrix must be equal to the length of x\nHere col(C) = %i and length(x) = %i", colC, _nDim);
    RuntimeException::selfThrow(err);
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
  // XXX with a variable TimeStep, we're screwed
  //  unsigned int _nSteps = (_model->finalT() - model()->t0())/_timeDiscretisation->currentTimeStep();
  //  _dataPlot.reset(new SimpleMatrix(_nSteps+10, 1+rowC));
  _k = 0;
  // -> event
  _storedY.reset(new SimpleVector(rowC));
  //  (_data[_eSensor])["StoredY"] = storedY;
  // set the dimension of the output
  _YDim = rowC;
  _DSx = _DS->x();
  *_storedY = prod(*_matC, *_DSx);
}

void LinearSensor::capture()
{
  *_storedY = prod(*_matC, *_DSx);
  //  (*_dataPlot)(_k, 0) = _timeDiscretisation->currentTime();
  //  _dataPlot->setSubRow(_k, 1, _storedY);
  _k++;
}
AUTO_REGISTER_SENSOR(100, LinearSensor);
