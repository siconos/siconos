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
#include "SensorPosition.hpp"
#include "SensorFactory.hpp"
#include "ioMatrix.hpp"
#include "DynamicalSystem.hpp"
#include "Model.hpp"
#include "TimeDiscretisation.hpp"
#include "NonSmoothDynamicalSystem.hpp"

using namespace std;
using namespace SensorFactory;

SensorPosition::SensorPosition(int name, SP::TimeDiscretisation t, SP::Model m): Sensor(name, t, m), _nSteps(2000)
{}

SensorPosition::~SensorPosition()
{
  ioMatrix io("resultSensor.dat", "ascii");
  io.write(*_dataPlot, "noDim");
}

void SensorPosition::initialize()
{
  // Call initialize of base class
  Sensor::initialize();

  // --- Get the values to be plotted ---
  // -> saved in a matrix dataPlot
  unsigned int outputSize = 3;
  _dataPlot.reset(new SimpleMatrix(_nSteps, outputSize));
  _k = 0;
}

void SensorPosition::capture()
{
  (*_dataPlot)(_k, 0) = _timeDiscretisation->currentTime();
  (*_dataPlot)(_k, 1) = (*_model->nonSmoothDynamicalSystem()->dynamicalSystemNumber(0)->x())(0);
  (*_dataPlot)(_k, 2) = (*_model->nonSmoothDynamicalSystem()->dynamicalSystemNumber(0)->x())(3);
  _k++;
}

SensorPosition* SensorPosition::convert(Sensor* s)
{
  SensorPosition* sp = dynamic_cast<SensorPosition*>(s);
  return sp;
}


AUTO_REGISTER_SENSOR(1, SensorPosition);


