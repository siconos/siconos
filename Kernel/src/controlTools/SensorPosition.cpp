/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "SensorPosition.h"
#include "SensorFactory.h"
#include "ioMatrix.h"
#include "DynamicalSystem.h"
#include "Model.h"
#include "TimeDiscretisation.h"
#include "NonSmoothDynamicalSystem.h"

using namespace std;
using namespace SensorFactory;

SensorPosition::SensorPosition(int name, TimeDiscretisation* t): Sensor(name, t), nSteps(2000)
{}

SensorPosition::~SensorPosition()
{
  ioMatrix io("resultSensor.dat", "ascii");
  io.write(*dataPlot, "noDim");
  delete dataPlot;
}

void SensorPosition::initialize()
{
  // Call initialize of base class
  Sensor::initialize();

  // --- Get the values to be plotted ---
  // -> saved in a matrix dataPlot
  unsigned int outputSize = 3;
  dataPlot = new SimpleMatrix(nSteps, outputSize);
  k = 0;
}

void SensorPosition::capture()
{
  (*dataPlot)(k, 0) = timeDiscretisation->getCurrentTime();
  (*dataPlot)(k, 1) = (*model->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0)->getXPtr())(0);
  (*dataPlot)(k, 2) = (*model->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0)->getXPtr())(3);
  k++;
}

SensorPosition* SensorPosition::convert(Sensor* s)
{
  SensorPosition* sp = dynamic_cast<SensorPosition*>(s);
  return sp;
}


AUTO_REGISTER_SENSOR(1, SensorPosition);


