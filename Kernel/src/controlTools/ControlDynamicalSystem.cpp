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
#include "ControlDynamicalSystem.hpp"

using namespace std;

ControlDynamicalSystem::ControlDynamicalSystem(double t0, double T, double h):
  _t0(t0), _T(T), _h(h), _theta(0.5)
{
}

void ControlDynamicalSystem::setTheta(unsigned int newTheta)
{
  _theta = newTheta;
}

void ControlDynamicalSystem::initialize(SP::SiconosVector x0)
{
  _x0 = x0;
  _nDim = _x0->size();

  // Simulation part
  _model.reset(new Model(_t0, _T));
  _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(_processDS);
  _processTD.reset(new TimeDiscretisation(_t0, _h));
  _processSimulation.reset(new TimeStepping(_processTD, 0));
  _processSimulation->setName("plant simulation");
  _processIntegrator.reset(new Moreau(_processDS, _theta));
  _processSimulation->insertIntegrator(_processIntegrator);
  _model->initialize(_processSimulation);

  // Control part
  _CM.reset(new ControlManager(_model));

  // Output
  _N = ceil((_T - _t0) / _h) + 10; // Number of time steps
  _dataM.reset(new SimpleMatrix(_N, 2 * _nDim + 1)); // we save all the states and the control
}


void ControlDynamicalSystem::addSensorPtr(SP::Sensor newSensor)
{
  _CM->addSensorPtr(newSensor);
}

void ControlDynamicalSystem::addActuatorPtr(SP::Actuator newActuator)
{
  _CM->addActuatorPtr(newActuator);
}

void ControlDynamicalSystem::run()
{
  SP::SiconosVector xProc = _processDS->x();
  SP::SiconosVector uProc = _processDS->z();
  (*_dataM)(0, 0) = _t0;
  for (unsigned int i = 0; i < _nDim; i++)
  {
    (*_dataM)(0, i + 1) = (*xProc)(i);
    (*_dataM)(0, _nDim + i + 1) = (*uProc)(i);
  }

  SP::EventsManager eventsManager = _processSimulation->eventsManager();
  unsigned int k = 0;
  boost::progress_display show_progress(_N);
  boost::timer time;
  time.restart();
  SP::Event nextEvent;

  while (_processSimulation->nextTime() < _T)
  {
    _processSimulation->computeOneStep();
    nextEvent = eventsManager->followingEvent(eventsManager->currentEvent());
    if (nextEvent->getType() == 1)
    {
      k++;
      (*_dataM)(k, 0) = _processSimulation->nextTime();
      for (unsigned int i = 0; i < _nDim; i++)
      {
        (*_dataM)(k, i + 1) = (*xProc)(i);
        (*_dataM)(k, _nDim + i + 1) = (*uProc)(i);
      }
      ++show_progress;
    }
    _processSimulation->nextStep();
  }

  _elapsedTime = time.elapsed();
  _dataM->resize(k, 2 * _nDim + 1);
}
