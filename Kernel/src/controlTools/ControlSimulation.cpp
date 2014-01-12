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
#include "TimeDiscretisation.hpp"
#include "ModelingTools.hpp"
#include "SimulationTools.hpp"
#include "ControlManager.hpp"
#include "Sensor.hpp"
#include "Actuator.hpp"
#include "ControlSimulation.hpp"
#include <boost/progress.hpp>
#include <boost/timer.hpp>



ControlSimulation::ControlSimulation(double t0, double T, double h, SP::SiconosVector x0):
  _t0(t0), _T(T), _h(h), _theta(0.5), _elapsedTime(0.0), _N(0), _nDim(0), _x0(x0)
{
}

void ControlSimulation::setTheta(unsigned int newTheta)
{
  _theta = newTheta;
}

void ControlSimulation::initialize()
{
  _nDim = _x0->size();

  // Simulation part
  _model.reset(new Model(_t0, _T));
  _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(_processDS);
  _processTD.reset(new TimeDiscretisation(_t0, _h));
  _processSimulation.reset(new TimeStepping(_processTD, 0));
  _processSimulation->setName("plant simulation");
  _processIntegrator.reset(new ZeroOrderHoldOSI(_processDS));
  _processSimulation->insertIntegrator(_processIntegrator);
  _model->initialize(_processSimulation);

  // Control part
  _CM.reset(new ControlManager(_processSimulation));

  // Output
  _N = ceil((_T - _t0) / _h) + 10; // Number of time steps
  _dataM.reset(new SimpleMatrix(_N, _nDim + 1)); // we save the system state 
}


void ControlSimulation::addSensorPtr(SP::Sensor newSensor, SP::TimeDiscretisation td)
{
  _CM->addSensorPtr(newSensor, td);
}

void ControlSimulation::addActuatorPtr(SP::Actuator newActuator, SP::TimeDiscretisation td)
{
  _CM->addActuatorPtr(newActuator, td);
}

void ControlSimulation::run()
{
  SP::SiconosVector xProc = _processDS->x();
  (*_dataM)(0, 0) = _t0;
  for (unsigned int i = 0; i < _nDim; i++)
  {
    (*_dataM)(0, i + 1) = (*xProc)(i);
  }

  SP::EventsManager eventsManager = _processSimulation->eventsManager();
  unsigned int k = 0;
  boost::progress_display show_progress(_N);
  boost::timer time;
  time.restart();
  SP::Event nextEvent;

  while (_processSimulation->hasNextEvent())
  {
    nextEvent = eventsManager->nextEvent();
    if (nextEvent->getType() == TD_EVENT)
    {
      _processSimulation->computeOneStep();
      k++;
      (*_dataM)(k, 0) = _processSimulation->nextTime();
      for (unsigned int i = 0; i < _nDim; i++)
      {
        (*_dataM)(k, i + 1) = (*xProc)(i);
      }
      ++show_progress;
    }
    _processSimulation->nextStep();
  }

  _elapsedTime = time.elapsed();
  _dataM->resize(k, _nDim + 1);
}
