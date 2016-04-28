/* Siconos-Kernel, Copyright INRIA 2005-2015
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
#include "Observer.hpp"
#include "ControlSimulation.hpp"
#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include "Model.hpp"

#include "ControlSimulation_impl.hpp"

ControlSimulation::ControlSimulation(double t0, double T, double h):
  _t0(t0), _T(T), _h(h), _theta(0.5), _elapsedTime(0.0), _N(0), _saveOnlyMainSimulation(false), _silent(false)
{
  _model.reset(new Model(_t0, _T));
  _processTD.reset(new TimeDiscretisation(_t0, _h));

}

void ControlSimulation::initialize()
{
  std::pair<unsigned, std::string> res;
  _dataLegend = "time";

  // Simulation part
  _model->setSimulation(_processSimulation);
  _model->initialize();
  // Control part
  _CM->initialize(*_model);

  // Output
  _N = (unsigned)ceil((_T - _t0) / _h) + 10; // Number of time steps
  DynamicalSystemsGraph& DSG0 = *_model->nonSmoothDynamicalSystem()->topology()->dSG(0);
  InteractionsGraph& IG0 = *_model->nonSmoothDynamicalSystem()->topology()->indexSet0();
  res = getNumberOfStates(DSG0, IG0);
  _nDim = res.first;
  _dataLegend += res.second;
  if (!_saveOnlyMainSimulation)
  {
    // iter over controller and observer
    const Actuators& allActuators = _CM->getActuators();
    for (ActuatorsIterator it = allActuators.begin(); it != allActuators.end(); ++it)
    {
      if ((*it)->getInternalModel())
      {
        Topology& topo = *(*it)->getInternalModel()->nonSmoothDynamicalSystem()->topology();
        res = getNumberOfStates(*topo.dSG(0), *topo.indexSet0());
        _nDim += res.first;
        _dataLegend += res.second;
      }
    }
    const Observers& allObservers = _CM->getObservers();
    for (ObserversIterator it = allObservers.begin(); it != allObservers.end(); ++it)
    {
      if ((*it)->getInternalModel())
      {
        Topology& topo = *(*it)->getInternalModel()->nonSmoothDynamicalSystem()->topology();
        res = getNumberOfStates(*topo.dSG(0), *topo.indexSet0());
        _nDim += res.first;
        _dataLegend += res.second;
      }
    }
  }
  _dataM.reset(new SimpleMatrix(_N, _nDim + 1)); // we save the system state 
}

void ControlSimulation::setTheta(unsigned int newTheta)
{
  _theta = newTheta;
}

void ControlSimulation::addDynamicalSystem(SP::DynamicalSystem ds, const std::string& name)
{
  _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(ds);
  _model->nonSmoothDynamicalSystem()->topology()->setOSI(ds, _processIntegrator);

  if (!name.empty())
  {
    _model->nonSmoothDynamicalSystem()->setName(ds, name);
  }
}

void ControlSimulation::addSensor(SP::Sensor sensor, const double h)
{
  SP::TimeDiscretisation td(new TimeDiscretisation(_t0, h));
  _CM->addSensorPtr(sensor, td);
}

void ControlSimulation::addActuator(SP::Actuator actuator, const double h)
{
  SP::TimeDiscretisation td(new TimeDiscretisation(_t0, h));
  _CM->addActuatorPtr(actuator, td);
}

void ControlSimulation::addObserver(SP::Observer observer, const double h)
{
  SP::TimeDiscretisation td(new TimeDiscretisation(_t0, h));
  _CM->addObserverPtr(observer, td);
}

void ControlSimulation::storeData(unsigned indx)
{
  unsigned startingColumn = 1;
  startingColumn = storeAllStates(indx, startingColumn, *_DSG0, *_IG0, *_dataM);
  if (!_saveOnlyMainSimulation)
  {
    // iter over controller and observer
    const Actuators& allActuators = _CM->getActuators();
    for (ActuatorsIterator it = allActuators.begin(); it != allActuators.end(); ++it)
    {
      if ((*it)->getInternalModel())
      {
        Topology& topo = *(*it)->getInternalModel()->nonSmoothDynamicalSystem()->topology();
        startingColumn = storeAllStates(indx, startingColumn, *topo.dSG(0), *topo.indexSet0(), *_dataM);
      }
    }
    const Observers& allObservers = _CM->getObservers();
    for (ObserversIterator it = allObservers.begin(); it != allObservers.end(); ++it)
    {
      if ((*it)->getInternalModel())
      {
        Topology& topo = *(*it)->getInternalModel()->nonSmoothDynamicalSystem()->topology();
        startingColumn = storeAllStates(indx, startingColumn, *topo.dSG(0), *topo.indexSet0(), *_dataM);
      }
    }
  }
}


