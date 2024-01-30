/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "ControlSimulation.hpp"

#include "Actuator.hpp"
#include "ControlManager.hpp"
#include "DynamicalSystem.hpp"
#include "Interaction.hpp"
#include "Observer.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "Simulation.hpp"
#include "TimeDiscretisation.hpp"
#include "Topology.hpp"

namespace siconos::control::internal {
static inline std::pair<unsigned, std::string> getNumberOfStates(
    siconos::graphs::DynamicalSystemsGraph& DSG0, siconos::graphs::InteractionsGraph& IG0)
{
  std::string legend;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsvi, dsvdend;
  unsigned nb = 0;
  unsigned counter = 0;
  for (std::tie(dsvi, dsvdend) = DSG0.vertices(); dsvi != dsvdend; ++dsvi) {
    auto& x = *DSG0.bundle(*dsvi)->x();
    nb += x.size();

    std::string nameDS;
    if (DSG0.name.hasKey(*dsvi)) {
      nameDS = DSG0.name[*dsvi];
    }
    else {
      nameDS = "unknownDS" + std::to_string(counter);
      ++counter;
    }

    for (decltype(x.size()) i = 0; i < x.size(); ++i) {
      legend.append(" " + nameDS + "_" + std::to_string(i));
    }

    if (DSG0.u.hasKey(*dsvi)) {
      auto sizeU = DSG0.u[*dsvi]->size();
      nb += sizeU;
      for (decltype(sizeU) i = 0; i < sizeU; ++i) {
        legend.append(" " + nameDS + "_u_" + std::to_string(i));
      }
    }

    if (DSG0.e.hasKey(*dsvi)) {
      auto sizeE = DSG0.e[*dsvi]->size();
      for (decltype(sizeE) i = 0; i < sizeE; ++i) {
        legend.append(" " + nameDS + "_e_" + std::to_string(i));
      }
      nb += DSG0.e[*dsvi]->size();
    }
  }

  siconos::graphs::InteractionsGraph::VIterator ivi, ivdend;
  counter = 0;
  for (std::tie(ivi, ivdend) = IG0.vertices(); ivi != ivdend; ++ivi) {
    std::string nameInter;
    if (IG0.name.hasKey(*ivi)) {
      nameInter = IG0.name[*ivi];
    }
    else {
      nameInter = "unknownInteraction" + std::to_string(counter);
      ++counter;
    }
    auto& y = *IG0.bundle(*ivi)->y(0);
    nb += y.size();
    for (decltype(y.size()) i = 0; i < y.size(); ++i) {
      legend.append(" " + nameInter + "_y_" + std::to_string(i));
    }

    auto& lambda = *IG0.bundle(*ivi)->lambda(0);
    nb += lambda.size();
    for (decltype(lambda.size()) i = 0; i < lambda.size(); ++i) {
      legend.append(" " + nameInter + "_lambda_" + std::to_string(i));
    }
  }

  return std::make_pair(nb, legend);
}

/** store all the states of the graph in a matrix
 * \param indx row index in the matrix
 * \param startColumn the starting column
 * \param DSG0 the graph of DynamicalSystem
 * \param IG0 the graph of Interaction
 * \param data the matrix where to save the data
 * \return the last written column
 */
static inline unsigned storeAllStates(unsigned indx, unsigned startColumn,
                                      siconos::graphs::DynamicalSystemsGraph& DSG0,
                                      siconos::graphs::InteractionsGraph& IG0,
                                      siconos::algebra::SimpleMatrix& data)
{
  siconos::graphs::DynamicalSystemsGraph::VIterator dsvi, dsvdend;
  auto column = startColumn;
  for (std::tie(dsvi, dsvdend) = DSG0.vertices(); dsvi != dsvdend; ++dsvi) {
    auto i = column;
    auto& x = *DSG0.bundle(*dsvi)->x();
    for (decltype(x.size()) j = 0; j < x.size(); ++i, ++j) {
      data(indx, i) = x(j);
    }
    column += x.size();

    if (DSG0.u.hasKey(*dsvi)) {
      auto& u = *DSG0.u[*dsvi];
      for (decltype(u.size()) j = 0; j < u.size(); ++i, ++j) {
        data(indx, i) = u(j);
      }
      column += u.size();
    }

    if (DSG0.e.hasKey(*dsvi)) {
      auto& e = *DSG0.e[*dsvi];
      for (decltype(e.size()) j = 0; j < e.size(); ++i, ++j) {
        data(indx, i) = e(j);
      }
      column += e.size();
    }
  }

  siconos::graphs::InteractionsGraph::VIterator ivi, ivdend;
  for (std::tie(ivi, ivdend) = IG0.vertices(); ivi != ivdend; ++ivi) {
    auto i = column;
    auto& y = *IG0.bundle(*ivi)->y(0);
    for (decltype(y.size()) j = 0; j < y.size(); ++i, ++j) {
      data(indx, i) = y(j);
    }
    column += y.size();

    auto& lambda = *IG0.bundle(*ivi)->lambda(0);
    for (decltype(lambda.size()) j = 0; j < lambda.size(); ++i, ++j) {
      data(indx, i) = lambda(j);
    }
    column += lambda.size();
  }

  return column;
}
}  // namespace siconos::control::internal

siconos::control::ControlSimulation::ControlSimulation(double t0, double T, double h)
    : _t0(t0), _T(T), _h(h)
{
  _nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(_t0, _T);
  _processTD = std::make_shared<siconos::simulation::TimeDiscretisation>(_t0, _h);
}

void siconos::control::ControlSimulation::initialize()
{
  std::pair<unsigned, std::string> res;
  _dataLegend = "time";

  // Simulation part
  _processSimulation->setNonSmoothDynamicalSystemPtr(_nsds);
  // _model->setSimulation(_processSimulation);
  // _model->initialize();
  // Control part
  _CM->initialize(*_nsds);

  // Output
  _N = (unsigned)ceil((_T - _t0) / _h) + 10;  // Number of time steps
  auto& DSG0 = *_nsds->topology()->dSG(0);
  auto& IG0 = *_nsds->topology()->indexSet0();
  res = siconos::control::internal::getNumberOfStates(DSG0, IG0);
  _nDim = res.first;
  _dataLegend += res.second;
  if (!_saveOnlyMainSimulation) {
    // iter over controller and observer
    const auto& allActuators = _CM->getActuators();
    for (auto it : allActuators) {
      if (it->getInternalNSDS()) {
        auto& topo = *(it)->getInternalNSDS()->topology();
        res = siconos::control::internal::getNumberOfStates(*topo.dSG(0), *topo.indexSet0());
        _nDim += res.first;
        _dataLegend += res.second;
      }
    }
    const auto& allObservers = _CM->getObservers();
    for (auto it : allObservers) {
      if (it->getInternalNSDS()) {
        auto& topo = *(it)->getInternalNSDS()->topology();
        res = siconos::control::internal::getNumberOfStates(*topo.dSG(0), *topo.indexSet0());
        _nDim += res.first;
        _dataLegend += res.second;
      }
    }
  }
  _dataM = std::make_shared<siconos::algebra::SimpleMatrix>(_N, _nDim + 1);
  // we save the system state
}

void siconos::control::ControlSimulation::setTheta(unsigned int newTheta)
{
  _theta = newTheta;
}

void siconos::control::ControlSimulation::addDynamicalSystem(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds, const std::string& name)
{
  _nsds->insertDynamicalSystem(ds);
  _processSimulation->associate(_processIntegrator, ds);

  if (!name.empty()) {
    _nsds->setName(ds, name);
  }
}

void siconos::control::ControlSimulation::addSensor(std::shared_ptr<Sensor> sensor,
                                                    const double h)
{
  if (h < _h) {
    THROW_EXCEPTION(
        "The timestep for a sensor cannot be smaller than the one for the simulation");
  }

  auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(_t0, h);
  _CM->addSensorPtr(sensor, td);
}

void siconos::control::ControlSimulation::addActuator(std::shared_ptr<Actuator> actuator,
                                                      const double h)
{
  if (h < _h) {
    THROW_EXCEPTION(
        "The timestep for an actuator cannot be smaller than the one for the simulation");
  }

  auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(_t0, h);
  _CM->addActuatorPtr(actuator, td);
}

void siconos::control::ControlSimulation::addObserver(std::shared_ptr<Observer> observer,
                                                      const double h)
{
  if (h < _h) {
    THROW_EXCEPTION(
        "The timestep for an observer cannot be smaller than the one for the simulation");
  }

  auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(_t0, h);
  _CM->addObserverPtr(observer, td);
}

void siconos::control::ControlSimulation::storeData(unsigned indx)
{
  unsigned startingColumn = 1;
  startingColumn =
      siconos::control::internal::storeAllStates(indx, startingColumn, *_DSG0, *_IG0, *_dataM);

  if (!_saveOnlyMainSimulation) {
    // iter over controller and observer
    const auto& allActuators = _CM->getActuators();
    for (auto it : allActuators) {
      if (it->getInternalNSDS()) {
        auto& topo = *(it)->getInternalNSDS()->topology();
        startingColumn = siconos::control::internal::storeAllStates(
            indx, startingColumn, *topo.dSG(0), *topo.indexSet0(), *_dataM);
      }
    }
    const auto& allObservers = _CM->getObservers();
    for (auto it : allObservers) {
      if (it->getInternalNSDS()) {
        auto& topo = *(it)->getInternalNSDS()->topology();
        startingColumn = siconos::control::internal::storeAllStates(
            indx, startingColumn, *topo.dSG(0), *topo.indexSet0(), *_dataM);
      }
    }
  }
}
