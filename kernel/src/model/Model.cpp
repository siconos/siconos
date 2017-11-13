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

#include "Model.hpp"
#include "SiconosConfig.h"
#include "NonSmoothDynamicalSystem.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "EventDriven.hpp"
#include "Topology.hpp"
#include "EventsManager.hpp"
#include "SiconosGraph.hpp" // For setOfGraph
#include "SimulationTypeDef.hpp"

#include "TypeName.hpp"


// --- CONSTRUCTORS ---

// --- Default (public) constructor ---
Model::Model(): _t(0.0), _t0(0.0), _T(0.0), _title("none"), _author("nobody"), _description("none"),
  _date("none")
{}

// --- From a minimum set of data ---
Model::Model(double newT0, double newT, const std::string& newTitle,
             const std::string& newAuthor, const std::string& newDescription,
             const std::string& newDate):
  _t(newT0), _t0(newT0), _T(-1), _title(newTitle),
  _author(newAuthor), _description(newDescription), _date(newDate)
{
  if (newT > _t0) _T = newT;
  else if (newT > 0 && newT <= _t0)
    RuntimeException::selfThrow
    ("Model::constructor from min data: Warning, final T lower than t0");

  /* empty */
  _nsds.reset(new NonSmoothDynamicalSystem());
  // else no T in the model!
}

Model::~Model()
{
  if (_simulation)
    _simulation->clear();
}

void Model::setSimulation(SP::Simulation newPtr)
{
  // Warning: this function may be used carefully because of the links
  // between Model and TimeDiscretisation The model of the simulation
  // input MUST be the current model.
  _simulation = newPtr;
}

void Model::initialize()
{

  assert(_simulation && "Model::initialize() error - The simulation object of this model is null.");

  // === topology init (computes Interaction sets, relative degrees ...) ===
  _nsds->topology()->initialize();

  // === Simulation init ===
  _simulation->initialize(shared_from_this());

  // symmetry in indexSets
  _nsds->topology()->setProperties();

}



// --- OTHER FUNCTIONS ---

void Model::display() const
{
  std::cout << " =========> Model named " << _title << ", written by " << _author << " (" << _date << ")." <<std::endl;
  std::cout << " ----- Description: " << _description <<std::endl;
  std::cout <<std::endl;
  std::cout << " Time runs from " << _t0 << " to " << _T <<std::endl;
  std::cout << " Current time is " << _t <<std::endl;
  std::cout <<std::endl;
  if (!_nsds) std::cout << "No NSDS linked to the Model" <<std::endl;
  if (_simulation) std::cout << "The simulation (name: " << _simulation->name() << ") is a " << Type::name(*_simulation) << "." <<std::endl;
  else std::cout << "No simulation attached to this model." <<std::endl;
  std::cout <<std::endl;
  std::cout << " ============================" <<std::endl;
}

void Model::setT(double newValue)
{
  _T = newValue;
  _simulation->updateT(newValue);
}
