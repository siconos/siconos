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
  if (_strat)
    _strat->clear();
}

void Model::setSimulationPtr(SP::Simulation newPtr)
{
  // Warning: this function may be used carefully because of the links
  // between Model and TimeDiscretisation The model of the simulation
  // input MUST be the current model.
  _strat = newPtr;
}

void Model::setNonSmoothDynamicalSystemPtr(SP::NonSmoothDynamicalSystem newPtr)
{
  _nsds = newPtr;
}

void Model::initialize(SP::Simulation simulation)
{
  _strat = simulation;

  assert(_strat && "Model::initialize() error - The simulation object of this model is null.");

  // === topology init (computes Interaction sets, relative degrees ...) ===
  _nsds->topology()->initialize();

  // === Simulation init ===
  _strat->initialize(shared_from_this());

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
  if (_strat) std::cout << "The simulation (name: " << _strat->name() << ") is a " << Type::name(*_strat) << "." <<std::endl;
  else std::cout << "No simulation attached to this model." <<std::endl;
  std::cout <<std::endl;
  std::cout << " ============================" <<std::endl;
}

void Model::setT(double newValue)
{
  _T = newValue;
  _strat->updateT(newValue);
}
