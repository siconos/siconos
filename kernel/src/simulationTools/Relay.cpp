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
#include "Relay.hpp"
#include <iostream>
#include <assert.h>
#include "Tools.hpp"
#include "Simulation.hpp"
#include "RelayNSL.hpp"
#include "OSNSMatrix.hpp"

// --- Numerics headers ---
#include "NonSmoothDrivers.h"
#include <Relay_Solvers.h>

#include <limits>

using namespace RELATION;


Relay::Relay(int numericsSolverId):
  LinearOSNS(numericsSolverId)
{
  _numerics_problem.reset(new RelayProblem);

  relay_setDefaultSolverOptions(NULL, &*_numerics_solver_options, numericsSolverId);


}
/* nslaw dispatch on bounds */

struct Relay::_BoundsNSLEffect : public SiconosVisitor
{

  using SiconosVisitor::visit;

  Relay* _parent;
  SP::Interaction _inter;
  unsigned int _pos;


  _BoundsNSLEffect(Relay *p, SP::Interaction inter, unsigned int pos) :
    _parent(p), _inter(inter), _pos(pos) {};

  void visit(const RelayNSL& nslaw)
  {

    for (unsigned i = 0; i <  _inter->nonSmoothLaw()->size(); ++i)
    {
      (*(_parent->lb()))(_pos + i) = nslaw.lb();
      (*(_parent->ub()))(_pos + i) = nslaw.ub();
    }
  }

  void visit(const ComplementarityConditionNSL& nslaw)
  {
    for (unsigned i = 0; i <  _inter->nonSmoothLaw()->size(); ++i)
    {
      (*(_parent->lb()))(_pos + i) = 0.0;
      (*(_parent->ub()))(_pos + i) = std::numeric_limits<double>::infinity();
    }
  }

};


void Relay::initialize(SP::Simulation sim)
{
  LinearOSNS::initialize(sim);
  //cout << "Relay::initialize" <<std::endl;


  // initialize memory for _lb and _ub
  if (! _lb)
    _lb.reset(new SiconosVector(maxSize()));
  else
  {
    if (_lb->size() != maxSize())
      _lb->resize(maxSize());
  }
  if (! _ub)
    _ub.reset(new SiconosVector(maxSize()));
  else
  {
    if (_ub->size() != maxSize())
      _ub->resize(maxSize());
  }
}


int Relay::compute(double time)
{
  int info = 0;
  // --- Prepare data for Relay computing ---
  bool cont = preCompute(time);
  if (!cont)
    return info;

  // fill _lb and _ub wiht the value of the NonSmooth Law

  InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());

  //cout << " _sizeOutput =" <<_sizeOutput <<std::endl;
  if (_lb->size() != _sizeOutput)
  {
    _lb->resize(_sizeOutput, false);
    _lb->zero();
  }
  if (_ub->size() != _sizeOutput)
  {
    _ub->resize(_sizeOutput, false);
    _ub->zero();
  }

  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet.bundle(*ui);

    // Compute q, this depends on the type of non smooth problem, on
    // the relation type and on the non smooth law
    unsigned int pos = indexSet.properties(*ui).absolute_position;
    SP::SiconosVisitor NSLEffect(new _BoundsNSLEffect(this, inter, pos));
    inter->nonSmoothLaw()->accept(*NSLEffect);
  }

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if (_sizeOutput != 0)
  {
    // The Relay in Numerics format
    RelayProblem numerics_problem;
    numerics_problem.M = &*_M->numericsMatrix();
    numerics_problem.q = _q->getArray();
    numerics_problem.lb = _lb->getArray();
    numerics_problem.ub = _ub->getArray();
    numerics_problem.size = _sizeOutput;

    //int nbSolvers = 1;
    // Call Relay Driver

    //      Relay_display(&numerics_problem);

    info = relay_driver(&numerics_problem, _z->getArray() , _w->getArray() ,
                        &*_numerics_solver_options);

    if (info != 0)
    {
      std::cout << "Warning : Problem in Relay resolution" <<std::endl;
    }

    // --- Recovering of the desired variables from Relay output ---
    postCompute();
  }

  return info;
}

void Relay::setSolverId(int solverId)
{
  // clear previous Solveroptions
  solver_options_delete(_numerics_solver_options.get());
  relay_setDefaultSolverOptions(NULL, _numerics_solver_options.get(), solverId);
}

void Relay::display() const
{
  std::cout << "======= Relay of size " << _sizeOutput << " with: " <<std::endl;
  LinearOSNS::display();
}

Relay::~Relay()
{
  solver_options_delete(&*_numerics_solver_options);
}

