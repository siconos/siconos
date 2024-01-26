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
#include "Relay.hpp"

#include "NumericsSolversNamespace.h"  // solver_options stuff
#include "SiconosVisitor.hpp"

// #include <iostream>
// #include <assert.h>
#include "Interaction.hpp"
#include "OSNSMatrix.hpp"
#include "RelayNSL.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "Tools.hpp"

// // --- Numerics headers ---
// #include "NonSmoothDrivers.h"
// #include <Relay_Solvers.h>

// #include <limits>

siconos::nonsmooth_formulations::Relay::Relay(int numericsSolverId)
    : Relay(std::shared_ptr<SolverOptions>(
          solver_options_create(numericsSolverId),
          solver_options_delete))
{
}

siconos::nonsmooth_formulations::Relay::Relay(std::shared_ptr<SolverOptions> options)
    : LinearOSNS(
          options)  //,  _numerics_problem(std::make_shared<RelayProblem>())
{
}

/* nslaw dispatch on bounds */

struct siconos::nonsmooth_formulations::Relay::_BoundsNSLEffect
    : public siconos::internal::SiconosVisitor {
  using siconos::internal::SiconosVisitor::visit;

  Relay* _parent{nullptr};
  std::shared_ptr<siconos::modeling::Interaction> _inter{nullptr};
  unsigned int _pos;

  _BoundsNSLEffect(Relay* p, std::shared_ptr<siconos::modeling::Interaction> inter,
                   unsigned int pos)
      : _parent(p), _inter(inter), _pos(pos){};

  void visit(const siconos::modeling::RelayNSL& nslaw) const override
  {
    for (unsigned i = 0; i < _inter->nonSmoothLaw()->size(); ++i) {
      (*(_parent->lb()))(_pos + i) = nslaw.lb();
      (*(_parent->ub()))(_pos + i) = nslaw.ub();
    }
  }

  void visit(const siconos::modeling::ComplementarityConditionNSL& nslaw) const override
  {
    for (unsigned i = 0; i < _inter->nonSmoothLaw()->size(); ++i) {
      (*(_parent->lb()))(_pos + i) = 0.0;
      (*(_parent->ub()))(_pos + i) = std::numeric_limits<double>::infinity();
    }
  }
};

void siconos::nonsmooth_formulations::Relay::initialize(
    std::shared_ptr<siconos::simulation::Simulation> sim)
{
  LinearOSNS::initialize(sim);
  // cout << "siconos::nonsmooth_formulations::Relay::initialize" <<std::endl;

  // initialize memory for _lb and _ub
  if (!_lb)
    _lb = std::make_shared<siconos::algebra::SiconosVector>(maxSize());
  else {
    if (_lb->size() != maxSize()) _lb->resize(maxSize());
  }
  if (!_ub)
    _ub = std::make_shared<siconos::algebra::SiconosVector>(maxSize());
  else {
    if (_ub->size() != maxSize()) _ub->resize(maxSize());
  }
}
bool siconos::nonsmooth_formulations::Relay::checkCompatibleNSLaw(siconos::modeling::NonSmoothLaw& nslaw)
{
  float type_number = static_cast<float>(siconos::types::type_value(nslaw));
  _nslawtype.insert(type_number);

  if (not(siconos::types::type_value(nslaw) ==
              siconos::modeling::Type::ComplementarityConditionNSL ||
          siconos::types::type_value(nslaw) == siconos::modeling::Type::RelayNSL)) {
    THROW_EXCEPTION(
        "\nsiconos::nonsmooth_formulations::Relay::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with Relay one step nonsmooth problem. \n \
                      Compatible siconos::modeling::NonSmoothLaw are: ComplementarityConditionNSL or RelayNSL\n");
    return false;
  }

  return true;
}

int siconos::nonsmooth_formulations::Relay::compute(double time)
{
  int info = 0;
  // --- Prepare data for Relay computing ---
  bool cont = preCompute(time);
  if (!cont) return info;

  // fill _lb and _ub wiht the value of the NonSmooth Law

  auto& indexSet = *simulation()->indexSet(indexSetLevel());

  // cout << " _sizeOutput =" <<_sizeOutput <<std::endl;
  if (_lb->size() != _sizeOutput) {
    _lb->resize(_sizeOutput, false);
    _lb->zero();
  }
  if (_ub->size() != _sizeOutput) {
    _ub->resize(_sizeOutput, false);
    _ub->zero();
  }

  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui) {
    auto inter = indexSet.bundle(*ui);

    // Compute q, this depends on the type of non smooth problem, on
    // the relation type and on the non smooth law
    auto pos = indexSet.properties(*ui).absolute_position;
    auto NSLEffect = std::make_shared<_BoundsNSLEffect>(this, inter, pos);
    inter->nonSmoothLaw()->accept(*NSLEffect);
  }

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if (_sizeOutput != 0) {
    // The Relay in Numerics format
    RelayProblem numerics_problem;
    numerics_problem.M = &*_M->numericsMatrix();
    numerics_problem.q = _q->getArray();
    numerics_problem.lb = _lb->getArray();
    numerics_problem.ub = _ub->getArray();
    numerics_problem.size = _sizeOutput;

    // int nbSolvers = 1;
    //  Call Relay Driver

    //      Relay_display(&numerics_problem);

    info = relay_driver(&numerics_problem, _z->getArray(), _w->getArray(),
                        &*_numerics_solver_options);

    if (info != 0) {
      std::cout << "Warning : Problem in Relay resolution" << std::endl;
    }

    // --- Recovering of the desired variables from Relay output ---
    postCompute();
  }

  return info;
}

void siconos::nonsmooth_formulations::Relay::display() const
{
  std::cout << "======= Relay of size " << _sizeOutput << " with: " << std::endl;
  LinearOSNS::display();
  std::cout << "lower bound : (_lb)" << std::endl;
  _lb->display();
  std::cout << "upper bound : (_ub)" << std::endl;
  _ub->display();
}
