/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "AVI.hpp"

#include "Interaction.hpp"
#include "NonSmoothLaw.hpp"
#include "NormalConeNSL.hpp"
#include "NumericsSolversNamespace.h"  // solver_options stuff
#include "NumericsToolsNamespace.h"    // polyhedron, NM_create ...
#include "OSNSMatrix.hpp"
#include "SiconosVector.hpp"
#include "SiconosVisitor.hpp"
#include "SimpleMatrix.hpp"
#include "Simulation.hpp"
#include "TypeName.hpp"  // check nslaw type, should be replaced by dynamic_cast or variant ?
// START visitor for nslaw
//  NOT USED FOR THE MOMENT
//  #if 0
//  struct siconos::nonsmooth_formulations::AVI::_BoundsNSLEffect : public
//  siconos::internal::SiconosVisitor
//  {

//   AVI* _parent{nullptr};
//   std::shared_ptr<siconos::modeling::Interaction> _inter{nullptr};
//   unsigned int _pos{0};

//   _BoundsNSLEffect(AVI *p, std::shared_ptr<siconos::modeling::Interaction> inter, unsigned
//   int pos) :
//     _parent(p), _inter(inter), _pos(pos) {};

//   void visit(const siconos::modeling::NormalConeNSL& nslaw) const override
//   {
//     if(_pos > 0)
//     {
//       S
//     }
//     // take the
//     auto& K = nslaw.K();
//     auto& H = nslaw.H();
//     _numerics_problem->size = nslaw.size();
//     _numerics_problem->d = nullptr;
//     _numerics_problem->poly->id = SICONOS_SET_POLYHEDRON;
//     _numerics_problem->poly->size_ineq = K.size();
//     _numerics_problem->poly->size_eq = 0;
//     _numerics_problem->poly->H = H.getArray();
//     _numerics_problem->poly->K = K.getArray();
//     _numerics_problem->poly->Heq = nullptr;
//     _numerics_problem->poly->Keq = nullptr;
//   }

//   void visit(const siconos::modeling::RelayNSL& nslaw) const override
//   {
//   }

//   void visit(const siconos::modeling::ComplementarityConditionNSL& nslaw) const override
//   {
//   }

// };
// #endif
// /*****************************************************
//  * END visitor for nslaw
// */

siconos::nonsmooth_formulations::AVI::AVI(int numericsSolverId)
    : AVI(std::shared_ptr<SolverOptions>(
          solver_options_create(numericsSolverId),
          solver_options_delete))
{
}

siconos::nonsmooth_formulations::AVI::AVI(
    std::shared_ptr<SolverOptions> options)
    : LinearOSNS(options),
      _numerics_problem(std::make_shared<AffineVariationalInequalities>())
{
  _numerics_problem->poly.split = new polyhedron;
}

siconos::nonsmooth_formulations::AVI::~AVI() noexcept { delete _numerics_problem->poly.split; }

void siconos::nonsmooth_formulations::AVI::initialize(
    std::shared_ptr<siconos::simulation::Simulation> sim)
{
  LinearOSNS::initialize(sim);

  // right now we support only one (1) NonsmoothLaw associated with this AVI
  // It is not clear whether having multiple NonsmoothLaw would be beneficial given the
  // exponential complexity of most solvers
  // TODO We should support RelayNSL with generic rectangles -- xhub
  auto& indexSet = *simulation()->indexSet(indexSetLevel());
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  unsigned nbInter = 0;
  for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui) {
    auto& nc =
        static_cast<siconos::modeling::NormalConeNSL&>(*indexSet.bundle(*ui)->nonSmoothLaw());
    assert(siconos::types::type_value(nc) == siconos::modeling::Type::NormalConeNSL &&
           "siconos::nonsmooth_formulations::AVI::initialize :: found a NonSmoothLaw that is "
           "not of the "
           "NormalConeNSL type! This is currently not supported");
    auto& K = nc.K();
    auto& H = nc.H();
    _numerics_problem->size = nc.size();
    _numerics_problem->d = nullptr;
    _numerics_problem->poly.split->id = SICONOS_SET_POLYHEDRON;
    _numerics_problem->poly.split->size_ineq = K.size();
    _numerics_problem->poly.split->size_eq = 0;
    _numerics_problem->poly.split->H = NM_create_from_data(
        NM_DENSE, K.size(), nc.size(), H.getArray());
    _numerics_problem->poly.split->K = K.getArray();
    _numerics_problem->poly.split->Heq = nullptr;
    _numerics_problem->poly.split->Keq = nullptr;

    // we do not support more than one interaction
    if (!(nbInter++ == 0))
      THROW_EXCEPTION(
          "siconos::nonsmooth_formulations::AVI::initialize :: more than one Interactions for "
          "this "
          "OneStepNSProblem is not support ATM!");
  }
}
bool siconos::nonsmooth_formulations::AVI::checkCompatibleNSLaw(
    siconos::modeling::NonSmoothLaw& nslaw)
{
  float type_number = static_cast<float>(siconos::types::type_value(nslaw));
  _nslawtype.insert(type_number);

  if (not(siconos::types::type_value(nslaw) == siconos::modeling::Type::NormalConeNSL)) {
    THROW_EXCEPTION(
        "\nsiconos::nonsmooth_formulations::AVI::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with AVVI one step nonsmooth problem. \n \
                      Compatible NonSmoothLaw are: NormalConeNSL\n");
    return false;
  }
  return true;
}

int siconos::nonsmooth_formulations::AVI::compute(double time)
{
  int info = 0;
  // --- Prepare data for AVI computing ---
  bool cont = preCompute(time);
  if (!cont) return info;

  if (_numerics_problem->size != _sizeOutput) {
    THROW_EXCEPTION(
        "siconos::nonsmooth_formulations::AVI::compute - size mismatch between AVI size and "
        "and the "
        "current size");
  }

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if (_sizeOutput != 0) {
    // The AVI in Numerics format
    _numerics_problem->M = _M->numericsMatrix().get();
    _numerics_problem->q = _q->getArray();

    info = avi_driver(_numerics_problem.get(), _z->getArray(),
                                         _w->getArray(), _numerics_solver_options.get());

    if (info != 0) {
      std::cout << "Warning : Problem in AVI resolution\n";
    }

    // --- Recovering of the desired variables from AVI output ---
    postCompute();
  }

  return info;
}

void siconos::nonsmooth_formulations::AVI::display() const
{
  std::cout << "======= AVI of size " << _sizeOutput << " with: \n";
  LinearOSNS::display();
}
